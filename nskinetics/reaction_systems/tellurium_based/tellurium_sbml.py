# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np
import pandas as pd
import tellurium as te

from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from ...exceptions import KineticSimulationError, EventCompilationError
from .events import _regenerate_and_resync

__all__ = ('TelluriumReactionSystem',)

_TIME_FACTORS = {'h': 1.0, 'hr': 1.0, 'min': 60.0, 'm': 60.0, 's': 3600.0, 'sec': 3600.0}
# 'kg/m3' (== 'g/L') is a *mass* concentration; it must live only in
# _MASS_CONC_UNITS. material_indexer checks the molar tuple first, so listing it
# in both would mis-resolve a mass unit to 'imol'.
_MOLAR_CONC_UNITS = ('M', 'mol/L')
_MASS_CONC_UNITS = ('g/L', 'kg/m3', 'kg/m^3')

# Pretty LaTeX renderings for known concentration units, used for the y-axis
# label of plot_kinetic_results. Unknown units fall back to the raw string.
_CONC_UNIT_LABELS = {
    'g/L': r'\mathrm{g} \cdot \mathrm{L}^{-1}',
    'kg/m3': r'\mathrm{kg} \cdot \mathrm{m}^{-3}',
    'kg/m^3': r'\mathrm{kg} \cdot \mathrm{m}^{-3}',
    'M': r'\mathrm{mol} \cdot \mathrm{L}^{-1}',
    'mol/L': r'\mathrm{mol} \cdot \mathrm{L}^{-1}',
}


class TelluriumReactionSystem():
    """
    Wrapper around a Tellurium extended RoadRunner model, adding a Python event
    API, unit-aware value access, and convenience methods to simulate the model
    and plot the results.

    Parameters
    ----------
    te : tellurium.roadrunner.extended_roadrunner.ExtendedRoadRunner
        A Tellurium extended RoadRunner object.
    units : dict, optional
        Mapping with keys ``'time'`` and ``'conc'``. Defaults to
        ``{'time': 'min', 'conc': 'M'}``.
    reset : callable, optional
        Custom reset function ``reset(trs, **kwargs)`` invoked by
        :meth:`reset`. When given it fully replaces the default reset (it should
        call ``trs._te.reset()`` itself if a base reset is wanted). Defaults to
        ``None`` (reset the underlying RoadRunner model). Also settable after
        construction via :attr:`reset_func`. Distinct from the boolean ``reset``
        argument of :meth:`simulate`.

    Notes
    -----
    :meth:`simulate` integrates the model and stores the full trajectory on
    this object (:attr:`results`, :attr:`results_array`, :attr:`results_dict`,
    :attr:`results_col_names`); :meth:`plot_simulation_results` (aliases
    :meth:`plot_time_course`, :meth:`plot_trajectory`) plots the most recent
    run. The reaction system is the source of truth for these full-trajectory
    results: drivers such as :class:`nskinetics.units.NSKBatchReactor` call
    :meth:`simulate` and read the results back off this object rather than
    storing their own copy. The underlying RoadRunner object remains directly
    accessible as :attr:`_te` for lower-level use.
    """
    def __init__(self, te, units=None, reset=None):
        self._te = te
        if units is not None:
            self._units = units
        else:
            self._units = {'time': 'min', 'conc': 'M'}
        self._reset_func = reset
        self.events = []
        self._events_compiled = False
        # Most-recent simulation results (set by simulate); the TRS is the
        # source of truth for full-trajectory kinetic results.
        self.results = None
        self.results_array = None
        self.results_dict = None
        self.results_col_names = None

    @classmethod
    def from_sbml(cls, filepath):
        """Load an SBML file into a Tellurium RoadRunner and wrap it.

        Parameters
        ----------
        filepath : str
            Path to an SBML (``.xml``) file.

        Returns
        -------
        TelluriumReactionSystem
            Wrapper with default units ``{'time': 'min', 'conc': 'M'}``; set
            :attr:`_units` or pass ``units`` to :class:`TelluriumReactionSystem`
            directly if the model uses other units.
        """
        return cls(te.loadSBMLModel(filepath))

    # --- value access -------------------------------------------------------
    def get_value(self, selection):
        """Return the current value of a roadrunner selection (plain name or
        ``'[conc]'``)."""
        return self._te[selection]

    def set_value(self, selection, value):
        """Set a roadrunner selection to ``value``."""
        self._te[selection] = value

    def set_concentration_from_stream(self, var, amount, volume):
        """Set the bracket-stripped model variable ``var`` to ``amount/volume``
        and return the value. Replicates the legacy concentration assignment."""
        plain = var.replace('[', '').replace(']', '')
        value = amount / volume
        self._te[plain] = value
        return value

    # --- units --------------------------------------------------------------
    @property
    def time_factor(self):
        """Hours per model time unit."""
        time = self._units.get('time', '')
        try:
            return _TIME_FACTORS[time.lower()]
        except (AttributeError, KeyError):
            raise KineticSimulationError(
                f'Unrecognized time units {time!r}; '
                f'expected one of {sorted(_TIME_FACTORS)}.')

    @property
    def material_indexer(self):
        """biosteam stream indexer matching the model's concentration units."""
        conc = self._units.get('conc')
        if conc in _MOLAR_CONC_UNITS:
            return 'imol'
        elif conc in _MASS_CONC_UNITS:
            return 'imass'
        raise KineticSimulationError(f'Unrecognized concentration units {conc!r}.')

    def validate_units(self):
        """Raise ``KineticSimulationError`` if time/conc units are unrecognized."""
        time_u = self._units.get('time', '').lower()
        if time_u not in _TIME_FACTORS:
            raise KineticSimulationError(
                f'Unrecognized time units {self._units.get("time")!r}; '
                f'expected one of {sorted(_TIME_FACTORS)}.')
        conc_u = self._units.get('conc')
        if conc_u not in _MOLAR_CONC_UNITS and conc_u not in _MASS_CONC_UNITS:
            raise KineticSimulationError(
                f'Unrecognized concentration units {conc_u!r}.')

    # --- events -------------------------------------------------------------
    def add_event(self, event):
        """Append an :class:`Event` to the pending list (compiled by
        :meth:`compile_events`). Returns the event."""
        self.events.append(event)
        self._events_compiled = False
        return event

    def remove_event(self, name):
        """Remove a pending event by name."""
        self.events = [e for e in self.events if e.name != name]

    def clear_events(self):
        """Drop all pending events."""
        self.events = []

    def compile_events(self):
        """Inject all pending events into the model, then regenerate once.

        Guarded so repeated calls do not double-inject. Raises
        :class:`EventCompilationError` with context on failure.

        Notes
        -----
        Compiling events regenerates the model and resets all model state to
        its SBML-defined origin values as a side effect. Set any custom
        initial conditions or parameter overrides AFTER compiling events, not
        before, or they will be discarded.
        """
        if self._events_compiled:
            return
        r = self._te
        try:
            for event in self.events:
                event.compile(r, force_regenerate=False)
            _regenerate_and_resync(r)
        except EventCompilationError:
            raise
        except Exception as e:
            raise EventCompilationError(
                f'Failed to compile {len(self.events)} event(s): {e}') from e
        self._events_compiled = True

    # --- reset --------------------------------------------------------------
    @property
    def reset_func(self):
        """The custom reset callable ``reset(trs, **kwargs)``, or ``None`` for
        the default (plain ``self._te.reset()``)."""
        return self._reset_func

    @reset_func.setter
    def reset_func(self, func):
        self._reset_func = func

    def reset(self, **kwargs):
        """Reset model state (not structure).

        If a custom reset function was supplied (constructor ``reset=`` or the
        :attr:`reset_func` property), it is called as ``reset(self, **kwargs)``
        and is *fully responsible* for restoring state (typically by calling
        ``self._te.reset()`` itself); the default reset is not additionally
        applied. Otherwise the underlying RoadRunner model is reset. Extra
        kwargs are forwarded to the custom function and ignored by the default,
        so callers may pass bookkeeping flags (e.g. ``reset_spike_cap``).

        Not to be confused with the boolean ``reset`` parameter of
        :meth:`simulate`, which selects *whether* to call this method before
        integrating.
        """
        if self._reset_func is not None:
            self._reset_func(self, **kwargs)
        else:
            self._te.reset()

    # --- simulation ---------------------------------------------------------
    def simulate(self, t0=0.0, t_end=None, n_steps=1000, selections=None,
                 reset=False):
        """Integrate the model and store the results on this object.

        Runs the underlying RoadRunner engine over ``[t0, t_end]`` and records
        the trajectory as :attr:`results` (DataFrame), :attr:`results_array`
        (2-D ``numpy.ndarray``), :attr:`results_dict` (``{column: values}``),
        and :attr:`results_col_names` (list of column names). This object is
        the source of truth for full-trajectory kinetic results; drivers such
        as :class:`nskinetics.units.NSKBatchReactor` call this and read the
        results back off the reaction system.

        Parameters
        ----------
        t0 : float
            Start time, in the model's own time units. Defaults to ``0.0``.
        t_end : float
            End time, in the model's own time units. Required.
        n_steps : int
            Number of output time points. Defaults to ``1000``.
        selections : sequence of str, optional
            RoadRunner selections to record (e.g. ``'time'`` and concentration
            selections like ``'[S]'``). If ``None``, defaults to ``'time'``
            followed by the concentration selection of every floating species.
        reset : bool
            If ``True``, reset model state to its origin before simulating.
            Defaults to ``False`` so callers that set initial conditions or
            parameter overrides just before simulating are not clobbered.

        Returns
        -------
        numpy.ndarray
            The 2-D results array (also stored as :attr:`results_array`).

        Raises
        ------
        KineticSimulationError
            If ``t_end`` is not provided, or the underlying integration fails.
        """
        if t_end is None:
            raise KineticSimulationError(
                'simulate() requires t_end (in the model time units).')
        if reset:
            self.reset()
        if selections is None:
            selections = ['time'] + [f'[{s}]'
                                     for s in self._te.getFloatingSpeciesIds()]
        cols = list(selections)
        try:
            raw = np.array(self._te.simulate(t0, t_end, n_steps, cols))
        except Exception as e:
            raise KineticSimulationError(f'Kinetic simulation failed: {e}') from e
        self.results_col_names = cols
        self.results_array = raw
        self.results_dict = {cols[i]: raw[:, i] for i in range(len(cols))}
        self.results = pd.DataFrame(raw, columns=cols)
        return raw

    # --- plotting -----------------------------------------------------------
    def plot_simulation_results(self, variables=None, labels=None,
                                xlim=None, ylim=None, markers=None,
                                save_fig=False, filename=None, figwidth=3.9,
                                **kwargs):
        """Plot concentration vs. time for the most recent simulation.

        Also available under the aliases :meth:`plot_time_course` and
        :meth:`plot_trajectory`.

        Parameters
        ----------
        variables : sequence of str, optional
            Result columns (selections) to plot. If ``None``, every result
            column whose name contains ``'['`` and ``']'`` (the concentration
            selections) is plotted.
        labels : sequence of str, optional
            Legend labels matching ``variables`` in order. Defaults to the raw
            selection names.
        xlim, ylim : tuple, optional
            Axis limits. Default x is ``(0, t_max)``; default y is
            ``(0, 1.1 * max)`` over the plotted concentrations.
        markers : sequence, optional
            Vertical reference lines. Each item is an x-position, or a
            ``(x, label)`` pair; drawn as dashed gray lines. Defaults to none.
        save_fig : bool
            If ``True``, save the figure to ``filename`` (600 dpi, tight,
            white background). Defaults to ``False``.
        filename : str, optional
            Output path used when ``save_fig`` is ``True``.
        figwidth : float
            Figure width in inches. Defaults to ``3.9``.
        **kwargs
            Named vertical event lines, as ``label=time`` pairs. Each draws a
            labeled dashed line at ``time`` (in the plotted time units) and adds
            it to the legend; underscores in the label become spaces. A value of
            ``None`` is skipped, so a caller may pass an unavailable time
            through harmlessly. Drivers such as
            :class:`nskinetics.units.NSKBatchReactor` use this to mark the
            aeration- and fermentation-end times.

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axes.Axes)

        Raises
        ------
        KineticSimulationError
            If no simulation results are stored yet.
        """
        d = self.results_dict
        if d is None:
            raise KineticSimulationError(
                'No simulation results; call simulate() first.')
        if variables is None:
            variables = [k for k in self.results_col_names
                         if '[' in k and ']' in k]
        if not variables:
            raise KineticSimulationError(
                'No concentration variables to plot; pass `variables` '
                'explicitly.')
        if labels is None:
            labels = list(variables)

        plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
        plt.rcParams['font.size'] = '12'
        n_minor_ticks = 4

        fig = plt.figure()
        fig.set_figwidth(figwidth)
        ax = plt.subplot(111)

        time = d['time']
        for var, label in zip(variables, labels):
            ax.plot(time, d[var], label=label)

        time_u = self._units.get('time', '')
        conc_u = self._units.get('conc', '')
        conc_label = _CONC_UNIT_LABELS.get(conc_u, conc_u)
        ax.set_xlabel(r"$\bfTime$" + f' [{time_u}]')
        ax.set_ylabel(r"$\bfConcentration$" + ' [' + f'${conc_label}$' + ']')

        ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks + 1))
        ax.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks + 1))
        ax.tick_params(axis='x', which='both', direction='inout', top=True,
                       width=0.65)
        ax.tick_params(axis='y', which='both', direction='inout', right=True,
                       width=0.65)

        if xlim is not None:
            ax.set_xlim(xlim)
        else:
            ax.set_xlim((0.0, float(time.max())))
        if ylim is not None:
            ax.set_ylim(ylim)
        else:
            ymax = max(float(np.max(d[var])) for var in variables)
            ax.set_ylim((0.0, 1.1 * ymax))

        # Draw vertical lines within the finalized axes box, then reassert the
        # limits so adding the line collections cannot rescale the view.
        xlim_final = ax.get_xlim()
        ylim_final = ax.get_ylim()
        y0, y1 = ylim_final
        if markers:
            for m in markers:
                x = m[0] if isinstance(m, (tuple, list)) else m
                ax.vlines(x=[x], ymin=y0, ymax=y1, linestyles='dashed',
                          linewidth=1.0, color='gray')
        # Named event lines (e.g. aeration-/fermentation-end times injected by
        # NSKBatchReactor); None values are skipped.
        _event_colors = ('dimgray', 'black', 'saddlebrown', 'teal')
        ci = 0
        for name, x in kwargs.items():
            if x is None:
                continue
            ax.vlines(x=[float(x)], ymin=y0, ymax=y1, linestyles='dashed',
                      linewidth=1.2, color=_event_colors[ci % len(_event_colors)],
                      label=name.replace('_', ' '))
            ci += 1
        ax.set_xlim(xlim_final)
        ax.set_ylim(ylim_final)

        ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left',
                  edgecolor='white')

        if save_fig:
            plt.savefig(filename, transparent=False, facecolor='white',
                        bbox_inches='tight', dpi=600)
        return fig, ax

    # Aliases so both kinetic modelers and TEA practitioners find the plot
    # under a familiar name.
    plot_time_course = plot_simulation_results
    plot_trajectory = plot_simulation_results
