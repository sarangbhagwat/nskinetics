# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Shared time-course plotting core (private).

Used by :meth:`nskinetics.KineticModel.plot_simulation_results` (the most
recent run stored on a model) and :meth:`nskinetics.SimulationResult.plot`
(a stored trajectory), so both draw with the same styling. Kept in its own
module so :mod:`nskinetics.engine.conditions` need not import
:mod:`nskinetics.engine.kinetic_model`.
"""
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from ..exceptions import KineticSimulationError

__all__ = ()

# Pretty LaTeX renderings for known concentration units, used for the y-axis
# label of plot_kinetic_results. Unknown units fall back to the raw string.
_CONC_UNIT_LABELS = {
    'g/L': r'\mathrm{g} \cdot \mathrm{L}^{-1}',
    'kg/m3': r'\mathrm{kg} \cdot \mathrm{m}^{-3}',
    'kg/m^3': r'\mathrm{kg} \cdot \mathrm{m}^{-3}',
    'M': r'\mathrm{mol} \cdot \mathrm{L}^{-1}',
    'mol/L': r'\mathrm{mol} \cdot \mathrm{L}^{-1}',
}

# Renderings for the non-concentration quantities a trajectory may also carry
# (working volume, most commonly), reachable via plot_simulation_results'
# `ylabel_unit`. Merged with the concentration units for lookup.
_OTHER_UNIT_LABELS = {
    'L': r'\mathrm{L}',
    'mL': r'\mathrm{mL}',
    'm3': r'\mathrm{m}^{3}',
    'm^3': r'\mathrm{m}^{3}',
}

_UNIT_LABELS = {**_CONC_UNIT_LABELS, **_OTHER_UNIT_LABELS}

# Colors cycled through for named vertical event lines.
_EVENT_COLORS = ('dimgray', 'black', 'saddlebrown', 'teal')


def _bf(text):
    """Render ``text`` as a bold mathtext run, escaping spaces.

    Mathtext collapses literal spaces, so a multi-word axis label such as
    ``'Working volume'`` needs its spaces escaped to survive ``\\mathbf``.
    (The bare ``\\bf`` switch is not accepted by matplotlib >= 3.11's mathtext
    parser when letters follow it directly, so the braced form is required.)
    """
    return r'$\mathbf{' + str(text).replace(' ', r'\ ') + '}$'


def _style_axes(ax, n_minor_ticks=4):
    """Apply the package's tick styling (minor ticks, in-out ticks on all
    four sides) to ``ax``."""
    ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks + 1))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks + 1))
    ax.tick_params(axis='x', which='both', direction='inout', top=True,
                   width=0.65)
    ax.tick_params(axis='y', which='both', direction='inout', right=True,
                   width=0.65)


def _set_font():
    """Set the package's default font family and size (rcParams)."""
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['font.size'] = '12'


def _plot_time_course(data, col_names, units, variables=None, labels=None,
                      xlim=None, ylim=None, markers=None, save_fig=False,
                      filename=None, figwidth=3.9, ylabel=None,
                      ylabel_unit=None, **kwargs):
    """Plot trajectory columns against ``data['time']``.

    The plotting core shared by
    :meth:`nskinetics.KineticModel.plot_simulation_results` and
    :meth:`nskinetics.SimulationResult.plot`; see the former for the meaning
    of every keyword argument.

    Parameters
    ----------
    data : mapping
        ``{column: 1-D array}`` including ``'time'`` (a dict or a DataFrame).
    col_names : sequence of str
        Column names in order; with ``variables=None`` every name containing
        ``'['`` and ``']'`` (the concentration selections) is plotted.
    units : dict
        ``{'time': ..., 'conc': ...}`` used for the axis labels.

    Returns
    -------
    (matplotlib.figure.Figure, matplotlib.axes.Axes)

    Raises
    ------
    KineticSimulationError
        If there is no concentration variable to plot and ``variables`` is
        not given.
    """
    if variables is None:
        variables = [k for k in col_names if '[' in k and ']' in k]
    if not variables:
        raise KineticSimulationError(
            'No concentration variables to plot; pass `variables` '
            'explicitly.')
    if labels is None:
        labels = list(variables)

    _set_font()

    fig = plt.figure()
    fig.set_figwidth(figwidth)
    ax = plt.subplot(111)

    time = np.asarray(data['time'])
    for var, label in zip(variables, labels):
        ax.plot(time, np.asarray(data[var]), label=label)

    time_u = units.get('time', '')
    y_quantity = 'Concentration' if ylabel is None else ylabel
    y_u = units.get('conc', '') if ylabel_unit is None else ylabel_unit
    y_u_label = _UNIT_LABELS.get(y_u, y_u)
    ax.set_xlabel(_bf('Time') + f' [{time_u}]')
    ax.set_ylabel(_bf(y_quantity) + ' [' + f'${y_u_label}$' + ']')

    _style_axes(ax)

    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim((0.0, float(time.max())))
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ymax = max(float(np.max(np.asarray(data[var]))) for var in variables)
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
    ci = 0
    for name, x in kwargs.items():
        if x is None:
            continue
        ax.vlines(x=[float(x)], ymin=y0, ymax=y1, linestyles='dashed',
                  linewidth=1.2, color=_EVENT_COLORS[ci % len(_EVENT_COLORS)],
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
