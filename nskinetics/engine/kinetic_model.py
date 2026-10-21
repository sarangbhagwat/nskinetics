# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import copy as _copy
import os
import pickle
from collections.abc import Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field, fields

import libsbml
import numpy as np
import pandas as pd
import roadrunner.roadrunner as _rr_module
import tellurium as te

from ..exceptions import KineticSimulationError, EventCompilationError
from ..utils import _sbml_to_table as _sbt
from .events import _regenerate_and_resync
# The plotting core and its unit-label tables live in a private module shared
# with SimulationResult.plot(); the names stay importable from here.
from ._plotting import (_bf, _CONC_UNIT_LABELS, _OTHER_UNIT_LABELS,  # noqa: F401
                        _UNIT_LABELS, _plot_time_course)
from .conditions import Condition, SimulationResult, _thaw

__all__ = ('KineticModel',)

# Reset modes accepted by KineticModel.reset (None resolves to 'custom' when a
# reset_func is set, else to 'baseline').
_RESET_MODES = (None, 'baseline', 'roadrunner', 'custom')

# Instance attributes every KineticModel owns; anything else public found on an
# instance (e.g. the shipped model's ``default_max_n_glu_spikes``) is an "extra"
# attribute that copy() / _ModelBlueprint carry over.
_OWN_PUBLIC_ATTRS = frozenset({'events', 'results_df', 'results_array',
                               'results_dict', 'results_col_names'})

_TIME_FACTORS = {'h': 1.0, 'hr': 1.0, 'min': 60.0, 'm': 60.0, 's': 3600.0, 'sec': 3600.0}
# 'kg/m3' (== 'g/L') is a *mass* concentration; it must live only in
# _MASS_CONC_UNITS. material_indexer checks the molar tuple first, so listing it
# in both would mis-resolve a mass unit to 'imol'.
_MOLAR_CONC_UNITS = ('M', 'mol/L')
_MASS_CONC_UNITS = ('g/L', 'kg/m3', 'kg/m^3')

# RoadRunner selections that can be read but never set: reaction rates (ids of
# reactions), rates of change ("S'") and MCA coefficients (prefixes below).
_MCA_PREFIXES = ('uec(', 'ec(', 'ucc(', 'cc(')

# File suffixes that mark a one-line from_antimony() argument as a path.
_ANTIMONY_SUFFIXES = ('.ant', '.antimony', '.txt', '.sb')


@contextmanager
def _keep_roadrunner_properties():
    """Keep ``r.<id>`` attribute access working on other live models.

    libroadrunner installs one Python property per model id (``r.k``,
    ``r.S``, ...) on the ``roadrunner.roadrunner.RoadRunner`` *class*, and
    every load (``RoadRunner._makeProperties``) first deletes the previous
    model's. A load nskinetics performs on the user's behalf --
    :meth:`KineticModel.copy`, and hence every :func:`nskinetics.run_many`
    and :func:`nskinetics.scan` -- would otherwise silently turn
    ``r.k = ...`` on every other live model into a plain instance attribute
    the model never sees (e.g. the shipped model's reset function). Wrap such
    loads in this context manager: properties the newly loaded model does
    not define itself are put back afterwards (the new model's own win).
    Explicit user loads keep libroadrunner's behaviour; item access
    (``r['k']``) is always safe.
    """
    rr = getattr(_rr_module, 'RoadRunner', None)
    names = list(getattr(rr, '_properties', None) or ())
    saved = {n: rr.__dict__[n] for n in names if n in rr.__dict__}
    try:
        yield
    finally:
        for name, prop in saved.items():
            if name not in rr.__dict__:
                setattr(rr, name, prop)
                rr._properties.append(name)


def _init_base(selection):
    """Return ``x`` for an ``'init(x)'`` or ``'init([x])'`` selection, else
    ``None``."""
    if not (selection.startswith('init(') and selection.endswith(')')):
        return None
    base = selection[5:-1]
    if base.startswith('[') and base.endswith(']'):
        base = base[1:-1]
    return base


def _preserve_ids(preserve):
    """Normalize a :meth:`KineticModel.reset` ``preserve`` argument to a tuple
    (a bare string is one id) and reject ``init(...)`` selections.

    A reset never changes ``init(...)`` values, and writing one back would
    re-initialize the model, reverting every parameter edit; touches no
    state.

    Raises
    ------
    ValueError
        If an entry is an ``init(...)`` selection.
    """
    if isinstance(preserve, str):
        preserve = (preserve,)
    preserve = tuple(preserve)
    inits = [s for s in preserve
             if isinstance(s, str) and _init_base(s) is not None]
    if inits:
        raise ValueError(
            f'preserve takes current-value selections, not {inits}: init(...) '
            'values are not changed by a reset (writing one back would '
            're-initialize the model and revert parameter edits).')
    return preserve


class KineticModel():
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
        Custom reset function ``reset(model, **kwargs)`` invoked by
        :meth:`reset`. When given it fully replaces the default reset (it should
        call ``model.reset(mode='baseline')`` or ``model._te.reset()`` itself
        if a base reset is wanted; see :meth:`reset`). Defaults to
        ``None`` (restore the model's baseline; see :meth:`reset`). Also
        settable after construction via :attr:`reset_func`. Distinct from the
        boolean ``reset`` argument of :meth:`simulate`.

    Notes
    -----
    :meth:`simulate` integrates the model and stores the full trajectory on
    this object (:attr:`results_df`, :attr:`results_array`, :attr:`results_dict`,
    :attr:`results_col_names`); :meth:`plot_simulation_results` (aliases
    :meth:`plot_time_course`, :meth:`plot_trajectory`) plots the most recent
    run. the kinetic model is the source of truth for these full-trajectory
    results: drivers such as :class:`nskinetics.units.NSKBatchReactor` call
    :meth:`simulate` and read the results back off this object rather than
    storing their own copy. The underlying RoadRunner object remains directly
    accessible as :attr:`_te` for lower-level use.

    Without a custom reset function, :meth:`reset` restores the model's
    *baseline*: every floating species, compartment, rate-rule variable and
    event-assigned variable goes back to its ``init(x)`` value, while global
    parameter values the user set (e.g. ``km._te.k = 0.4``) persist. Use
    :meth:`copy` for an independent model rebuilt from the origin SBML, and
    :meth:`overrides` for temporary value changes.

    For reproducible runs, describe each run as a :class:`Condition` and call
    :meth:`run`, which resets, applies the condition, simulates, restores the
    overridden parameters and returns a :class:`SimulationResult`;
    :func:`nskinetics.run_many` and :func:`nskinetics.scan` run many
    conditions on a private copy.
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
        # The origin SBML is read live from the RoadRunner (getSBML(); see
        # _get_origin_sbml); this is only a fallback for wrapped objects that
        # lack getSBML (e.g. test stubs).
        self._origin_sbml = None
        # ((event ids, SBML hash), event-assigned ids) parsed from the SBML;
        # see _sbml_event_targets.
        self._event_targets_cache = None
        # (structure signature, frozenset of RoadRunner keys); see _model_keys.
        self._model_keys_cache = None
        # Most-recent simulation results (set by simulate); the TRS is the
        # source of truth for full-trajectory kinetic results.
        self.results_df = None
        self.results_array = None
        self.results_dict = None
        self.results_col_names = None
        # Values the run that produced the results_* slots applied and then
        # restored ({id: value}; set by run(), cleared by simulate()), so that
        # compute_flux_summary can replay the trajectory with them.
        self._results_overrides = {}

    # --- construction -------------------------------------------------------
    @classmethod
    def from_antimony(cls, text_or_path, units=None, reset=None):
        """Load an Antimony model (text or file) and wrap it.

        Parameters
        ----------
        text_or_path : str or os.PathLike
            Antimony model text, or the path of a file containing it. A string
            without line breaks that names an existing file is read as a path.
        units : dict, optional
            Mapping with keys ``'time'`` and ``'conc'``. Defaults to
            ``{'time': 'min', 'conc': 'M'}``.
        reset : callable, optional
            Custom reset function; see :class:`KineticModel`.

        Returns
        -------
        KineticModel
            Wrapper around a new RoadRunner.

        Raises
        ------
        TypeError
            If ``text_or_path`` is neither a string nor a path.
        FileNotFoundError
            If a path-like argument names no file, or a one-line string that
            looks like a path (no ``;`` or ``=``, and a ``.ant``,
            ``.antimony``, ``.txt`` or ``.sb`` suffix or a path separator)
            names no file.
        Exception
            Raised by Tellurium/Antimony if the text is not valid Antimony.
        """
        if isinstance(text_or_path, os.PathLike):
            path = os.fspath(text_or_path)      # open() raises if missing
        elif isinstance(text_or_path, str):
            path = None
            if '\n' not in text_or_path:
                if os.path.isfile(text_or_path):
                    path = text_or_path
                elif (';' not in text_or_path and '=' not in text_or_path
                      and (text_or_path.lower().endswith(_ANTIMONY_SUFFIXES)
                           or os.sep in text_or_path
                           or '/' in text_or_path)):
                    # Otherwise Antimony reports an unrelated parse error.
                    raise FileNotFoundError(
                        f'No such Antimony file: {text_or_path!r}')
        else:
            raise TypeError('from_antimony() expects Antimony text or a path; '
                            f'got {type(text_or_path).__name__}.')
        if path is not None:
            with open(path, encoding='utf-8') as f:
                text = f.read()
        else:
            text = text_or_path
        return cls(te.loadAntimonyModel(text), units=units, reset=reset)

    @classmethod
    def from_sbml(cls, filepath, units=None, reset=None):
        """Load an SBML file into a Tellurium RoadRunner and wrap it.

        Parameters
        ----------
        filepath : str or os.PathLike
            Path to an SBML (``.xml``) file (an SBML string is also accepted,
            as by ``tellurium.loadSBMLModel``).
        units : dict, optional
            Mapping with keys ``'time'`` and ``'conc'``. Defaults to
            ``{'time': 'min', 'conc': 'M'}``; pass the model's own units if it
            uses others.
        reset : callable, optional
            Custom reset function; see :class:`KineticModel`.

        Returns
        -------
        KineticModel
            Wrapper around a new RoadRunner.

        Raises
        ------
        RuntimeError
            Raised by RoadRunner/Tellurium if the file cannot be opened or is
            not valid SBML.
        """
        r = te.loadSBMLModel(os.fspath(filepath)
                             if isinstance(filepath, os.PathLike) else filepath)
        return cls(r, units=units, reset=reset)

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

    def _model_keys(self):
        """``frozenset(r.model.keys())``, cached on the model's structure.

        A large model has many thousand keys (every MCA coefficient
        combination), which take milliseconds to list, while :meth:`run`
        validates ids on every call. The cache is keyed on the ids of the
        global parameters, compartments, species and reactions, so it
        refreshes if the wrapped RoadRunner is replaced or regenerated with a
        different structure.
        """
        r = self._te
        signature = (tuple(r.getGlobalParameterIds()),
                     tuple(r.getCompartmentIds()),
                     tuple(r.getFloatingSpeciesIds()),
                     tuple(r.getBoundarySpeciesIds()),
                     tuple(r.getReactionIds()))
        cached = getattr(self, '_model_keys_cache', None)
        if cached is not None and cached[0] == signature:
            return cached[1]
        keys = frozenset(r.model.keys())
        self._model_keys_cache = (signature, keys)
        return keys

    def _validate_ids(self, ids, *, settable=True):
        """Raise ``ValueError`` if any id is not a selection of this model.

        Touches no model state, so callers can validate every id before
        mutating anything.

        Parameters
        ----------
        ids : iterable of str
            RoadRunner selections (plain ids such as ``'k'``, concentration
            selections such as ``'[S]'``, or ``'init(x)'``).
        settable : bool
            If ``True`` (default), also reject ids defined by an assignment
            rule (also as ``'[x]'``, ``'init(x)'`` or ``'init([x])'``),
            which RoadRunner cannot set independently; read-only
            selections: reaction ids (rates), rates of change (``"S'"``) and
            MCA coefficients (``'ec(J, S)'`` and the like); and the
            ``init(...)`` writes that permanently rewrite the model (see
            Notes). If ``False``, the ids are checked as *output* selections:
            any model key or ``'time'`` is accepted.

        Raises
        ------
        ValueError
            Listing every unknown (and, if ``settable``, rule-defined,
            read-only or model-rewriting) id.

        Notes
        -----
        In RoadRunner (2.10) writing an ``init(x)`` selection edits the SBML
        document itself, and writing the old value back does not undo the
        edit; later runs, :meth:`copy` and every batch row would see it. Two
        forms are therefore rejected (the others, e.g. ``init(k)`` of a
        plain parameter or ``init(env)`` of a compartment, leave the document
        unchanged):

        * ``init(x)`` / ``init([x])`` of an id defined by an *initial
          assignment*: the write deletes the assignment. Set ``init()`` of
          its inputs instead (``init(k)`` for ``S = 20*k``).
        * ``init()`` of a species in the other form than the SBML declares
          its initial value: ``init(S)`` (an amount) for a species declared
          with an initial concentration (the usual Antimony ``S = 50``)
          redeclares it as an initial amount, and ``init([S])`` for one
          declared with an initial amount (e.g. ``substanceOnly``) the other
          way round, which changes later runs that set a volume. Use the
          declared form (``init([S])`` and ``init(S)`` respectively).
        """
        ids = list(ids)
        if not ids:
            return
        r = self._te
        known = self._model_keys()
        unknown = [i for i in ids
                   if i not in known and (settable or i != 'time')]
        rule_defined, read_only, initial_assigned, wrong_form = [], [], [], []
        if settable:
            assigned = set(r.getAssignmentRuleIds())
            reactions = set(r.getReactionIds())
            with_ia, species, declared = set(), set(), None
            if any(i in known and _init_base(i) is not None for i in ids):
                # Only runs that write init(...) pay for these lookups.
                try:
                    with_ia = set(r.getInitialAssignmentIds())
                    species = (set(r.getFloatingSpeciesIds())
                               | set(r.getBoundarySpeciesIds()))
                except AttributeError:          # stub objects
                    pass
            for i in ids:
                if i not in known:
                    continue
                # 'init(x)', '[x]' and 'init([x])' all name x; setting any of
                # them fails if x is defined by an assignment rule.
                base = i
                if base.startswith('init(') and base.endswith(')'):
                    base = base[5:-1]
                if base.startswith('[') and base.endswith(']'):
                    base = base[1:-1]
                if base in assigned:
                    rule_defined.append(i)
                elif (i in reactions or i.endswith("'")
                      or i.startswith(_MCA_PREFIXES)):
                    read_only.append(i)
                elif _init_base(i) is not None:
                    # Writes that would rewrite the SBML document for good
                    # (see Notes).
                    if base in with_ia:
                        initial_assigned.append(i)
                    elif base in species:
                        if declared is None:
                            declared = self._species_initial_forms()
                        as_conc = i[5:-1].startswith('[')
                        form = declared.get(base)
                        if form == 'concentration' and not as_conc:
                            wrong_form.append(
                                f"{i!r} ({base} is declared with an initial "
                                f"concentration; use 'init([{base}])')")
                        elif form == 'amount' and as_conc:
                            wrong_form.append(
                                f"{i!r} ({base} is declared with an initial "
                                f"amount; use 'init({base})')")
        msgs = []
        if unknown:
            msgs.append(f'unknown model id(s) {unknown}')
        if rule_defined:
            msgs.append(f'id(s) defined by assignment rules (cannot be set) '
                        f'{rule_defined}')
        if read_only:
            msgs.append(f'read-only selection(s) (reaction rates, rates of '
                        f'change or MCA coefficients; cannot be set) '
                        f'{read_only}')
        if initial_assigned:
            msgs.append(f'init() of id(s) defined by an initial assignment '
                        f'{initial_assigned} (writing init() deletes the '
                        f'assignment from the model; set init() of its '
                        f'inputs instead, e.g. init(k))')
        if wrong_form:
            msgs.append('init() of species in the form their SBML does not '
                        'declare (writing it would redeclare the initial '
                        'value for good): ' + ', '.join(wrong_form))
        if msgs:
            raise ValueError('; '.join(msgs) + '.')

    def _species_initial_forms(self):
        """``{species id: 'amount' | 'concentration' | None}``: how the origin
        SBML declares each species' initial value (``None`` if neither is
        set, e.g. for one defined by an initial assignment).

        Cached on the SBML text; touches no model state. Used by
        :meth:`_validate_ids` to reject ``init()`` writes in the other form.
        """
        sbml = self._get_origin_sbml()
        cached = getattr(self, '_species_forms_cache', None)
        if cached is not None and cached[0] == sbml:
            return cached[1]
        forms = {}
        model = (libsbml.readSBMLFromString(sbml).getModel()
                 if sbml else None)
        if model is not None:
            for s in model.getListOfSpecies():
                forms[s.getId()] = (
                    'amount' if s.isSetInitialAmount()
                    else 'concentration' if s.isSetInitialConcentration()
                    else None)
        self._species_forms_cache = (sbml, forms)
        return forms

    @contextmanager
    def overrides(self, mapping=None, **values):
        """Temporarily set model values; restore the previous ones on exit.

        Parameters
        ----------
        mapping : mapping of str to float, optional
            ``{selection: value}`` pairs, e.g. ``{'k': 0.4, '[S]': 20.0}``.
            Values are in the model's own units.
        **values : float
            More ``selection=value`` pairs (plain ids only); merged over
            ``mapping``.

        Yields
        ------
        KineticModel
            This model, with the overrides applied.

        Raises
        ------
        ValueError
            If any id is unknown, defined by an assignment rule or read-only
            (a reaction rate, rate of change or MCA coefficient), or is an
            ``init()`` form that would rewrite the model's SBML for good
            (see Notes); raised on entry, before any value is changed.

        Notes
        -----
        The previous values are read before anything is set and written back
        in a ``finally`` block, so they are restored even if the body raises.
        What is restored is the underlying quantity: a species given as a
        concentration (``'[S]'``) gets its previous *amount* back, so the
        restore does not depend on volumes the body changed. Values the body
        itself changes (e.g. by simulating) are not otherwise touched.

        ``init(x)`` entries are set first. RoadRunner re-initializes the
        whole model whenever an ``init(...)`` value is written: time, state
        and every global parameter return to their initial values, with
        initial assignments re-evaluated. The model's other settings (global
        parameters and boundary species not defined by an initial
        assignment) keep their current values, and on exit every setting is
        restored, as are the integrator's settings (the re-initialization
        resets some, e.g. rk45's ``epsilon``); the state, however, is
        re-initialized on entry and again on exit. As for :meth:`run`, ``init()`` of an id defined by an initial
        assignment (set ``init()`` of its inputs instead) and ``init()`` of a
        species in the other form than its SBML declares (``init([S])`` for
        the usual concentration-declared species) are rejected: RoadRunner
        would rewrite the SBML, and restoring the value would not undo it.

        Examples
        --------
        >>> with km.overrides(k=0.4):          # doctest: +SKIP
        ...     km.simulate(0, 48, 481, reset=True)
        """
        changes = dict(mapping) if mapping is not None else {}
        changes.update(values)
        self._validate_ids(changes)
        r = self._te
        targets, settings = self._restore_targets(changes)
        previous = {t: r[t] for t in targets}
        integrator = (self._integrator_snapshot()
                      if any(_init_base(t) is not None for t in targets)
                      else None)
        try:
            self._write_values(changes.items(),
                               {p: previous[p] for p in settings})
            yield self
        finally:
            # init(...) targets come first: writing one re-initializes the
            # model, so everything else is written back after it.
            for t in targets:
                r[t] = previous[t]
            self._reapply_integrator(integrator)

    def _restore_targets(self, selections):
        """What to read before, and write back after, temporarily setting
        ``selections``; touches no state.

        The underlying quantities rather than the selections as given: a
        species concentration ``'[x]'`` maps to the amount ``x`` (a
        concentration written back after the volume changed would change the
        amount), and an ``init(...)`` selection also brings every setting
        (:meth:`_user_parameter_values`), since writing it re-initializes the
        model (see :meth:`_write_values`).

        Returns
        -------
        (list of str, tuple of str)
            The targets in write-back order (``init(...)`` selections first),
            and the setting ids among them that were added because of an
            ``init(...)`` selection (empty without one).
        """
        r = self._te
        try:
            species = (set(r.getFloatingSpeciesIds())
                       | set(r.getBoundarySpeciesIds()))
        except AttributeError:          # stub objects
            species = set()
        inits, others = [], []
        for sel in selections:
            if _init_base(sel) is not None:
                inits.append(sel)
            elif (sel.startswith('[') and sel.endswith(']')
                  and sel[1:-1] in species):
                others.append(sel[1:-1])
            else:
                others.append(sel)
        settings = tuple(self._user_parameter_values()) if inits else ()
        targets = list(dict.fromkeys([*inits, *others, *settings]))
        return targets, settings

    def _write_values(self, items, settings=None):
        """Write ``(selection, value)`` pairs, ``init(...)`` selections first.

        RoadRunner (2.10) re-initializes the whole model whenever an
        ``init(x)`` selection is written: time, state *and every global
        parameter* return to their initial values, and initial assignments
        are re-evaluated. So the model's settings are written back right after
        the ``init(...)`` selections -- except the ``init(...)`` targets
        themselves and ids defined by an initial assignment, which take their
        new, re-evaluated values -- and so are the integrator's settings,
        which the re-initialization resets for some integrators (see
        :meth:`_integrator_snapshot`). The other selections are then written
        in the given order.

        Parameters
        ----------
        items : iterable of (str, float)
            Selections and values.
        settings : mapping of str to float, optional
            The setting values to write back after the ``init(...)``
            selections; read now (:meth:`_user_parameter_values`) if not
            given. Unused without an ``init(...)`` selection.
        """
        r = self._te
        items = list(items)
        inits = [(s, v) for s, v in items if _init_base(s) is not None]
        if inits:
            if settings is None:
                settings = self._user_parameter_values()
            integrator = self._integrator_snapshot()
            skip = {_init_base(s) for s, _ in inits}
            skip.update(r.getInitialAssignmentIds())
            for sel, value in inits:
                r[sel] = value
            for p, value in settings.items():
                if p not in skip:
                    r[p] = value
            self._reapply_integrator(integrator)
        for sel, value in items:
            if _init_base(sel) is None:
                r[sel] = value

    def _integrator_snapshot(self):
        """The integrator's name and settings, for :meth:`_reapply_integrator`.

        Re-initializing the model (on every ``init(x)`` write) resets the
        settings of some integrators to their defaults in RoadRunner 2.10
        (rk45's ``epsilon``, gillespie's ``seed``; CVODE's survive).

        Returns
        -------
        tuple or None
            ``(name, {setting: value})``, or ``None`` for an object without
            ``getIntegrator`` (stubs).
        """
        try:
            integrator = self._te.getIntegrator()
        except AttributeError:          # stub objects
            return None
        return (integrator.getName(),
                {k: integrator.getValue(k) for k in integrator.getSettings()})

    def _reapply_integrator(self, snapshot):
        """Restore an :meth:`_integrator_snapshot` (only differing values)."""
        if snapshot is not None:
            _configure_integrator(self._te, *snapshot)

    # --- units --------------------------------------------------------------
    @property
    def time_factor(self):
        """Model time units per hour (e.g. ``60.0`` for a model in minutes)."""
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

    def state_selections(self):
        """List the RoadRunner selections that fully capture integrator state.

        The returned selections are everything :func:`compute_flux_summary`
        must write back into the model to re-evaluate a stored trajectory:
        ``'time'``, each settable compartment, every floating species as a
        concentration selection (``[species]``), every rate-rule variable,
        and every variable assigned by an event -- both the :class:`Event`
        objects in ``events`` and events defined in the SBML/Antimony
        source itself. Variables defined by assignment rules are excluded
        because they cannot be set independently.

        Returns
        -------
        list of str
            Ordered, de-duplicated selections, ``'time'`` first.
        """
        r = self._te
        try:
            assigned = set(r.getAssignmentRuleIds())
        except (AttributeError, RuntimeError):
            # older roadrunner builds lack the accessor, and it raises on a
            # model that has not been loaded yet; either way, no assignment
            # rules are known and nothing needs excluding.
            assigned = set()
        floating = set(r.getFloatingSpeciesIds())
        sels = ['time']

        def _add(sel):
            if sel not in sels:
                sels.append(sel)

        for c in r.getCompartmentIds():
            if c not in assigned:
                _add(c)
        for s in r.getFloatingSpeciesIds():
            # roadrunner lists assignment-rule species here too; writing one
            # back raises RuntimeError ('Could not set value for NON conserved
            # moiety floating species'), so they must be filtered out as well.
            if s not in assigned:
                _add(f'[{s}]')
        for v in r.getRateRuleIds():
            if v not in assigned:
                _add(v)
        event_targets = [var for event in self.events for var in event.do]
        event_targets.extend(self._sbml_event_targets())
        for var in event_targets:
            if var in assigned:
                continue
            _add(f'[{var}]' if var in floating else var)
        return sels

    def _sbml_event_targets(self):
        """Ids assigned by events present in the loaded SBML, cached.

        Covers events written in the Antimony/SBML source as well as compiled
        :class:`Event` objects. The cache is keyed on the event ids and the
        SBML text (``RoadRunner.getSBML()``), so it follows every structural
        edit: events added (e.g. by :meth:`compile_events`) or removed, and
        assignments added to or removed from an existing event (e.g.
        ``r.addEventAssignment(...)``).

        Returns
        -------
        tuple of str
            Event-assignment variables in document order, de-duplicated; empty
            for a model without events or an object that is not a loaded
            RoadRunner.
        """
        r = self._te
        try:
            ids = tuple(r.model.getEventIds())
        except AttributeError:
            return ()
        if not ids:
            # Models without events never pay for reading the SBML.
            return ()
        sbml = r.getSBML()
        key = (ids, hash(sbml))
        # getattr: a model unpickled from an older version lacks the slot.
        cached = getattr(self, '_event_targets_cache', None)
        if cached is not None and cached[0] == key:
            return cached[1]
        targets = []
        model = libsbml.readSBMLFromString(sbml).getModel()
        if model is not None:
            for event in model.getListOfEvents():
                for ea in event.getListOfEventAssignments():
                    var = ea.getVariable()
                    if var and var not in targets:
                        targets.append(var)
        targets = tuple(targets)
        self._event_targets_cache = (key, targets)
        return targets

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

        From then on the model's origin SBML (``RoadRunner.getSBML()``, which
        :meth:`copy` rebuilds from) includes the compiled events.
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
        self._event_targets_cache = None

    def _get_origin_sbml(self):
        """Return the model's origin SBML: the live ``RoadRunner.getSBML()``.

        ``getSBML()`` is the structural document RoadRunner keeps and edits
        (``addEvent`` and the other structural methods): the loaded model plus
        every structural edit -- events compiled by :meth:`compile_events` or
        directly by ``Event.compile`` / ``FeedSpike.compile`` -- and values
        written to ``init(x)`` selections, with initial assignments intact.
        Unlike ``getCurrentSBML()`` it never bakes run state or parameter
        edits (``r.k = 0.3``) into the document. It is read afresh on every
        call (a cached copy could go stale, e.g. after a structural edit or a
        replaced RoadRunner), so :meth:`copy` and the :attr:`parameters` /
        :attr:`species` tables always describe the model :meth:`run` runs.
        Reading it touches no model state.

        Returns
        -------
        str or None
            The SBML text; for a wrapped object without ``getSBML`` (e.g. a
            test stub), the ``_origin_sbml`` fallback (``None`` unless set).
        """
        get_sbml = getattr(self._te, 'getSBML', None)
        if get_sbml is None:
            return getattr(self, '_origin_sbml', None)
        return get_sbml()

    # --- reset --------------------------------------------------------------
    @property
    def reset_func(self):
        """The custom reset callable ``reset(model, **kwargs)``, or ``None`` for
        the default (baseline restore; see :meth:`reset`)."""
        return self._reset_func

    @reset_func.setter
    def reset_func(self, func):
        self._reset_func = func

    @property
    def baseline_variables(self):
        """Non-species variables a baseline reset restores from ``init(x)``.

        These are :meth:`state_selections` minus ``'time'`` and the
        floating-species concentration selections (which
        ``RoadRunner.reset()`` already restores): settable compartments,
        rate-rule variables and event-assigned non-species variables.

        Returns
        -------
        tuple of str
            Empty when the wrapped object does not expose RoadRunner's
            accessors (e.g. a test stub).
        """
        r = self._te
        if not all(hasattr(r, accessor) for accessor in (
                'getCompartmentIds', 'getFloatingSpeciesIds',
                'getRateRuleIds')):
            return ()
        # Unguarded: any other error (e.g. a broken Event object) must surface
        # rather than silently degrading the reset to the plain one.
        sels = self.state_selections()
        return tuple(s for s in sels
                     if s != 'time' and not (s.startswith('[') and s.endswith(']')))

    def _resolve_reset_mode(self, mode):
        """Validate a reset ``mode`` and resolve ``None``; touches no state."""
        if mode not in _RESET_MODES:
            raise ValueError(f'Unknown reset mode {mode!r}; expected one of '
                             f'{_RESET_MODES}.')
        if mode is None:
            return 'custom' if self._reset_func is not None else 'baseline'
        if mode == 'custom' and self._reset_func is None:
            raise ValueError("reset(mode='custom') requires a reset_func.")
        return mode

    def _check_reset_kwargs(self, kwargs):
        """Validate :meth:`reset` keyword arguments without touching state.

        Checks ``mode`` (as :meth:`reset` does) and that every ``preserve``
        id is a settable model id and not an ``init(...)`` selection; other
        keywords belong to a custom reset function and are not checked.

        Raises
        ------
        ValueError
            For an unknown or unusable ``mode`` or a bad ``preserve`` id.
        """
        self._resolve_reset_mode(kwargs.get('mode'))
        preserve = _preserve_ids(kwargs.get('preserve', ()))
        self._validate_ids(preserve)

    def reset(self, mode=None, preserve=(), **kwargs):
        """Reset model state (not structure).

        Parameters
        ----------
        mode : {None, 'baseline', 'roadrunner', 'custom'}
            How to reset:

            * ``None`` (default): ``'custom'`` when a :attr:`reset_func` is
              set (behaviour unchanged), otherwise ``'baseline'``.
            * ``'baseline'``: every :attr:`baseline_variables` entry ``v`` is
              set back to ``init(v)``, then ``self._te.reset()`` restores
              floating species, rate-rule variables and time (now against the
              origin compartment volumes), then the baseline variables are
              set to ``init(v)`` once more. Compartments and event-assigned
              variables (spike counters, working volume, flags) therefore start
              every run from their origin values, while global parameter values
              the user set (e.g. ``km._te.k = 0.4``) persist, since roadrunner
              leaves ``init(k)`` untouched when ``k`` is set.
            * ``'roadrunner'``: the plain ``self._te.reset()`` (the pre-baseline
              default, which leaves compartments and event-assigned parameters
              at their post-run values).
            * ``'custom'``: ``reset_func(self, **kwargs)``, which is *fully
              responsible* for restoring state: typically it calls
              ``self.reset(mode='baseline')`` and then does its extra work,
              or calls ``self._te.reset()`` itself, setting compartment
              volumes and event-assigned variables before that call as well
              as after it (RoadRunner's reset converts initial
              concentrations to amounts at the *current* volumes).
        preserve : sequence of str
            Selections read just before the reset and written back right after
            it, e.g. ``('n_max',)`` to keep a modified spike cap. Current
            values only: ``init(...)`` selections are rejected (a reset never
            changes them, and writing one back would re-initialize the model,
            reverting parameter edits). Never forwarded to ``reset_func``.
        **kwargs
            Forwarded to a custom reset function and ignored otherwise, so
            callers may pass bookkeeping flags (e.g. ``reset_spike_cap``).

        Raises
        ------
        ValueError
            If ``mode`` is unknown, ``mode='custom'`` without a
            :attr:`reset_func`, or a ``preserve`` id is an ``init(...)``
            selection or not a settable model selection (unknown, defined by
            an assignment rule, or read-only); raised before any state is
            changed.

        Notes
        -----
        If the wrapped object lacks RoadRunner's accessors (e.g. a test stub),
        ``'baseline'`` degrades to the plain ``self._te.reset()``.

        Not to be confused with the boolean ``reset`` parameter of
        :meth:`simulate`, which selects *whether* to call this method before
        integrating.
        """
        mode = self._resolve_reset_mode(mode)
        preserve = _preserve_ids(preserve)
        r = self._te
        if preserve and hasattr(r, 'model'):
            # Unknown, rule-defined or read-only ids could be read but not
            # written back after the reset; reject them up front.
            self._validate_ids(preserve)
        kept = {}
        for sel in preserve:
            try:
                kept[sel] = r[sel]
            except RuntimeError as e:
                raise ValueError(f'Cannot preserve {sel!r}: {e}') from e
        try:
            if mode == 'custom':
                self._reset_func(self, **kwargs)
            elif mode == 'roadrunner':
                r.reset()
            else:
                variables = self.baseline_variables
                # Restore the baseline variables BEFORE RoadRunner's reset as
                # well as after it: reset() converts species initial values to
                # amounts with the *current* compartment volumes, so resetting
                # under a post-run volume gives wrong initial concentrations
                # (an initial assignment is off by the volume ratio; a plain
                # initial concentration by an ulp, which CVODE amplifies to
                # ~1e-4 relative across events). Restoring first makes the
                # reset bit-identical to a freshly loaded model.
                for v in variables:
                    r[v] = r[f'init({v})']
                r.reset()
                for v in variables:
                    r[v] = r[f'init({v})']
        finally:
            for sel, value in kept.items():
                r[sel] = value

    # --- simulation ---------------------------------------------------------
    @staticmethod
    def _check_times(times):
        """Validate an explicit output grid; return it as a list of floats.

        Raises
        ------
        TypeError
            If ``times`` is a string or mapping, or holds values that are not
            numbers (as :class:`Condition` does).
        ValueError
            If ``times`` is not a 1-D sequence of at least two finite, strictly
            increasing numbers.
        """
        if isinstance(times, (str, bytes, Mapping)):
            raise TypeError('times must be a sequence of numbers.')
        try:
            grid = np.asarray(times, dtype=float)
        except (TypeError, ValueError) as e:
            raise TypeError(f'times must be a sequence of numbers: {e}') from e
        if grid.ndim != 1 or grid.size < 2:
            raise ValueError('times must be a 1-D sequence of at least two '
                             'output times.')
        if not np.all(np.isfinite(grid)):
            raise ValueError('times must be finite.')
        if not np.all(np.diff(grid) > 0):
            raise ValueError('times must be strictly increasing.')
        return [float(t) for t in grid]

    def simulate(self, t0=None, t_end=None, n_steps=1000, selections=None,
                 reset=False, times=None):
        """Integrate the model and store the results on this object.

        Runs the underlying RoadRunner engine over ``[t0, t_end]`` (or over an
        explicit ``times`` grid) and records the trajectory as
        :attr:`results_df` (DataFrame), :attr:`results_array` (2-D
        ``numpy.ndarray``), :attr:`results_dict` (``{column: values}``), and
        :attr:`results_col_names` (list of column names). This object is the
        source of truth for full-trajectory kinetic results; drivers such as
        :class:`nskinetics.units.NSKBatchReactor` call this and read the
        results back off the kinetic model.

        Parameters
        ----------
        t0 : float, optional
            Start time, in the model's own time units. ``None`` (default)
            means ``0.0`` with ``t_end`` and ``times[0]`` with ``times``.
            Integration with ``times`` starts at ``times[0]``, so an explicit
            ``t0`` must then equal ``times[0]`` (include ``t0`` in ``times``
            to start earlier).
        t_end : float
            End time, in the model's own time units. Required unless ``times``
            is given; must not be combined with ``times``.
        n_steps : int
            Number of output time points. Defaults to ``1000``. Ignored with
            ``times``.
        selections : sequence of str, optional
            RoadRunner selections to record (e.g. ``'time'`` and concentration
            selections like ``'[S]'``). If ``None``, defaults to ``'time'``
            followed by the concentration selection of every floating species.
        reset : bool
            If ``True``, call :meth:`reset` (default mode) before simulating.
            Defaults to ``False`` so callers that set initial conditions or
            parameter overrides just before simulating are not clobbered.
        times : sequence of float, optional
            Explicit output grid in the model's own time units; strictly
            increasing, at least two points. Rows are reported exactly at
            these times (events still fire in between).

        Returns
        -------
        numpy.ndarray
            The 2-D results array (also stored as :attr:`results_array`).

        Raises
        ------
        TypeError
            If ``times`` is not a sequence of numbers.
        ValueError
            If ``times`` is invalid (not 1-D, fewer than two points, not
            finite or not strictly increasing) or combined with ``t_end`` or
            a ``t0`` other than ``times[0]``; raised before any model state is
            touched.
        KineticSimulationError
            If neither ``t_end`` nor ``times`` is provided, or the underlying
            integration fails.
        """
        grid = None
        if times is not None:
            if t_end is not None:
                raise ValueError('Pass either t_end or times to simulate(), '
                                 'not both.')
            grid = self._check_times(times)
            if t0 is not None and t0 != grid[0]:
                raise ValueError(f't0={t0!r} conflicts with times[0]='
                                 f'{grid[0]!r}; integration starts at '
                                 'times[0] (include t0 in times to start '
                                 'there).')
        elif t_end is None:
            raise KineticSimulationError(
                'simulate() requires t_end (in the model time units).')
        elif t0 is None:
            t0 = 0.0
        if reset:
            self.reset()
        if selections is None:
            selections = ['time'] + [f'[{s}]'
                                     for s in self._te.getFloatingSpeciesIds()]
        cols = list(selections)
        try:
            if grid is not None:
                raw = np.array(self._te.simulate(times=grid, selections=cols))
            else:
                raw = np.array(self._te.simulate(t0, t_end, n_steps, cols))
        except Exception as e:
            raise KineticSimulationError(f'Kinetic simulation failed: {e}') from e
        self.results_col_names = cols
        self.results_array = raw
        self.results_dict = {cols[i]: raw[:, i] for i in range(len(cols))}
        self.results_df = pd.DataFrame(raw, columns=cols)
        # New results, run with the model's current values (run() records its
        # temporary ones after this returns).
        self._results_overrides = {}
        return raw

    def _default_selections(self):
        """``'time'`` plus the concentration of every floating species (the
        columns :meth:`simulate` records by default)."""
        return ['time'] + [f'[{s}]' for s in self._te.getFloatingSpeciesIds()]

    def _state_ids(self):
        """Ids whose values a baseline reset restores: :meth:`state_selections`
        (minus ``'time'``) plus the plain (amount) ids of floating species.

        Anything else a run sets (global parameters, boundary species,
        ``init(x)``...) would otherwise outlive the run.
        """
        ids = set(self.state_selections())
        ids.discard('time')
        ids.update(self._te.getFloatingSpeciesIds())
        return ids

    def _run_selections(self, condition, extra_selections=()):
        """Columns to record for ``condition``: its selections (or the
        defaults) with ``'time'`` first, then any ``extra_selections`` not
        already present."""
        sels = (list(condition.selections) if condition.selections is not None
                else self._default_selections())
        if 'time' not in sels:
            sels.insert(0, 'time')
        for s in extra_selections:
            if s not in sels:
                sels.append(s)
        return sels

    def _check_condition(self, condition, extra_selections=()):
        """Validate ``condition`` against this model without touching state.

        Returns
        -------
        list of str
            The selections the run will record.

        Raises
        ------
        TypeError
            If ``condition`` is not a :class:`Condition`.
        ValueError
            Listing unknown, rule-defined or read-only ids among the
            condition's ``parameters`` / ``initial`` and unknown ids among the
            selections.
        """
        if not isinstance(condition, Condition):
            raise TypeError(f'Expected a Condition; got '
                            f'{type(condition).__name__}.')
        msgs = []
        try:
            self._validate_ids([*condition.parameters, *condition.initial])
        except ValueError as e:
            msgs.append(str(e))
        sels = self._run_selections(condition, extra_selections)
        try:
            self._validate_ids(sels, settable=False)
        except ValueError as e:
            msgs.append(f'selections: {e}')
        try:
            self._check_reset_kwargs(condition.reset_kwargs)
        except ValueError as e:
            msgs.append(f'reset_kwargs: {e}')
        if msgs:
            label = f' {condition.label!r}' if condition.label else ''
            raise ValueError(f'Condition{label} does not fit this model: '
                             + ' '.join(msgs))
        return sels

    def _initial_order(self, initial):
        """``initial`` (or ``parameters``) as ``[(selection, value), ...]`` in
        the order :meth:`run` applies it: ``init(...)`` selections (compartments
        first; :meth:`run` writes these before the reset), then compartment
        volumes, then everything else, each group in the given order.

        Setting a compartment keeps species amounts (diluting their
        concentrations) and setting ``init(x)`` also sets ``x``, so applying
        the entries in key order would make equal conditions run differently.
        This order lets every concentration and current value given hold at
        ``t0``.
        """
        try:
            compartments = set(self._te.getCompartmentIds())
        except AttributeError:          # stub objects
            compartments = set()

        def phase(sel):
            if sel.startswith('init(') and sel.endswith(')'):
                return 0 if sel[5:-1] in compartments else 1
            return 2 if sel in compartments else 3

        return sorted(initial.items(), key=lambda item: phase(item[0]))

    def run(self, condition, *, extra_selections=()):
        """Simulate one :class:`Condition` and return its trajectory.

        Steps, in order: validate every id of the condition (nothing is
        changed if one is invalid); set the ``init(x)`` selections of either
        mapping, if any (see Notes); reset the model (if ``condition.reset``,
        via :meth:`reset` with ``condition.reset_kwargs``); set the other
        ``condition.parameters``; set the other ``condition.initial`` entries
        (see Notes for the order); integrate over ``[t0, t_end]`` with
        ``n_steps`` points, or over ``condition.times``.

        Parameters
        ----------
        condition : Condition
            What to run. Times are in model time units, values in the model's
            own units.
        extra_selections : sequence of str, optional
            Additional columns to record after the condition's own (e.g. the
            selections a metric needs); ids already recorded are skipped.

        Returns
        -------
        SimulationResult
            The trajectory (a copy of :attr:`results_df`), this model's units,
            the condition and its label.

        Raises
        ------
        TypeError
            If ``condition`` is not a :class:`Condition`.
        ValueError
            If an id in ``parameters`` / ``initial`` is unknown, defined by an
            assignment rule or read-only, or an ``init()`` form that would
            rewrite the model (see Notes), a selection is unknown, or
            ``reset_kwargs`` are rejected by :meth:`reset`; raised before any
            model state changes.
        KineticSimulationError
            If the integration fails.

        Notes
        -----
        ``parameters`` and then ``initial`` are each applied in a fixed order
        that does not depend on the order of their keys, so equal conditions
        always run the same: first ``init(x)`` selections (compartments
        first), then compartment volumes, then everything else (``'[S]'``,
        amounts, counters...); each group keeps the given order. Setting a
        volume keeps species *amounts*, so a volume alone scales the
        starting concentrations, while every concentration given alongside
        it holds in that volume at ``t0``.

        ``init(x)`` selections, whether in ``initial`` or ``parameters``, are
        written first, before the reset. RoadRunner re-initializes the whole
        model whenever an ``init(...)`` value is written: time, state and
        every global parameter return to their initial values, and initial
        assignments are re-evaluated. An ``init(x)`` entry therefore changes
        the model's *definition* for the run: ``initial={'init(k)': v}``
        starts it as a model defined with ``k = v`` would start, including
        every initial assignment that reads ``k``. The model's other
        settings (global parameters and boundary species not defined by an
        initial assignment, e.g. an edited ``km._te['k3']``) are written back
        right after, and the reset, ``parameters`` and the rest of
        ``initial`` follow as usual. Two ``init()`` forms would rewrite the
        model's SBML for good and are rejected with ``ValueError``:
        ``init()`` of an id defined by an initial assignment (set ``init()``
        of its inputs instead, e.g. ``init(k)`` for ``S = 20*k``), and
        ``init()`` of a species in the other form than its SBML declares
        (``init([S])`` for a species declared with an initial concentration,
        the usual Antimony ``S = 50``; ``init(S)`` for one declared with an
        initial amount).

        ``parameters`` -- and any ``initial`` entry that a baseline reset
        would not restore by itself, such as a global parameter, a boundary
        species or an ``init(x)`` selection -- are restored to their pre-run
        values in a ``finally`` block, so a condition never leaks into later
        runs. The underlying quantity is restored: a boundary species given
        as a concentration (``'[B]'``) gets its pre-run *amount* back, and
        the model's settings -- and the integrator's, which re-initializing
        the model resets for some integrators (e.g. rk45's ``epsilon``) --
        are restored along with ``init(x)`` entries.
        Floating species, compartments and other state are left at their
        end-of-run values, as after :meth:`simulate` -- except after a run
        with ``init(x)`` entries: restoring them re-initializes the model
        (for a compartment, RoadRunner also re-derives the species values in
        it), so the model is then at its initial state, not its end-of-run
        state. Read the final state from the returned result (e.g.
        ``result.final``).

        Like :meth:`simulate`, the run also fills the legacy
        :attr:`results_df` / :attr:`results_array` / :attr:`results_dict` /
        :attr:`results_col_names` slots, so :meth:`plot_simulation_results`
        and :func:`nskinetics.compute_flux_summary` keep working. It also
        records the values it set and then restored (``parameters`` and
        non-state ``initial`` entries), and
        :func:`~nskinetics.compute_flux_summary` re-evaluates the rate laws
        with them; the next :meth:`simulate` clears them.

        Examples
        --------
        >>> cond = nsk.Condition(initial={'[S]': 50}, t_end=48, n_steps=481)  # doctest: +SKIP
        >>> res = km.run(cond)                         # doctest: +SKIP
        >>> res.final['[P]'], res.at(12.0)['[S]']      # doctest: +SKIP
        """
        sels = self._check_condition(condition, extra_selections)
        r = self._te
        # Both mappings in the fixed phase order (all parameters, then all
        # initial), so equal conditions run the same whatever their key order.
        initial = self._initial_order(condition.initial)
        entries = [*self._initial_order(condition.parameters), *initial]
        # Writing an init(...) value re-initializes the model, so those go
        # first (compartments first), before the reset; see Notes.
        inits = self._initial_order({s: v for s, v in entries
                                     if _init_base(s) is not None})
        state = self._state_ids() if initial else ()
        targets, settings = self._restore_targets(
            [*condition.parameters,
             *(sel for sel, _ in initial if sel not in state)])
        # Read everything before the first write.
        previous = {t: r[t] for t in targets}
        integrator = self._integrator_snapshot() if inits else None
        try:
            self._write_values(inits, {p: previous[p] for p in settings})
            if condition.reset:
                # Thawed: a custom reset_func gets plain dicts/lists, not the
                # condition's read-only mappings/tuples.
                self.reset(**_thaw(condition.reset_kwargs))
            for sel, value in entries:
                if _init_base(sel) is None:
                    r[sel] = value
            # The values the rate laws see, for compute_flux_summary's replay
            # after the restore (recorded below only if the run succeeds).
            applied = {t: float(r[t]) for t in targets
                       if _init_base(t) is None and t not in state}
            if condition.times is not None:
                self.simulate(selections=sels, times=condition.times)
            else:
                self.simulate(condition.t0, condition.t_end,
                              condition.n_steps, sels)
        finally:
            # init(...) targets come first: writing one re-initializes the
            # model, so everything else is written back after it (and the
            # integrator's settings, which that can reset, after that).
            for t in targets:
                r[t] = previous[t]
            self._reapply_integrator(integrator)
        self._results_overrides = applied
        return SimulationResult(df=self.results_df.copy(),
                                units=dict(self._units), condition=condition,
                                label=condition.label)

    # --- copies -------------------------------------------------------------
    def _user_parameter_values(self):
        """Current values of the settings a user may have edited.

        Every global parameter and boundary species that is neither an
        assignment-rule, rate-rule nor event target; those are state
        (restored by a baseline reset) or derived, not settings. Neither kind
        is restored by ``RoadRunner.reset()``, so an edit persists across
        resets and a copy must carry it to run the same.

        Only ids the SBML document defines are listed: RoadRunner adds
        globals of its own (e.g. the conserved-moiety totals ``_CSUM0``...
        under ``conservedMoietyAnalysis``) that a model rebuilt from the SBML
        could not be given.

        Returns
        -------
        dict
            ``{id: float}``: parameter values, and boundary species as
            amounts (plain ids), in the model's own units.
        """
        r = self._te
        excluded = set(r.getAssignmentRuleIds()) | set(r.getRateRuleIds())
        excluded.update(var for event in self.events for var in event.do)
        excluded.update(self._sbml_event_targets())
        ids = [*r.getGlobalParameterIds(), *r.getBoundarySpeciesIds()]
        defined = _sbml_parameter_and_species_ids(self._get_origin_sbml())
        return {i: float(r[i]) for i in ids
                if i not in excluded and i in defined}

    def _extra_attributes(self):
        """Public instance attributes beyond the ones every model owns
        (e.g. ``default_max_n_glu_spikes`` on the shipped model)."""
        return {k: v for k, v in vars(self).items()
                if not k.startswith('_') and k not in _OWN_PUBLIC_ATTRS}

    def _blueprint(self):
        """Return a picklable :class:`_ModelBlueprint` of this model."""
        return _ModelBlueprint.from_model(self)

    def copy(self):
        """Return an independent model rebuilt from this model's origin SBML.

        Returns
        -------
        KineticModel
            A new model on its own RoadRunner, at its baseline state (as after
            ``reset(mode='baseline')``; the custom :attr:`reset_func` is never
            called), carrying: the user-edited values of global parameters
            and boundary species (every one that is not an assignment-rule,
            rate-rule or event target; reset never restores these, so edits
            persist), the ``events`` list (shallow; already compiled into
            the origin SBML when compiled here), units (copied),
            :attr:`reset_func`, the integrator name and the settings that
            differ from its defaults (tolerances, step limits, seed...),
            RoadRunner's ``conservedMoietyAnalysis`` flag, and extra public
            instance attributes (e.g. ``default_max_n_glu_spikes``;
            ``list``/``dict``/``set`` values are shallow-copied, others
            shared). The legacy ``results_*`` slots are not copied.

        Raises
        ------
        Exception
            Raised by Tellurium/RoadRunner if the origin SBML cannot be loaded
            (e.g. a wrapped object that exposes no SBML).

        Notes
        -----
        The copy is rebuilt from the *origin* SBML, ``RoadRunner.getSBML()``
        read at the time of the call: the loaded document plus every
        structural edit made on the RoadRunner (events compiled by
        :meth:`compile_events` or directly by ``Event.compile`` /
        ``FeedSpike.compile``, ``r.addEvent(...)``...) and values written to
        ``init(x)`` selections, with origin values and initial assignments
        intact. It is not exported with ``getCurrentSBML()``, which bakes the
        post-run state into the initial values and drops initial
        assignments. The run state and the legacy ``results_*`` slots are not
        carried. Mutating the copy never affects this model, and vice versa.

        Loading the copy keeps ``r.<id>`` attribute access working on other
        live models (libroadrunner would otherwise remove their class-level
        properties).
        """
        return self._blueprint().build()

    # --- introspection ------------------------------------------------------
    def _value_or_nan(self, selection):
        try:
            return float(self._te[selection])
        except RuntimeError:
            return np.nan

    def _introspect(self):
        """Build the (parameters, species) tables from the origin SBML."""
        model = libsbml.readSBMLFromString(self._get_origin_sbml()).getModel()
        if model is None:
            raise KineticSimulationError('The model has no readable SBML.')
        kinds = {}
        for var in _column(_sbt.extract_initial_assignments(model), 'symbol'):
            kinds[var] = 'initial_assignment'
        rules = _sbt.extract_rules(model)
        for rule_type, var in zip(_column(rules, 'rule_type'),
                                  _column(rules, 'variable')):
            # Rules take precedence over an initial assignment of the same id.
            if rule_type == 'assignment':
                kinds[var] = 'assignment_rule'
            elif rule_type == 'rate':
                kinds[var] = 'rate_rule'
        event_targets = set(_column(
            _sbt.extract_events(model)['event_assignments'], 'variable'))
        event_targets.update(var for event in self.events for var in event.do)
        used = _reaction_usage(model)

        param_rows = []
        params = _sbt.extract_parameters(model)
        for pid, sbml_value, units, constant in zip(
                _column(params, 'id'), _column(params, 'value'),
                _column(params, 'units'), _column(params, 'constant')):
            kind = kinds.get(pid, 'parameter')
            if kind == 'assignment_rule':
                # init(x) of a rule-defined id evaluates the rule on the
                # current state; it is not an origin value.
                initial = np.nan
            else:
                initial = self._value_or_nan(f'init({pid})')
                if np.isnan(initial) and sbml_value != '':
                    initial = float(sbml_value)
            param_rows.append({
                'id': pid,
                'value': self._value_or_nan(pid),
                'initial_value': initial,
                'units': units or '',
                # roadrunner's addEvent leaves an event target's SBML
                # constant='true' in place, so event targets are reported as
                # non-constant explicitly.
                'constant': bool(constant) and pid not in event_targets,
                'kind': kind,
                'assigned_by_events': pid in event_targets,
                'used_in_reactions': tuple(used.get(pid, ())),
            })
        species_rows = []
        species = _sbt.extract_species(model)
        for sid, units, constant, comp, boundary in zip(
                _column(species, 'id'), _column(species, 'substance_units'),
                _column(species, 'constant'), _column(species, 'compartment'),
                _column(species, 'boundary_condition')):
            kind = kinds.get(sid, 'species')
            species_rows.append({
                'id': sid,
                'value': self._value_or_nan(f'[{sid}]'),
                'initial_value': (np.nan if kind == 'assignment_rule' else
                                  self._value_or_nan(f'init([{sid}])')),
                'units': units or '',
                'constant': bool(constant) and sid not in event_targets,
                'kind': kind,
                'assigned_by_events': sid in event_targets,
                'used_in_reactions': tuple(used.get(sid, ())),
                'compartment': comp or '',
                'boundary': bool(boundary),
            })
        return (_table(param_rows, _PARAMETER_COLUMNS),
                _table(species_rows, _SPECIES_COLUMNS))

    @property
    def parameters(self):
        """Global parameters of the model, as a table.

        Returns
        -------
        pandas.DataFrame
            Indexed by parameter id, with columns ``value`` (current value),
            ``initial_value`` (``init(id)``: for a plain parameter, a
            rate-rule or an event target, the origin value, which edits of
            the id itself leave unchanged; for kind ``'initial_assignment'``,
            the assignment evaluated with the *current* values of its inputs,
            so it follows edits of those; NaN for assignment-rule ids),
            ``units`` (SBML unit id or ``''``), ``constant`` (SBML attribute,
            and ``False`` for event targets), ``kind`` (``'parameter'``,
            ``'assignment_rule'``, ``'rate_rule'`` or
            ``'initial_assignment'``), ``assigned_by_events`` (bool) and
            ``used_in_reactions`` (tuple of the reaction ids whose kinetic law
            names the parameter directly). Values are in the model's own
            units.

        Raises
        ------
        KineticSimulationError
            If the origin SBML holds no readable model.

        Notes
        -----
        Structure is read from the origin SBML with libSBML (rules, initial
        assignments, events and kinetic-law ASTs), not inferred from id
        prefixes; compiled :class:`Event` objects count as events. A kinetic
        law's local parameters shadow globals of the same id.
        """
        return self._introspect()[0]

    @property
    def species(self):
        """Species of the model, as a table.

        Returns
        -------
        pandas.DataFrame
            Indexed by species id, with the same columns as
            :attr:`parameters` (``value`` and ``initial_value`` are
            concentrations, i.e. the ``[id]`` and ``init([id])`` selections,
            in the model's concentration units; as for parameters, the
            ``initial_value`` of a species defined by an initial assignment
            is that assignment evaluated with the current values of its
            inputs; ``units`` is the SBML substance-units id;
            ``kind`` is ``'species'`` unless a rule or initial assignment
            defines the species; ``used_in_reactions`` lists reactions in
            which it participates or which name it in their kinetic law), plus
            ``compartment`` and ``boundary`` (SBML boundary condition).

        Raises
        ------
        KineticSimulationError
            If the origin SBML holds no readable model.
        """
        return self._introspect()[1]

    # --- plotting -----------------------------------------------------------
    def plot_simulation_results(self, variables=None, labels=None,
                                xlim=None, ylim=None, markers=None,
                                save_fig=False, filename=None, figwidth=3.9,
                                ylabel=None, ylabel_unit=None, **kwargs):
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
        ylabel : str, optional
            Quantity name for the y-axis, rendered in bold. Defaults to
            ``'Concentration'``. Pass e.g. ``'Volume'`` when plotting a
            non-concentration selection such as the working-volume compartment.
        ylabel_unit : str, optional
            Unit shown in brackets after ``ylabel``. Defaults to this system's
            ``'conc'`` unit. Known units (``'g/L'``, ``'M'``, ``'L'``, ``'m3'``,
            …) are rendered as LaTeX; anything else is shown verbatim. Pass e.g.
            ``'L'`` alongside ``ylabel='Volume'``.
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
        return _plot_time_course(
            d, self.results_col_names, self._units, variables=variables,
            labels=labels, xlim=xlim, ylim=ylim, markers=markers,
            save_fig=save_fig, filename=filename, figwidth=figwidth,
            ylabel=ylabel, ylabel_unit=ylabel_unit, **kwargs)

    # Aliases so both kinetic modelers and TEA practitioners find the plot
    # under a familiar name.
    plot_time_course = plot_simulation_results
    plot_trajectory = plot_simulation_results


# --- introspection helpers ----------------------------------------------------
_PARAMETER_COLUMNS = ('value', 'initial_value', 'units', 'constant', 'kind',
                      'assigned_by_events', 'used_in_reactions')
_SPECIES_COLUMNS = _PARAMETER_COLUMNS + ('compartment', 'boundary')


def _column(df, name):
    """Return column ``name`` of an ``_sbml_to_table`` frame as a list
    (empty when the frame has no rows, and hence no columns)."""
    return list(df[name]) if name in df.columns else []


def _table(rows, columns):
    """DataFrame indexed by ``'id'`` with ``columns``, also when empty."""
    df = pd.DataFrame(rows, columns=('id',) + tuple(columns))
    return df.set_index('id')


def _sbml_parameter_and_species_ids(sbml):
    """Ids of the global parameters and species an SBML document defines
    (an empty set for ``None`` or unreadable text)."""
    model = (libsbml.readSBMLFromString(sbml).getModel()
             if sbml is not None else None)
    if model is None:
        return set()
    return ({p.getId() for p in model.getListOfParameters()}
            | {s.getId() for s in model.getListOfSpecies()})


def _ast_names(ast):
    """All identifier names referenced in a libSBML math AST."""
    names = set()
    stack = [ast] if ast is not None else []
    while stack:
        node = stack.pop()
        if node.isName():
            names.add(node.getName())
        for i in range(node.getNumChildren()):
            stack.append(node.getChild(i))
    return names


def _reaction_usage(model):
    """Map each id to the reactions that use it, in document order.

    A reaction uses an id if its kinetic law names it (local parameters
    shadow globals) or, for species, if the species is a reactant, product or
    modifier.
    """
    used = {}
    for reaction in model.getListOfReactions():
        names = set()
        if reaction.isSetKineticLaw():
            law = reaction.getKineticLaw()
            names = _ast_names(law.getMath())
            names -= {p.getId() for p in law.getListOfLocalParameters()}
            names -= {p.getId() for p in law.getListOfParameters()}
        for refs in (reaction.getListOfReactants(),
                     reaction.getListOfProducts(),
                     reaction.getListOfModifiers()):
            names.update(ref.getSpecies() for ref in refs)
        rid = reaction.getId()
        for name in names:
            used.setdefault(name, []).append(rid)
    return used


# --- blueprint ----------------------------------------------------------------
def _detach(value):
    """Shallow-copy builtin containers so a copy's attribute is not aliased
    to the original's; other objects are shared."""
    if isinstance(value, (list, dict, set, bytearray)):
        return _copy.copy(value)
    return value


def _configure_integrator(r, name, settings):
    """Select RoadRunner ``r``'s integrator ``name`` and apply ``settings``.

    Only settings that differ from the integrator's current values are
    written: some defaults do not survive a ``setValue`` round trip
    (gillespie's ``maximum_num_steps=0`` is stored as the wrong type, and
    every later simulate fails with a variant error).
    """
    if r.getIntegrator().getName() != name:
        r.setIntegrator(name)
    integrator = r.getIntegrator()
    for key, value in settings.items():
        if integrator.getValue(key) != value:
            integrator.setValue(key, value)


@dataclass
class _ModelBlueprint:
    """Picklable recipe that rebuilds a :class:`KineticModel` (private).

    :meth:`KineticModel.copy` is ``_ModelBlueprint.from_model(km).build()``.
    Process-pool workers receive a blueprint (it pickles as long as its
    ``reset_func`` and extra attributes do; see :meth:`check_picklable`) and
    build one model each.

    Attributes
    ----------
    sbml : str
        The model's origin SBML, ``RoadRunner.getSBML()`` as read when the
        blueprint was made: the loaded document plus structural edits
        (compiled events...) and ``init(x)`` values, with origin values.
    units : dict
        ``{'time': ..., 'conc': ...}``.
    events : list of Event
        The model's :attr:`KineticModel.events` (shallow).
    events_compiled : bool
        Whether ``events`` are already part of ``sbml``.
    reset_func : callable or None
        Custom reset function (must be importable, e.g. module level, to
        pickle).
    integrator_name : str
        RoadRunner integrator name (e.g. ``'cvode'``).
    integrator_settings : dict
        ``{setting: value}`` from the integrator (tolerances, step limits,
        seed...); :meth:`build` writes only those that differ from the new
        integrator's defaults.
    parameter_values : dict
        Current values of user-settable global parameters and boundary
        species (amounts; see :meth:`KineticModel._user_parameter_values`).
    attributes : dict
        Extra public instance attributes.
    model_class : type
        Class to instantiate (``KineticModel`` or a subclass; built through
        ``KineticModel.__init__``).
    conserved_moiety_analysis : bool
        RoadRunner's ``conservedMoietyAnalysis`` flag.
    """
    sbml: str
    units: dict
    events: list = field(default_factory=list)
    events_compiled: bool = False
    reset_func: object = None
    integrator_name: str = 'cvode'
    integrator_settings: dict = field(default_factory=dict)
    parameter_values: dict = field(default_factory=dict)
    attributes: dict = field(default_factory=dict)
    model_class: type = KineticModel
    conserved_moiety_analysis: bool = False

    @classmethod
    def from_model(cls, model):
        """Capture ``model`` (a :class:`KineticModel`) without mutating it."""
        r = model._te
        integrator = r.getIntegrator()
        return cls(
            sbml=model._get_origin_sbml(),
            units=dict(model._units),
            events=list(model.events),
            events_compiled=bool(model._events_compiled),
            reset_func=model._reset_func,
            integrator_name=integrator.getName(),
            integrator_settings={k: integrator.getValue(k)
                                 for k in integrator.getSettings()},
            parameter_values=model._user_parameter_values(),
            attributes={k: _detach(v)
                        for k, v in model._extra_attributes().items()},
            model_class=type(model),
            conserved_moiety_analysis=bool(
                getattr(r, 'conservedMoietyAnalysis', False)),
        )

    def build(self):
        """Build a fresh :class:`KineticModel` on its own RoadRunner.

        Returns
        -------
        KineticModel
            At its baseline state: the origin SBML is loaded (keeping other
            live models' ``r.<id>`` properties; see
            :func:`_keep_roadrunner_properties`), conserved-moiety analysis
            switched on if it was, the integrator configured, parameter
            values applied, then ``reset(mode='baseline')`` (never the custom
            ``reset_func``) re-evaluates initial assignments against those
            values and restores compartments and event-assigned variables.
        """
        with _keep_roadrunner_properties():
            r = te.loadSBMLModel(self.sbml)
            if self.conserved_moiety_analysis:
                r.conservedMoietyAnalysis = True
        _configure_integrator(r, self.integrator_name,
                              self.integrator_settings)
        for pid, value in self.parameter_values.items():
            r[pid] = value
        km = self.model_class.__new__(self.model_class)
        KineticModel.__init__(km, r, units=dict(self.units),
                              reset=self.reset_func)
        km.events = list(self.events)
        km._events_compiled = self.events_compiled
        for name, value in self.attributes.items():
            setattr(km, name, _detach(value))
        # Start where a baseline reset of the original starts (also with a
        # compartment initial assignment on a carried parameter).
        km.reset(mode='baseline')
        return km

    def check_picklable(self):
        """Raise ``ValueError`` naming the first field that cannot be pickled.

        Raises
        ------
        ValueError
            E.g. for a ``reset_func`` that is a lambda or closure, or an extra
            attribute holding an unpicklable object; such a model cannot be
            sent to process workers.
        """
        candidates = [(f.name, getattr(self, f.name)) for f in fields(self)
                      if f.name != 'attributes']
        candidates += [(f'attributes[{k!r}]', v)
                       for k, v in self.attributes.items()]
        for name, value in candidates:
            try:
                pickle.dumps(value)
            except (pickle.PicklingError, TypeError, AttributeError) as e:
                hint = (' Define the reset function at module level (not as a '
                        'lambda or closure) to use process workers.'
                        if name == 'reset_func' else '')
                raise ValueError(f'The model cannot be sent to worker '
                                 f'processes: {name} is not picklable '
                                 f'({e}).{hint}') from e
