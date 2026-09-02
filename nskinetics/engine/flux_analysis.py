# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Flux and inhibition analysis of a simulated kinetic model.

Re-evaluates a model's own rate laws along its stored trajectory to obtain
per-reaction cumulative flux, and per-(reaction, inhibitor) flux-weighted
fraction of potential flux lost, by zeroing inhibition coefficients one group
at a time. Everything is derived from the model itself; no rate law is
re-implemented here.
"""

import re
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..exceptions import KineticSimulationError

__all__ = ('FluxSummary', 'compute_flux_summary')

#: Columns of the long-format CSV written by :meth:`FluxSummary.to_csv`.
_CSV_COLUMNS = ('reaction', 'quantity', 'key', 'value')

#: The SBML csymbol ``definitionURL`` fragment marking a reference to
#: simulation time inside a MathML expression.
_TIME_CSYMBOL = 'symbols/time'

#: Matches one ``<kineticLaw>`` element (and its contents) in an SBML document.
_KINETIC_LAW_RE = re.compile(r'<kineticLaw\b.*?</kineticLaw>', re.S)

#: Matches one ``<assignmentRule>``: its target variable and its math body.
_ASSIGNMENT_RULE_RE = re.compile(
    r'<assignmentRule\b[^>]*\bvariable="([^"]+)"(.*?)</assignmentRule>', re.S)

#: Matches a MathML identifier reference; the capture is the whole name, so
#: membership tests compare complete tokens rather than substrings.
_CI_RE = re.compile(r'<ci\b[^>]*>\s*([^<]*?)\s*</ci>')


def _read_float(val):
    """Parse one CSV ``value`` cell, or ``None`` if it is blank/unparseable.

    Parameters
    ----------
    val : str or float
        Raw cell contents as read back from the CSV.

    Returns
    -------
    float or None
    """
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


@dataclass
class FluxSummary:
    """Cumulative flux and inhibition strengths for one simulated run.

    Parameters
    ----------
    label : str or None
        Human-readable scenario label (used as a panel header).
    reaction_ids : list of str
        Reactions summarized, in the order requested.
    cumulative_mass : dict
        ``{reaction_id: grams processed}`` over the run (extensive integral).
    cumulative_flux : dict
        ``{reaction_id: g/L of final broth}`` (``cumulative_mass`` / final
        volume). Each value is grams of *that step's own substrate* consumed,
        so two reactions are directly comparable where they compete for the
        same substrate (a branch point), but not along a chain: the
        stoichiometry, and hence the mass carried per unit of flux, differs
        from step to step.
    fraction_lost : dict
        ``{reaction_id: {inhibitor: fraction}}``; ``1 - actual/uninhibited``.
        Negative for enhancement (e.g. product-accelerated decay). The
        "uninhibited" counterfactual is *not* a re-simulation: it is this
        model's own rate law re-evaluated along the actual (inhibited)
        trajectory with that inhibitor's coefficients set to zero, so it
        measures the flux lost at the states the system really visited.
        Reactions with no entry in ``inhibition_map`` are absent from this
        dict (use ``.get()``), and a fraction is ``0.0`` when the uninhibited
        integral is exactly zero (a dormant reaction).
    fraction_lost_all : dict
        ``{reaction_id: fraction}`` with every mapped coefficient zeroed at
        once — same counterfactual convention and same omissions as
        ``fraction_lost``.
    final_volume : float
        Broth volume at the last trajectory row (model volume units).
    final_concentrations : dict
        ``{species_id: concentration}`` at the last row (bare ids, no brackets).
    t_end : float
        Final time (model time units).
    inhibitors : list of str
        Distinct inhibitor names present in the mapping, in first-seen order.
    """
    label: str | None
    reaction_ids: list
    cumulative_mass: dict
    cumulative_flux: dict
    fraction_lost: dict
    fraction_lost_all: dict
    final_volume: float
    final_concentrations: dict
    t_end: float
    inhibitors: list = field(default_factory=list)

    def to_csv(self, path):
        """Write this summary to a long-format CSV.

        The file has the columns ``reaction, quantity, key, value``, where
        ``quantity`` is one of ``cumulative_mass``, ``cumulative_flux``,
        ``fraction_lost``, ``fraction_lost_all``, ``final_concentration`` or
        ``meta``, and ``key`` is the inhibitor/species name (or the metadata
        name for ``meta`` rows) and is empty where it does not apply. Rows
        that are not per-reaction (``final_concentration``, ``meta``) carry an
        empty ``reaction``. Reactions absent from ``fraction_lost`` /
        ``fraction_lost_all`` (i.e. with no entry in the inhibition map) get
        no rows for those quantities, so :meth:`from_csv` reproduces the same
        omissions rather than inventing entries. ``inhibitors`` is written as a
        single ``|``-joined field, so an inhibitor name must not itself contain
        ``|``; an empty ``label`` is written as an empty field and reads back
        as ``None``.

        Parameters
        ----------
        path : str or path-like
            Destination CSV path; an existing file is overwritten.

        Returns
        -------
        None
        """
        rows = []
        for rid in self.reaction_ids:
            rows.append((rid, 'cumulative_mass', '', self.cumulative_mass[rid]))
            rows.append((rid, 'cumulative_flux', '', self.cumulative_flux[rid]))
            if rid in self.fraction_lost_all:
                rows.append((rid, 'fraction_lost_all', '',
                             self.fraction_lost_all[rid]))
            for inh, v in self.fraction_lost.get(rid, {}).items():
                rows.append((rid, 'fraction_lost', inh, v))
        for sp, v in self.final_concentrations.items():
            rows.append(('', 'final_concentration', sp, v))
        rows.append(('', 'meta', 'final_volume', self.final_volume))
        rows.append(('', 'meta', 't_end', self.t_end))
        rows.append(('', 'meta', 'label', self.label if self.label else ''))
        rows.append(('', 'meta', 'inhibitors', '|'.join(self.inhibitors)))
        pd.DataFrame(rows, columns=list(_CSV_COLUMNS)).to_csv(path, index=False)

    @classmethod
    def from_csv(cls, path):
        """Reconstruct a :class:`FluxSummary` written by :meth:`to_csv`.

        Every column is read as text (``dtype=str``, ``keep_default_na=False``)
        so blanks stay blank and floats are converted from their full
        ``repr``, i.e. they round-trip exactly. A ``fraction_lost_all`` row
        whose value is blank or ``nan`` is skipped, so a summary written by an
        older version (which emitted such a row for unmapped reactions) still
        reloads with that reaction absent from the dict. Every field is
        recovered exactly except an empty-string ``label``, which reads back as
        ``None`` (the two are not distinguished on disk); ``inhibitors`` is
        split on ``|``, so a name containing that character does not survive.

        Parameters
        ----------
        path : str or path-like
            CSV written by :meth:`to_csv`.

        Returns
        -------
        FluxSummary
        """
        df = pd.read_csv(path, dtype={c: str for c in _CSV_COLUMNS},
                         keep_default_na=False)
        reaction_ids, cmass, cflux, floss, fall, fconc = [], {}, {}, {}, {}, {}
        label, final_volume, t_end, inhibitors = None, 1.0, 0.0, []
        for _, row in df.iterrows():
            q, rid, key, val = row['quantity'], row['reaction'], row['key'], \
                row['value']
            if q == 'cumulative_mass':
                if rid not in reaction_ids:
                    reaction_ids.append(rid)
                cmass[rid] = float(val)
            elif q == 'cumulative_flux':
                cflux[rid] = float(val)
            elif q == 'fraction_lost_all':
                v = _read_float(val)
                if v is not None and not np.isnan(v):
                    fall[rid] = v
            elif q == 'fraction_lost':
                floss.setdefault(rid, {})[key] = float(val)
            elif q == 'final_concentration':
                fconc[key] = float(val)
            elif q == 'meta':
                if key == 'final_volume':
                    final_volume = float(val)
                elif key == 't_end':
                    t_end = float(val)
                elif key == 'label':
                    label = val or None
                elif key == 'inhibitors':
                    inhibitors = [x for x in str(val).split('|') if x]
        return cls(label=label, reaction_ids=reaction_ids,
                   cumulative_mass=cmass, cumulative_flux=cflux,
                   fraction_lost=floss, fraction_lost_all=fall,
                   final_volume=final_volume, final_concentrations=fconc,
                   t_end=t_end, inhibitors=inhibitors)


def _resolve_model(model_or_unit):
    km = getattr(model_or_unit, 'nsk_kinetic_model', None)
    if km is not None:
        return km
    if getattr(model_or_unit, 'results_df', None) is not None or \
            hasattr(model_or_unit, 'state_selections'):
        return model_or_unit
    raise TypeError(
        'compute_flux_summary expects a KineticModel or an NSKBatchReactor '
        f'(got {type(model_or_unit).__name__}).')


def _end_index(df, t_end, unit):
    """Index of the last trajectory row inside the ``0..t_end`` window.

    Parameters
    ----------
    df : pandas.DataFrame
        Trajectory, with a ``time`` column.
    t_end : float
        Inclusive upper bound of the integration window.
    unit : NSKBatchReactor or None
        The reactor the trajectory came from, if any. When ``t_end`` is that
        reactor's own harvest time *and* the row it points at really is at
        that time, its stored ``tau_index`` is used verbatim, so the last
        included row is exactly the row the reactor reports as
        ``nsk_results_specific_tau``. A ``tau_index`` that is stale with
        respect to the stored trajectory falls through to the time scan below.

    Returns
    -------
    int

    Raises
    ------
    ValueError
        If no row satisfies ``time <= t_end``.
    """
    if unit is not None and getattr(unit, 'tau', None) == t_end:
        idx = getattr(unit, 'tau_index', None)
        if idx is not None and 0 <= idx < len(df) \
                and df['time'].iat[int(idx)] == t_end:
            return int(idx)
    t = df['time'].to_numpy()
    keep = np.flatnonzero(t <= t_end)
    if not len(keep):
        raise ValueError(
            f'No trajectory row at or before t_end={t_end} '
            f'(first recorded time is {float(t[0])}).')
    return int(keep[-1])


def _write_order(selections):
    """Order selections so a write-back lands at the intended values.

    Non-concentration state (``time``, compartments, parameters, rate-rule
    variables) is written first and concentration selections (``[species]``)
    last: writing a compartment volume rescales the concentrations derived
    from the held amounts, so concentrations must be set afterwards. Every
    write-back in this module (trajectory rows and the state restore) goes
    through this order, so neither depends on the order
    :meth:`~KineticModel.state_selections` happens to emit.

    Parameters
    ----------
    selections : sequence of str
        RoadRunner selections to reorder.

    Returns
    -------
    list of str
    """
    return ([c for c in selections if not c.startswith('[')]
            + [c for c in selections if c.startswith('[')])


def _time_dependent_variables(sbml):
    """Assignment-rule variables whose value depends on simulation time.

    Seeded with the rules whose math contains the SBML time csymbol, then
    closed transitively: a rule that references (as a whole ``<ci>`` token) a
    variable already known to be time-dependent is itself time-dependent.

    Parameters
    ----------
    sbml : str
        The model's current SBML document.

    Returns
    -------
    set of str
        Variable ids; empty when no assignment rule reads time.
    """
    rules = {var: body for var, body in _ASSIGNMENT_RULE_RE.findall(sbml)}
    refs = {var: set(_CI_RE.findall(body)) for var, body in rules.items()}
    time_vars = {var for var, body in rules.items() if _TIME_CSYMBOL in body}
    changed = bool(time_vars)
    while changed:
        changed = False
        for var, names in refs.items():
            if var not in time_vars and names & time_vars:
                time_vars.add(var)
                changed = True
    return time_vars


def _kinetic_laws_use_time(r):
    """Whether any of the model's rate laws reads simulation time.

    A rate law reads time either directly, through the SBML time csymbol, or
    indirectly, by referencing an assignment-rule variable that (transitively)
    depends on time -- e.g. ``temp := T0 + ramp*time`` used in an Arrhenius
    rate law. Both count: :func:`compute_flux_summary` re-evaluates rate laws
    at written-back states, so an undetected indirect dependence would
    evaluate every rate at one frozen clock value and silently corrupt every
    integral.

    A reference to time anywhere else does not count -- an assignment rule no
    rate law reads (the shipped model's ``prod_EtOH := s_EtOH/time``), or an
    event trigger -- because the rates are read directly rather than
    integrated. A model whose current SBML cannot be retrieved is treated as
    time-dependent, i.e. the conservative answer.

    Parameters
    ----------
    r : roadrunner.RoadRunner
        The loaded model.

    Returns
    -------
    bool
    """
    try:
        sbml = r.getCurrentSBML()
    except (AttributeError, RuntimeError):
        return True
    time_vars = _time_dependent_variables(sbml)
    for law in _KINETIC_LAW_RE.findall(sbml):
        if _TIME_CSYMBOL in law:
            return True
        if time_vars and time_vars & set(_CI_RE.findall(law)):
            return True
    return False


def _apply_row(r, arrs, cols, i):
    """Write trajectory row ``i`` back into the model, in ``cols`` order.

    Parameters
    ----------
    r : roadrunner.RoadRunner
        The loaded model.
    arrs : dict
        ``{selection: numpy array over the trajectory}`` (hoisted out of the
        caller's loop, so no per-row DataFrame lookup happens here).
    cols : sequence of str
        Selections to write, already in :func:`_write_order` order.
    i : int
        Trajectory row index.
    """
    for c in cols:
        r[c] = arrs[c][i]


def _restore(r, snap, ordered_cols):
    for c in ordered_cols:
        r[c] = snap[c]


def _rates_along(r, arrs, cols, idx_of, n):
    """Re-evaluate the requested reactions' rates along a stored trajectory.

    This is the module's single write-back loop: every base and counterfactual
    integral goes through it, so the state written per row is identical in
    both.

    Parameters
    ----------
    r : roadrunner.RoadRunner
        The loaded model.
    arrs : dict
        ``{selection: numpy array}`` as for :func:`_apply_row`.
    cols : sequence of str
        Selections written back on each row.
    idx_of : dict
        ``{reaction_id: index into getReactionRates()}``.
    n : int
        Number of trajectory rows.

    Returns
    -------
    dict
        ``{reaction_id: numpy array of rates}``.
    """
    items = tuple(idx_of.items())
    out = {rid: np.empty(n) for rid, _j in items}
    for i in range(n):
        _apply_row(r, arrs, cols, i)
        rates = r.getReactionRates()
        for rid, j in items:
            out[rid][i] = rates[j]
    return out


def _cum(r, arrs, cols, idx_of, n, t):
    """Trapezoidal integral of each requested reaction's re-evaluated rate.

    Parameters
    ----------
    r, arrs, cols, idx_of, n
        As for :func:`_rates_along`.
    t : numpy.ndarray
        Trajectory times.

    Returns
    -------
    dict
        ``{reaction_id: cumulative extensive flux}``.
    """
    return {rid: float(np.trapezoid(v, t))
            for rid, v in _rates_along(r, arrs, cols, idx_of, n).items()}


def _frac(base, counterfactual):
    if counterfactual == 0.0:
        return 0.0
    return 1.0 - base / counterfactual


def compute_flux_summary(model_or_unit, inhibition_map, reactions=None,
                         label=None, t_end=None, write_time=None):
    """Compute cumulative flux and inhibition strengths for a simulated run.

    For a reactor the window is 0..harvest time (``unit.tau``); the reactor
    simulates to ``tau_max`` but harvests at ``tau``, so integrating the whole
    stored trajectory would count flux the process never sees.

    Parameters
    ----------
    model_or_unit : KineticModel or NSKBatchReactor
        A model whose :meth:`~KineticModel.simulate` has run, or a reactor
        whose ``nsk_kinetic_model`` has (via the reactor's own simulation).
    inhibition_map : dict
        ``{parameter_id: (reaction_id, inhibitor_name)}``. Each parameter is a
        coefficient whose removal (set to zero) lifts that inhibitor's effect
        on that reaction. Denominator constants (e.g. ``K_6e``) are valid.
    reactions : sequence of str, optional
        Reactions to summarize. Defaults to every reaction in the model.
        Mapping entries naming a reaction outside this list are validated but
        otherwise ignored.
    label : str, optional
        Scenario label carried on the summary.
    t_end : float, optional
        Integrate only rows with ``time <= t_end``, truncating the trajectory
        before any integration; ``FluxSummary.t_end``, ``final_volume`` and
        ``final_concentrations`` then come from the last included row.
        Defaults to the reactor's harvest time ``unit.tau`` when
        ``model_or_unit`` is a reactor, and to the whole trajectory for a bare
        :class:`~nskinetics.KineticModel`.
    write_time : bool, optional
        Whether to write ``'time'`` back on every replayed row. ``None`` (the
        default) decides automatically from the model's rate laws; ``True`` or
        ``False`` forces it. See the Notes.

    Returns
    -------
    FluxSummary

    Raises
    ------
    ValueError
        If a mapping/`reactions` entry names a reaction or parameter absent
        from the model, or if no trajectory row satisfies ``time <= t_end``.
        Raised before any model state is written.
    KineticSimulationError
        If the model has no results, if the trajectory lacks a required state
        column (see :meth:`~KineticModel.state_selections`), or if the final
        broth volume is zero or non-finite.

    Notes
    -----
    The trajectory is replayed by writing each row's state back into the
    model and reading its reaction rates. ``'time'`` is written back only when
    a rate law actually reads it -- directly through the SBML time csymbol, or
    indirectly through an assignment-rule variable that depends on time (see
    :func:`_kinetic_laws_use_time`). A model with compiled native events (as
    the shipped one has) would otherwise have its time-triggered events
    re-armed by the backwards time jump at the start of every counterfactual
    pass. ``'time'`` is snapshotted and restored either way, and is always
    recorded as a column, so the restore contract is unchanged.

    ``write_time`` is the escape hatch from that detection, which reads the
    model's SBML text: pass ``True`` to force the write-back for a model whose
    time dependence the parser cannot see, or ``False`` to suppress it.
    """
    km = _resolve_model(model_or_unit)
    unit = model_or_unit if km is not model_or_unit else None
    df = km.results_df
    if df is None:
        raise KineticSimulationError(
            'No simulation results found; call simulate() first.')
    r = km._te
    rxn_ids = list(r.getReactionIds())
    if reactions is None:
        reactions = list(rxn_ids)
    else:
        reactions = list(reactions)

    # --- validate before touching state ---
    bad_rxn = [x for x in reactions if x not in rxn_ids]
    bad_rxn += [rid for (_p, (rid, _i)) in inhibition_map.items()
                if rid not in rxn_ids]
    if bad_rxn:
        raise ValueError(f'Unknown reaction id(s): {sorted(set(bad_rxn))}')
    params = set(r.getGlobalParameterIds())
    bad_par = [p for p in inhibition_map if p not in params]
    if bad_par:
        raise ValueError(f'Unknown parameter id(s): {sorted(bad_par)}')

    state_cols = km.state_selections()
    missing = [c for c in state_cols if c not in df.columns]
    if missing:
        raise KineticSimulationError(
            'Trajectory is missing state columns required for re-evaluation: '
            f'{missing}. Record KineticModel.state_selections() (an '
            'NSKBatchReactor does this automatically).')

    # Truncate to the integration window BEFORE integrating: a reactor
    # simulates to tau_max but harvests at tau, so everything downstream
    # (integrals, final volume, final concentrations, t_end) must see only the
    # rows up to the harvest row.
    if t_end is None and unit is not None:
        t_end = getattr(unit, 'tau', None)
    if t_end is not None:
        df = df.iloc[:_end_index(df, t_end, unit) + 1]

    # `time` is part of the snapshot/restore contract regardless, but it is
    # only WRITTEN back per row when a rate law actually reads it: the shipped
    # model carries compiled native SBML events whose triggers depend on time,
    # and jumping time backwards outside an integration step re-arms them.
    ordered_cols = _write_order(state_cols)
    if write_time is None:
        write_time = _kinetic_laws_use_time(r)
    write_cols = (ordered_cols if write_time
                  else [c for c in ordered_cols if c != 'time'])
    t = df['time'].to_numpy()
    n = len(df)
    # Hoist the per-column arrays out of the write-back loops: one numpy
    # lookup per row per column instead of a DataFrame scalar access.
    arrs = {c: df[c].to_numpy() for c in write_cols}
    idx_of = {rid: rxn_ids.index(rid) for rid in reactions}

    # inhibitor order (first-seen)
    inhibitors = []
    for _p, (_rid, inh) in inhibition_map.items():
        if inh not in inhibitors:
            inhibitors.append(inh)

    # group params by reaction then inhibitor
    by_rxn = defaultdict(lambda: defaultdict(list))
    for p, (rid, inh) in inhibition_map.items():
        by_rxn[rid][inh].append(p)

    snap = {c: r[c] for c in ordered_cols}
    snap_par = {p: r[p] for p in inhibition_map}

    def _restore_model():
        _restore(r, snap, ordered_cols)
        for p, v in snap_par.items():
            r[p] = v

    try:
        cumulative_mass = _cum(r, arrs, write_cols, idx_of, n, t)
        fraction_lost = {}
        fraction_lost_all = {}
        for rid in reactions:
            inhdict = by_rxn.get(rid)
            if not inhdict:
                continue
            fl = {}
            all_params = []
            for inh, plist in inhdict.items():
                all_params += plist
                saved = {p: r[p] for p in plist}
                for p in plist:
                    r[p] = 0.0
                cf = _cum(r, arrs, write_cols, {rid: idx_of[rid]}, n, t)[rid]
                for p, v in saved.items():
                    r[p] = v
                fl[inh] = _frac(cumulative_mass[rid], cf)
            fraction_lost[rid] = fl
            saved = {p: r[p] for p in all_params}
            for p in all_params:
                r[p] = 0.0
            cf_all = _cum(r, arrs, write_cols, {rid: idx_of[rid]}, n, t)[rid]
            for p, v in saved.items():
                r[p] = v
            fraction_lost_all[rid] = _frac(cumulative_mass[rid], cf_all)
    except BaseException as exc:
        # The model is restored on the way out either way, but a failure in
        # the restore must chain from the original error rather than mask it.
        try:
            _restore_model()
        except BaseException as restore_exc:
            raise restore_exc from exc
        raise
    _restore_model()

    # env is the model's broth-volume compartment; fall back to 1.0
    final_volume = float(df['env'].iloc[-1]) if 'env' in df.columns else 1.0
    if not np.isfinite(final_volume) or final_volume == 0.:
        raise KineticSimulationError(
            f'Final broth volume is {final_volume}; cumulative flux per litre '
            'is undefined. Check the trajectory\'s "env" column.')
    cumulative_flux = {rid: cumulative_mass[rid] / final_volume
                       for rid in reactions}
    final_conc = {}
    for s in r.getFloatingSpeciesIds():
        col = f'[{s}]'
        if col in df.columns:
            final_conc[s] = float(df[col].iloc[-1])
    return FluxSummary(
        label=label, reaction_ids=reactions,
        cumulative_mass=cumulative_mass, cumulative_flux=cumulative_flux,
        fraction_lost=fraction_lost, fraction_lost_all=fraction_lost_all,
        final_volume=final_volume, final_concentrations=final_conc,
        t_end=float(t[-1]), inhibitors=inhibitors)
