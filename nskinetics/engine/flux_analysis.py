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

from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

from ..exceptions import KineticSimulationError

__all__ = ('FluxSummary', 'compute_flux_summary')


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
        ``{reaction_id: g/L of final broth}`` (``cumulative_mass`` / final volume).
    fraction_lost : dict
        ``{reaction_id: {inhibitor: fraction}}``; ``1 - actual/uninhibited``.
        Negative for enhancement (e.g. product-accelerated decay).
    fraction_lost_all : dict
        ``{reaction_id: fraction}`` with every mapped coefficient zeroed at once.
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


def _apply_row(r, df, state_cols, i):
    # Set non-concentration state (compartments, params, rate-rule vars) first,
    # then concentrations: writing a compartment volume rescales the
    # concentrations derived from held amounts, so concentrations must be set
    # last to land at their recorded values.
    for c in state_cols:
        if not c.startswith('['):
            r[c] = df[c].iat[i]
    for c in state_cols:
        if c.startswith('['):
            r[c] = df[c].iat[i]


def _rates_along(r, df, state_cols, idx_of, n):
    out = {rid: np.empty(n) for rid in idx_of}
    for i in range(n):
        _apply_row(r, df, state_cols, i)
        rates = r.getReactionRates()
        for rid, j in idx_of.items():
            out[rid][i] = rates[j]
    return out


def _cum_one(r, df, state_cols, rid_index, n, t):
    vals = np.empty(n)
    for i in range(n):
        _apply_row(r, df, state_cols, i)
        vals[i] = r.getReactionRates()[rid_index]
    return float(np.trapezoid(vals, t))


def _frac(base, counterfactual):
    if counterfactual == 0.0:
        return 0.0
    return 1.0 - base / counterfactual


def compute_flux_summary(model_or_unit, inhibition_map, reactions=None,
                         label=None):
    """Compute cumulative flux and inhibition strengths for a simulated run.

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

    Returns
    -------
    FluxSummary

    Raises
    ------
    ValueError
        If a mapping/`reactions` entry names a reaction or parameter absent
        from the model. Raised before any model state is written.
    KineticSimulationError
        If the model has no results, or the trajectory lacks a required state
        column (see :meth:`~KineticModel.state_selections`).
    """
    km = _resolve_model(model_or_unit)
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

    state_no_time = [c for c in state_cols if c != 'time']
    t = df['time'].to_numpy()
    n = len(df)
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

    snap = {c: r[c] for c in state_no_time}
    snap_par = {p: r[p] for p in inhibition_map}
    try:
        base = _rates_along(r, df, state_cols, idx_of, n)
        cumulative_mass = {rid: float(np.trapezoid(base[rid], t))
                           for rid in reactions}
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
                cf = _cum_one(r, df, state_cols, idx_of[rid], n, t)
                for p, v in saved.items():
                    r[p] = v
                fl[inh] = _frac(cumulative_mass[rid], cf)
            fraction_lost[rid] = fl
            saved = {p: r[p] for p in all_params}
            for p in all_params:
                r[p] = 0.0
            cf_all = _cum_one(r, df, state_cols, idx_of[rid], n, t)
            for p, v in saved.items():
                r[p] = v
            fraction_lost_all[rid] = _frac(cumulative_mass[rid], cf_all)
    finally:
        for c, v in snap.items():
            r[c] = v
        for p, v in snap_par.items():
            r[p] = v

    # env is the model's broth-volume compartment; fall back to 1.0
    final_volume = float(df['env'].iloc[-1]) if 'env' in df.columns else 1.0
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
