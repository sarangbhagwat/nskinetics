# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Kinetic-parameter classification for the shipped *S. cerevisiae*
ethanol/isobutanol model, and a helper that names a parameter change.

Every kinetic parameter has a *role* (what kind of term it is in its rate
law) and belongs to one or more pathway *modules* (through the reactions it
appears in). A change to a set of parameters is then described as direction
+ role + effector + modules, e.g. ``stronger isobutanol inhibition of
glycolysis/fermentation, overflow/acetate and growth``.

Operation parameters -- fed-batch feeding, aeration staging, dilution -- are
not kinetics and are excluded from descriptions; they are enumerated in
:data:`OPERATION_PARAMETERS` only so the classification can be proven
complete against the model. RoadRunner also lists assignment-rule outputs
(``qO2``, the ``y_*`` yields, the ``curr_*`` mirrors) as global parameters;
those are derived quantities and are neither listed nor described.

This module is pure Python and imports nothing from the package, so its
tables can be read without loading the model.
"""

import math
from dataclasses import dataclass

__all__ = ('ROLES', 'MODULES', 'MODULE_LABELS', 'REACTION_MODULES',
           'ParameterInfo', 'KINETIC_PARAMETERS', 'OPERATION_PARAMETERS',
           'snapshot_parameters', 'diff_parameters')

#: Parameter roles, in the order their clauses appear in a description.
ROLES = ('capacity', 'affinity', 'substrate_regulation',
         'product_inhibition', 'product_self_inhibition',
         'lethality', 'lethality_threshold', 'initial_state')

#: Pathway modules, in display order.
MODULES = ('glycolysis_fermentation', 'respiration', 'overflow_acetate',
           'growth', 'physiological_state', 'ehrlich')

#: Display strings for :data:`MODULES`.
MODULE_LABELS = {
    'glycolysis_fermentation': 'glycolysis/fermentation',
    'respiration': 'respiration',
    'overflow_acetate': 'overflow/acetate',
    'growth': 'growth',
    'physiological_state': 'physiological state',
    'ehrlich': 'Ehrlich branch',
}

#: ``{reaction_id: module}`` for the kinetic reactions. The continuous-mode
#: dilution reactions (``s_glu_in``, ``*_out``; D = 0 in fed-batch use) are
#: deliberately absent.
REACTION_MODULES = {
    'r1': 'glycolysis_fermentation',   # glycolysis
    'r3': 'glycolysis_fermentation',   # pyruvate decarboxylase
    'r6': 'glycolysis_fermentation',   # alcohol dehydrogenase
    'r2': 'respiration',               # pyruvate -> TCA
    'r5': 'respiration',               # acetate -> TCA
    'r4': 'overflow_acetate',          # acetaldehyde -> acetate
    'r7': 'growth',                    # growth on glucose
    'r8': 'growth',                    # growth on acetate
    'r9': 'physiological_state',       # AcDH synthesis
    'r10': 'physiological_state',      # active-biomass decay
    'r11': 'physiological_state',      # AcDH degradation
    'r13': 'ehrlich',                  # acetolactate synthase
    'r14': 'ehrlich',                  # KARI
    'r15': 'ehrlich',                  # DHAD
    'r16': 'ehrlich',                  # KDC + ADH
}


@dataclass(frozen=True)
class ParameterInfo:
    """Classification of one kinetic parameter.

    Parameters
    ----------
    role : str
        One of :data:`ROLES`.
    reactions : tuple of str
        Reaction ids the parameter appears in (empty for initial-state
        parameters, which sit in rate rules rather than rate laws).
    modules : tuple of str
        Pathway modules, derived from ``reactions`` through
        :data:`REACTION_MODULES` unless given explicitly.
    effector : str or None
        The product or substrate the parameter couples the reaction to:
        ``'ethanol'``, ``'acetate'``, ``'isobutanol'`` for product terms,
        ``'glucose'`` or ``'acetaldehyde'`` for substrate regulation,
        ``None`` for plain capacities, affinities and initial state.
    """
    role: str
    reactions: tuple
    modules: tuple
    effector: str = None


def _p(role, *reactions, effector=None, modules=None):
    if modules is None:
        seen = []
        for rxn in reactions:
            m = REACTION_MODULES[rxn]
            if m not in seen:
                seen.append(m)
        modules = tuple(seen)
    return ParameterInfo(role, tuple(reactions), tuple(modules), effector)


#: ``{parameter_id: ParameterInfo}`` -- the source of truth. Insertion order
#: is the order members are listed inside a verbose clause.
KINETIC_PARAMETERS = {
    # --- capacity: multiplies the whole rate term (Vmax-like) ---
    'k_1h': _p('capacity', 'r1'),
    'k_1l': _p('capacity', 'r1'),
    'k_1e': _p('capacity', 'r1'),
    'k_2': _p('capacity', 'r2'),
    'k_3': _p('capacity', 'r3'),
    'k_4': _p('capacity', 'r4'),
    'k_5': _p('capacity', 'r5'),
    'k_5e': _p('capacity', 'r5'),
    'k_6': _p('capacity', 'r6'),
    'k_7': _p('capacity', 'r7'),
    'k_8': _p('capacity', 'r8'),
    'k_9': _p('capacity', 'r9'),
    'k_9e': _p('capacity', 'r9'),
    'k_9c': _p('capacity', 'r9'),
    'k_10': _p('capacity', 'r10'),
    'k_11': _p('capacity', 'r11'),
    'k_13': _p('capacity', 'r13'),
    'k_14': _p('capacity', 'r14'),
    'k_15': _p('capacity', 'r15'),
    'k_16': _p('capacity', 'r16'),
    # --- affinity: S/(S + K) saturation constant ---
    'K_1h': _p('affinity', 'r1'),
    'K_1l': _p('affinity', 'r1'),
    'K_1e': _p('affinity', 'r1'),
    'K_2': _p('affinity', 'r2'),
    'K_3': _p('affinity', 'r3'),
    'K_4': _p('affinity', 'r4'),
    'K_5': _p('affinity', 'r5'),
    'K_5e': _p('affinity', 'r5', 'r8'),
    'K_6': _p('affinity', 'r6'),
    'K_7': _p('affinity', 'r7'),
    'K_9': _p('affinity', 'r9'),
    'K_9e': _p('affinity', 'r9'),
    'K_13': _p('affinity', 'r13'),
    'K_14': _p('affinity', 'r14'),
    'K_15': _p('affinity', 'r15'),
    'K_16': _p('affinity', 'r16'),
    # --- substrate regulation: 1/(1 + K*S) repression or signal saturation ---
    'K_1i': _p('substrate_regulation', 'r1', effector='acetaldehyde'),
    'K_2i': _p('substrate_regulation', 'r2', effector='glucose'),
    'K_5i': _p('substrate_regulation', 'r5', 'r8', effector='glucose'),
    'K_9i': _p('substrate_regulation', 'r9', effector='glucose'),
    # --- product inhibition (sublethal): exp(-k*P) ---
    'k_1ie': _p('product_inhibition', 'r1', effector='ethanol'),
    'k_1ia': _p('product_inhibition', 'r1', effector='acetate'),
    'k_1ii': _p('product_inhibition', 'r1', effector='isobutanol'),
    'k_4ie': _p('product_inhibition', 'r4', effector='ethanol'),
    'k_4ia': _p('product_inhibition', 'r4', effector='acetate'),
    'k_4ii': _p('product_inhibition', 'r4', effector='isobutanol'),
    'k_6ia': _p('product_inhibition', 'r6', effector='acetate'),
    'k_6ii': _p('product_inhibition', 'r6', effector='isobutanol'),
    'k_7ie': _p('product_inhibition', 'r7', effector='ethanol'),
    'k_7ia': _p('product_inhibition', 'r7', effector='acetate'),
    'k_7ii': _p('product_inhibition', 'r7', effector='isobutanol'),
    'k_16ie': _p('product_inhibition', 'r16', effector='ethanol'),
    'k_16ia': _p('product_inhibition', 'r16', effector='acetate'),
    # --- product self-inhibition: K_*e*P denominator and k_*r*P reverse term ---
    'K_6e': _p('product_self_inhibition', 'r6', effector='ethanol'),
    'k_6r': _p('product_self_inhibition', 'r6', effector='ethanol'),
    'K_16i': _p('product_self_inhibition', 'r16', effector='isobutanol'),
    'k_16r': _p('product_self_inhibition', 'r16', effector='isobutanol'),
    # --- lethality: steepness of the threshold-gated exp on r10 ---
    'k_10ie': _p('lethality', 'r10', effector='ethanol'),
    'k_10ia': _p('lethality', 'r10', effector='acetate'),
    'k_10ii': _p('lethality', 'r10', effector='isobutanol'),
    # --- lethality threshold: onset concentration of the r10 enhancement ---
    'P_10e': _p('lethality_threshold', 'r10', effector='ethanol'),
    'P_10a': _p('lethality_threshold', 'r10', effector='acetate'),
    'P_10i': _p('lethality_threshold', 'r10', effector='isobutanol'),
    # --- initial state: rate-rule variables with explicit initial values ---
    'X_a': _p('initial_state', modules=('physiological_state',)),
    'X_AcDH': _p('initial_state', modules=('physiological_state',)),
}

#: Parameters that describe how the reactor is run, not the kinetics.
#: Skipped by the diff/description helpers; listed so the coverage test can
#: assert every model parameter is accounted for.
OPERATION_PARAMETERS = frozenset({
    'D', 'is_aerobic', 'stage_1_max_time', 'stage_1_max_x',
    'anaerobic_growth_mult',
    'max_n_glu_spikes', 'n_glu_spikes',
    'threshold_conc_glu_spike', 'target_conc_glu_spike',
    'conc_glu_feed_spike',
    'last_vol_glu_feed_added', 'tot_vol_glu_feed_added',
    'glucose_feed_spikeDelay', 'glucose_feed_spike_dissolveDelay',
})


# --- helpers ----------------------------------------------------------------

def snapshot_parameters(model):
    """Read the current value of every kinetic parameter off a model.

    Parameters
    ----------
    model : KineticModel or roadrunner.RoadRunner
        The shipped *S. cerevisiae* ethanol/isobutanol model (a
        :class:`~nskinetics.KineticModel` is unwrapped via ``_te``).

    Returns
    -------
    dict
        ``{parameter_id: float}`` for every id in :data:`KINETIC_PARAMETERS`.

    Raises
    ------
    ValueError
        If the model lacks any kinetic parameter, i.e. it is not the shipped
        model.

    Notes
    -----
    ``X_a`` and ``X_AcDH`` are rate-rule variables, so their values evolve
    during a simulation. Take the baseline and the comparison snapshot at the
    same point of a run (typically right after :meth:`~nskinetics.KineticModel.reset`).
    """
    r = getattr(model, '_te', model)
    ids = set(r.getGlobalParameterIds())
    missing = [p for p in KINETIC_PARAMETERS if p not in ids]
    if missing:
        raise ValueError(
            f'{missing[0]!r} is not a parameter of this model; '
            'snapshot_parameters expects the shipped S. cerevisiae '
            'ethanol/isobutanol model.')
    return {p: float(r[p]) for p in KINETIC_PARAMETERS}


def _as_values(obj):
    """Return ``obj`` as a plain ``{id: value}`` dict, snapshotting a model."""
    if isinstance(obj, dict):
        return dict(obj)
    return snapshot_parameters(obj)


def _check_known(values, side):
    unknown = [k for k in values
               if k not in KINETIC_PARAMETERS and k not in OPERATION_PARAMETERS]
    if unknown:
        raise ValueError(
            f'Unknown parameter(s) in {side}: {unknown}; expected kinetic or '
            'operation parameters of the shipped model.')


def diff_parameters(baseline, current):
    """Return the kinetic parameters whose value differs between two states.

    Parameters
    ----------
    baseline, current : dict or KineticModel or roadrunner.RoadRunner
        ``{parameter_id: value}`` mappings; a model is snapshotted with
        :func:`snapshot_parameters` first. ``current`` may be a partial
        override dict (e.g. ``SCENARIO_B_EHRLICH``): kinetic keys absent from
        it are treated as unchanged. Keys in :data:`OPERATION_PARAMETERS` are
        skipped on either side.

    Returns
    -------
    dict
        ``{parameter_id: (old, new)}`` in :data:`KINETIC_PARAMETERS` order.
        Values within ``math.isclose(rel_tol=1e-12)`` are unchanged.

    Raises
    ------
    ValueError
        If either side carries a key that is neither a kinetic nor an
        operation parameter, or if ``current`` carries a kinetic key that
        ``baseline`` lacks.
    """
    base = _as_values(baseline)
    cur = _as_values(current)
    _check_known(base, 'baseline')
    _check_known(cur, 'current')
    missing = [k for k in cur if k in KINETIC_PARAMETERS and k not in base]
    if missing:
        raise ValueError(
            f'{missing} present in current but not in baseline; the baseline '
            'must carry every parameter being compared.')
    out = {}
    for p in KINETIC_PARAMETERS:
        if p in base and p in cur:
            old, new = float(base[p]), float(cur[p])
            if not math.isclose(old, new, rel_tol=1e-12, abs_tol=0.0):
                out[p] = (old, new)
    return out
