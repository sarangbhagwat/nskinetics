# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Kinetic-only scenario presets for the shipped *S. cerevisiae*
ethanol/isobutanol model.

The two configurations differ only in the engineered Ehrlich pathway (r13-r16):
scenario A leaves it off (all five rate constants -- k_13, k_14, k_15, k_16,
k_16r -- zero, so no isobutanol is made); scenario B turns it on. Every product-inhibition coefficient is already at its
scenario-B value in the shipped antimony (it has no effect in A because
isobutanol stays zero), so only the r13-r16 rate constants change here.

Values mirror ``parameter-distributions_corn_IBO_EtOH_B.xlsx`` in the
(read-only) isobutanol biorefinery, which is the source of truth. The
fed-batch feeding strategy (spike count, thresholds) differs between scenarios
too, but that is a caller concern and is NOT set here.
"""

__all__ = ('apply_scenario_A', 'apply_scenario_B',
           'SCENARIO_B_EHRLICH', 'SCENARIO_A_EHRLICH')

# r13-r16 rate constants (K_* saturation/inhibition constants are already at
# their scenario-B values in the shipped antimony and are unchanged).
SCENARIO_B_EHRLICH = {
    'k_13': 5.81, 'k_14': 4.8, 'k_15': 4.8, 'k_16': 2.82, 'k_16r': 0.0125,
}
SCENARIO_A_EHRLICH = {
    'k_13': 0.0, 'k_14': 0.0, 'k_15': 0.0, 'k_16': 0.0, 'k_16r': 0.0,
}


def _apply(model, values):
    r = getattr(model, '_te', model)
    for name, val in values.items():
        r[name] = val
    return model


def apply_scenario_A(model):
    """Set the shipped model to scenario A (Ehrlich branch off).

    Parameters
    ----------
    model : KineticModel
        The model to mutate (typically the shipped ``te_r``).

    Returns
    -------
    KineticModel
        The same model, mutated.
    """
    return _apply(model, SCENARIO_A_EHRLICH)


def apply_scenario_B(model):
    """Set the shipped model to scenario B (Ehrlich branch on).

    Parameters
    ----------
    model : KineticModel
        The model to mutate (typically the shipped ``te_r``).

    Returns
    -------
    KineticModel
        The same model, mutated.
    """
    return _apply(model, SCENARIO_B_EHRLICH)
