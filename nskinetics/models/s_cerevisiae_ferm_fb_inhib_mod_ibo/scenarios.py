# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Kinetic-only scenario presets for the shipped *S. cerevisiae*
ethanol/isobutanol model.

The two configurations differ only in the engineered Ehrlich pathway
(r13-r17): scenario A leaves it off (all four rate constants -- k_13, k_14,
k_15, k_16 -- zero, so no isobutanol is made); scenario B turns it on. Every
product-inhibition coefficient is already at its scenario-B value in the
shipped antimony (it has no effect in A because isobutanol stays zero), so
only the r13-r16 rate constants change here. r16 is irreversible
Michaelis-Menten in KIV (since 2026-09-13), so there is no reverse term to
switch.

``k_17``, the Adh6 capacity of the reductase r17 that the 2026-09-15 split
broke out of the lumped r16, is deliberately NOT a scenario key. Adh6 is
native and constitutive, so it keeps its nonzero default in both scenarios;
the branch is gated at its KIV entry by ``k_16``, and with ``k_16 = 0`` no
isobutyraldehyde is ever made for r17 to reduce. This is also what lets a
caller that knows only the four historical capacities (the isobutanol
workbooks) still switch the branch on correctly.

Values mirror ``parameter-distributions_corn_IBO_EtOH_B.xlsx`` in the
(read-only) isobutanol biorefinery, which is the source of truth. The
fed-batch feeding strategy (spike count, thresholds) differs between scenarios
too, but that is a caller concern and is NOT set here.
"""

__all__ = ('apply_scenario_A', 'apply_scenario_B',
           'SCENARIO_B_EHRLICH', 'SCENARIO_A_EHRLICH')

# r13-r16 rate constants (K_* saturation constants are already at their
# scenario-B values in the shipped antimony and are unchanged; k_17/K_17 and
# the rest of the Adh6 set are constitutive, not scenario-switched).
SCENARIO_B_EHRLICH = {
    'k_13': 5.81, 'k_14': 4.8, 'k_15': 4.8, 'k_16': 2.82,
}
SCENARIO_A_EHRLICH = {
    'k_13': 0.0, 'k_14': 0.0, 'k_15': 0.0, 'k_16': 0.0,
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
