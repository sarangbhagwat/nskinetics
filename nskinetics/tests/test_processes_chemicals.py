# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Fast tests for the sugar-prep + fermentation factory's own chemical set.

The set is deliberately factory-specific (future factories in
``nskinetics.processes`` will each ship their own ``create_<factory>_chemicals``
beside them), and deliberately light: no ``te_r`` kinetic model and no
``biorefineries`` import, so this module runs in the default
``-m "not slow"`` suite and in CI.
"""

import sys

from nskinetics.processes import create_sugar_prep_and_fermentation_chemicals

EXPECTED_IDS = {'Water', 'Ethanol', 'Isobutanol', 'AceticAcid', 'Glucose',
                'Sucrose', 'CO2', 'O2', 'N2', 'NH3', 'Yeast'}


def test_chemicals_compile_and_contents():
    chems = create_sugar_prep_and_fermentation_chemicals()
    assert EXPECTED_IDS <= set(chems.IDs)
    Yeast = chems.Yeast
    assert Yeast.phase_ref == 's'
    assert Yeast.atoms == {'C': 1.0, 'H': 1.61, 'O': 0.56}
    # Synonyms used by the factory's units and the biorefinery convention.
    assert chems.H2O is chems.Water
    assert chems.DryYeast is chems.Yeast
    # Yeast thermo mirrors the cane/corn-biorefinery definition: Hf tied to
    # glucose per unit mass (growth heats deliberately ignored).
    Glucose = chems.Glucose
    assert abs(Yeast.Hf - Glucose.Hf / Glucose.MW * Yeast.MW) < 1e-6


def test_chemicals_do_not_pull_heavy_imports():
    assert not [m for m in sys.modules
                if 's_cerevisiae_ferm_fb_inhib_mod_ibo' in m]
    assert 'biorefineries' not in sys.modules
