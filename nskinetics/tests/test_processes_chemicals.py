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

import json

from nskinetics.processes import create_sugar_prep_and_fermentation_chemicals

from .isobutanol._support import assert_ran_ok, run_isolated

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


_PROBE_SENTINEL = 'NSK_MODULES '


def test_chemicals_do_not_pull_heavy_imports():
    """Importing the chemicals factory AND building the chemical set loads
    neither the shipped kinetic model, nor the upstream ``biorefineries``
    package, nor the isobutanol biorefinery. Probed in a fresh interpreter
    (run from the repository root), because this pytest process's
    ``sys.modules`` holds whatever any other test imported. (``hensmith`` is
    deliberately not checked: biosteam >= 2.54 imports it whenever it is
    installed.)"""
    code = ('from nskinetics.processes import create_sugar_prep_and_fermentation_chemicals\n'
            'chems = create_sugar_prep_and_fermentation_chemicals()\n'
            "assert 'Yeast' in chems.IDs\n"
            'import json, sys\n'
            f'print({_PROBE_SENTINEL!r} + json.dumps(sorted(sys.modules)))\n')
    proc = run_isolated(['-c', code], timeout=600, simulation=False)
    out = assert_ran_ok(proc, _PROBE_SENTINEL)
    line = [ln for ln in out.splitlines() if ln.startswith(_PROBE_SENTINEL)][-1]
    modules = json.loads(line[len(_PROBE_SENTINEL):])
    assert 'nskinetics.processes' in modules
    assert not [m for m in modules if 's_cerevisiae_ferm_fb_inhib_mod_ibo' in m]
    assert not [m for m in modules
                if m == 'biorefineries' or m.startswith('biorefineries.')]
    assert not [m for m in modules if m.startswith('nskinetics.biorefineries.isobutanol')]
