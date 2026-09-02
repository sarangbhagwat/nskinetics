# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Flux-map tests exercising a fully simulated ``NSKBatchReactor``.

Marked ``slow``: the module-scoped fixture builds and simulates the shipped
sugar-prep + fermentation system, which imports the shipped ``te_r`` kinetic
model and the corn-biorefinery chemical set. Every heavy import therefore
lives inside the fixture/test bodies -- ``nskinetics.tests`` is imported by
``nskinetics/__init__.py``, so a module-level ``biosteam``/``biorefineries``
import here would pollute ``sys.modules`` and break the import-guard tests.
"""

import os

import pytest

pytestmark = pytest.mark.slow


@pytest.fixture(scope='module')
def simulated_V406():
    pytest.importorskip('biorefineries')
    import biosteam as bst
    import thermosteam as tmo
    from biorefineries import corn
    chems = [c for c in corn.chemicals.create_chemicals()]
    chems.append(tmo.Chemical('Isobutanol'))
    chems.append(tmo.Chemical('AceticAcid'))
    tmo.settings.set_thermo(chems)
    bst.main_flowsheet.set_flowsheet('test_flux_map')
    from nskinetics.processes import create_sugar_prep_and_fermentation_system
    system = create_sugar_prep_and_fermentation_system()
    system.simulate()
    return system.flowsheet.unit.V406


def test_reactor_records_state_selection_columns(simulated_V406):
    V406 = simulated_V406
    km = V406.nsk_kinetic_model
    columns = list(V406.nsk_results_df.columns)
    cols = set(columns)
    for sel in km.state_selections():
        assert sel in cols, f'missing state column {sel}'
    # The state selections are appended, so the pre-existing columns must
    # survive unmoved: 'time' stays first and the mapped species remain.
    assert columns[0] == 'time'
    assert '[s_EtOH]' in cols
    # Every column carries a full trajectory, not just an endpoint.
    assert len(V406.nsk_results_df) > 1


def test_scenario_presets_set_exact_values():
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        te_r, apply_scenario_A, apply_scenario_B)
    r = te_r._te
    try:
        apply_scenario_B(te_r)
        assert r.k_13 == 5.81 and r.k_14 == 4.8 and r.k_15 == 4.8
        assert r.k_16 == 2.82 and r.k_16r == 0.0125
        # inhibition coefficients already at B in the shipped model
        assert r.k_7ii == 0.15 and r.k_1ii == 0.075
        apply_scenario_A(te_r)
        assert r.k_13 == 0.0 and r.k_14 == 0.0 and r.k_15 == 0.0
        assert r.k_16 == 0.0 and r.k_16r == 0.0
    finally:
        apply_scenario_A(te_r)
        te_r.reset()


def test_shipped_spec_maps_only_model_reactions_and_params():
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        te_r, FLUX_MAP_SPEC)
    r = te_r._te
    rxn = set(r.getReactionIds())
    par = set(r.getGlobalParameterIds())
    for rid in FLUX_MAP_SPEC.edges:
        assert rid in rxn, f'edge reaction {rid} not in model'
    for p, (rid, _inh) in FLUX_MAP_SPEC.inhibition_map.items():
        assert p in par, f'inhibition param {p} not in model'
        assert rid in rxn, f'inhibition reaction {rid} not in model'


def test_end_to_end_both_scenarios(tmp_path, simulated_V406):
    import matplotlib.pyplot as plt
    from nskinetics import compute_flux_summary
    from nskinetics.visualization import draw_flux_map
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        te_r, FLUX_MAP_SPEC, apply_scenario_A, apply_scenario_B)
    V406 = simulated_V406
    imap = FLUX_MAP_SPEC.inhibition_map
    rxns = list(FLUX_MAP_SPEC.edges)
    try:
        apply_scenario_A(te_r)
        V406.system.simulate()
        sa = compute_flux_summary(V406, imap, reactions=rxns, label='Scenario A')
        apply_scenario_B(te_r)
        V406.system.simulate()
        sb = compute_flux_summary(V406, imap, reactions=rxns, label='Scenario B')
    finally:
        apply_scenario_A(te_r)
        V406.system.simulate()
    fig, axes = draw_flux_map([sa, sb], FLUX_MAP_SPEC, save_dir=str(tmp_path))
    assert len(axes) == 2
    assert os.path.exists(os.path.join(tmp_path, 'flux_map.pdf'))
    assert sa.cumulative_flux['r6'] > 0        # ethanol pathway active in A
    assert sb.cumulative_flux['r16'] > 0       # isobutanol pathway active in B
    assert sa.cumulative_flux['r16'] == 0.0    # Ehrlich off in A
    plt.close(fig)
