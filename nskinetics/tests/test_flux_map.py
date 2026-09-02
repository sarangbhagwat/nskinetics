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
model and the corn-biorefinery chemical set. The module is deliberately not
registered in ``nskinetics/tests/__init__.py``, so a plain ``import
nskinetics`` never pulls any of that in. Every heavy import nonetheless lives
inside the fixture/test bodies: pytest imports this module during collection
even for a run that deselects it, and a module-level
``biosteam``/``biorefineries`` import at that point would pollute
``sys.modules`` and break the import-guard tests collected alongside it.
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
        apply_scenario_A(te_r)
        assert r.k_13 == 0.0 and r.k_14 == 0.0 and r.k_15 == 0.0
        assert r.k_16 == 0.0 and r.k_16r == 0.0
    finally:
        apply_scenario_A(te_r)
        te_r.reset()


def test_shipped_model_needs_no_time_write_back():
    # `prod_EtOH := s_EtOH/time` is the model's only time-dependent assignment
    # rule and no rate law reads it, so the replay must leave `time` alone --
    # the model's compiled native events trigger on it.
    from nskinetics.engine import flux_analysis as fa
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r
    r = te_r._te
    assert 'prod_EtOH' in fa._time_dependent_variables(r.getCurrentSBML())
    assert fa._kinetic_laws_use_time(r) is False


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
    # the computed reaction list is the drawn edges plus mapped-but-undrawn
    # reactions (r10, biomass decay), and every one of them is a real reaction
    reactions = FLUX_MAP_SPEC.reactions
    assert reactions[:len(FLUX_MAP_SPEC.edges)] == list(FLUX_MAP_SPEC.edges)
    assert 'r10' in reactions
    for rid in reactions:
        assert rid in rxn, f'reaction {rid} not in model'


def _stoichiometry(r, species, reaction):
    """Stoichiometric coefficient of ``species`` in ``reaction``, from the model."""
    import numpy as np
    m = r.getFullStoichiometryMatrix()
    rows = list(getattr(m, 'rownames', None) or r.getFloatingSpeciesIds())
    cols = list(getattr(m, 'colnames', None) or r.getReactionIds())
    return float(np.asarray(m)[rows.index(species), cols.index(reaction)])


def _balance_error(V406, summary, species, reaction):
    """Relative gap between a product's cumulative flux and its accumulation.

    ``cumulative_mass[reaction] * stoichiometry`` is everything the model made
    of ``species`` inside the integration window; with no outflow reaction for
    it, that must equal the change in its total amount (concentration times
    broth volume) over the same window.
    """
    r = V406.nsk_kinetic_model._te
    df = V406.nsk_results_df.iloc[:V406.tau_index + 1]
    col = f'[{species}]'
    made = summary.cumulative_mass[reaction] * _stoichiometry(r, species,
                                                              reaction)
    accumulated = (df[col].iat[-1] * df['env'].iat[-1]
                   - df[col].iat[0] * df['env'].iat[0])
    return abs(made - accumulated) / abs(accumulated)


def test_product_balances_close_in_both_scenarios(simulated_V406):
    # The strongest available check that the re-evaluated rates really are the
    # model's own: what r6 (and, in scenario B, r16) made must equal what
    # accumulated in the broth, since nothing else produces or consumes them.
    from nskinetics import compute_flux_summary
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        FLUX_MAP_SPEC, apply_scenario_A, apply_scenario_B)
    V406 = simulated_V406
    model = V406.nsk_kinetic_model
    imap, reactions = FLUX_MAP_SPEC.inhibition_map, FLUX_MAP_SPEC.reactions
    try:
        apply_scenario_A(model)
        V406.system.simulate()
        sa = compute_flux_summary(V406, imap, reactions=reactions)
        assert _balance_error(V406, sa, 's_EtOH', 'r6') < 1e-3
        apply_scenario_B(model)
        V406.system.simulate()
        sb = compute_flux_summary(V406, imap, reactions=reactions)
        assert _balance_error(V406, sb, 's_EtOH', 'r6') < 1e-3
        assert _balance_error(V406, sb, 's_IBO', 'r16') < 1e-3
    finally:
        apply_scenario_A(model)
        V406.system.simulate()


def test_compute_flux_summary_restores_the_shipped_model(simulated_V406):
    # Every selection that defines integrator state, and every coefficient the
    # counterfactuals zero, must come back exactly as it was found.
    from nskinetics import compute_flux_summary
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        FLUX_MAP_SPEC)
    V406 = simulated_V406
    km = V406.nsk_kinetic_model
    r = km._te
    before = {c: r[c] for c in km.state_selections()}
    before_params = {p: r[p] for p in FLUX_MAP_SPEC.inhibition_map}
    assert 'time' in before and before_params
    compute_flux_summary(V406, FLUX_MAP_SPEC.inhibition_map,
                         reactions=FLUX_MAP_SPEC.reactions)
    for c, v in before.items():
        assert r[c] == v, f'{c} not restored'
    for p, v in before_params.items():
        assert r[p] == v, f'{p} not restored'


def test_helper_rejects_a_reactor_carrying_a_foreign_model():
    # The presets are applied to the reactor's OWN kinetic model, so a reactor
    # running something else must be refused rather than silently mutating the
    # shipped te_r and drawing a figure of the wrong model. The stand-in model
    # is a stub, not a second tellurium model: loading another Antimony model
    # into this process perturbs the shipped model's own integration (it
    # harvests several hours earlier afterwards), which would poison every
    # test that runs after this one.
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        draw_scenario_flux_map)

    class _ForeignRoadRunner:
        @staticmethod
        def getGlobalParameterIds():
            return ['k1', 'K1']

    class _ForeignModel:
        _te = _ForeignRoadRunner()

    class _FakeUnit:
        nsk_kinetic_model = _ForeignModel()

    with pytest.raises(ValueError, match='k_13'):
        draw_scenario_flux_map(_FakeUnit())


def test_end_to_end_both_scenarios(tmp_path, simulated_V406):
    import matplotlib.pyplot as plt
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        draw_scenario_flux_map)
    V406 = simulated_V406
    fig, axes, (sa, sb) = draw_scenario_flux_map(V406, save_dir=str(tmp_path))
    try:
        assert len(axes) == 2
        assert os.path.exists(os.path.join(tmp_path, 'flux_map.pdf'))
        assert sa.cumulative_flux['r6'] > 0     # ethanol pathway active in A
        assert sb.cumulative_flux['r16'] > 0    # isobutanol pathway active in B
        assert sa.cumulative_flux['r16'] == 0.0  # Ehrlich off in A
        # the integration window ends at the harvest row, not at tau_max
        assert sa.t_end == V406.tau
        assert sa.t_end < V406.nsk_results_df['time'].iat[-1]
        assert sa.final_concentrations['s_EtOH'] == pytest.approx(
            V406.nsk_results_specific_tau_dict['[s_EtOH]'], rel=1e-6)
        # r10 has no edge but is mapped, so the shipped path still computes it;
        # ethanol exceeds P_10e, so decay is ENHANCED (negative fraction lost)
        assert 'r10' in sa.fraction_lost
        assert sa.fraction_lost['r10']['ethanol'] < 0
    finally:
        plt.close(fig)
