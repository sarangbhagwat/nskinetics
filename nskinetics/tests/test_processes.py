# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Structure tests for the ``nskinetics.processes`` sugar-prep + fermentation
system factory.

Marked ``slow``: building the factory with defaults imports the shipped
``te_r`` kinetic model (which configures and simulates itself on import) and
needs the corn-biorefinery chemical set, so these tests are excluded from the
default ``-m "not slow"`` run and from CI.
"""

import pytest

pytestmark = pytest.mark.slow

EXPECTED_UNIT_IDS = {
    'S301',
    'F301', 'F301_P0', 'F301_P1', 'M301', 'H301',
    'F302', 'F302_P0', 'F302_P1', 'M302', 'H302',
    'V406',
    'K330', 'V330',
}

EXPECTED_MAP = {
    '[s_glu]': 'Glucose',
    '[x]': 'Yeast',
    '[s_EtOH]': 'Ethanol',
    '[s_IBO]': 'Isobutanol',
    '[s_acetate]': 'AceticAcid',
}

EXPECTED_TRACK_VARS = [
    'y_EtOH_glu_added',
    'y_EtOH_glu_consumed',
    'y_IBO_glu_added',
    'y_IBO_glu_consumed',
    'y_EtOH_IBO_glu_added',
    'curr_n_glu_spikes',
    'curr_a',
    'prod_EtOH',
    'curr_tot_vol_glu_feed_added',
    'curr_env',
]


@pytest.fixture(scope='module')
def factory_system():
    """Build the factory once, with defaults, against the same chemical set
    the isobutanol biorefinery uses (corn chemicals + Isobutanol + AceticAcid)."""
    pytest.importorskip('biorefineries')
    import biosteam as bst
    import thermosteam as tmo
    from biorefineries import corn

    chems = [c for c in corn.chemicals.create_chemicals()]
    chems.append(tmo.Chemical('Isobutanol'))
    chems.append(tmo.Chemical('AceticAcid'))
    tmo.settings.set_thermo(chems)
    bst.main_flowsheet.set_flowsheet('test_processes_factory')

    from nskinetics.processes import create_sugar_prep_and_fermentation_system
    return create_sugar_prep_and_fermentation_system()


def _units_by_id(system):
    return {u.ID: u for u in system.units}


def test_unit_ids_and_types(factory_system):
    import biosteam as bst
    from nskinetics.units import FermentationSaccharomycesEthanolIsobutanol
    u = _units_by_id(factory_system)
    assert set(u) == EXPECTED_UNIT_IDS
    assert isinstance(u['S301'], bst.Splitter)
    assert isinstance(u['F301'], bst.MultiEffectEvaporator)
    assert isinstance(u['F302'], bst.MultiEffectEvaporator)
    for pump in ('F301_P0', 'F301_P1', 'F302_P0', 'F302_P1'):
        assert isinstance(u[pump], bst.units.Pump)
    assert isinstance(u['M301'], bst.units.Mixer)
    assert isinstance(u['M302'], bst.units.Mixer)
    assert isinstance(u['H301'], bst.units.HXutility)
    assert isinstance(u['H302'], bst.units.HXutility)
    assert isinstance(u['V406'], FermentationSaccharomycesEthanolIsobutanol)
    assert isinstance(u['K330'], bst.units.IsothermalCompressor)
    assert isinstance(u['V330'], bst.units.IsenthalpicValve)
    assert u['V330'].line == 'Valve'


def test_connectivity(factory_system):
    u = _units_by_id(factory_system)
    # Factory inlets: saccharified slurry -> S301; seed -> V406 ins[1].
    assert u['S301'].ins[0] is factory_system.ins[0]
    assert u['V406'].ins[1] is factory_system.ins[1]
    # Initial-feed train: S301-0 -> F301 -> P0 -> M301 -> H301 -> V406 ins[0].
    assert u['F301'].ins[0] is u['S301'].outs[0]
    assert u['F301_P0'].ins[0] is u['F301'].outs[0]
    assert u['F301_P1'].ins[0] is u['F301'].outs[1]
    assert u['M301'].ins[0] is u['F301_P0'].outs[0]
    assert u['H301'].ins[0] is u['M301'].outs[0]
    assert u['V406'].ins[0] is u['H301'].outs[0]
    # Spike train: S301-1 -> F302 -> P0 -> M302 -> H302 -> V406 ins[2].
    assert u['F302'].ins[0] is u['S301'].outs[1]
    assert u['F302_P0'].ins[0] is u['F302'].outs[0]
    assert u['F302_P1'].ins[0] is u['F302'].outs[1]
    assert u['M302'].ins[0] is u['F302_P0'].outs[0]
    assert u['H302'].ins[0] is u['M302'].outs[0]
    assert u['V406'].ins[2] is u['H302'].outs[0]
    # Aeration loop: K330 -> V330 -> V406 ins[3].
    assert u['V330'].ins[0] is u['K330'].outs[0]
    assert u['V406'].ins[3] is u['V330'].outs[0]
    # Factory outlets: vent, effluent, and both evaporator condensates.
    assert factory_system.outs[0] is u['V406'].outs[0]
    assert factory_system.outs[1] is u['V406'].outs[1]
    assert factory_system.outs[2] is u['F301_P1'].outs[0]
    assert factory_system.outs[3] is u['F302_P1'].outs[0]


def test_specifications_attached(factory_system):
    u = _units_by_id(factory_system)
    for ID in ('F301', 'F302', 'M301', 'M302', 'H301', 'H302', 'K330', 'V330'):
        assert len(u[ID].specifications) == 1, f'{ID} missing its specification'
    # The corn-coupled feed-flow spec is deliberately NOT built by the factory.
    assert len(u['V406'].specifications) == 0


def test_defaults_match_inline_system(factory_system):
    u = _units_by_id(factory_system)
    S301, F301, F302, M301, M302 = (u[i] for i in
                                    ('S301', 'F301', 'F302', 'M301', 'M302'))
    H301, H302, V406, K330 = (u[i] for i in ('H301', 'H302', 'V406', 'K330'))
    assert (S301.split == 0.8).all()
    for F in (F301, F302):
        assert F.V == 0.1
        assert tuple(F.P) == (101325, 73581, 50892, 32777, 20000)
    assert M301.water_to_sugar_mol_ratio == 100.
    assert M302.water_to_sugar_mol_ratio == 100.
    assert H301.T == 32 + 273.15
    assert H302.T == 32 + 273.15
    assert K330.P == 3e7
    assert K330.eta == 0.6
    assert V406.tau == 3 * 24
    assert V406.tau_max == 3 * 24
    assert V406.tau_update_policy == ('max', '[s_EtOH]')
    assert V406.n_decimal_places_for_tau_update_policy == 0
    assert V406.n_simulation_steps == 1000
    assert V406.sugar_IDs == ('Glucose',)
    assert V406.spike_feed_index == 2
    assert V406.perform_hydrolysis is False
    assert V406.map_species_to_chemicals == EXPECTED_MAP
    assert V406.track_vars == EXPECTED_TRACK_VARS
    assert V406.stage_1_max_x == 5.0
    assert V406.stage_1_max_time == 25.0
    # Default kinetic model is the shipped te_r singleton.
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r
    assert V406.nsk_kinetic_model is te_r


def test_fed_batch_strategy_specification(factory_system):
    """The factory builds a FedBatchStrategySpecification wired to its own
    units, with the inline isobutanol system's values, and attaches it to the
    fermentor under both the descriptive name and the short alias."""
    from nskinetics.units import (FedBatchStrategySpecification,
                                  ConcentrationActuator)
    u = _units_by_id(factory_system)
    V406 = u['V406']
    spec = V406.fed_batch_strategy_specification
    assert isinstance(spec, FedBatchStrategySpecification)
    assert V406.fbs_spec is spec
    assert spec.fermentation_reactor is V406
    assert spec.splitter is u['S301']
    for actuator, unit_ID, attr, lb, ub in (
            (spec.feed_concentrator, 'F301', 'V', 0.0, 0.8),
            (spec.feed_diluter, 'M301', 'water_to_sugar_mol_ratio', 0.0, 100_000),
            (spec.spike_concentrator, 'F302', 'V', 0.0, 0.8),
            (spec.spike_diluter, 'M302', 'water_to_sugar_mol_ratio', 0.0, 100_000),
            ):
        assert isinstance(actuator, ConcentrationActuator)
        assert actuator.unit is u[unit_ID]
        assert actuator.attr == attr
        assert actuator.lb == lb
        assert actuator.ub == ub
    assert spec.feed_units_sequential == [
        u[i] for i in ('F301', 'F301_P0', 'F301_P1', 'M301', 'H301')]
    assert spec.spike_units_sequential == [
        u[i] for i in ('F302', 'F302_P0', 'F302_P1', 'M302', 'H302')]
    assert spec.species_IDs == ['Glucose']
    assert spec.solvent_ID == 'Water'
    # Spec concentrations: the factory-kwarg defaults; tau_max: reused from
    # the factory's own tau_max default (3 * 24 h).
    assert spec.target_conc == 220.0
    assert spec.threshold_conc == 210.0
    assert spec.spike_conc == 600.0
    assert spec.tau_max == 72.0
    cv = spec.control_variables
    assert cv.spike_conc_var == 'conc_glu_feed_spike'
    assert cv.target_conc_var == 'target_conc_glu_spike'
    assert cv.threshold_conc_var == 'threshold_conc_glu_spike'
    assert cv.volume_col is None
    assert cv.feed_volume_added_col is None
    assert spec.baseline_specifications == {'target_conc': 221.25,
                                            'threshold_conc': 217.125,
                                            'spike_conc': 600.0,
                                            'tau_max': 120.0}


def test_set_thermo_builds_standalone():
    """set_thermo=True activates the factory's own shipped chemical set and
    the system builds without any biorefineries-provided thermo. Runs last in
    this module on purpose: it replaces the global thermo."""
    import biosteam as bst
    import thermosteam as tmo
    from nskinetics.processes import create_sugar_prep_and_fermentation_system
    bst.main_flowsheet.set_flowsheet('test_processes_set_thermo')
    system = create_sugar_prep_and_fermentation_system(
        ID='set_thermo_sys', set_thermo=True)
    chems = tmo.settings.chemicals
    for ID in ('Water', 'Glucose', 'Yeast', 'Ethanol', 'Isobutanol',
               'AceticAcid', 'Sucrose', 'CO2', 'O2', 'N2', 'NH3'):
        assert ID in chems.IDs
    assert {u.ID for u in system.units} == EXPECTED_UNIT_IDS
