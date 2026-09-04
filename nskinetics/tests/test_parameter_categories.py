# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests for the shipped model's kinetic-parameter classification and the
parameter-change description helper.

Marked ``slow`` as a whole and deliberately not registered in
``nskinetics/tests/__init__.py``: the shipped subpackage builds ``te_r`` on
import, and ``test_processes_contract.test_import_is_lightweight`` asserts
that nothing named for the shipped model is in ``sys.modules`` during the
``-m "not slow"`` run. Even the dict-only tests below import the
classification from that subpackage, so every such import lives inside a
test body and the module runs only under ``-m ""`` (same precedent as
``test_flux_map.py``).
"""

import pytest

pytestmark = pytest.mark.slow


def _pc():
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        parameter_categories)
    return parameter_categories


# --- structure ------------------------------------------------------------

def test_roles_and_modules_are_closed_sets():
    pc = _pc()
    assert set(pc.REACTION_MODULES.values()) <= set(pc.MODULES)
    assert set(pc.MODULE_LABELS) == set(pc.MODULES)
    for pid, info in pc.KINETIC_PARAMETERS.items():
        assert info.role in pc.ROLES, pid
        assert set(info.modules) <= set(pc.MODULES), pid
        for rxn in info.reactions:
            assert rxn in pc.REACTION_MODULES, (pid, rxn)


def test_modules_derive_from_reactions():
    pc = _pc()
    for pid, info in pc.KINETIC_PARAMETERS.items():
        if info.reactions:
            derived = []
            for rxn in info.reactions:
                m = pc.REACTION_MODULES[rxn]
                if m not in derived:
                    derived.append(m)
            assert info.modules == tuple(derived), pid
        else:
            assert info.role == 'initial_state', pid
            assert info.modules == ('physiological_state',), pid


def test_two_module_parameters():
    pc = _pc()
    assert pc.KINETIC_PARAMETERS['K_5e'].modules == ('respiration', 'growth')
    assert pc.KINETIC_PARAMETERS['K_5i'].modules == ('respiration', 'growth')
    assert pc.KINETIC_PARAMETERS['K_5i'].effector == 'glucose'
    assert pc.KINETIC_PARAMETERS['K_1i'].effector == 'acetaldehyde'


def test_every_role_assignment_matches_the_spec_table():
    # Part 1 of
    # docs/superpowers/specs/2026-09-02-parameter-categories-design.md
    pc = _pc()
    expected = {
        'capacity': frozenset(
            'k_1h k_1l k_1e k_2 k_3 k_4 k_5 k_5e k_6 k_7 k_8 k_9 k_9e k_9c '
            'k_10 k_11 k_13 k_14 k_15 k_16'.split()),
        'affinity': frozenset(
            'K_1h K_1l K_1e K_2 K_3 K_4 K_5 K_5e K_6 K_7 K_9 K_9e K_13 K_14 '
            'K_15 K_16'.split()),
        'substrate_regulation': frozenset('K_1i K_2i K_5i K_9i'.split()),
        'product_inhibition': frozenset(
            'k_1ie k_1ia k_1ii k_4ie k_4ia k_4ii k_6ia k_6ii k_7ie k_7ia '
            'k_7ii k_16ie k_16ia'.split()),
        'product_self_inhibition': frozenset('K_6e k_6r K_16i k_16r'.split()),
        'lethality': frozenset('k_10ie k_10ia k_10ii'.split()),
        'lethality_threshold': frozenset('P_10e P_10a P_10i'.split()),
        'initial_state': frozenset('X_a X_AcDH'.split()),
    }
    assert set(expected) == set(pc.ROLES)
    actual = {role: frozenset(p for p, info in pc.KINETIC_PARAMETERS.items()
                              if info.role == role)
              for role in pc.ROLES}
    assert actual == expected


def test_kinetic_and_operation_sets_are_disjoint():
    pc = _pc()
    assert not set(pc.KINETIC_PARAMETERS) & pc.OPERATION_PARAMETERS
    assert len(pc.KINETIC_PARAMETERS) == 65
    assert len(pc.OPERATION_PARAMETERS) == 16


def test_flux_map_inhibition_map_is_consistent():
    pc = _pc()
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import FLUX_MAP_SPEC
    for pid, (rxn, inhibitor) in FLUX_MAP_SPEC.inhibition_map.items():
        info = pc.KINETIC_PARAMETERS[pid]
        assert rxn in info.reactions, pid
        assert info.effector == inhibitor, pid


# --- model coverage -------------------------------------------------------

def test_every_model_parameter_is_classified():
    pc = _pc()
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r
    r = te_r._te
    assigned = set(r.getAssignmentRuleIds())
    kinetic = set(pc.KINETIC_PARAMETERS)
    operation = set(pc.OPERATION_PARAMETERS)
    for pid in r.getGlobalParameterIds():
        hits = (pid in kinetic) + (pid in operation) + (pid in assigned)
        assert hits == 1, f'{pid} classified {hits} times'
    ids = set(r.getGlobalParameterIds())
    assert kinetic <= ids
    assert operation <= ids
    assert set(pc.REACTION_MODULES) <= set(r.getReactionIds())


# --- snapshot / diff ------------------------------------------------------

def test_diff_reports_changed_kinetic_parameters_in_source_order():
    pc = _pc()
    base = {'k_7ii': 0.04, 'k_1ii': 0.04, 'k_4ii': 0.04, 'k_7': 1.203}
    cur = {'k_7ii': 0.15, 'k_1ii': 0.15, 'k_4ii': 0.15, 'k_7': 1.203}
    d = pc.diff_parameters(base, cur)
    assert list(d) == ['k_1ii', 'k_4ii', 'k_7ii']
    assert d['k_7ii'] == (0.04, 0.15)


def test_diff_partial_current_treats_missing_as_unchanged():
    pc = _pc()
    base = {'k_13': 0.0, 'k_14': 0.0, 'k_7': 1.203}
    assert pc.diff_parameters(base, {'k_13': 5.81}) == {'k_13': (0.0, 5.81)}


def test_diff_skips_operation_parameters():
    pc = _pc()
    base = {'k_7': 1.203, 'max_n_glu_spikes': 5}
    cur = {'k_7': 1.203, 'max_n_glu_spikes': 13}
    assert pc.diff_parameters(base, cur) == {}


def test_diff_unknown_key_raises():
    pc = _pc()
    with pytest.raises(ValueError, match='not_a_param'):
        pc.diff_parameters({'k_7': 1.0}, {'k_7': 1.0, 'not_a_param': 1.0})
    with pytest.raises(ValueError, match='qO2'):
        pc.diff_parameters({'qO2': 1.0}, {})


def test_diff_key_only_in_current_raises():
    pc = _pc()
    with pytest.raises(ValueError, match='k_13'):
        pc.diff_parameters({'k_7': 1.0}, {'k_13': 5.81})


def test_diff_tolerates_float_noise():
    pc = _pc()
    assert pc.diff_parameters({'k_7': 1.0}, {'k_7': 1.0 + 1e-14}) == {}


def test_snapshot_reads_every_kinetic_parameter_from_te_r():
    pc = _pc()
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r
    snap = pc.snapshot_parameters(te_r)
    assert set(snap) == set(pc.KINETIC_PARAMETERS)
    assert all(isinstance(v, float) for v in snap.values())
    assert snap['k_7'] == pytest.approx(1.203)
    # a raw roadrunner works too, and gives the same values
    assert pc.snapshot_parameters(te_r._te) == snap
    # the subpackage re-exports the module's names, not copies of them
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        describe_parameter_change, KINETIC_PARAMETERS)
    assert describe_parameter_change is pc.describe_parameter_change
    assert KINETIC_PARAMETERS is pc.KINETIC_PARAMETERS


def test_snapshot_rejects_a_foreign_model():
    pc = _pc()
    import tellurium as te
    r = te.loadAntimonyModel(
        'model toy() S = 1; k = 1; J: S => ; k*S; end')
    with pytest.raises(ValueError, match='k_1h'):
        pc.snapshot_parameters(r)


# --- describe -------------------------------------------------------------

def test_describe_isobutanol_inhibition_raised():
    pc = _pc()
    base = {'k_7ii': 0.04, 'k_1ii': 0.04, 'k_4ii': 0.04}
    cur = {'k_7ii': 0.15, 'k_1ii': 0.15, 'k_4ii': 0.15}
    assert pc.describe_parameter_change(base, cur) == (
        'stronger isobutanol inhibition of glycolysis/fermentation, '
        'overflow/acetate and growth')


def test_describe_verbose_lists_members_in_source_order():
    pc = _pc()
    base = {'k_7ii': 0.04, 'k_1ii': 0.04, 'k_4ii': 0.04}
    cur = {'k_7ii': 0.15, 'k_1ii': 0.15, 'k_4ii': 0.15}
    assert pc.describe_parameter_change(base, cur, verbose=True) == (
        'stronger isobutanol inhibition of glycolysis/fermentation, '
        'overflow/acetate and growth '
        '(k_1ii 0.04 -> 0.15, k_4ii 0.04 -> 0.15, k_7ii 0.04 -> 0.15)')
    # the exact phrase quoted in tutorial 05
    assert pc.describe_parameter_change(
        {'k_1ie': 0.02, 'k_7ie': 0.04},
        {'k_1ie': 0.04, 'k_7ie': 0.08}, verbose=True) == (
        'stronger ethanol inhibition of glycolysis/fermentation and growth '
        '(k_1ie 0.02 -> 0.04, k_7ie 0.04 -> 0.08)')


def test_describe_ehrlich_on_from_partial_override():
    pc = _pc()
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        SCENARIO_A_EHRLICH, SCENARIO_B_EHRLICH)
    assert pc.describe_parameter_change(SCENARIO_A_EHRLICH, SCENARIO_B_EHRLICH) == (
        'Ehrlich branch on; isobutanol self-inhibition of Ehrlich branch on')
    assert pc.describe_parameter_change(SCENARIO_B_EHRLICH, SCENARIO_A_EHRLICH) == (
        'Ehrlich branch off; isobutanol self-inhibition of Ehrlich branch off')


def test_describe_threshold_and_lethality_words():
    pc = _pc()
    assert pc.describe_parameter_change({'P_10e': 100.}, {'P_10e': 120.}) == \
        'later ethanol death onset'
    assert pc.describe_parameter_change({'P_10a': 5.}, {'P_10a': 3.}) == \
        'earlier acetate death onset'
    assert pc.describe_parameter_change({'k_10ie': 0.04}, {'k_10ie': 0.08}) == \
        'steeper ethanol lethality'


def test_describe_affinity_and_regulation_words():
    pc = _pc()
    assert pc.describe_parameter_change({'K_7': 0.0101}, {'K_7': 0.02}) == \
        'lower substrate affinity in growth'
    assert pc.describe_parameter_change({'K_2i': 0.101}, {'K_2i': 0.05}) == \
        'weaker glucose regulation of respiration'
    assert pc.describe_parameter_change({'K_5i': 440.}, {'K_5i': 500.}) == \
        'stronger glucose regulation of respiration and growth'


def test_describe_self_inhibition_and_initial_state():
    pc = _pc()
    assert pc.describe_parameter_change({'K_6e': 0.057}, {'K_6e': 0.1}) == \
        'stronger ethanol self-inhibition of glycolysis/fermentation'
    assert pc.describe_parameter_change({'X_a': 0.1, 'X_AcDH': 0.0075},
                                        {'X_a': 0.2, 'X_AcDH': 0.005}) == \
        'higher initial X_a; lower initial X_AcDH'


def test_describe_disjoint_capacity_changes_stay_separate():
    pc = _pc()
    base = {'k_7': 1.203, 'k_1h': 0.584}
    cur = {'k_7': 1.5, 'k_1h': 0.5}
    assert pc.describe_parameter_change(base, cur) == \
        'faster growth; slower glycolysis/fermentation'


def test_describe_overlapping_opposite_changes_collapse_to_retuned():
    pc = _pc()
    base = {'k_1h': 0.584, 'k_1l': 1.43}
    cur = {'k_1h': 0.7, 'k_1l': 1.0}
    assert pc.describe_parameter_change(base, cur) == \
        'retuned glycolysis/fermentation capacity'
    base = {'k_7ie': 0.04, 'k_7ia': 0.12}
    cur = {'k_7ie': 0.08, 'k_7ia': 0.06}
    # different effectors never collapse
    assert pc.describe_parameter_change(base, cur) == \
        'stronger ethanol inhibition of growth; weaker acetate inhibition of growth'


def test_describe_mixed_clause_outside_capacity():
    pc = _pc()
    assert pc.describe_parameter_change({'K_1h': 1., 'K_1l': 2.},
                                        {'K_1h': 2., 'K_1l': 1.}) == (
        'retuned substrate affinity in glycolysis/fermentation')


def test_describe_non_switchable_zero_uses_up_down_words():
    pc = _pc()
    assert pc.describe_parameter_change({'K_7': 0.0}, {'K_7': 0.01}) == \
        'lower substrate affinity in growth'


def test_every_parameter_describes_in_every_direction():
    pc = _pc()
    for p in pc.KINETIC_PARAMETERS:
        for old, new in ((1., 2.), (2., 1.), (0., 1.), (1., 0.)):
            text = pc.describe_parameter_change({p: old}, {p: new})
            assert isinstance(text, str), (p, old, new)
            assert text, (p, old, new)
            assert text != 'no kinetic parameter changes', (p, old, new)


def test_describe_empty_diff():
    pc = _pc()
    assert pc.describe_parameter_change({'k_7': 1.0}, {'k_7': 1.0}) == \
        'no kinetic parameter changes'


def test_describe_against_te_r_snapshot():
    pc = _pc()
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        te_r, SCENARIO_B_EHRLICH)
    snap = pc.snapshot_parameters(te_r)
    assert pc.describe_parameter_change(snap, SCENARIO_B_EHRLICH) == (
        'Ehrlich branch on; isobutanol self-inhibition of Ehrlich branch on')


# --- categorize -----------------------------------------------------------

def test_categorize_groups_by_role_and_module():
    pc = _pc()
    out = pc.categorize(['k_7', 'K_5e', 'k_7ii'])
    assert out == {
        ('capacity', 'growth'): ['k_7'],
        ('affinity', 'respiration'): ['K_5e'],
        ('affinity', 'growth'): ['K_5e'],
        ('product_inhibition', 'growth'): ['k_7ii'],
    }


def test_categorize_rejects_operation_and_unknown():
    pc = _pc()
    with pytest.raises(ValueError, match='operation'):
        pc.categorize(['D'])
    with pytest.raises(ValueError, match='nope'):
        pc.categorize(['nope'])
