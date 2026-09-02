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


def test_kinetic_and_operation_sets_are_disjoint():
    pc = _pc()
    assert not set(pc.KINETIC_PARAMETERS) & pc.OPERATION_PARAMETERS
    assert len(pc.KINETIC_PARAMETERS) == 65
    assert len(pc.OPERATION_PARAMETERS) == 14


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
