# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""KineticModel.parameters / .species introspection tables (F1 item 10)."""
import libsbml
import numpy as np
import tellurium as te

import nskinetics as nsk

TOY = """
model introspect()
  compartment env = 1; species S in env, P in env, E in env; species $B in env;
  S = 50; P = 0; B = 2; E = 2*k2;
  k = 0.5; k2 = 0.25; Km = 3; threshold = 5; n_spk = 0; z = 1; q = 4*k;
  k2 has per_h;
  curr_env := env;
  z' = 0.1;
  J: S => P; k*S*env;
  J2: S + E => P + E; k2*S*E/(Km + S)*B;
  J3: P => ; k*P;
  E1: at (S < threshold): S = S + 10, n_spk = n_spk + 1;
end
"""


def _model():
    # Give J3 a local parameter `k` that shadows the global one.
    doc = libsbml.readSBMLFromString(te.antimonyToSBML(TOY))
    local = doc.getModel().getReaction('J3').getKineticLaw()
    local = local.createLocalParameter()
    local.setId('k')
    local.setValue(7.0)
    km = nsk.KineticModel.from_sbml(libsbml.writeSBMLToString(doc),
                                    units={'time': 'h', 'conc': 'g/L'})
    km.add_event(nsk.Event(when='time >= 3', do={'threshold': '1'},
                           name='lower_threshold'))
    km.compile_events()
    return km


def test_parameter_kinds_and_flags():
    km = _model()
    p = km.parameters
    assert list(p.columns) == ['value', 'initial_value', 'units', 'constant',
                               'kind', 'assigned_by_events',
                               'used_in_reactions']
    assert p.index.name == 'id'
    assert p.loc['curr_env', 'kind'] == 'assignment_rule'
    assert p.loc['z', 'kind'] == 'rate_rule'
    assert p.loc['q', 'kind'] == 'initial_assignment'
    assert p.loc['Km', 'kind'] == 'parameter'
    # event targets: a native Antimony event and a compiled nsk.Event
    assert p.loc['n_spk', 'assigned_by_events']
    assert p.loc['threshold', 'assigned_by_events']
    assert not p.loc['n_spk', 'constant']
    assert not p.loc['k', 'assigned_by_events'] and p.loc['k', 'constant']
    assert set(p.index[p['assigned_by_events']]) == {'n_spk', 'threshold'}
    assert p.loc['k2', 'units'] == 'per_h' and p.loc['k', 'units'] == ''
    assert np.isnan(p.loc['curr_env', 'initial_value'])


def test_used_in_reactions_walks_kinetic_law_asts():
    p = _model().parameters
    # J3's local k shadows the global one, so only J uses the global k
    assert p.loc['k', 'used_in_reactions'] == ('J',)
    assert p.loc['k2', 'used_in_reactions'] == ('J2',)
    assert p.loc['Km', 'used_in_reactions'] == ('J2',)
    assert p.loc['threshold', 'used_in_reactions'] == ()


def test_values_track_user_edits_initial_values_do_not():
    km = _model()
    km._te.Km = 4.5
    p = km.parameters
    assert p.loc['Km', 'value'] == 4.5
    assert p.loc['Km', 'initial_value'] == 3.0
    assert p.loc['q', 'value'] == 2.0


def test_species_table():
    km = _model()
    s = km.species
    assert list(s.columns) == ['value', 'initial_value', 'units', 'constant',
                               'kind', 'assigned_by_events',
                               'used_in_reactions', 'compartment', 'boundary']
    assert set(s.index) == {'S', 'P', 'E', 'B'}
    assert s.loc['S', 'compartment'] == 'env'
    assert s.loc['S', 'kind'] == 'species'
    assert s.loc['E', 'kind'] == 'initial_assignment'
    assert s.loc['E', 'initial_value'] == 0.5
    assert s.loc['S', 'assigned_by_events'] and not s.loc['P', 'assigned_by_events']
    assert s.loc['S', 'used_in_reactions'] == ('J', 'J2')
    assert s.loc['P', 'used_in_reactions'] == ('J', 'J2', 'J3')
    assert s.loc['B', 'used_in_reactions'] == ('J2',)    # rate law only
    assert s.loc['B', 'boundary'] and not s.loc['S', 'boundary']
    km.simulate(0, 1, 11, reset=True)
    s = km.species
    assert s.loc['S', 'initial_value'] == 50.0
    assert s.loc['S', 'value'] < 50.0


def test_query_by_kind_replaces_name_prefix_filters():
    p = _model().parameters
    fittable = p.query("kind == 'parameter' and not assigned_by_events")
    assert set(fittable.index) == {'k', 'k2', 'Km'}


def test_introspection_does_not_touch_model_state():
    km = _model()
    km.simulate(0, 2, 5, reset=True)
    before = {sel: km._te[sel] for sel in km.state_selections()}
    km.parameters, km.species
    assert {sel: km._te[sel] for sel in km.state_selections()} == before
