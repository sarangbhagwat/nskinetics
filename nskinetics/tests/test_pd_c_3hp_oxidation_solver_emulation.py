# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of ``nskinetics.models.pd_c_3hp_oxidation.solver_emulation``: the
shipped Pd/C 3-HP oxidation models against the model nodes the authors
deposited for Fig. 9a, Fig. 9b and SI Fig. S12a (H2O2) of Kim et al.,
Appl. Catal. B: Environ. Energy 386 (2026) 126403.

The deposited curves are solver steps, reproduced step for step by the
method of MATLAB ``ode23s`` at default tolerances (emulated; MATLAB was not
run), and lie up to 2.0e-3 from the converged solution, so the tests
that solve the models accurately (``test_pd_c_3hp_oxidation.py``) can hold
them only to that error. Here the authors' solver run is repeated with the
rates of the shipped models, and the result is compared with the deposited
nodes one by one:

* the number of steps must be the deposited one;
* every node time must agree within ``TIME_LIMIT``: two units of the last
  digit of the Fig. 9 file, which prints five decimals (2e-5 min);
* every value must agree within ``VALUE_LIMIT`` = 1e-5, one unit of the
  last of the five decimals the files print. (Rounding alone accounts for
  5e-6; the measured largest differences are 4.4e-6 to 6.8e-6.)

The limit on the node times has a small margin over rounding: the measured
differences are 6e-6 to 8e-6 min, and random relative errors of 1e-15 in
the rates leave them within the limit, whereas errors of 1e-14 can exceed
it (measured outside the tests). A failure on another platform that shows
only in the node times therefore points to rounding, not to the model.

The limit on the values is 85 to 200 times smaller than the solver error
that the accurate-solver tests have to allow for (8.5e-4 and 2.0e-3). The
constants are the sets of ``O2_NODE_PARAMETER_SETS`` and
``H2O2_NODE_PARAMETER_SETS``, which were inferred from the nodes; with
Table 3 as printed the run comes within 4.0e-4 (O2) and 4.4e-4 (H2O2),
which is what the rounding of Table 3 costs.

The negative controls each change one thing that was inferred and show
that the limits then fail. They test the shipped constants together with
the shipped structure: none of them refits the constants of an alternative
(``NOTES.md`` says which statements rest on such refits, which are not in
the repository).

Each run takes some twenty steps and a few hundredths of a second, so this
module is not marked slow. It imports pytest at top level and is therefore
not registered in ``nskinetics/tests/__init__.py``.
"""
import os

import numpy as np
import pytest
import tellurium as te

import nskinetics as nsk
import nskinetics.models.pd_c_3hp_oxidation as pkg
from nskinetics.models.pd_c_3hp_oxidation import (
    H2O2_NODE_PARAMETER_SETS, NORMALIZED, O2_NODE_PARAMETER_SETS,
    emulate_authors_run, emulate_fig9a, emulate_fig9b, load_h2o2_model,
    load_o2_model, node_deviations, ode23s, simulate_at_times)

_PKG_DIR = os.path.dirname(os.path.abspath(pkg.__file__))

VALUE_LIMIT = 1e-5
TIME_LIMIT = 2e-5               # min

# Paper Table 3 as printed. k_1 to k_6 in 1/(g_cat min), k_7 in 1/min.
TABLE_3_H2O2 = {'k_1': 5.31, 'k_2': 3.11e1, 'k_3': 7.32e-4, 'k_4': 1.15e0,
                'k_5': 3.03e-3, 'k_6': 4.10e-9, 'k_7': 6.48e-3}
TABLE_3_O2 = {'k_1_ref': 7.14, 'k_2_ref': 1.39e2, 'k_3_ref': 3.16e-4,
              'k_4_ref': 5.86e-2, 'k_5_ref': 1.13e-2, 'k_6_ref': 5.52e-3}
M_CAT = 0.2                     # g
X_O2_REF = 0.005085             # 15*0.000339

# Measured for each deposited curve: number of solver steps, and the rms
# and largest difference between the emulated and the deposited values.
MEASURED = {
    'fig9b': (19, 2.38e-6, 4.97e-6),
    'fig9a': (23, 3.02e-6, 6.77e-6),
    'figS12a_h2o2': (24, 2.56e-6, 4.43e-6),
    }


def _reference(curve):
    """(deposited nodes, end of the run in min, compared columns) of a
    deposited curve."""
    if curve == 'fig9b':
        tables, nodes = pkg.read_fig9b_o2(), 'model'
    else:
        tables = pkg.read_fig9a_h2o2()
        nodes = 'model' if curve == 'fig9a' else 'h2o2_trace'
    columns = ['eq_H2O2'] if nodes == 'h2o2_trace' else list(NORMALIZED)
    return (tables[nodes], 1.1*tables['experiment']['time'].iloc[-1],
            columns)


def _load(curve):
    """The shipped model with the constant set of a deposited curve."""
    if curve == 'fig9b':
        return load_o2_model(parameter_set='fig9b_nodes')
    return load_h2o2_model(curve + '_nodes')


def _compare(km, curve):
    nodes, t_end, columns = _reference(curve)
    return node_deviations(emulate_authors_run(km, t_end), nodes, columns)


def _reproduces(deviations, curve):
    """Whether an emulated run is the deposited curve to its print rounding."""
    return bool(deviations['steps'] == MEASURED[curve][0]
                and deviations['time'] <= TIME_LIMIT
                and deviations['max'] <= VALUE_LIMIT)


def _variant(antimony_filename, edits=(), **values):
    """A model built from an edited Antimony text, with values set on it.

    ``edits`` holds (old, new, count) replacements; ``count`` is the number
    of places the old text must occur in.
    """
    with open(os.path.join(_PKG_DIR, antimony_filename),
              encoding='utf-8') as f:
        text = f.read()
    for old, new, count in edits:
        assert text.count(old) == count, old
        text = text.replace(old, new)
    km = nsk.KineticModel(te.loada(text), units={'time': 'min', 'conc': 'M'})
    for name, value in values.items():
        km.set_value(name, value)
    return km


_OPA_IN_THE_SUM = ('/(HP3 + MA + AAc + OXA + D)',
                   '/(HP3 + OPA + MA + AAc + OXA + D)', 5)


# --- 1. the solver ------------------------------------------------------------

def test_ode23s_step_rules_on_a_linear_problem():
    # y' = -y from 1: J f = y, so the first step is
    # 1/(1.25*sqrt(0.5)/rtol**(1/3)).
    t, y, info = ode23s(lambda y: -y, 5., [1.])
    assert info['first_step'] == pytest.approx(
        1e-3**(1./3.)/(1.25*np.sqrt(0.5)), rel=1e-6)
    assert t[0] == 0. and t[1] == pytest.approx(info['first_step'])
    assert t[-1] == 5.
    assert (np.diff(t) > 0.).all() and np.diff(t).max() <= 0.5 + 1e-12
    assert info['failed'] == [] and max(info['error']) <= 1.
    assert len(info['error']) == len(t) - 1 == len(y) - 1
    # third-order steps at a relative tolerance of 1e-3
    assert np.abs(y[:, 0] - np.exp(-t)).max() < 2e-3
    # a step never grows beyond max_step
    t, _, _ = ode23s(lambda y: -y, 5., [1.], max_step=0.2)
    assert np.diff(t).max() <= 0.2 + 1e-12


def test_ode23s_lets_a_step_grow_at_most_fivefold():
    # After a step without a rejection the next one is the step divided by
    # 1.25*(error/rtol)**(1/3), at most five times it: the limit applies
    # when the error is below (0.2/1.25)**3 = 0.0041 of the tolerance. In
    # the three deposited runs that is the case in the last two steps
    # only, where the next step, if there is one, is cut short at t_end, so
    # the limit is tested here on y' = -100 y with y(0) = 1, integrated to
    # t = 100. Once y is far below atol/rtol the error estimate collapses,
    # and the steps grow from 0.04 to max_step = 10.
    t, _, info = ode23s(lambda y: -100.*y, 100., [1.])
    steps, errors = np.diff(t), np.array(info['error'])
    assert info['failed'] == [] and len(steps) == 56
    assert errors[42:46] == pytest.approx(
        [7.29e-3, 3.78e-3, 9.57e-4, 5.20e-5], rel=2e-3)
    # error above the limit: the step grows by 1/(1.25*error**(1/3)) = 4.13
    assert steps[43]/steps[42] == pytest.approx(
        1./(1.25*errors[42]**(1./3.)), rel=1e-12)
    assert steps[43]/steps[42] == pytest.approx(4.1265, abs=1e-3)
    # error below it, twice: exactly five times the step, where the
    # unlimited rule would give 5.14 and 8.12 times
    assert steps[44]/steps[43] == pytest.approx(5., rel=1e-12)
    assert steps[45]/steps[44] == pytest.approx(5., rel=1e-12)
    assert 1./(1.25*errors[43:45]**(1./3.)) == pytest.approx([5.14, 8.12],
                                                             abs=0.01)
    # then max_step takes over
    assert steps[45] == pytest.approx(4.3827, abs=1e-4)
    assert steps[46:55] == pytest.approx(10., rel=1e-12)
    # every step follows the rule
    growth = np.where(errors > 0.16**3, 1./(1.25*errors**(1./3.)), 5.)
    assert steps[1:55] == pytest.approx(
        np.minimum(10., growth[:54]*steps[:54]), rel=1e-12)


def test_ode23s_stretches_a_step_that_would_end_just_before_t_end():
    # A step that would end within 10 % of its own length before t_end is
    # stretched to end there, even beyond max_step; one that would end
    # further from it is followed by a short last step. The last step of
    # each deposited run is cut short at t_end, so the stretching is tested
    # here on y' = -y with max_step = 0.2: after the first step (0.1131)
    # every step is 0.2.
    def run(t_end):
        t, _, info = ode23s(lambda y: -y, t_end, [1.], max_step=0.2)
        assert info['failed'] == [] and t[-1] == t_end
        return t, np.diff(t)

    t, steps = run(5.)
    assert steps[0] == pytest.approx(0.113137, abs=1e-6)
    assert steps[1:-1] == pytest.approx(0.2, rel=1e-12)
    # ... so a run that ends 1.08 steps after t[5] stretches its last step
    _, stretched = run(t[5] + 1.08*0.2)
    assert len(stretched) == 6
    assert stretched[:5] == pytest.approx(steps[:5], rel=1e-12)
    assert stretched[5] == pytest.approx(0.216, rel=1e-9)
    # ... and one that ends 1.12 steps after it takes a full step and a
    # short one
    _, two = run(t[5] + 1.12*0.2)
    assert len(two) == 7
    assert two[5:] == pytest.approx([0.2, 0.024], rel=1e-9)


def test_ode23s_weighs_the_error_with_the_larger_of_the_two_states():
    # The error of a step is taken relative to max(|y|, |y_new|, atol/rtol).
    # The deposited runs do not test the |y_new| term (without it all three
    # take the same steps), so it is tested here on a growing solution,
    # y' = y from 1, with the first step written out for a scalar.
    rtol = 1e-3
    h = rtol**(1./3.)/(1.25*np.sqrt(0.5))
    d, e32 = 1./(2. + np.sqrt(2.)), 6. + np.sqrt(2.)
    W = 1. - h*d
    k1 = 1./W
    f1 = 1. + 0.5*h*k1
    k2 = (f1 - k1)/W + k1
    y_new = 1. + h*k2
    k3 = (y_new - e32*(k2 - f1) - 2.*(k1 - 1.))/W
    error = h/6.*abs(k1 - 2.*k2 + k3)/y_new/rtol
    assert (y_new, error) == pytest.approx((1.1198504, 0.05785923), rel=1e-7)

    t, y, info = ode23s(lambda y: y, 10., [1.])
    assert t[1] == pytest.approx(h, rel=1e-12)
    assert y[1, 0] == pytest.approx(y_new, rel=1e-12)
    assert info['error'][0] == pytest.approx(error, rel=1e-9)
    # relative to |y| alone the estimate would be 12 % larger, and the
    # second step, which follows from it, 0.2253 instead of 0.2340
    assert error*y_new == pytest.approx(0.06479368, rel=1e-7)
    assert t[2] - t[1] == pytest.approx(h/(1.25*error**(1./3.)), rel=1e-9)
    assert t[2] - t[1] == pytest.approx(0.23401, abs=1e-5)
    assert info['failed'] == [] and len(t) - 1 == 44


def test_ode23s_rejected_step_rules_on_a_quadrature():
    # None of the three deposited runs has a rejected step, so the rules for
    # one are tested here on y' = cos(10 t), written with the time as a
    # first state. J f is zero at the start, so the first step tried is
    # max_step, 0.5. It fails by a factor of 600 and is retried with a tenth
    # of it (the lower limit of the factor), and after that twice with half
    # the step. The step taken, 0.0125, is kept for the next one.
    t, y, info = ode23s(lambda y: np.array([1., np.cos(10.*y[0])]), 5.,
                        [0., 0.])
    failed, steps = info['failed'], np.diff(t)
    assert info['first_step'] == 0.5
    assert [time for time, _, _ in failed[:3]] == [0., 0., 0.]
    assert [step for _, step, _ in failed[:3]] == pytest.approx(
        [0.5, 0.05, 0.025], rel=1e-12)
    assert [error for _, _, error in failed[:3]] == pytest.approx(
        [600.38, 10.363, 2.6008], rel=1e-4)
    assert steps[:2] == pytest.approx([0.0125, 0.0125], rel=1e-9)
    # The next rejection, by a factor of 1.84 at the tenth step: the step
    # is retried with 0.8*(rtol/error)**(1/3) = 0.653 of it and kept.
    time, step, error = failed[3]
    assert time == t[9] and error == pytest.approx(1.8428, rel=1e-4)
    assert steps[9] == pytest.approx(0.8*step/error**(1./3.), rel=1e-9)
    assert steps[10] == pytest.approx(steps[9], rel=1e-9)
    # the solution is sin(10 t)/10
    assert np.abs(y[:, 1] - 0.1*np.sin(10.*t)).max() < 1e-3


def test_ode23s_tolerances_other_than_the_defaults_miss_the_nodes():
    nodes, t_end, columns = _reference('fig9b')
    km = _load('fig9b')
    for options in (dict(rtol=5e-4), dict(rtol=2e-3), dict(atol=5e-7),
                    dict(atol=2e-6)):
        deviations = node_deviations(
            emulate_authors_run(km, t_end, **options), nodes, columns)
        assert deviations['steps'] != 19 or deviations['time'] > 10.


# --- 2. the deposited nodes, to their print rounding --------------------------

@pytest.mark.parametrize('curve', list(MEASURED))
def test_emulation_returns_the_deposited_nodes(curve):
    nodes, t_end, columns = _reference(curve)
    emulated = {'fig9a': emulate_fig9a,
                'figS12a_h2o2': lambda: emulate_fig9a('figS12a_h2o2_nodes'),
                'fig9b': emulate_fig9b}[curve]()
    assert emulated['time'].iloc[-1] == t_end
    deviations = node_deviations(emulated, nodes, columns)
    steps, rms, largest = MEASURED[curve]
    assert deviations['steps'] == steps
    if curve == 'figS12a_h2o2':
        # the file holds every step of the trace but the last one
        assert deviations['deposited_steps'] == steps - 1
        assert emulated['time'].iloc[-2] == pytest.approx(1537.08028,
                                                          abs=TIME_LIMIT)
    else:
        assert deviations['deposited_steps'] == steps
    assert deviations['time'] <= TIME_LIMIT, deviations
    assert deviations['max'] <= VALUE_LIMIT, deviations
    assert deviations['rms'] == pytest.approx(rms, rel=0.05)
    assert deviations['max'] == pytest.approx(largest, rel=0.05)
    # no step was rejected, and none was accepted narrowly, so that rounding
    # differences between machines cannot change the sequence of steps
    info = emulated.attrs['ode23s']
    assert max(info['error']) < 0.6
    assert info['failed'] == []
    # Two step-size rules are not tested by these runs. An error small
    # enough for a step to grow fivefold occurs in the last two steps only,
    # and the last step is cut short at t_end, not stretched to it.
    errors, steps = np.array(info['error']), np.diff(emulated['time'])
    assert (errors[:-2] > 0.16**3).all()
    growth = (1./(1.25*errors[-2]**(1./3.)) if errors[-2] > 0.16**3 else 5.)
    assert steps[-1] < 0.95*min(0.1*t_end, growth*steps[-2])
    assert _reproduces(deviations, curve)


@pytest.mark.parametrize('curve', list(MEASURED))
def test_emulation_is_the_loaded_model(curve):
    # The same run through the loader and through the convenience function,
    # and an accurate solution of the same model at the node times: the
    # emulated steps carry the solver error of the deposited curve.
    nodes, t_end, columns = _reference(curve)
    km = _load(curve)
    emulated = emulate_authors_run(km, t_end)
    assert _reproduces(node_deviations(emulated, nodes, columns), curve)
    accurate = simulate_at_times(km, emulated['time'], columns)
    solver_error = np.abs(emulated[columns].to_numpy()
                          - accurate[columns].to_numpy()).max()
    expected = {'fig9b': 2.00e-3, 'fig9a': 8.54e-4,
                'figS12a_h2o2': 0.0115}[curve]
    assert solver_error == pytest.approx(expected, rel=0.03)
    assert solver_error > 15.*VALUE_LIMIT


@pytest.mark.parametrize('curve', list(MEASURED))
def test_script_states_do_not_change_a_step(curve):
    # H2O, H2 and CO2 of the posted script enter no rate; integrated
    # along or not, the steps are the same.
    _, t_end, _ = _reference(curve)
    km = _load(curve)
    with_states = emulate_authors_run(km, t_end)
    without = emulate_authors_run(km, t_end, script_states=False)
    names = ['time'] + list(NORMALIZED)
    assert np.abs(with_states[names].to_numpy()
                  - without[names].to_numpy()).max() < 1e-12
    extra = set(with_states.columns) - set(without.columns)
    assert extra == {'fig9b': {'state_H2O', 'state_H2', 'state_CO2'}}.get(
        curve, set())
    # their rows of the stoichiometric matrix of the posted O2 script
    # (columns r1 to r6): 1 0 0 0 0 0, 0 0 0 1 0 0 and 0 0 1 1 0 0
    assert pkg.solver_emulation._SCRIPT_ROWS == {
        'H2O': {'R1': 1.}, 'H2': {'R4': 1.}, 'CO2': {'R3': 1., 'R4': 1.}}


def test_emulated_run_starts_from_the_authors_initial_state():
    emulated = emulate_fig9b()
    first = emulated.iloc[0]
    assert first['state_HP3'] == pytest.approx(0.000339, rel=1e-14)
    assert first['state_H2O'] == pytest.approx(1. - 0.000339, rel=1e-14)
    assert first['state_yO2'] == 1.
    assert emulated.columns[-10:].tolist() == [
        'state_' + name for name in ('HP3', 'OPA', 'MA', 'AAc', 'OXA', 'D',
                                     'H2O', 'H2', 'yO2', 'CO2')]
    first = emulate_fig9a().iloc[0]
    assert first['state_H2O2'] == pytest.approx(15*0.000339, rel=1e-14)
    assert first['eq_H2O2'] == pytest.approx(15., rel=1e-14)


def test_opa_inset_of_fig_s12b_is_the_same_solver_run():
    # The 13 vertices of the OPA curve in the SI Fig. S12b inset (digitized)
    # are the first 13 nodes of the emulated Fig. 9b run (t = 0 and its
    # first 12 steps).
    inset = pkg.read_figS12b_opa_inset_digitized()
    steps = emulate_fig9b().iloc[:len(inset)]
    assert len(inset) == 13
    assert (np.abs(steps['time'].to_numpy() - inset['time'].to_numpy())
            <= 3.*inset['time_unc'].to_numpy() + 1e-9).all()
    difference = np.abs(steps['x_OPA'].to_numpy() - inset['x_OPA'].to_numpy())
    assert difference.max() == pytest.approx(1.9e-7, rel=0.1)
    assert difference.max() < 0.015*steps['x_OPA'].max()


def test_emulation_leaves_the_model_reset_and_keeps_set_values():
    km = load_o2_model()
    km.set_value('m_cat', 0.4)
    before = simulate_at_times(km, [0., 30., 600.])
    emulate_authors_run(km, 1848.)
    assert km.get_value('m_cat') == 0.4
    assert km.get_value('[HP3]') == 0.06 and km.get_value('yO2') == 1.
    assert simulate_at_times(km, [0., 30., 600.]).equals(before)


def test_emulation_refuses_a_model_without_a_deposited_curve():
    other = nsk.KineticModel(te.loada('S1 -> S2; k*S1; k = 1; S1 = 1'),
                             units={'time': 'min', 'conc': 'M'})
    with pytest.raises(ValueError, match='no deposited curve'):
        emulate_authors_run(other, 100.)


# --- 3. Table 3 as printed, and the constant sets of the nodes ----------------

def test_table_3_as_printed_comes_close_but_not_to_the_print_rounding():
    # What rounding Table 3 to three digits costs, measured.
    o2 = node_deviations(emulate_fig9b('table3'),
                         pkg.read_fig9b_o2()['model'])
    assert o2['steps'] == 19
    assert (o2['time'], o2['rms'], o2['max']) == pytest.approx(
        (0.0405, 1.37e-4, 4.04e-4), rel=0.02)
    h2o2 = node_deviations(emulate_fig9a('table3'),
                           pkg.read_fig9a_h2o2()['model'])
    assert h2o2['steps'] == 23
    assert (h2o2['time'], h2o2['rms'], h2o2['max']) == pytest.approx(
        (0.471, 1.96e-4, 4.40e-4), rel=0.02)
    assert not _reproduces(o2, 'fig9b') and not _reproduces(h2o2, 'fig9a')


def test_node_sets_relative_to_table_3():
    o2 = O2_NODE_PARAMETER_SETS['fig9b_nodes']
    assert [o2[name]/TABLE_3_O2[name] for name in TABLE_3_O2] == pytest.approx(
        [1.00007, 0.99770, 1.00028, 1.00013, 0.99743, 1.00023], abs=6e-6)
    # ... every one of them prints as Table 3 prints it
    assert {name: float('%.3g' % value) for name, value in o2.items()} \
        == TABLE_3_O2
    fig9a = H2O2_NODE_PARAMETER_SETS['fig9a_nodes']
    assert [fig9a[name]/TABLE_3_H2O2[name] for name in TABLE_3_H2O2] \
        == pytest.approx([1.00242, 0.99399, 1.00522, 0.99756, 0.99827, 1.,
                          1.00129], abs=6e-6)
    # ... of which only k_4 (and k_6, kept as printed) prints as Table 3
    assert [name for name, value in fig9a.items()
            if float('%.3g' % value) == TABLE_3_H2O2[name]] == ['k_4', 'k_6']
    trace = H2O2_NODE_PARAMETER_SETS['figS12a_h2o2_nodes']
    assert [trace[name]/TABLE_3_H2O2[name] for name in TABLE_3_H2O2] \
        == pytest.approx([1.04354, 1.05228, 0.91569, 0.96122, 0.78103, 1.,
                          1.07045], abs=6e-6)


def test_default_models_are_table_3_and_the_loaders_check_the_set_name():
    for km in (load_o2_model(), load_o2_model(parameter_set='table3')):
        assert {name: km.get_value(name) for name in TABLE_3_O2} == TABLE_3_O2
    for km in (load_h2o2_model(), load_h2o2_model('table3')):
        assert {name: km.get_value(name)
                for name in TABLE_3_H2O2} == TABLE_3_H2O2
    with pytest.raises(ValueError, match="'table3', 'fig9b_nodes'"):
        load_o2_model(parameter_set='fig9a_nodes')
    for name in ('fig9b_nodes', 'Table3', None, ['fig9a_nodes']):
        with pytest.raises(ValueError, match="'table3', 'fig9a_nodes', "
                                             "'figS12a_h2o2_nodes'"):
            load_h2o2_model(name)
    # the temperature and the constant set can be given together
    km = load_o2_model(T=348.15, parameter_set='fig9b_nodes')
    assert km.get_value('T') == 348.15
    assert km.get_value('k_1_ref') == 7.1405267


@pytest.mark.parametrize('curve', ['fig9b', 'fig9a', 'figS12a_h2o2'])
def test_node_set_survives_compile_events_and_reset(curve):
    km = _load(curve)
    selected = dict(O2_NODE_PARAMETER_SETS['fig9b_nodes'] if curve == 'fig9b'
                    else H2O2_NODE_PARAMETER_SETS[curve + '_nodes'])
    before = _compare(km, curve)
    # never fires; an event may only assign to a variable of the model
    km.add_event(nsk.Event(when='time < 0', do={'D': 'D'}))
    km.compile_events()
    km.reset()
    for name, value in selected.items():
        assert km.get_value(name) == value, name
    assert _compare(km, curve) == before


# --- 4. negative controls: O2 model (Fig. 9b) ---------------------------------

def test_o2_factor_of_0p0051_misses_the_fig9b_nodes():
    # Only the products k_i*m_cat*x_O2_ref are determined, so the test is
    # made against Table 3: the products of the node set, over Table 3 as
    # printed times 0.2 g, are the oxidant factor that Table 3 implies.
    node_set = O2_NODE_PARAMETER_SETS['fig9b_nodes']
    implied = {name: node_set[name]*X_O2_REF/TABLE_3_O2[name]
               for name in ('k_1_ref', 'k_2_ref', 'k_4_ref')}
    assert list(implied.values()) == pytest.approx(
        [0.0050854, 0.0050733, 0.0050857], abs=5e-8)
    # A printed value v stands for v - 0.5 to v + 0.5 units of its third
    # digit, which bounds the factor from both sides:
    half_unit = {'k_1_ref': 0.005, 'k_2_ref': 0.5, 'k_4_ref': 0.00005}
    lowest = max(implied[name]*TABLE_3_O2[name]
                 / (TABLE_3_O2[name] + half_unit[name]) for name in implied)
    highest = min(implied[name]*TABLE_3_O2[name]
                  / (TABLE_3_O2[name] - half_unit[name]) for name in implied)
    assert (lowest, highest) == pytest.approx((0.0050818, 0.0050889),
                                              abs=5e-8)
    assert lowest < X_O2_REF < highest < 0.0051

    km = _load('fig9b')
    assert _reproduces(_compare(km, 'fig9b'), 'fig9b')
    # 0.0051 with the constants of the node set ...
    km.set_value('x_O2_ref', 0.0051)
    deviations = _compare(km, 'fig9b')
    assert not _reproduces(deviations, 'fig9b')
    assert (deviations['time'], deviations['max']) == pytest.approx(
        (1.80, 9.4e-4), rel=0.03)
    # ... and with the constants of reactions 1, 2 and 4 at the edge of
    # what Table 3 can stand for that favours it most.
    for name in implied:
        km.set_value(name, TABLE_3_O2[name] - half_unit[name])
    deviations = _compare(km, 'fig9b')
    assert not _reproduces(deviations, 'fig9b')
    assert (deviations['time'], deviations['max']) == pytest.approx(
        (1.36, 7.0e-4), rel=0.03)
    assert deviations['max'] > 50.*VALUE_LIMIT
    # The other three constants were those of the node set; with Table 3
    # as printed for them it is 5.8e-4.
    for name in ('k_3_ref', 'k_5_ref', 'k_6_ref'):
        km.set_value(name, TABLE_3_O2[name])
    deviations = _compare(km, 'fig9b')
    assert (deviations['time'], deviations['max']) == pytest.approx(
        (1.54, 5.8e-4), rel=0.03)
    assert deviations['max'] > 50.*VALUE_LIMIT


def test_o2_fraction_held_constant_misses_the_fig9b_nodes():
    # The posted script integrates the O2 state and uses yO2 = 1.0 in its
    # rates. With the constants of the node set that misses the nodes.
    km = _variant(pkg.o2_antimony_filename,
                  [('x_O2_ref*yO2*', 'x_O2_ref*', 3)],
                  **O2_NODE_PARAMETER_SETS['fig9b_nodes'])
    deviations = _compare(km, 'fig9b')
    assert deviations['steps'] == 19
    assert not _reproduces(deviations, 'fig9b')
    assert (deviations['time'], deviations['max']) == pytest.approx(
        (0.149, 7.7e-5), rel=0.03)
    # unedited, the variant is the loaded model
    km = _variant(pkg.o2_antimony_filename,
                  **O2_NODE_PARAMETER_SETS['fig9b_nodes'])
    assert _compare(km, 'fig9b') == _compare(_load('fig9b'), 'fig9b')


@pytest.mark.parametrize('curve, largest', [
    ('fig9b', 0.038), ('fig9a', 0.094)])
def test_opa_in_the_normalization_misses_the_nodes(curve, largest):
    filename = {'fig9b': pkg.o2_antimony_filename,
                'fig9a': pkg.h2o2_antimony_filename}[curve]
    constants = {'fig9b': O2_NODE_PARAMETER_SETS['fig9b_nodes'],
                 'fig9a': H2O2_NODE_PARAMETER_SETS['fig9a_nodes']}[curve]
    deviations = _compare(_variant(filename, [_OPA_IN_THE_SUM], **constants),
                          curve)
    # the steps are the same; only the plotted quantity differs
    assert deviations['steps'] == MEASURED[curve][0]
    assert deviations['time'] <= TIME_LIMIT
    assert deviations['max'] == pytest.approx(largest, rel=0.03)
    assert not _reproduces(deviations, curve)


# --- 5. negative controls: H2O2 model (Fig. 9a, SI Fig. S12a) -----------------

def test_k7_multiplied_by_the_catalyst_mass_misses_the_fig9a_nodes():
    # Table 3 gives k_7 in 1/min (footnote b); k_1 to k_6 are per gram.
    constants = H2O2_NODE_PARAMETER_SETS['fig9a_nodes']
    km = _variant(pkg.h2o2_antimony_filename,
                  [('k_7*H2O2*env', 'k_7*m_cat*H2O2*env', 1)], **constants)
    deviations = _compare(km, 'fig9a')
    assert deviations['steps'] == 16
    assert deviations['max'] == pytest.approx(0.395, rel=0.03)
    assert not _reproduces(deviations, 'fig9a')
    km = load_h2o2_model('fig9a_nodes')
    km.set_value('k_7', M_CAT*constants['k_7'])
    assert _compare(km, 'fig9a') == deviations


def test_h2o2_trace_is_another_run_than_the_fig9a_curve():
    trace = pkg.read_fig9a_h2o2()['h2o2_trace']
    for parameter_set in ('fig9a_nodes', 'table3'):
        deviations = node_deviations(emulate_fig9a(parameter_set), trace,
                                     ['eq_H2O2'])
        # one step fewer than the trace has, and node times far off
        assert deviations['steps'] == 23
        assert deviations['time'] > 50.
    curve = pkg.read_fig9a_h2o2()['model']
    deviations = node_deviations(emulate_fig9a('figS12a_h2o2_nodes'), curve)
    assert deviations['steps'] == 24
    assert deviations['max'] > 0.05


# --- 6. how tightly the nodes pin each constant -------------------------------

# Relative change of each constant that the limits detect: the smallest
# power of ten for which the run then misses a limit at least twofold.
# None: not detected by a change of 10 %.
DETECTED_CHANGE = {
    # k_1_ref to k_6_ref
    'fig9b': (1e-6, 1e-6, 1e-3, 1e-3, 1e-6, 1e-4),
    # k_1 to k_7; k_6 is not identifiable
    'fig9a': (1e-6, 1e-6, 1e-3, 1e-4, 1e-5, None, 1e-6),
    'figS12a_h2o2': (1e-6, 1e-6, 1e-3, 1e-4, 1e-4, None, 1e-6),
    }


def _constants(curve):
    if curve == 'fig9b':
        return O2_NODE_PARAMETER_SETS['fig9b_nodes']
    return H2O2_NODE_PARAMETER_SETS[curve + '_nodes']


@pytest.mark.parametrize('curve', list(DETECTED_CHANGE))
def test_how_tightly_the_nodes_pin_each_constant(curve):
    km = _load(curve)
    constants = _constants(curve)
    assert len(constants) == len(DETECTED_CHANGE[curve])
    for (name, value), change in zip(constants.items(),
                                     DETECTED_CHANGE[curve]):
        if change is None:
            # not detected: a change of 10 % leaves the nodes as they are
            km.set_value(name, 1.1*value)
            assert _reproduces(_compare(km, curve), curve), name
        else:
            km.set_value(name, (1. + change)*value)
            assert not _reproduces(_compare(km, curve), curve), name
        km.set_value(name, value)
    assert _reproduces(_compare(km, curve), curve)
