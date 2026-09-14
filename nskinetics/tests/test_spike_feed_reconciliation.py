# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Reactor-boundary conservation of the fed-batch spike feed.

``NSKBatchReactor`` initialises the kinetic model from the initial charge
only and lets the model's ``FeedSpike`` events add spike volume; the glucose
the flowsheet actually delivers through the spike inlet is otherwise ignored.
These tests pin that the reactor now reconciles the two inside its own run
(the ``spike_feed_reconciler`` hook the process factory wires to the
``FedBatchStrategySpecification``), exposes the residual, and raises only
when reconciliation cannot close the balance. Design and evidence:
``docs/reports/fed-batch-spike-feed-reconciliation.md``.

Marked ``slow`` as a whole and deliberately not registered in
``nskinetics/tests/__init__.py``: it builds the quickstart system on the
shipped ``te_r`` model (imported at call time), so every import lives inside
a fixture or test body. Never load a second Antimony model in this process:
doing so perturbs the shipped model.
"""

import numpy as np
import pytest

pytestmark = pytest.mark.slow

# Diagnostic references at the baseline strategy (CLAUDE.md, 2026-09-13).
REF_TAU = 61.8378
REF_SPLIT = 0.81898
REF_N_SPIKES = 3
REF_ETHANOL_OUT = 17441.6     # kg/hr
REF_MASS_CLOSURE = -0.0034    # (out - in)/in; NH3 zeroing + rounding


@pytest.fixture(scope='module')
def system():
    """Build the quickstart system once, WITHOUT simulating it: the first
    test below exercises the reactor before any strategy is imposed, the
    second is the first to call ``system.simulate()``. Tests run in file
    order and every later test leaves the system at the simulated
    baseline."""
    import os
    os.environ.setdefault('MPLBACKEND', 'Agg')
    import biosteam as bst
    from nskinetics.processes import create_sugar_prep_and_fermentation_system
    bst.main_flowsheet.set_flowsheet('test_spike_feed_reconciliation')
    return create_sugar_prep_and_fermentation_system(
        ID='spike_feed_reconciliation_sys', set_thermo=True)


def _units(system):
    u = system.flowsheet.unit
    return u.V406, u.S301, u.V406.fbs_spec


def _split(S301):
    return float(np.atleast_1d(S301.split)[0])


def _mass_closure(V406):
    m_in = sum(i.F_mass for i in V406.ins)
    m_out = sum(o.F_mass for o in V406.outs)
    return (m_out - m_in) / m_in


def _independent_residual(V406, spec):
    """The report's own residual, computed without the reactor's attributes."""
    d = V406.nsk_results_specific_tau_dict
    env, added = d['curr_env'], d['curr_tot_vol_glu_feed_added']
    V0 = V406.ins[0].ivol['Water'] + V406.ins[1].ivol['Water']
    implied = added / (env - added) * spec.spike_conc * V0
    delivered = V406.ins[2].imass['Glucose']
    return (implied - delivered) / delivered


def test_factory_wires_the_specification_as_reconciler(system):
    V406, S301, spec = _units(system)
    assert V406.spike_feed_reconciler is spec


def test_raw_actuators_are_reconciled_before_any_strategy_is_imposed(system):
    """Simulate the feed trains and the reactor on the freshly built system,
    before load_specifications ever ran: the evaporators, dilution mixers and
    splitter sit at their raw construction values, so the initial feed is far
    from the nominal target concentration (~79 g/L vs 220). This is what the
    isobutanol biorefinery's initialising ``system.simulate()`` does (it
    builds the factory with ``mockup=True``, so no system specification runs
    first). The reactor must still close the balance — deriving the split
    from the measured balance, not from the nominal concentrations, whose
    fixed point never converges here."""
    V406, S301, spec = _units(system)
    spec._simulate_upstream_units()
    assert spec.get_feed_conc() < 0.5 * spec.target_conc  # raw actuators
    split_raw = _split(S301)
    V406.simulate()
    assert abs(V406.spike_feed_residual) <= V406.spike_feed_reconciliation_tol
    assert abs(_independent_residual(V406, spec)) <= 1e-3
    assert 1 <= V406.n_spike_feed_reconciliation_passes <=         V406.spike_feed_max_reconciliation_passes
    assert _split(S301) != pytest.approx(split_raw, abs=1e-2)
    # Only the controlled species is reconciled here: the raw spike stream is
    # far more dilute than the spike concentration the model ran with, so its
    # excess water is not in the effluent and total mass does NOT close on
    # this run (it does once load_desired_concs has imposed the strategy;
    # see test_baseline_residual_pin).
    assert _mass_closure(V406) < -0.1


def test_baseline_residual_pin(system):
    """At the baseline strategy the split and the run already agree, so the
    reactor exits on pass 0 with the residual below tolerance and the
    diagnostic references unchanged. First call of ``system.simulate()`` on
    the module's system."""
    system.simulate()
    V406, S301, spec = _units(system)
    d = V406.nsk_results_specific_tau_dict
    assert abs(V406.spike_feed_residual) <= 1e-3
    assert V406.spike_feed_residual == pytest.approx(
        _independent_residual(V406, spec), abs=1e-6)
    assert V406.n_spike_feed_reconciliation_passes == 0
    assert V406.spike_feed_delivered == pytest.approx(
        V406.ins[2].imass['Glucose'])
    assert d['curr_n_glu_spikes'] == REF_N_SPIKES
    assert d['time'] == pytest.approx(REF_TAU, abs=1e-3)
    assert _split(S301) == pytest.approx(REF_SPLIT, abs=1e-4)
    assert _mass_closure(V406) == pytest.approx(REF_MASS_CLOSURE, abs=5e-4)


def test_forced_split_is_reconciled_inside_the_reactor_run(system):
    """Force the splitter away from its consistent value and re-run only the
    feed trains and the reactor (no load_specifications): the reactor must
    move the split back, close the residual, and reproduce the baseline
    effluent — instead of emitting an effluent that lost 11.6 % of its mass."""
    V406, S301, spec = _units(system)
    split0 = np.array(S301.split, copy=True)
    try:
        S301.split = 0.5
        spec._simulate_upstream_units()
        V406.simulate()
        assert V406.n_spike_feed_reconciliation_passes >= 1
        assert abs(V406.spike_feed_residual) <= V406.spike_feed_reconciliation_tol
        assert abs(_independent_residual(V406, spec)) <= 1e-3
        assert _split(S301) == pytest.approx(REF_SPLIT, abs=1e-3)
        assert V406.nsk_results_specific_tau_dict['curr_n_glu_spikes'] == REF_N_SPIKES
        assert V406.outs[1].imass['Ethanol'] == pytest.approx(
            REF_ETHANOL_OUT, rel=5e-3)
        assert _mass_closure(V406) == pytest.approx(REF_MASS_CLOSURE, abs=5e-4)
        # The reconciliation passes ran with the spike count frozen at the
        # count the first pass reached (the model cap is left at that count).
        assert V406.nsk_kinetic_model.get_value('max_n_glu_spikes') == REF_N_SPIKES
    finally:
        S301.split = split0
        spec._simulate_upstream_units()
        V406.simulate()


def test_unreconcilable_spike_feed_raises_after_the_passes(system):
    """With a reconciler that cannot move the split, the reactor exhausts its
    passes and raises MassBalanceError naming implied and delivered."""
    from nskinetics.exceptions import MassBalanceError

    class _Stuck:
        def spike_feed_balance(self, reactor, minimal_feed, spike_feed):
            return spec.spike_feed_balance(reactor, minimal_feed, spike_feed)

        def reconcile_spike_feed(self, reactor):
            pass

    V406, S301, spec = _units(system)
    split0 = np.array(S301.split, copy=True)
    V406.spike_feed_reconciler = _Stuck()
    try:
        S301.split = 0.5
        spec._simulate_upstream_units()
        with pytest.raises(MassBalanceError) as excinfo:
            V406.simulate()
        msg = str(excinfo.value)
        assert f'{V406.spike_feed_implied:.4g}' in msg
        assert f'{V406.spike_feed_delivered:.4g}' in msg
        assert V406.n_spike_feed_reconciliation_passes == \
            V406.spike_feed_max_reconciliation_passes
    finally:
        V406.spike_feed_reconciler = spec
        S301.split = split0
        spec._simulate_upstream_units()
        V406.simulate()


def test_reconciler_off_leaves_the_residual_unreported(system):
    """Without a reconciler the reactor runs as before: no reconciliation
    passes and the residual attributes are None."""
    V406, S301, spec = _units(system)
    V406.spike_feed_reconciler = None
    try:
        V406.simulate()
        assert V406.spike_feed_residual is None
        assert V406.n_spike_feed_reconciliation_passes == 0
    finally:
        V406.spike_feed_reconciler = spec
        V406.simulate()
