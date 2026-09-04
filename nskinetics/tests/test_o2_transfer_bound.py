# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests for the kLa-bounded oxygen transfer in the shipped *S. cerevisiae*
ethanol/isobutanol model.

The model caps the respiratory + growth-associated O2 demand of its aerobic
fluxes at ``OTR_max = kLa*C_O2_sat`` through the assignment-rule factor
``f_O2``, which replaces the binary ``is_aerobic`` gate on r2, r5, r8 and the
aerobic share of r7 (see
``docs/superpowers/specs/2026-09-04-o2-transfer-bound-design.md``).

Marked ``slow`` as a whole and deliberately not registered in
``nskinetics/tests/__init__.py``: the shipped subpackage builds ``te_r`` on
import, and ``test_processes_contract.test_import_is_lightweight`` asserts
that nothing named for the shipped model is in ``sys.modules`` during the
``-m "not slow"`` run. Every import of it therefore lives inside a test body
(same precedent as ``test_parameter_categories.py``). Never load a second
Antimony model in this process: doing so perturbs the shipped model.
"""

import contextlib
import math

import numpy as np
import pytest

pytestmark = pytest.mark.slow

# Selections recorded by every run. Bracketed ids are concentrations (g/L);
# reaction ids are rates (g/h); the rest are parameters/assignment rules.
COLS = ['time', '[s_glu]', '[x]', '[s_EtOH]', 'is_aerobic',
        'f_O2', 'OTR_max', 'OUR_demand', 'OUR_gated', 'OUR_committed',
        'qO2_TCA_growth_only', 'r2', 'r5', 'r7', 'r8', 'v2_0']


def _model():
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        te_r, reset_nsk_kinetic_model)
    return te_r, reset_nsk_kinetic_model


@contextlib.contextmanager
def _batch(*, kLa=None, aerobic=True, stage_1_max_x=math.inf):
    """Reset the shipped model to a 100 g/L glucose batch (the model's default
    10 -> 100 g/L glucose refills stay active) and yield the raw RoadRunner
    object; restores ``kLa`` and the stage-1 cutoffs afterwards so later
    tests (and the process factory, which sets them itself) see the SBML
    defaults."""
    te_r, reset = _model()
    r = te_r._te
    old_kLa = r.kLa
    try:
        reset(te_r)
        r.s_glu = 100.
        r.x = 1.
        r.stage_1_max_time = math.inf
        r.stage_1_max_x = stage_1_max_x
        r.is_aerobic = 1 if aerobic else 0
        if kLa is not None:
            r.kLa = kLa
        yield r
    finally:
        r.kLa = old_kLa
        r.stage_1_max_time = math.inf
        r.stage_1_max_x = math.inf
        reset(te_r)


def _run(r, t_end=60., n=601):
    res = r.simulate(0, t_end, n, COLS)
    return {c: np.asarray(res[:, i], dtype=float) for i, c in enumerate(COLS)}


def _realized_volumetric_uptake(d):
    # qO2_TCA_growth_only is per gram biomass; times [x] it is mmol O2/L/h.
    return d['qO2_TCA_growth_only'] * d['[x]']


def test_defaults_and_derived_cap():
    te_r, _ = _model()
    r = te_r._te
    assert r.kLa == pytest.approx(200.0)
    assert r.C_O2_sat == pytest.approx(0.20)
    assert r.OTR_max == pytest.approx(40.0)
    assert set(r.getAssignmentRuleIds()) >= {
        'OTR_max', 'v2_0', 'v5_0', 'v7_0', 'v8_0',
        'OUR_gated', 'OUR_committed', 'OUR_demand', 'f_O2'}


def test_default_cap_is_non_binding_in_a_staged_batch():
    """A 100 g/L batch that stops aeration at x = 5 g/L, i.e. the same
    ``stage_1_max_x`` cutoff style the shipped process uses -- not the shipped
    V406 fed-batch itself. What is pinned here is the local claim: under that
    staging the peak aerated demand stays well below the default cap, so
    ``f_O2`` is identically 1 while aerated and the cap changes nothing. The
    corresponding claim for the shipped fed-batch is carried separately, by the
    in-repo diagnostic's ``min_f_O2_aerobic`` probe."""
    # Peak demand while aerated is ~13.5 mmol/L/h here, comfortably below the
    # 40 mmol/L/h default cap.
    with _batch(stage_1_max_x=5.0) as r:
        d = _run(r)
    aerobic = d['is_aerobic'] == 1.0
    assert aerobic.any() and (~aerobic).any()          # the stage switch fired
    assert np.all(d['f_O2'][aerobic] == 1.0)
    assert np.all(d['f_O2'][~aerobic] == 0.0)
    assert d['OUR_demand'][aerobic].max() < d['OTR_max'][0]
    with _batch(stage_1_max_x=5.0, kLa=1e9) as r:
        d_unbounded = _run(r)
    # The gate itself is bit-identical: raising kLa cannot change a trajectory
    # whose f_O2 is already 1 throughout stage 1 and 0 (via is_aerobic) after.
    assert np.array_equal(d['f_O2'], d_unbounded['f_O2'])
    # The states agree to ~1e-6 rather than exactly. After the aeration switch
    # the *demand* keeps climbing and crosses OTR_max at t ~ 7 h; f_O2 stays 0
    # there (is_aerobic = 0), but the piecewise condition is still a root the
    # integrator stops on, so the default-kLa run takes a slightly different
    # step sequence than one whose cap is never approached. Same solver-path
    # noise CLAUDE.md documents for the in-repo diagnostic.
    # (Measured max rel diff 7.9e-7; 2e-6 leaves headroom for solver noise.)
    for col in ('[s_glu]', '[x]', '[s_EtOH]'):
        np.testing.assert_allclose(d[col], d_unbounded[col], rtol=2e-6, atol=1e-9)


def test_default_cap_binds_in_fully_aerobic_batch():
    # Aerobic throughout, biomass climbs to ~14 g/L and the unbounded demand
    # peaks near 69 mmol/L/h -- above the 40 mmol/L/h default cap.
    with _batch() as r:
        d = _run(r)
    f = d['f_O2']
    assert f.min() < 1.0
    assert np.all((f >= 0.0) & (f <= 1.0))
    realized = _realized_volumetric_uptake(d)
    otr = d['OTR_max']
    assert np.all(realized <= otr * (1 + 1e-6))
    binding = f < 1.0
    # Exact form: where the cap binds, realized uptake equals the cap.
    np.testing.assert_allclose(realized[binding], otr[binding], rtol=1e-6)
    # The gate scales the gated rate laws, term for term.
    np.testing.assert_allclose(d['r2'], f * d['v2_0'], rtol=1e-9, atol=1e-12)
    # Against an uncapped run the bound costs respiration and, through it,
    # biomass: r8 (growth on acetate) is gated, and r7 is not throttled at the
    # default anaerobic_growth_mult = 1. It does *not* push carbon into
    # ethanol -- r3 (pyruvate -> acetaldehyde) is a Hill term saturated far
    # above K_3 = 5e-7, so throttling r2 leaves pyruvate standing rather than
    # raising the ethanol branch; less active biomass means slightly less
    # ethanol too (110.215 vs 110.363 g/L at t = 60 h).
    with _batch(kLa=1e9) as r:
        d_unbounded = _run(r)
    assert np.all(d_unbounded['f_O2'] == 1.0)
    assert (_realized_volumetric_uptake(d).max()
            < _realized_volumetric_uptake(d_unbounded).max())
    assert d['[x]'][-1] < d_unbounded['[x]'][-1]


def test_lower_kla_lowers_the_cap():
    with _batch(kLa=100.0) as r:
        d = _run(r)
    assert d['OTR_max'][0] == pytest.approx(20.0)
    realized = _realized_volumetric_uptake(d)
    assert np.all(realized <= 20.0 * (1 + 1e-6))
    assert d['f_O2'].min() < 1.0


def test_committed_growth_alone_can_saturate_transfer():
    # The piecewise's second branch: OUR_committed >= OTR_max, i.e. the
    # ungated growth share of r7 alone needs more O2 than the vessel can
    # transfer, so nothing is left for the gated fluxes and f_O2 -> 0.
    # Measured on this fully aerobic 100 g/L batch: OUR_committed peaks at
    # ~6.10 mmol/L/h, so kLa = 1 (OTR_max = 0.2 mmol/L/h) puts the branch in
    # reach immediately -- it holds at 251 of the 601 sampled steps, from
    # t = 0 h. (kLa = 2/5/10 also reach it, later; 1 is the clearest.)
    with _batch(kLa=1.0) as r:
        d = _run(r)
    assert d['OTR_max'][0] == pytest.approx(0.2)
    aerobic = d['is_aerobic'] == 1.0
    assert aerobic.all()
    saturated = aerobic & (d['OUR_committed'] >= d['OTR_max'])
    assert saturated.any()
    assert np.all(d['f_O2'][saturated] == 0.0)
    # With f_O2 = 0 the gated fluxes are off, so the only O2 actually taken up
    # is the committed growth share -- never more.
    off = d['f_O2'] == 0.0
    realized = _realized_volumetric_uptake(d)
    assert np.all(realized[off] <= d['OUR_committed'][off] + 1e-9)


def test_anaerobic_stage_forces_zero():
    with _batch(aerobic=False) as r:
        d = _run(r, t_end=20., n=201)
    assert np.all(d['f_O2'] == 0.0)
    for rxn in ('r2', 'r5', 'r8'):
        assert np.all(d[rxn] == 0.0), rxn
    assert np.all(d['r7'] > 0.0)                        # growth continues


def test_zero_gated_demand_has_no_nan():
    # At t = 0 of a batch there is no pyruvate or acetate yet, so the gated
    # demand is exactly 0 while growth's committed share is not: the
    # piecewise must take its first branch without dividing by OUR_gated.
    # (x = 0 itself is not a usable state of this model: the pre-existing
    # X_a' rate rule divides by x*env.)
    with _batch() as r:
        d = _run(r, t_end=5., n=51)
    assert d['OUR_gated'][0] == 0.0
    assert d['OUR_committed'][0] > 0.0
    assert d['f_O2'][0] == 1.0
    for col, arr in d.items():
        assert not np.isnan(arr).any(), col
