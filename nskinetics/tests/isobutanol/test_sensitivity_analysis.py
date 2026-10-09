# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Offline tests of ``sensitivity_analysis.py`` (the Sobol' / Shapley analysis
of the metabolic_split_12d campaign space).

No ``load()``, no simulation: ``sa``, ``ko`` and ``eb`` are loaded BY FILE
PATH (never through the package), so these tests are sim-safe on any
numba-cache state. Converted from the script
``analyses/test_sensitivity_analysis_offline.py``; each test is one of its 13
numbered checks, with every assertion and tolerance kept.

The campaign-space tests need the B parameter workbook (``openpyxl``); the
surrogate tests need ``scikit-learn``. The two surrogate-fitting tests that
take ~20 s and ~100 s (9 and 10) are also marked ``slow``.
"""
import itertools
import json
import os
import types

import numpy as np
import pytest

from ._support import PKG_DIR, load_by_path, skip_unless_modules

pytestmark = pytest.mark.isobutanol

needs_workbooks = skip_unless_modules('pandas', 'openpyxl')
needs_sklearn = skip_unless_modules('sklearn')

#: Feeding baseline of the campaign (scenario A).
BASE_FEED = dict(target_conc=221.25, threshold_conc=217.125, spike_conc=600.0)
TOL = 0.02


@pytest.fixture(scope='module')
def sa():
    return load_by_path('sa', os.path.join(PKG_DIR, 'sensitivity_analysis.py'))


@pytest.fixture(scope='module')
def ko():
    return load_by_path('ko', os.path.join(PKG_DIR, 'kinetic_optimization.py'))


@pytest.fixture(scope='module')
def eb():
    return load_by_path('eb', os.path.join(PKG_DIR, 'enzyme_burden.py'))


@pytest.fixture(scope='module')
def campaign(sa, ko, eb):
    """The real split_12d space, built offline as the driver builds it.

    Scenario-A-like live baselines: the B workbook's rows with the Ehrlich
    rates zeroed and k_17 at the corrected antimony value (what A carries
    live).
    """
    kb = {**ko.workbook_kinetic_baselines('B'),
          'k_13': 0.0, 'k_14': 0.0, 'k_15': 0.0, 'k_16': 0.0, 'k_17': 0.1077}
    preset, engine_kwargs = sa.campaign_engine_kwargs(
        ko, 'ethanol_isobutanol', 'metabolic_split_12d')
    space, _ = ko.build_search_space(kb, **engine_kwargs)
    bm = eb.BurdenModel.from_reference(kb)
    scalar_feasible = ko.feasibility_predicate(
        burden_on=True, volume_on=True, burden_model=bm,
        parameter_groups=engine_kwargs['parameter_groups'], kinetic_baselines=kb,
        baseline_model_kwargs=BASE_FEED, baseline_max_n_spikes=16,
        volume_cap=20.0, group_references=engine_kwargs['group_references'])
    design = sa.design_record(
        search_space=space, parameter_groups=engine_kwargs['parameter_groups'],
        group_references=engine_kwargs['group_references'], kinetic_baselines=kb,
        burden_model=bm, eb=eb, baseline_model_kwargs=BASE_FEED,
        baseline_max_n_spikes=16, volume_cap=20.0,
        target_conc_max=ko.TARGET_CONC_MAX, meta={'seed': 7})
    design_rt = json.loads(json.dumps(design))
    vf = sa.VectorizedFeasibility.from_design(design_rt)     # from the JSON round trip
    return types.SimpleNamespace(kb=kb, preset=preset, engine_kwargs=engine_kwargs,
                                 space=space, bm=bm, scalar_feasible=scalar_feasible,
                                 design=design, design_rt=design_rt, vf=vf)


def _box(U):
    return np.ones(len(U), dtype=bool)


@needs_workbooks
def test_campaign_space_and_design_record(campaign):
    """1. campaign_engine_kwargs reproduces the 12-variable split_12d space;
    design_record is JSON round-trippable."""
    space, preset, engine_kwargs = campaign.space, campaign.preset, campaign.engine_kwargs
    assert list(space) == ['k_3', 'k_6', 'k_13', 'k_17', 'glycolysis',
                           'ehrlich_downstream', 'inhib_ethanol',
                           'inhib_isobutanol', 'inhib_acetate', 'threshold_conc',
                           'target_delta', 'max_n_spikes'], list(space)
    assert preset['scenario'] == 'A'
    # the preset's A-anchored override must win over the B-workbook-derived band
    for name, band in (preset['param_bounds_override'] or {}).items():
        assert engine_kwargs['param_bounds_override'][name] == band, name
    design_rt = json.loads(json.dumps(campaign.design))
    assert list(design_rt['search_space']) == list(space)
    assert design_rt['meta']['seed'] == 7


@needs_workbooks
def test_unit_to_values_matches_ko(sa, ko, campaign):
    """2. unit_to_values == ko.unit_to_external, elementwise, on 500 random
    points (log, linear, int)."""
    space = campaign.space
    rng = np.random.default_rng(0)
    U = rng.random((500, len(space)))
    vals = sa.unit_to_values(U, space)
    for r in range(len(U)):
        ext = ko.unit_to_external(U[r], space)
        for name in space:
            assert np.isclose(vals[name][r], ext[name], rtol=1e-12, atol=0.0), (r, name)
    assert vals['max_n_spikes'].dtype.kind == 'i'


@needs_workbooks
def test_feasible_sobol_stream_and_vectorized_feasibility(sa, ko, campaign):
    """3. Feasible Sobol' stream: deterministic, nested, feasible; the
    vectorized feasibility equals ko.feasibility_predicate on 3000 points and
    its phi_M equals BurdenModel.evaluate."""
    space, scalar_feasible, vf = campaign.space, campaign.scalar_feasible, campaign.vf
    engine_kwargs, kb, bm = campaign.engine_kwargs, campaign.kb, campaign.bm

    def take(n, seed):
        return list(itertools.islice(sa.feasible_sobol_stream(
            space, scalar_feasible, seed, ko.unit_to_external), n))
    a, b = take(40, 11), take(40, 11)
    assert a == b                                   # deterministic
    assert take(15, 11) == a[:15]                   # nested prefix
    assert take(40, 12) != a                        # seed matters
    idx = [c for c, _ in a]
    assert idx == sorted(idx) and len(set(idx)) == len(idx)
    assert all(scalar_feasible(v) for _, v in a)
    assert all(isinstance(v['max_n_spikes'], int) for _, v in a)
    assert vf.d == 12 and vf.names == list(space)
    U = np.random.default_rng(1).random((3000, vf.d))
    vec = vf(U)
    sca = np.array([scalar_feasible(ko.unit_to_external(u, space)) for u in U])
    assert (vec == sca).all(), f'{(vec != sca).sum()} disagreements'
    assert 0.05 < vec.mean() < 0.95, vec.mean()     # a real, non-trivial constraint
    phi = vf.phi_M(U[:50])
    for r in range(50):
        ext = ko.unit_to_external(U[r], space)
        res = bm.evaluate(ko.expand_grouped_values(
            ext, engine_kwargs['parameter_groups'], kb, engine_kwargs['group_references']))
        assert np.isclose(phi[r], res.Phi_M, rtol=1e-10), r
    print(f'   feasible fraction of the campaign box: {vec.mean():.3f}')


def test_ishigami_on_a_box(sa):
    """4. Unconstrained box: analytic Ishigami indices (first / total /
    closed / Shapley) within 0.02; Shapley sums to 1."""
    def ishigami(U):
        x = -np.pi + 2*np.pi*U
        return np.sin(x[:, 0]) + 7*np.sin(x[:, 1])**2 + 0.1*x[:, 2]**4*np.sin(x[:, 0])
    S4, fb4 = sa.all_closed_indices({'y': ishigami}, _box, 3, n_base=2**15,
                                    n_replicates=3, seed=0)
    m4 = S4['y'].mean(axis=0)
    assert S4['y'].shape == (3, 8) and m4[0] == 0.0 and m4[7] == 1.0
    assert fb4.max() == 0.0                                  # a box never needs the fallback
    assert np.allclose(sa.first_order(m4, 3), [0.3139, 0.4424, 0.0], atol=TOL)
    assert np.allclose(sa.total_order(m4, 3), [0.5576, 0.4424, 0.2437], atol=TOL)
    assert abs(m4[0b101] - 0.5576) < TOL                     # closed {x1, x3} = 1 - ST_2
    sh4 = sa.shapley_effects(m4, 3)
    assert np.allclose(sh4, [0.3139 + 0.2437/2, 0.4424, 0.2437/2], atol=TOL)
    assert abs(sh4.sum() - 1.0) < 1e-12


def test_additive_function_on_a_box(sa):
    """5. Additive function on a box: Shapley == first-order == total ==
    analytic 9:4:1."""
    def additive(U):
        return 3*U[:, 0] + 2*U[:, 1] + U[:, 2]
    S5, _ = sa.all_closed_indices({'y': additive}, _box, 3, n_base=2**14,
                                  n_replicates=2, seed=1)
    m5 = S5['y'].mean(axis=0)
    expect5 = np.array([9, 4, 1])/14
    assert np.allclose(sa.first_order(m5, 3), expect5, atol=TOL)
    assert np.allclose(sa.shapley_effects(m5, 3), expect5, atol=TOL)
    assert np.allclose(sa.total_order(m5, 3), expect5, atol=TOL)


def test_constrained_triangle_domain(sa):
    """6. CONSTRAINED domain (triangle x1 + x2 < 1): uniform-feasible
    sampling, conditional partners, S_{1}(X1) = 1, S_{1}(X1+X2) = 1/4
    analytic, Shapley symmetric."""
    def tri(U):
        return U[:, 0] + U[:, 1] < 1.0
    X6 = sa.sample_feasible(20000, 2, tri, np.random.default_rng(2))
    assert tri(X6).all() and abs(X6[:, 0].mean() - 1/3) < 0.01      # uniform on the triangle
    Xp6, nfb6 = sa.conditional_partners(X6, 0b01, tri, np.random.default_rng(3))
    assert (Xp6[:, 0] == X6[:, 0]).all() and tri(Xp6).all()          # x1 frozen, partner feasible
    assert nfb6 < 20                                                 # fallback only at x1 ~ 1
    S6, _ = sa.all_closed_indices({'x1': lambda U: U[:, 0],
                                   'sum': lambda U: U[:, 0] + U[:, 1]},
                                  tri, 2, n_base=2**15, n_replicates=3, seed=4)
    assert abs(S6['x1'].mean(axis=0)[0b01] - 1.0) < 1e-9             # Y = X1: S_{1} = 1 exactly
    # Y = X1 + X2 on the triangle: E[Y|X1] = (1 + X1)/2, Var(X1) = Var(Y) = 1/18 -> S_{1} = 1/4
    msum = S6['sum'].mean(axis=0)
    assert abs(msum[0b01] - 0.25) < TOL and abs(msum[0b10] - 0.25) < TOL
    assert np.allclose(sa.shapley_effects(msum, 2), [0.5, 0.5], atol=TOL)


def test_subset_search_on_a_hand_built_table(sa):
    """7. best_subsets / smallest_subset_reaching / mask_names on a
    hand-built index table."""
    names7 = ['a', 'b', 'c']
    S7 = np.zeros(8)
    S7[0b001], S7[0b010], S7[0b100] = 0.10, 0.30, 0.05
    S7[0b011], S7[0b101], S7[0b110] = 0.85, 0.20, 0.40
    S7[0b111] = 1.0
    best7 = sa.best_subsets(S7, 3, sizes=(1, 2, 3))
    assert best7[1] == (0b010, 0.30) and best7[2] == (0b011, 0.85) and best7[3] == (0b111, 1.0)
    assert sa.smallest_subset_reaching(S7, 3, 0.8) == (0b011, 0.85)
    assert sa.smallest_subset_reaching(S7, 3, 0.9) == (0b111, 1.0)
    assert sa.mask_names(0b011, names7) == ('a', 'b') and sa.mask_names(0b100, names7) == ('c',)


def test_tail_variance_share(sa):
    """8. tail_variance_share: 1 - within-variance of the non-tail group /
    Var(y), against a hand calculation."""
    y8 = np.array([-10.0, -8.0, 1.0, 2.0, 3.0, 4.0])
    # exact: 1 - (variance carried by the y >= 0 group alone, within-group) / Var(y)
    p_hi = 4/6
    share8 = 1.0 - p_hi*np.var(y8[2:])/np.var(y8)
    assert abs(sa.tail_variance_share(y8, 0.0) - share8) < 1e-12
    assert sa.tail_variance_share(np.array([1.0, 2.0, 3.0]), 0.0) == 0.0


@pytest.mark.slow
@needs_sklearn
def test_surrogate_selection(sa):
    """9. Surrogate selection: GP chosen for a smooth response, trees fit a
    step, pure noise flagged unreliable."""
    rng9 = np.random.default_rng(5)
    U9 = rng9.random((600, 3))
    smooth = np.sin(2*np.pi*U9[:, 0]) + U9[:, 1]**2
    sur_smooth = sa.fit_surrogates(U9, smooth, seed=0)
    assert sur_smooth.name == 'gp' and sur_smooth.q2['gp'] > 0.95 and sur_smooth.reliable
    step = 5.0*((U9[:, 0] > 0.5) & (U9[:, 1] > 0.5))
    sur_step = sa.fit_surrogates(U9, step, seed=0)
    assert sur_step.q2['hgb'] > 0.9 and sur_step.reliable
    assert sur_step.name == max(sur_step.q2, key=sur_step.q2.get)   # chosen = best CV Q2
    noise = rng9.standard_normal(600)
    sur_noise = sa.fit_surrogates(U9, noise, seed=0)
    assert not sur_noise.reliable and max(sur_noise.q2.values()) < sa.RELIABLE_Q2
    assert sur_smooth.predict(U9[:7]).shape == (7,)


@pytest.mark.slow
@needs_workbooks
@needs_sklearn
def test_end_to_end_on_the_real_feasible_domain(sa, campaign):
    """10. End to end on a synthetic response over the REAL feasible
    split_12d domain: surrogate + all 4094 subsets recover the planted
    k_13 x ehrlich_downstream pair as the best 2-subset and the top-2
    Shapley effects."""
    vf = campaign.vf
    j13, jed = vf.names.index('k_13'), vf.names.index('ehrlich_downstream')

    def truth(U):
        return 4.0*U[:, j13]*U[:, jed] + 0.3*U[:, 0]        # conjunctive, like the Ehrlich pathway
    U10 = sa.sample_feasible(1500, vf.d, vf, np.random.default_rng(6))
    sur10 = sa.fit_surrogates(U10, truth(U10), seed=0)
    assert sur10.reliable, sur10.q2
    # d = 12 -> 4094 subsets: keep the offline check light (small base sample)
    S10, fb10 = sa.all_closed_indices({'y': sur10.predict}, vf, vf.d, n_base=256,
                                      n_replicates=1, seed=7)
    m10 = S10['y'][0]
    best10 = sa.best_subsets(m10, vf.d, sizes=(2,))
    assert set(sa.mask_names(best10[2][0], vf.names)) == {'k_13', 'ehrlich_downstream'}, best10
    sh10 = sa.shapley_effects(m10, vf.d)
    assert abs(sh10.sum() - 1.0) < 1e-9
    assert set(np.argsort(sh10)[-2:]) == {j13, jed}
    print(f'   max fallback fraction over subsets: {fb10.max():.4f}')


@needs_sklearn
def test_surrogate_candidates(sa):
    """11. Candidate selection: GP skipped unless asked for (the --gp-metrics
    cost saver); an empty or unknown candidate raises; the default still
    cross-validates both at unchanged Q2."""
    rng11 = np.random.default_rng(8)
    U11 = rng11.random((300, 3))
    y11 = np.sin(2*np.pi*U11[:, 0]) + U11[:, 1]**2
    sur_hgb = sa.fit_surrogates(U11, y11, seed=0, candidates=('hgb',))
    assert sur_hgb.name == 'hgb' and set(sur_hgb.q2) == {'hgb'}, sur_hgb.q2
    assert np.isfinite(sur_hgb.predict(U11[:7])).all()
    sur_gp = sa.fit_surrogates(U11, y11, seed=0, candidates=('gp',))
    assert sur_gp.name == 'gp' and set(sur_gp.q2) == {'gp'}, sur_gp.q2
    for bad in ((), ('rf',), ('hgb', 'rf')):
        with pytest.raises(ValueError):
            sa.fit_surrogates(U11, y11, seed=0, candidates=bad)
    sur_both = sa.fit_surrogates(U11, y11, seed=0)
    assert set(sur_both.q2) == {'gp', 'hgb'}, sur_both.q2
    assert sur_both.q2['gp'] == sur_gp.q2['gp'] and sur_both.q2['hgb'] == sur_hgb.q2['hgb']


@needs_workbooks
def test_screening_search_space(sa, ko, eb, campaign):
    """12. Screening measure: six axes re-specified (linear [0, ceiling]
    pathway axes, 0.1x-4x native band), order and other entries unchanged,
    both unit maps agree and invert, the feasible co-production corner mass
    rises from ~0 to > 3 %, three error paths."""
    space, kb, vf, bm = campaign.space, campaign.kb, campaign.vf, campaign.bm
    engine_kwargs = campaign.engine_kwargs
    scr = sa.screening_search_space(space, kb)
    assert list(scr) == list(space), 'insertion order changed'
    for n in sa.SCREENING_LINEAR_AXES:
        assert scr[n] == {'low': 0.0, 'high': space[n]['high'], 'log': False}, (n, scr[n])
    for n in ('k_3', 'k_6'):
        m_lo, m_hi = sa.SCREENING_NATIVE_BAND
        assert scr[n] == {'low': m_lo*kb[n], 'high': m_hi*kb[n], 'log': True}, (n, scr[n])
    assert scr['glycolysis'] == {'low': 0.1, 'high': 4.0, 'log': True}, scr['glycolysis']
    for n in space:
        if n not in sa.SCREENING_LINEAR_AXES + sa.SCREENING_NATIVE_AXES:
            assert scr[n] == space[n], n
    assert space['k_13']['log'] is True and space['k_13']['low'] > 0, 'the campaign space was mutated'
    # The measure is carried by the existing linear branch of both maps: vectorized
    # == scalar, external_to_unit inverts it, and u = 0 lands EXACTLY on 0.
    U12 = np.random.default_rng(12).random((50, len(scr)))
    U12[0] = 0.0
    U12[1] = 1.0
    vals12 = sa.unit_to_values(U12, scr)
    assert all(vals12[n][0] == 0.0 for n in sa.SCREENING_LINEAR_AXES)
    for i in (0, 1, 7):
        ext = ko.unit_to_external(U12[i], scr)
        for n in scr:
            assert np.isclose(float(ext[n]), float(vals12[n][i]), rtol=1e-12, atol=1e-12), (n, ext[n], vals12[n][i])
        back = ko.external_to_unit(ext, scr)
        for j, n in enumerate(scr):
            if not scr[n].get('int'):
                assert abs(back[j] - U12[i, j]) < 1e-9, (n, back[j], U12[i, j])
    # The point of the measure: the co-production corner gets real feasible mass.
    design_scr = json.loads(json.dumps(sa.design_record(
        search_space=scr, parameter_groups=engine_kwargs['parameter_groups'],
        group_references=engine_kwargs['group_references'], kinetic_baselines=kb,
        burden_model=bm, eb=eb, baseline_model_kwargs=BASE_FEED,
        baseline_max_n_spikes=16, volume_cap=20.0,
        target_conc_max=ko.TARGET_CONC_MAX, meta={'seed': 7})))
    vf_scr = sa.VectorizedFeasibility.from_design(design_scr)
    assert vf_scr.names == vf.names
    U12b = np.random.default_rng(13).random((20000, vf.d))

    def _corner(v):
        return ((v['k_13'] >= 2.0) & (v['k_17'] >= 1.0) & (v['ehrlich_downstream'] >= 0.4)
                & (v['k_3'] >= 0.1*kb['k_3']) & (v['k_6'] >= 0.1*kb['k_6']))
    mass_campaign = float((_corner(sa.unit_to_values(U12b, space)) & vf(U12b)).mean())
    mass_screening = float((_corner(sa.unit_to_values(U12b, scr)) & vf_scr(U12b)).mean())
    print(f'   feasible co-production corner mass: campaign {mass_campaign:.4%}, screening {mass_screening:.2%}')
    assert mass_screening > 0.03 and mass_screening > 20*max(mass_campaign, 1e-5), (mass_campaign, mass_screening)
    assert 0.15 < vf_scr(U12b).mean() < 0.45, vf_scr(U12b).mean()     # burden + volume feasible fraction
    # error paths: a missing axis, a nonpositive native baseline, an int-typed axis
    for bad_space, bad_kb, exc in (
            ({k: v for k, v in space.items() if k != 'k_17'}, kb, KeyError),
            (space, {**kb, 'k_3': 0.0}, ValueError),
            ({**space, 'k_13': {**space['k_13'], 'int': True}}, kb, ValueError)):
        with pytest.raises(exc):
            sa.screening_search_space(bad_space, bad_kb)


@needs_workbooks
def test_custom_search_space(sa, ko, eb, campaign):
    """13. Custom measure: every axis from a bands dict (absolute units;
    order-free, int axis kept), the shipped bands file matches, content-hash
    tag, both unit maps agree, six error paths."""
    space, kb, vf, bm = campaign.space, campaign.kb, campaign.vf, campaign.bm
    engine_kwargs = campaign.engine_kwargs
    bands13 = {n: {'low': 0.0, 'high': 4.0*(kb[n] if n in ('k_3', 'k_6') else 1.0), 'log': False}
               for n in ('k_3', 'k_6', 'glycolysis')}
    bands13.update({n: {'low': 0.0, 'high': space[n]['high'], 'log': False}
                    for n in ('k_13', 'k_17', 'ehrlich_downstream')})
    bands13.update({n: {'low': 0.75, 'high': 1.5, 'log': False}
                    for n in ('inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate')})
    bands13.update({n: {'low': space[n]['low'], 'high': space[n]['high'], 'log': space[n]['log']}
                    for n in ('threshold_conc', 'target_delta', 'max_n_spikes')})
    bands13 = {n: bands13[n] for n in reversed(list(space))}     # order must not matter
    cus = sa.custom_search_space(space, bands13)
    assert list(cus) == list(space), 'insertion order changed'
    assert cus['max_n_spikes'] == {'low': 0, 'high': 50, 'log': False, 'int': True}, cus['max_n_spikes']
    assert cus['k_3'] == {'low': 0.0, 'high': 4.0*kb['k_3'], 'log': False}, cus['k_3']
    assert space['k_3']['log'] is True, 'the campaign space was mutated'
    # The shipped bands file (package data, data/sobol_bands/) is this same set
    # (to the file's 4 significant figures).
    with open(os.path.join(PKG_DIR, 'data', 'sobol_bands',
                           'linear_0-4x_2026-09-22.json')) as fh:
        shipped = {k: v for k, v in json.load(fh).items() if not k.startswith('_')}
    shipped_space = sa.custom_search_space(space, shipped)
    for n in space:
        for key in ('low', 'high'):
            assert np.isclose(shipped_space[n][key], cus[n][key], rtol=1e-3, atol=0), (n, key)
        assert shipped_space[n]['log'] == cus[n]['log'], n
    # Hash: order-free, 8 hex, changes with any band.
    tag = sa.search_space_tag(bands13)
    assert len(tag) == 8 and tag == sa.search_space_tag(dict(reversed(list(bands13.items()))))
    assert tag != sa.search_space_tag({**bands13, 'k_3': {**bands13['k_3'], 'high': 20.0}})
    # Both unit maps agree on the all-linear space and u = 0 lands on 0.
    U13 = np.random.default_rng(14).random((50, len(cus)))
    U13[0] = 0.0
    vals13 = sa.unit_to_values(U13, cus)
    assert all(vals13[n][0] == 0.0 for n in ('k_3', 'k_6', 'glycolysis', 'k_13', 'k_17'))
    for i in (0, 3, 9):
        ext = ko.unit_to_external(U13[i], cus)
        for n in cus:
            assert np.isclose(float(ext[n]), float(vals13[n][i]), rtol=1e-12, atol=1e-12), (n, ext[n], vals13[n][i])
    design_cus = json.loads(json.dumps(sa.design_record(
        search_space=shipped_space, parameter_groups=engine_kwargs['parameter_groups'],
        group_references=engine_kwargs['group_references'], kinetic_baselines=kb,
        burden_model=bm, eb=eb, baseline_model_kwargs=BASE_FEED,
        baseline_max_n_spikes=16, volume_cap=20.0,
        target_conc_max=ko.TARGET_CONC_MAX, meta={'seed': 7})))
    vf_cus = sa.VectorizedFeasibility.from_design(design_cus)
    U13b = np.random.default_rng(15).random((20000, vf.d))
    feas13 = float(vf_cus(U13b).mean())
    print(f'   shipped custom bands: burden + volume feasible fraction {feas13:.1%}')
    assert 0.01 < feas13 < 1.0, feas13
    # error paths: a missing axis, an unknown axis, low >= high, a log band from 0, a log / fractional int axis
    for bad, exc in (({k: v for k, v in bands13.items() if k != 'k_17'}, KeyError),
                     ({**bands13, 'k_99': {'low': 0, 'high': 1, 'log': False}}, KeyError),
                     ({**bands13, 'k_3': {'low': 2.0, 'high': 2.0, 'log': False}}, ValueError),
                     ({**bands13, 'k_3': {'low': 0.0, 'high': 2.0, 'log': True}}, ValueError),
                     ({**bands13, 'max_n_spikes': {'low': 1, 'high': 50, 'log': True}}, ValueError),
                     ({**bands13, 'max_n_spikes': {'low': 0, 'high': 49.5, 'log': False}}, ValueError)):
        with pytest.raises(exc):
            sa.custom_search_space(space, bad)
