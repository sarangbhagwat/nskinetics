# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of nskinetics.analysis.profile_likelihood and ProfileResult."""
import dataclasses
import math

import matplotlib
import numpy as np
import pandas as pd
import pytest
import tellurium as te
from scipy.optimize import brentq
from scipy.stats import chi2
from scipy.stats import t as student_t

import nskinetics as nsk
from nskinetics.analysis import (CalibrationProblem, Experiment, FitParameter,
                                 FitResult, Observable, ProfileResult,
                                 fit_parameters, profile_likelihood)
from nskinetics.analysis import calibration as cal

matplotlib.use('Agg')

UNITS = {'time': 'h', 'conc': 'g/L'}

MM_TOY = """
model mm()
  S = 10; P = 0; Vmax = 1; Km = 2;
  J: S => P; Vmax*S/(Km+S);
end
"""

DECAY_TOY = """
model decay()
  S = 10; P = 0; k = 0.5;
  J: S => P; k*S;
end
"""

# Only the product k1*k2 is identifiable.
PRODUCT_TOY = """
model prod()
  S = 10; k1 = 0.5; k2 = 1;
  J: S => ; k1*k2*S;
end
"""

T = np.linspace(0.5, 6.0, 12)


def _model(text):
    return nsk.KineticModel(te.loadAntimonyModel(text), units=UNITS)


def _decay_problem(sigma=0.3, params=None, **kwargs):
    """S = 10 e^(-0.5 t) with seeded noise of standard deviation 0.3."""
    S = 10 * np.exp(-0.5 * T)
    S = S + 0.3 * np.random.default_rng(0).standard_normal(S.size)
    data = pd.DataFrame({'observable': 'S', 'time': T, 'value': S})
    exp = Experiment('E', nsk.Condition(t_end=6), data)
    return CalibrationProblem(
        _model(DECAY_TOY), [exp], [Observable('S', '[S]', sigma=sigma)],
        params or [FitParameter('k', (0.01, 10), initial=0.2)], **kwargs)


@pytest.fixture(scope='module')
def mm_fit():
    """Michaelis-Menten substrate depletion (Vmax 1, Km 2) with seeded
    noise (sigma 0.1)."""
    km = _model(MM_TOY)
    t = np.linspace(0.5, 14, 28)
    res = km.copy().run(nsk.Condition(times=np.concatenate([[0.0], t]),
                                      selections=('[S]',)))
    S = res['[S]'][1:] + 0.1 * np.random.default_rng(0).standard_normal(
        t.size)
    data = pd.DataFrame({'observable': 'S', 'time': t, 'value': S})
    prob = CalibrationProblem(
        km, [Experiment('E', nsk.Condition(t_end=14), data)],
        [Observable('S', '[S]', sigma=0.1)],
        [FitParameter('Vmax', (0.1, 10), initial=0.5),
         FitParameter('Km', (0.1, 20), initial=5)])
    return fit_parameters(prob, n_starts=2, seed=0)


@pytest.fixture(scope='module')
def mm_profiles(mm_fit):
    return profile_likelihood(mm_fit, ['Vmax', 'Km'], n_points=11)


@pytest.fixture(scope='module')
def decay_fit():
    return fit_parameters(_decay_problem(), n_starts=2, seed=0)


def _crossing(profile, side, threshold):
    """The threshold crossing on one side, by linear interpolation of
    sqrt(delta) against the fitting-space value between the bracketing grid
    values (exact for a quadratic profile; computed here from the profile
    table independently of ProfileResult)."""
    center = profile[profile['side'] == 'center'].iloc[0]
    rows = profile[profile['side'] == side]
    rows = rows.iloc[::-1] if side == 'lower' else rows
    x0, r0 = center['x'], math.sqrt(max(center['delta'], 0.0))
    for _, row in rows.iterrows():
        r = math.sqrt(row['delta'])
        if row['delta'] >= threshold:
            return x0 + (math.sqrt(threshold) - r0) / (r - r0) * (row['x']
                                                                 - x0)
        x0, r0 = row['x'], r
    return None


# --- (c) a well-determined pair: profile and Wald agree -----------------------
def test_threshold_is_the_chi2_quantile():
    assert cal.PROFILE_THRESHOLD == chi2.ppf(0.95, 1)
    assert math.isclose(cal.PROFILE_THRESHOLD, 3.841458820694124)


def test_michaelis_menten_profile_matches_wald(mm_fit, mm_profiles):
    truth = {'Vmax': 1.0, 'Km': 2.0}
    assert list(mm_profiles) == ['Vmax', 'Km']
    for name, prof in mm_profiles.items():
        lo, hi = prof.confidence_interval
        wald = mm_fit.ci.loc[name].to_numpy()
        # both 95% intervals contain the truth
        assert lo < truth[name] < hi
        assert wald[0] < truth[name] < wald[1]
        assert prof.is_identifiable and prof.n_failed == 0
        assert math.isclose(prof.estimate, mm_fit.estimates[name])
        # the profile interval is within 20% of the Wald interval
        assert np.allclose([lo, hi], wald, rtol=0.2)
        assert 0.8 < (hi - lo) / (wald[1] - wald[0]) < 1.2
        # 11 grid values, sorted, the optimum in the middle
        frame = prof.profile
        assert len(frame) == 11
        assert list(frame['side']) == ['lower'] * 5 + ['center'] + \
            ['upper'] * 5
        assert np.all(np.diff(frame['value']) > 0)
        center = frame.iloc[5]
        assert center['x'] == mm_fit.x[mm_fit.names.index(name)]
        assert center['cost'] == mm_fit.cost
        assert np.allclose(frame['value'], np.exp(frame['x']), rtol=1e-15)
        # the grid spans span = 3 SE in the fitting space
        se = mm_fit.se[name]
        assert np.isclose(frame['x'].iloc[-1] - frame['x'].iloc[0], 6 * se)
        # delta = 2 (cost - cost_min) >= 0; the other parameter moves
        assert prof.reference_cost <= mm_fit.cost
        assert np.all(frame['delta'] >= 0)
        assert np.allclose(frame['delta'], 2 * (frame['cost']
                                                - prof.reference_cost))
        other = 'Km' if name == 'Vmax' else 'Vmax'
        assert frame[f'estimate:{other}'].nunique() == 11
        assert (frame['status'] == 'ok').all()
        # interpolated crossings (linear in sqrt(delta))
        for side, end in zip(('lower', 'upper'), (lo, hi)):
            assert math.isclose(math.exp(_crossing(frame, side,
                                                   prof.threshold)), end,
                                rel_tol=1e-12)


def test_profile_points_are_reoptimized(mm_fit, mm_profiles):
    """At each grid value the other parameter minimizes the objective:
    the problem's own cost there equals the profile cost and is not
    lowered by nudging the other parameter."""
    prob = mm_fit.problem
    row = mm_profiles['Vmax'].profile.iloc[2]
    theta = {'Vmax': row['value'], 'Km': row['estimate:Km']}
    assert math.isclose(prob.cost(theta), row['cost'], rel_tol=1e-9)
    for f in (0.99, 1.01):
        nudged = dict(theta, Km=theta['Km'] * f)
        assert prob.cost(nudged) > row['cost']


# --- (d) a non-identifiable rate does not close --------------------------------
def test_non_identifiable_product_profile_does_not_close():
    S = 10 * np.exp(-0.5 * T)
    data = pd.DataFrame({'observable': 'S', 'time': T, 'value': S})
    prob = CalibrationProblem(
        _model(PRODUCT_TOY), [Experiment('E', nsk.Condition(t_end=6), data)],
        [Observable('S', '[S]', sigma=0.1)],
        [FitParameter('k1', (0.05, 5), initial=0.3),
         FitParameter('k2', (0.1, 10), initial=1.2)])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    # only k1*k2 is determined: the estimates are anti-correlated
    assert fit.correlation.loc['k1', 'k2'] < -0.99
    assert 'high_correlation' in set(fit.flags['flag'])
    assert np.isinf(fit.se['k1'])
    prof = profile_likelihood(fit, 'k1', n_points=7)['k1']
    assert not prof.is_identifiable
    assert prof.confidence_interval == (0.0, math.inf)
    frame = prof.profile
    # an infinite SE: the grid runs to both bounds, and k2 compensates
    assert frame['value'].iloc[0] == pytest.approx(0.05, rel=1e-12)
    assert frame['value'].iloc[-1] == pytest.approx(5.0, rel=1e-12)
    assert frame['delta'].max() < 1e-6
    k = frame['value'] * frame['estimate:k2']
    assert np.allclose(k[(frame['value'] >= 0.06) & (frame['value'] <= 4)],
                       0.5, rtol=1e-6)
    # over two decades the log axis keeps its own (decade) tick labels
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    fig, ax = prof.plot()
    assert ax.get_xscale() == 'log'
    assert not isinstance(ax.xaxis.get_major_formatter(), ScalarFormatter)
    plt.close(fig)


def test_local_and_observable_parameters_are_fixed_where_they_apply():
    """Profiling a local initial value and a virtual observable scale: at
    every grid value the profile cost is the problem's own cost at the
    fixed value and the re-optimized others."""
    rows = []
    for name, S0 in (('E1', 8.0), ('E2', 12.0)):
        S = S0 * np.exp(-0.5 * T)
        rows += [(name, 'S', t, v) for t, v in zip(T, S)]
        if name == 'E1':
            rows += [(name, 'P', t, 2.0 * v) for t, v in zip(T, S0 - S)]
    data = nsk.analysis.Dataset(pd.DataFrame(
        rows, columns=['experiment', 'observable', 'time', 'value']))
    exps = [Experiment(n, nsk.Condition(t_end=6), data[n])
            for n in ('E1', 'E2')]
    prob = CalibrationProblem(
        _model(DECAY_TOY), exps,
        [Observable('S', '[S]', sigma=0.1),
         Observable('P', '[P]', scale='a', sigma=0.1)],
        [FitParameter('k', (0.01, 10), initial=0.3),
         FitParameter('[S]', (1, 20), scope='local'),
         FitParameter('a', (0.1, 10))])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert np.allclose(fit.estimates, [0.5, 8.0, 12.0, 2.0], rtol=1e-6)
    profiles = profile_likelihood(fit, ['[S]@E2', 'a'], n_points=5)
    for name, truth in (('[S]@E2', 12.0), ('a', 2.0)):
        prof = profiles[name]
        lo, hi = prof.confidence_interval
        assert prof.is_identifiable and lo < truth < hi
        for _, row in prof.profile.iterrows():
            theta = {n: row[f'estimate:{n}'] for n in fit.names if n != name}
            theta[name] = row['value']
            assert math.isclose(prob.cost(theta), row['cost'], rel_tol=1e-9)


# --- grid rules ---------------------------------------------------------------
def test_a_side_continues_past_the_span_until_it_crosses(decay_fit):
    full = profile_likelihood(decay_fit, 'k', n_points=11)['k']
    short = profile_likelihood(decay_fit, 'k', n_points=11, span=0.5)['k']
    assert len(full.profile) == 11
    # half an SE does not reach the threshold: doubling steps follow, and
    # the crossing's bracket is bisected down to the grid step
    assert len(short.profile) > 11
    step = 0.1 * decay_fit.se['k']
    for side, sign in (('lower', -1.0), ('upper', 1.0)):
        rows = short.profile[short.profile['side'] == side]
        d = sign * (rows['x'].to_numpy() - decay_fit.x[0])
        order = np.argsort(d)
        d, delta = d[order], rows['delta'].to_numpy()[order]
        assert np.allclose(d[:5], step * np.arange(1, 6))
        # the outermost value ends the doubling steps 2, 4, ..., 2^m
        n = (d[-1] - d[4]) / step + 2.0
        assert n >= 4 and n == pytest.approx(2.0 ** round(math.log2(n)))
        # every value from the first crossing outward crosses, and that
        # crossing is at most one grid step beyond the last value below
        first = int(np.argmax(delta >= short.threshold))
        assert first > 4 and (delta[first:] >= short.threshold).all()
        assert delta[first - 1] < short.threshold
        assert d[first] - d[first - 1] <= step * (1 + 1e-9)
    assert short.is_identifiable
    assert np.allclose(short.confidence_interval, full.confidence_interval,
                       rtol=1e-3)


def test_an_extended_side_is_bisected_to_an_accurate_end():
    # four points, k and S0 free: the profile of k is flat on the lower
    # side and needs three doublings there; S0 is linear in the model
    # given k, so the exact profile is closed-form
    t = np.array([0.5, 1.0, 1.5, 2.0])
    y = 10 * np.exp(-0.5 * t) + np.random.default_rng(5).standard_normal(4)
    data = pd.DataFrame({'observable': 'S', 'time': t, 'value': y})
    prob = CalibrationProblem(
        _model(DECAY_TOY), [Experiment('E', nsk.Condition(t_end=6), data)],
        [Observable('S', '[S]', sigma=1.0)],
        [FitParameter('k', (0.01, 10)), FitParameter('[S]', (1, 50))])
    fit = fit_parameters(prob, n_starts=4, seed=0)

    def ssr(k):
        e = np.exp(-k * t)
        s0 = np.clip(np.sum(y * e) / np.sum(e * e), 1, 50)
        return np.sum((s0 * e - y) ** 2)

    k_hat = fit.estimates['k']
    lower = brentq(lambda k: ssr(k) - ssr(k_hat) - cal.PROFILE_THRESHOLD,
                   0.01, k_hat)
    prof = profile_likelihood(fit, 'k', n_points=21)['k']
    assert len(prof.profile) > 21                      # bisection points
    assert prof.confidence_interval[0] == pytest.approx(lower, rel=0.01)
    assert np.allclose(prof.profile['chi2'],
                       [ssr(v) for v in prof.profile['value']], atol=1e-6)


def test_single_parameter_profile_is_the_chi2_curve(decay_fit):
    """With one free parameter nothing is re-optimized: delta is the
    weighted SSR difference itself, and its 95% interval is close to the
    Wald one for this nearly linear problem."""
    prof = profile_likelihood(decay_fit, 'k', n_points=7)['k']
    data = decay_fit.residuals['data'].to_numpy()
    chi = [np.sum(((10 * np.exp(-k * T) - data) / 0.3) ** 2)
           for k in prof.profile['value']]
    assert np.allclose(prof.profile['chi2'], chi, rtol=1e-7)
    assert np.allclose(prof.profile['delta'], prof.profile['chi2']
                       - decay_fit.chi2, atol=1e-6)
    assert list(prof.profile.columns) == [
        'side', 'x', 'value', 'cost', 'chi2', 'delta', 'status', 'success',
        'nfev', 'n_failed_evals', 'message']
    assert np.allclose(prof.confidence_interval, decay_fit.ci.loc['k'],
                       rtol=2e-3)


def test_relative_sigma_scales_the_profile_by_reduced_chi2():
    fit = fit_parameters(_decay_problem(sigma=1.0, absolute_sigma=False),
                         n_starts=1, seed=0)
    prof = profile_likelihood(fit, 'k', n_points=7)['k']
    assert prof.scale == fit.reduced_chi2
    assert np.allclose(prof.profile['delta'],
                       (prof.profile['chi2'] - fit.chi2) / fit.reduced_chi2,
                       atol=1e-9)
    assert np.allclose(prof.confidence_interval, fit.ci.loc['k'], rtol=2e-3)
    # the noise level is estimated: Student's t (Wald) and F(1, dof) = t**2
    # (profile) replace the normal and chi-square quantiles
    q = student_t.ppf(0.975, fit.dof)
    assert math.isclose(prof.threshold, q ** 2, rel_tol=1e-12)
    se = fit.se['k']
    assert np.allclose(fit.ci.loc['k'], np.exp([fit.x[0] - q * se,
                                                fit.x[0] + q * se]),
                       rtol=1e-12)


def _linear_fit(absolute_sigma, n=5):
    """y = a S(t) + b with S(t) = 10 e^(-0.5 t) fixed: linear in (a, b), so
    delta is exactly quadratic and the Wald interval exact."""
    t = np.linspace(0.5, 6.0, n)
    s = 10 * np.exp(-0.5 * t)
    y = 2.0 * s + 1.0 + 0.3 * np.random.default_rng(0).standard_normal(n)
    data = pd.DataFrame({'observable': 'y', 'time': t, 'value': y})
    prob = CalibrationProblem(
        _model(DECAY_TOY), [Experiment('E', nsk.Condition(t_end=6), data)],
        [Observable('y', '[S]', scale='a', offset='b', sigma=0.3)],
        [FitParameter('a', (0.1, 10), log=False),
         FitParameter('b', (-10, 10), log=False)],
        absolute_sigma=absolute_sigma)
    return fit_parameters(prob, n_starts=1, seed=0), s, y


@pytest.mark.parametrize('absolute_sigma', [True, False])
def test_linear_model_wald_and_profile_intervals_are_exact(absolute_sigma):
    fit, s, y = _linear_fit(absolute_sigma)
    # closed-form weighted least squares
    X = np.column_stack([s, np.ones_like(s)]) / 0.3
    beta, *_ = np.linalg.lstsq(X, y / 0.3, rcond=None)
    cov = np.linalg.inv(X.T @ X)
    dof = len(y) - 2
    if absolute_sigma:
        q = chi2.ppf(0.95, 1) ** 0.5
    else:
        cov = cov * np.sum((X @ beta - y / 0.3) ** 2) / dof
        q = student_t.ppf(0.975, dof)
    se = np.sqrt(np.diag(cov))
    assert np.allclose(fit.x, beta, rtol=1e-7)
    assert np.allclose(fit.se, se, rtol=1e-6)
    assert np.allclose(fit.ci['lower'], beta - q * se, rtol=1e-6)
    assert np.allclose(fit.ci['upper'], beta + q * se, rtol=1e-6)
    # the profile, interpolated in sqrt(delta), gives the same interval on
    # a coarse grid (linear interpolation in delta would be ~1% narrower)
    profiles = profile_likelihood(fit, ['a', 'b'], n_points=11)
    for name, prof in profiles.items():
        assert math.isclose(prof.threshold, q ** 2, rel_tol=1e-12)
        assert np.allclose(prof.confidence_interval, fit.ci.loc[name],
                           rtol=1e-6)


def test_estimate_at_a_bound_leaves_that_side_open():
    S = 10 * np.exp(-0.5 * T)
    data = pd.DataFrame({'observable': 'S', 'time': T, 'value': S})
    prob = CalibrationProblem(
        _model(DECAY_TOY), [Experiment('E', nsk.Condition(t_end=6), data)],
        [Observable('S', '[S]', sigma=0.1)], [FitParameter('k', (0.6, 2))])
    fit = fit_parameters(prob, n_starts=1, seed=0)
    assert fit.active_mask.tolist() == [-1]
    prof = profile_likelihood(fit, 'k', n_points=7)['k']
    assert list(prof.profile['side']) == ['center'] + ['upper'] * 3
    lo, hi = prof.confidence_interval
    assert lo == 0.0 and 0.6 < hi < 0.61
    assert not prof.is_identifiable


def test_a_profile_below_the_optimum_warns(tmp_path):
    """A fit stopped early: the profile finds lower costs, warns, and uses
    the lowest one as the reference."""
    fit = fit_parameters(_decay_problem(), n_starts=1, seed=0, max_nfev=1)
    assert not fit.success
    with pytest.warns(UserWarning, match='below the fit'):
        prof = profile_likelihood(fit, 'k', n_points=5)['k']
    assert prof.reference_cost == prof.profile['cost'].min() < fit.cost
    assert prof.profile['delta'].min() == 0.0
    # the interval is searched from the profile's minimum, not from the
    # estimate (whose delta is far above the threshold): it holds the
    # minimum, not the estimate, and is not identifiable by accident
    lo, hi = prof.confidence_interval
    vmin = prof.profile.loc[prof.profile['delta'].idxmin(), 'value']
    assert lo < vmin < hi
    assert not lo <= prof.estimate <= hi
    center = prof.profile[prof.profile['side'] == 'center']
    assert center['delta'].iloc[0] > prof.threshold
    path = tmp_path / 'prof.csv'
    prof.to_csv(path)
    back = ProfileResult.from_csv(path)
    assert back.confidence_interval == prof.confidence_interval
    assert back.is_identifiable == prof.is_identifiable


def test_two_workers_match_serial(mm_fit, mm_profiles):
    par = profile_likelihood(mm_fit, ['Vmax', 'Km'], n_points=11,
                             workers=2)
    for name, prof in mm_profiles.items():
        pd.testing.assert_frame_equal(par[name].profile, prof.profile)
        assert par[name].confidence_interval == prof.confidence_interval


# --- validation, CSV and plot --------------------------------------------------
def test_invalid_arguments_raise_before_any_run(mm_fit, monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(nsk.KineticModel, 'run',
                        lambda *a, **k: calls.append(1))
    path = tmp_path / 'fit.csv'
    mm_fit.to_csv(path)
    cases = [
        ((object(), 'Vmax'), {}, TypeError, 'FitResult'),
        ((FitResult.from_csv(path), 'Vmax'), {}, ValueError, 'no problem'),
        ((mm_fit, 'k'), {}, ValueError, 'free parameters'),
        ((mm_fit, []), {}, ValueError, 'free parameters'),
        ((mm_fit, ['Vmax', 'Vmax']), {}, ValueError, 'Duplicate'),
        ((mm_fit, 'Vmax'), {'n_points': 2}, ValueError, 'n_points'),
        ((mm_fit, 'Vmax'), {'span': 0}, ValueError, 'span'),
        ((mm_fit, 'Vmax'), {'span': math.nan}, ValueError, 'span'),
        ((mm_fit, 'Vmax'), {'workers': 0}, ValueError, 'workers'),
    ]
    for args, kwargs, error, match in cases:
        with pytest.raises(error, match=match):
            profile_likelihood(*args, **kwargs)
    no_dof = dataclasses.replace(mm_fit, absolute_sigma=False, chi2=0.0)
    with pytest.raises(ValueError, match='reduced chi-square'):
        profile_likelihood(no_dof, 'Vmax')
    # a robust loss with absolute_sigma=False (CalibrationProblem refuses it)
    robust_rel = dataclasses.replace(mm_fit, absolute_sigma=False,
                                     loss='soft_l1')
    with pytest.raises(ValueError, match='needs absolute sigmas'):
        profile_likelihood(robust_rel, 'Vmax')
    assert calls == []


def _assert_same_profile(a, b):
    for f in dataclasses.fields(ProfileResult):
        if f.name == 'profile':
            pd.testing.assert_frame_equal(a.profile, b.profile)
        else:
            assert getattr(a, f.name) == getattr(b, f.name), f.name
    assert a.confidence_interval == b.confidence_interval
    assert a.is_identifiable == b.is_identifiable
    assert a.reference_cost == b.reference_cost


def test_csv_round_trip_is_exact(mm_profiles, tmp_path):
    prof = mm_profiles['Km']
    path = tmp_path / 'profile.csv'
    prof.to_csv(path)
    back = ProfileResult.from_csv(path)
    _assert_same_profile(prof, back)
    again = tmp_path / 'again.csv'
    back.to_csv(again)
    assert path.read_text() == again.read_text()
    # an open profile with an infinite SE round-trips too
    open_prof = dataclasses.replace(prof, se=math.inf)
    open_prof.to_csv(path)
    _assert_same_profile(open_prof, ProfileResult.from_csv(path))
    # foreign and truncated files
    with pytest.raises(ValueError, match='not written by ProfileResult'):
        pd.DataFrame({'a': [1]}).to_csv(path, index=False)
        ProfileResult.from_csv(path)
    prof.to_csv(path)
    lines = path.read_text().splitlines()
    path.write_text('\n'.join(lines[:-1]) + '\n')
    with pytest.raises(ValueError, match='incomplete'):
        ProfileResult.from_csv(path)
    with pytest.raises(ValueError, match='lacks'):
        dataclasses.replace(prof, profile=prof.profile.drop(columns='cost'))


# dx/dt = k x^2 blows up at t = 1/k: for k > 1/0.9 a run to 0.9 fails.
BLOW_TOY = """
model blow()
  x = 1; k = 0.5;
  J: => x; k*x^2;
end
"""


def test_failed_grid_points_are_recorded_and_never_cross(tmp_path):
    import matplotlib.pyplot as plt
    t = np.linspace(0.1, 0.9, 9)
    data = pd.DataFrame({'observable': 'x', 'time': t,
                         'value': 1 / (1 - 0.5 * t)})
    prob = CalibrationProblem(
        _model(BLOW_TOY), [Experiment('E', nsk.Condition(t_end=0.9), data)],
        [Observable('x', 'x', sigma=3.0)],
        [FitParameter('k', (0.1, 3), initial=0.2)])
    fit = fit_parameters(prob, n_starts=2, seed=0)
    prof = profile_likelihood(fit, 'k', n_points=7)['k']
    frame = prof.profile
    failed = frame['status'] == 'failed'
    assert prof.n_failed == failed.sum() > 0
    # the runs past the blow-up fail, and their rise is unknown
    assert (frame.loc[failed, 'value'] > 1 / 0.9).all()
    assert frame.loc[failed, ['cost', 'delta']].isna().all().all()
    ok = frame.loc[~failed]
    assert np.isfinite(ok['delta']).all() and (ok['value'] < 1 / 0.9).all()
    # sigma 3: no ok point reaches the threshold on either side, and the
    # failed points are not taken as crossings, so both ends stay open
    assert (ok['delta'] < prof.threshold).all()
    assert frame['value'].iloc[0] == pytest.approx(0.1, rel=1e-12)
    assert prof.confidence_interval == (0.0, math.inf)
    assert not prof.is_identifiable
    path = tmp_path / 'profile.csv'
    prof.to_csv(path)
    back = ProfileResult.from_csv(path)
    _assert_same_profile(prof, back)
    again = tmp_path / 'again.csv'
    back.to_csv(again)
    assert path.read_text() == again.read_text()
    prof.plot()
    plt.close('all')


def test_plot_smoke(mm_profiles):
    import matplotlib.pyplot as plt
    prof = mm_profiles['Vmax']
    fig, ax = prof.plot()
    labels = [line.get_label() for line in ax.get_lines()]
    assert {'profile', 'Wald', 'estimate'} <= set(labels)
    assert ax.get_xscale() == 'log'
    profile_line = ax.get_lines()[labels.index('profile')]
    assert np.allclose(profile_line.get_ydata(), prof.profile['delta'])
    assert ax.get_ylabel() == r'$\Delta\chi^2$'
    fig0, ax0 = plt.subplots()
    assert prof.plot(ax=ax0) == (fig0, ax0)
    # a robust loss or absolute_sigma=False: not a chi-square difference
    for changes, label in (({'loss': 'soft_l1'}, r'$2\,\Delta\mathrm{cost}$'),
                           ({'scale': 2.0}, r'$2\,\Delta\mathrm{cost}/s^2$')):
        _, ax = dataclasses.replace(prof, **changes).plot()
        assert ax.get_ylabel() == label
    plt.close('all')


def test_plot_labels_plain_values_within_a_decade(mm_profiles):
    # A log axis spanning less than a decade has no decade ticks, and
    # matplotlib's log formatter would label crowded minor ticks such as
    # '3.99x10^-1'; the plot labels a few plain values instead.
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    prof = mm_profiles['Vmax']
    fig, ax = prof.plot()
    lo, hi = ax.get_xlim()
    assert ax.get_xscale() == 'log' and hi < 10 * lo
    assert isinstance(ax.xaxis.get_major_formatter(), ScalarFormatter)
    assert len(ax.xaxis.get_minorticklocs()) == 0
    fig.canvas.draw()
    ticks = [t for t in ax.xaxis.get_majorticklocs() if lo <= t <= hi]
    labels = [t.get_text() for t in ax.get_xticklabels() if t.get_text()]
    assert 2 <= len(ticks) <= 6
    assert not any('times' in text or '^' in text for text in labels)
    plt.close(fig)
