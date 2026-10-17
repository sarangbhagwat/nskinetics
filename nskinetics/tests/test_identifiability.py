# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Practical identifiability (F2 item 4; spec tests f, j, k), plus the
Cramer-Rao standard error against its closed form on a one-parameter
decay."""
import dataclasses
import math
import warnings

import matplotlib
import numpy as np
import pandas as pd
import pytest

import nskinetics as nsk
from nskinetics.analysis import (IdentifiabilityResult, identifiability,
                                 local_sensitivity)
from nskinetics.analysis import _fd
from nskinetics.analysis.identifiability import (NULL_LOADING,
                                                RELATIVE_FLOOR,
                                                _rank_test_doubtful,
                                                _svd_analysis)

matplotlib.use('Agg')

UNITS = {'time': 'h', 'conc': 'g/L'}
TIGHT = (1e-12, 1e-10)          # integrator (atol, rtol) well below rel_step
CVODE = (1e-12, 1e-6)           # CVODE's own default tolerances

# Only the product k1*k2 enters the rate: ln k1 + ln k2 is identifiable,
# ln k1 - ln k2 is not. k3 acts on a second step and is identifiable.
PRODUCT_TOY = """
model product()
  S = 10; P = 0; Q = 0; k1 = 0.5; k2 = 2; k3 = 0.3;
  J1: S => P; k1*k2*S;
  J2: P => Q; k3*P;
end
"""

# S = S0 e^(-kt): dS/dk = -t S and dS/dS0 = e^(-kt).
DECAY_TOY = """
model decay()
  S = 10; P = 0; k = 0.5;
  J: S => P; k*S;
end
"""

# The F2 event toy with an inert parameter z.
EVENT_TOY = """
model ev()
  S = 10; n = 0; k = 0.5; z = 1;
  J: S => ; k*S;
  E1: at S < 5: S = S + 5, n = n + 1;
end
"""


@pytest.fixture(scope='module')
def product_ls():
    km = nsk.KineticModel.from_antimony(PRODUCT_TOY, units=UNITS)
    return local_sensitivity(km, nsk.Condition(t_end=5, n_steps=21),
                             ['k1', 'k2', 'k3'], outputs=['[S]', '[P]'],
                             tolerances=TIGHT)


@pytest.fixture(scope='module')
def decay_ls():
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    return local_sensitivity(km, nsk.Condition(t_end=10, n_steps=11),
                             ['k', '[S]'], outputs=['[S]', '[P]'],
                             tolerances=TIGHT)


# --- (f) identifiability -----------------------------------------------------
def test_product_rate_is_rank_deficient_with_a_log_null_direction(product_ls):
    ident = identifiability(product_ls, sigma={'[S]': 0.5, '[P]': 0.5})
    assert ident.parameters == ['k1', 'k2', 'k3']
    assert ident.space == ['log', 'log', 'log']
    assert ident.numerical_rank == 2
    assert ident.condition_number > 1e8
    s = ident.singular_values
    assert s.shape == (3,) and np.all(np.diff(s) <= 0)
    assert s[2] < 1e-6 * s[0]
    # The FIM is S^T S; its eigenvalues are the squared singular values.
    assert np.allclose(np.linalg.eigvalsh(ident.fim.to_numpy())[::-1], s ** 2,
                       rtol=1e-9, atol=1e-9 * s[0] ** 2)
    # One null direction: +-(1, -1, 0)/sqrt(2) in (ln k1, ln k2, ln k3).
    null = ident.null_directions
    assert list(null.columns) == ['v3']
    assert list(null.index) == ident.parameters
    v = null['v3'].to_numpy()
    assert np.allclose(np.abs(v[:2]), 1 / math.sqrt(2), atol=1e-6)
    assert v[0] * v[1] < 0 and abs(v[2]) < 1e-6
    assert v[np.argmax(np.abs(v))] > 0              # deterministic sign
    # Cramer-Rao: infinite for k1 and k2, finite for k3; no crash.
    se = ident.se
    assert np.isinf(se['k1']) and np.isinf(se['k2'])
    assert np.isfinite(se['k3']) and se['k3'] > 0
    table = ident.table
    assert list(table['identifiable']) == [False, False, True]
    assert list(table['null_direction']) == [True, True, False]
    assert np.allclose(table.loc[['k1', 'k2'], 'null_loading'],
                       1 / math.sqrt(2), atol=1e-6)
    corr = ident.correlation
    assert np.isnan(corr.loc['k1', 'k2']) and corr.loc['k3', 'k3'] == 1.0
    # Collinearity: the (k1, k2) pair is (numerically) singular.
    assert ident.gamma(['k1', 'k2']) > 1e3
    col = ident.collinearity(subset_size=2)
    assert col.attrs['search'] == 'exhaustive' and len(col) == 3
    assert col['parameters'].iloc[0] == ('k1', 'k2')
    assert col['gamma'].iloc[0] > 1e3 and np.all(col['gamma'].iloc[1:] < 5)
    # The subset keeps exactly one of k1/k2 (equal importance: the first).
    best = ident.best_subset(gamma_max=10)
    assert sorted(best) == ['k1', 'k3']
    assert len({'k1', 'k2'} & set(best)) == 1
    # Without a collinearity cap the exactly collinear pair is still never
    # kept together (the subset must be full rank).
    assert sorted(ident.best_subset(gamma_max=math.inf)) == ['k1', 'k3']


def test_well_posed_two_parameter_toy_has_low_collinearity(decay_ls):
    ident = identifiability(decay_ls, sigma={'[S]': 0.5})
    assert ident.numerical_rank == 2 and ident.null_directions.shape == (2, 0)
    gamma = ident.gamma(['k', '[S]'])
    assert 1 <= gamma < 5
    col = ident.collinearity(2)
    assert col['parameters'].iloc[0] == ('k', '[S]')
    assert col['gamma'].iloc[0] == pytest.approx(gamma, rel=1e-12)
    assert np.all(np.isfinite(ident.se)) and ident.table['identifiable'].all()
    assert sorted(ident.best_subset()) == ['[S]', 'k']
    # A tighter cap than the pair's collinearity keeps the more important
    # parameter only.
    first = ident.delta_msqr.idxmax()
    assert ident.best_subset(gamma_max=1.0) == [first]


def test_cramer_rao_se_matches_closed_form_for_one_parameter_decay():
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    cond = nsk.Condition(t_end=10, n_steps=11)
    ls = local_sensitivity(km, cond, ['k'], outputs=['[S]'],
                           tolerances=TIGHT)
    t = np.linspace(0, 10, 11)
    S = 10 * np.exp(-0.5 * t)
    sigma = 0.5
    se_k = 1 / math.sqrt(np.sum((t * S / sigma) ** 2))
    ident = identifiability(ls, sigma={'[S]': sigma})
    table = ident.table
    assert table.loc['k', 'se_natural'] == pytest.approx(se_k, rel=1e-5)
    # In ln k: the relative standard error.
    assert ident.se['k'] == pytest.approx(se_k / 0.5, rel=1e-5)
    assert table.loc['k', 'rel_se'] == pytest.approx(se_k / 0.5, rel=1e-5)
    assert ident.fim.loc['k', 'k'] == pytest.approx(
        np.sum((0.5 * t * S / sigma) ** 2), rel=1e-5)
    assert ident.covariance.loc['k', 'k'] == pytest.approx(
        (se_k / 0.5) ** 2, rel=1e-5)
    assert ident.delta_msqr['k'] == pytest.approx(
        math.sqrt(np.mean((0.5 * t * S / sigma) ** 2)), rel=1e-5)
    # One sigma per time.
    sig = np.linspace(0.2, 1.2, 11)
    ident = identifiability(ls, sigma={'[S]': sig})
    assert ident.table.loc['k', 'se_natural'] == pytest.approx(
        1 / math.sqrt(np.sum((t * S / sig) ** 2)), rel=1e-5)
    assert ident.sigma == {'[S]': [float(v) for v in sig]}


def test_relative_sigma_uses_log_sensitivities_and_keeps_every_row(decay_ls):
    ident = identifiability(decay_ls, sigma=('relative', 0.1),
                            params=['k'])
    t = decay_ls.time
    # [S]: d ln S / d ln k = -k t. [P] = 0 at t0 is a row too (its sigma is
    # the floor): no row is dropped, and a zero row adds no information.
    assert ident.row_outputs == ['[S]'] * 11 + ['[P]'] * 11
    assert np.array_equal(ident.row_times, np.concatenate([t, t]))
    assert np.allclose(ident.sensitivity[:11, 0], -0.5 * t / 0.1, rtol=1e-5,
                       atol=1e-9)
    assert ident.sensitivity[11, 0] == 0.0
    assert ident.sigma == ('relative', 0.1)
    assert ident.relative_floor == 1e-3 == RELATIVE_FLOOR
    one = identifiability(decay_ls, sigma={'[S]': ('relative', 0.1)},
                          params=['k'])
    assert one.row_outputs == ['[S]'] * 11
    assert one.fim.loc['k', 'k'] == pytest.approx(
        np.sum((0.5 * t / 0.1) ** 2), rel=1e-5)
    # Absolute sigma: no floor is recorded.
    assert identifiability(decay_ls, {'[S]': 0.5}).relative_floor is None


def test_relative_sigma_is_bounded_below_by_the_floor():
    # S = S0 e^(-kt) falls below 1e-3 of its peak S0 after t = 13.8: there
    # sigma = rel * floor * S0, above that rel * S. In closed form,
    # S_i,k = -k t S / sigma and S_i,S0 = S / sigma (finite differences
    # at S / S0 ~ 5e-5 are accurate to ~1e-5).
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=20, n_steps=21),
                           ['k', '[S]'], outputs=['[S]', '[P]'],
                           tolerances=TIGHT)
    t = ls.time
    S = 10 * np.exp(-0.5 * t)
    P = 10 - S
    rel, floor = 0.1, 1e-3
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ident = identifiability(ls, ('relative', rel))
    sig_s = rel * np.maximum(S, floor * 10)
    sig_p = rel * np.maximum(P, floor * P.max())
    depleted = S < floor * 10
    assert depleted.sum() == 7 and not depleted[:14].any()
    s_rows = ident.sensitivity[:21]
    assert np.allclose(s_rows[:, 0], -0.5 * t * S / sig_s, rtol=1e-4,
                       atol=1e-12)
    assert np.allclose(s_rows[:, 1], S / sig_s, rtol=1e-4, atol=1e-12)
    # On the depleted rows: -k t S / (rel floor S0), not -k t / rel.
    assert np.allclose(s_rows[depleted, 0],
                       -0.5 * t[depleted] * S[depleted] / (rel * floor * 10),
                       rtol=1e-4)
    assert np.all(np.abs(s_rows[depleted, 0])
                  < 0.5 * t[depleted] / rel)
    # [P] starts at 0: its first rows sit on the floor of its own peak.
    p_rows = ident.sensitivity[21:]
    assert p_rows[0, 0] == 0.0 and np.abs(p_rows[0, 1]) < 1e-9
    assert np.allclose(p_rows[:, 0], 0.5 * t * S / sig_p, rtol=1e-4,
                       atol=1e-12)
    # The peak is over the whole grid: times= does not change a row.
    late = identifiability(ls, ('relative', rel), times=[14.0, 20.0])
    assert np.array_equal(late.sensitivity,
                          ident.sensitivity[[14, 20, 35, 41]])
    # A floor of 1 is a constant sigma = rel * peak: the absolute model.
    top = identifiability(ls, ('relative', rel), relative_floor=1.0)
    peak = np.abs(ls.nominal).max(axis=0)
    flat = identifiability(ls, {'[S]': rel * peak[0], '[P]': rel * peak[1]})
    assert top.relative_floor == 1.0 and flat.relative_floor is None
    assert np.allclose(top.sensitivity, flat.sensitivity, rtol=1e-12,
                       atol=0)


# Only k1 + k2 matters; run to t = 30, S is depleted to ~5e-12.
DEPLETED_TOY = """
model s()
  S = 10; P = 0; k1 = 0.5; k2 = 1;
  J: S => P; (k1 + k2)*S;
end
"""


def test_relative_floor_keeps_the_verdict_on_a_depleted_output():
    # Under a pure relative sigma, 0.1 |S|, the late rows of S are FD noise
    # over a tiny sigma and make the exact non-identifiability look full
    # rank with se ~ 1e-6.
    km = nsk.KineticModel.from_antimony(DEPLETED_TOY, units=UNITS)
    cond = nsk.Condition(t_end=30, n_steps=61)
    ls = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]', '[P]'])
    assert ls.nominal[-1, 0] < 1e-10
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ident = identifiability(ls, ('relative', 0.1))
        alone = identifiability(ls, {'[S]': ('relative', 0.1)})
    for result in (ident, alone):
        # The null singular value is ~3.6e-7 of the largest: a ~3x margin
        # to rank_rtol, about half of it O(h^2) truncation error (it
        # shrinks with a larger rel_step). A failure here on another
        # sundials/roadrunner build points to integrator noise rather than
        # to the floor logic.
        s = result.singular_values
        assert s[-1] / s[0] < result.rank_rtol, s[-1] / s[0]
        assert result.numerical_rank == 1
        assert np.isinf(result.se).all()
        assert result.gamma(['k1', 'k2']) > 1e3
        assert not result.table['identifiable'].any()
        assert result.best_subset() == ['k2']          # the more important
    assert len(ident.row_outputs) == 122                # no row dropped
    assert ident.rank_test_doubtful is False
    # Without an effective floor the noise rows break the verdict; the
    # error estimate, scaled by the realized weight amplification (1e12 at
    # the t0 row of P, which is 0), flags it.
    with pytest.warns(UserWarning, match='relative_floor=1e-15') as record:
        pure = identifiability(ls, ('relative', 0.1), relative_floor=1e-15)
    assert len(record) == 1
    assert 'Raise relative_floor to the default 0.001' in str(
        record[0].message)
    assert pure.numerical_rank == 2
    assert pure.rank_test_doubtful is True
    assert pure.fd_error_level == pytest.approx(ident.fd_error_level * 1e12,
                                                rel=1e-12)
    # Smaller floors than the default already lift the null (concepts.rst):
    # each is flagged.
    for floor in (1e-4, 1e-5):
        with pytest.warns(UserWarning, match=f'relative_floor={floor:g}, '
                                             'below the default'):
            low = identifiability(ls, ('relative', 0.1), relative_floor=floor)
        assert low.numerical_rank == 2 and low.rank_test_doubtful is True
    # Late times only: every row is below the floor, so the data carry no
    # usable information -- nothing identifiable, an empty subset.
    late = identifiability(ls, {'[S]': ('relative', 0.1)},
                           times=[25.0, 30.0])
    assert not late.table['identifiable'].any()
    assert late.best_subset() == []


def test_relative_sigma_keeps_well_scaled_rows():
    km = nsk.KineticModel.from_antimony("""
    model decay()
      S = 10; P = 0; k = 0.5; z = 1;
      J: S => P; k*S;
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=2, n_steps=5), ['k', 'z'],
                           outputs=['[S]', '[P]'])
    t = ls.time
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        one = identifiability(ls, {'[S]': ('relative', 0.1)}, params=['k'])
        both = identifiability(ls, ('relative', 0.1), params=['k'])
        inert = identifiability(ls, ('relative', 0.1), params=['z'])
    assert one.row_outputs == ['[S]'] * 5
    assert one.se['k'] == pytest.approx(
        1 / math.sqrt(np.sum((0.5 * t / 0.1) ** 2)), rel=1e-5)
    assert both.row_outputs == ['[S]'] * 5 + ['[P]'] * 5
    # an inert parameter: rank 0 and an infinite se, not an error
    assert inert.numerical_rank == 0 and np.isinf(inert.se['z'])
    assert len(inert.row_outputs) == 10
    # A single measurement time: d ln S / d ln k = -k t and
    # d ln P / d ln k = k t e^(-kt) / (1 - e^(-kt)) at t = 2.
    ls = local_sensitivity(km, nsk.Condition(t_end=6, n_steps=61), ['k'],
                           outputs=['[S]', '[P]'])
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ident = identifiability(ls, ('relative', 0.1), times=[2.0])
    e_p = math.exp(-1) / (1 - math.exp(-1))
    assert ident.row_outputs == ['[S]', '[P]']
    assert ident.se['k'] == pytest.approx(
        1 / math.sqrt((1 / 0.1) ** 2 + (e_p / 0.1) ** 2), rel=1e-3)


def test_growing_products_keep_every_row_on_a_fine_grid():
    # B and C grow from 0; on a fine grid their first rows lie far below
    # their peaks. They are all kept, with the floor's weight, and the
    # verdicts do not depend on the scale of A0.
    found = {}
    for A0, conc in ((1.0, 'g/L'), (1e-3, 'mol/L')):
        km = nsk.KineticModel.from_antimony(CHAIN_TOY.format(A0=A0),
                                            units={'time': 'h', 'conc': conc})
        ls = local_sensitivity(km, nsk.Condition(t_end=10, n_steps=1001),
                               ['k1', 'k2'], outputs=['[A]', '[B]', '[C]'])
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            ident = identifiability(ls, ('relative', 0.1))
            only_c = identifiability(ls, {'[C]': ('relative', 0.1)})
        c = ls.nominal[:, 2]
        assert (np.abs(c) < 1e-3 * c.max()).sum() > 5   # rows on the floor
        assert len(ident.row_outputs) == 3003
        assert len(only_c.row_outputs) == 1001
        for result in (ident, only_c):
            assert result.numerical_rank == 2
            assert result.table['identifiable'].all()
            assert sorted(result.best_subset()) == ['k1', 'k2']
        found[A0] = ident.singular_values, only_c.se.to_numpy()
    assert np.allclose(found[1.0][0], found[1e-3][0], rtol=1e-4)
    assert np.allclose(found[1.0][1], found[1e-3][1], rtol=1e-4)


def test_scaling_one_output_leaves_the_analysis_unchanged(decay_ls):
    # Output units: [P] in mg/L instead of g/L multiplies its values and
    # derivatives by 1000. A relative sigma (floor included) scales alike,
    # and so does an absolute sigma given in the new unit.
    nominal, raw = decay_ls.nominal.copy(), decay_ls.raw.copy()
    nominal[:, 1] *= 1000
    raw[:, 1, :] *= 1000
    milli = dataclasses.replace(decay_ls, nominal=nominal, raw=raw)
    pairs = [(('relative', 0.1), ('relative', 0.1)),
             ({'[P]': ('relative', 0.05)}, {'[P]': ('relative', 0.05)}),
             ({'[S]': 0.5, '[P]': 0.2}, {'[S]': 0.5, '[P]': 200.0})]
    for sigma, sigma_milli in pairs:
        a = identifiability(decay_ls, sigma, relative_floor=0.05)
        b = identifiability(milli, sigma_milli, relative_floor=0.05)
        assert np.allclose(a.sensitivity, b.sensitivity, rtol=1e-12, atol=0)
        assert np.allclose(a.singular_values, b.singular_values, rtol=1e-12,
                           atol=0)
        assert a.numerical_rank == b.numerical_rank
        pd.testing.assert_series_equal(a.table['identifiable'],
                                       b.table['identifiable'])
        assert a.best_subset() == b.best_subset()


def test_an_output_that_is_always_zero_has_no_relative_sigma():
    km = nsk.KineticModel.from_antimony("""
    model z()
      species Z;
      S = 10; P = 0; Z = 0; k = 0.5;
      J: S => P; k*S;
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=4, n_steps=9), ['k'],
                           outputs=['[S]', '[Z]'], tolerances=TIGHT)
    assert not ls.nominal[:, 1].any()
    with pytest.warns(UserWarning, match=r"\['\[Z\]'\] are 0 at every.*"
                      r"give them an absolute sigma in a mapping"):
        relative = identifiability(ls, ('relative', 0.1))
    assert relative.row_outputs == ['[S]'] * 9
    with pytest.raises(ValueError, match=r"No measured output.*\['\[Z\]'\]"):
        identifiability(ls, {'[Z]': ('relative', 0.1)})
    # An absolute sigma measures it (rows of zeros), also next to a relative
    # one for the others.
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ident = identifiability(ls, {'[S]': 0.5, '[Z]': 0.5})
        mixed = identifiability(ls, {'[S]': ('relative', 0.1), '[Z]': 0.5})
    assert ident.row_outputs == ['[S]'] * 9 + ['[Z]'] * 9
    assert not ident.sensitivity[9:].any()
    assert mixed.row_outputs == ['[S]'] * 9 + ['[Z]'] * 9
    assert np.array_equal(mixed.sensitivity[:9], relative.sensitivity)
    assert not mixed.sensitivity[9:].any()
    assert mixed.relative_floor == RELATIVE_FLOOR
    # When the only relative output is left out, every row is absolute: no
    # floor is recorded, and the matrix is the absolute one.
    with pytest.warns(UserWarning, match=r"\['\[Z\]'\] are 0 at every"):
        only_abs = identifiability(ls, {'[S]': 0.5, '[Z]': ('relative', 0.1)})
    assert only_abs.relative_floor is None
    assert np.array_equal(only_abs.sensitivity,
                          identifiability(ls, {'[S]': 0.5}).sensitivity)


def test_a_floor_without_relative_rows_does_not_gate_the_rank_test():
    # The error level (3 rtol / h = 3e-4 at rtol 1e-7) exceeds rank_rtol,
    # but s_min / s_max sits far above 100 x it. A lowered relative_floor
    # whose only relative output (Z, always 0) is left out shapes no row:
    # it is not recorded, does not scale the level, and changes nothing.
    km = nsk.KineticModel.from_antimony("""
    model z()
      S = 10; Z = 0; k = 0.5; c = 2; kz = 0;
      J: S => ; k*S;
      J2: => S; c;
      Jz: => Z; kz*S;
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=5, n_steps=21),
                           ['k', 'c'], outputs=['[S]', '[Z]'],
                           tolerances=(1e-12, 1e-7))
    sigma = {'[S]': 0.5, '[Z]': ('relative', 0.1)}
    found = {}
    for floor in (1e-3, 1e-4):
        with pytest.warns(UserWarning) as record:
            found[floor] = identifiability(ls, sigma, relative_floor=floor)
        messages = [str(w.message) for w in record]
        assert len(messages) == 1, messages
        assert "['[Z]'] are 0 at every" in messages[0]
        result = found[floor]
        assert result.row_outputs == ['[S]'] * 21
        assert result.relative_floor is None
        assert result.fd_error_level > result.rank_rtol
        assert result.rank_test_doubtful is False
    assert np.array_equal(found[1e-3].sensitivity, found[1e-4].sensitivity)


def test_a_sigma_mapping_can_mix_absolute_and_relative_entries(decay_ls,
                                                               tmp_path):
    # {'[S]': 0.5, '[P]': ('relative', 0.05)} is the absolute model on [S]
    # and the bounded relative one on [P]: the same matrix as the per-time
    # array rel * max(|P|, relative_floor * max_t |P|), whose peak is over
    # the whole grid whatever the selected times.
    k = decay_ls.outputs.index('[P]')
    P = decay_ls.nominal[:, k]
    for times in (None, [8, 2, 6]):
        idx = (np.arange(len(decay_ls.time)) if times is None
               else np.searchsorted(decay_ls.time, times))
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            mixed = identifiability(decay_ls, {'[S]': 0.5,
                                               '[P]': ('relative', 0.05)},
                                    times=times)
        array = 0.05 * np.maximum(np.abs(P[idx]), 1e-3 * np.max(np.abs(P)))
        manual = identifiability(decay_ls, {'[S]': 0.5, '[P]': array},
                                 times=times)
        assert np.allclose(mixed.sensitivity, manual.sensitivity,
                           rtol=1e-14, atol=0)
        assert mixed.row_outputs == manual.row_outputs
        assert mixed.sigma['[P]'] == ('relative', 0.05)
        assert mixed.relative_floor == RELATIVE_FLOOR
        assert manual.relative_floor is None
    # A list works like a tuple, and the floor applies to the relative entry.
    floored = identifiability(decay_ls,
                              {'[S]': 0.5, '[P]': ['relative', 0.05]},
                              relative_floor=0.5)
    assert floored.sigma == {'[S]': 0.5, '[P]': ('relative', 0.05)}
    assert floored.relative_floor == 0.5
    assert not np.array_equal(floored.sensitivity, mixed.sensitivity)
    # The record round-trips exactly through the CSV.
    floored.to_csv(tmp_path / 'mixed.csv')
    back = IdentifiabilityResult.from_csv(tmp_path / 'mixed.csv')
    _assert_same(floored, back)
    # A tag other than 'relative' is refused with the accepted forms.
    for bad in (('rel', 0.05), ('relative',), ('relative', 0.05, 1)):
        with pytest.raises(TypeError, match=r"sigma\['\[P\]'\] must be .*"
                                            r"\('relative', rel\)"):
            identifiability(decay_ls, {'[S]': 0.5, '[P]': bad})
    for bad in (0, -0.1, math.inf):
        with pytest.raises(ValueError, match=r"relative sigma of '\[P\]'"):
            identifiability(decay_ls, {'[P]': ('relative', bad)})


def test_a_sigma_that_is_not_finite_names_the_output(decay_ls):
    # An overflowing or underflowing sigma, or a NaN nominal output, is
    # reported on the output (without RuntimeWarnings), not blamed on the
    # parameters' sensitivities.
    nominal = decay_ls.nominal.copy()
    nominal[3, 1] = np.nan
    broken = dataclasses.replace(decay_ls, nominal=nominal)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        for rel in (1e308, 5e-324):
            with pytest.raises(ValueError, match=r"relative sigma of "
                                                 r"'\[S\]' is not finite"):
                identifiability(decay_ls, ('relative', rel))
        with pytest.raises(ValueError, match=r"nominal values of '\[P\]' "
                                             r"are not finite"):
            identifiability(broken, ('relative', 0.1))
        alone = identifiability(broken, {'[S]': ('relative', 0.1)})
        assert alone.numerical_rank == 2
        with pytest.raises(ValueError, match=r"sigma of '\[S\]' is too "
                                             r"small"):
            identifiability(decay_ls, {'[S]': 1e-310})


def test_relative_sigma_after_step_refinement_keeps_every_row():
    # The first firing lies 3e-5 (relative) after t = 1.4: k is refined
    # twice to a step of 1e-5.
    T = 1.4 * (1 + 3e-5)
    k = math.log(2) / T
    km = nsk.KineticModel.from_antimony(f"""
    model ev()
      S = 10; n = 0; k = {k!r};
      J: S => ; k*S;
      E1: at S < 5: S = S + 5, n = n + 1;
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=6, n_steps=61), ['k'],
                           outputs=['[S]'])
    flags = ls.flags.loc['k']
    assert flags['n_refinements'] == 2 and flags['reliable']
    assert flags['scheme'] == 'central'
    assert flags['step'] == pytest.approx(1e-5, rel=1e-12)
    # 3 rtol / h = 3e-5 > rank_rtol, but a single parameter has no
    # compensating combination to hide: no warning. d ln S / d ln k = -k t
    # exactly between firings.
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ident = identifiability(ls, ('relative', 0.05))
    assert ident.row_outputs == ['[S]'] * 61
    t = ls.time
    assert ident.se['k'] == pytest.approx(
        0.05 / math.sqrt(np.sum((k * t) ** 2)), rel=1e-3)


def test_refined_steps_are_reported_on_a_near_null_direction():
    # The refinement toy with k = k1 + k2: both columns are refined to a
    # step of 1e-5, so integration noise sits at 3 rtol / h = 3e-5 of
    # s_max, above rank_rtol. s_min / s_max (2e-9) is below rank_rtol:
    # rank 1, the conservative verdict, so no doubt and no warning.
    T = 1.4 * (1 + 3e-5)
    K = math.log(2) / T
    km = nsk.KineticModel.from_antimony(f"""
    model ev()
      S = 10; n = 0; k1 = {0.4 * K!r}; k2 = {0.6 * K!r};
      J: S => ; (k1 + k2)*S;
      E1: at S < 5: S = S + 5, n = n + 1;
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=6, n_steps=61),
                           ['k1', 'k2'], outputs=['[S]'])
    assert list(ls.flags['n_refinements']) == [2, 2]
    assert ls.flags['reliable'].all() and ls.flags['scheme'].eq(
        'central').all()
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ident = identifiability(ls, ('relative', 0.05))
    assert ident.numerical_rank == 1 and ident.rank_test_doubtful is False
    assert ident.fd_error_level == pytest.approx(3e-5, rel=1e-5)
    # Below 2e-9, rank_rtol keeps that value: one warning, with the rtol
    # that would bring the noise below rank_rtol (about rank_rtol h / 3,
    # not the default 1e-10 already in use).
    with pytest.warns(UserWarning, match='error level') as record:
        identifiability(ls, ('relative', 0.05), rank_rtol=1e-10)
    message = str(record[0].message)
    assert len(record) == 1 and '3e-05' in message
    assert 'tolerances=(1e-12, rtol), rtol below 2e-16' in message


def test_well_posed_problems_raise_no_error_level_warning():
    # The error level exceeds rank_rtol at a large step, with forward
    # differences or at a loose rtol, but s_min / s_max sits several times
    # above FD_ERROR_MARGIN times the level (about 6 x to 120 x here): no
    # warning. (At CVODE's own rtol 1e-6, the decay toy sits right at the
    # gate, as almost every multi-parameter problem does.)
    km = nsk.KineticModel.from_antimony("""
    model d()
      S = 10; P = 0; k = 0.5; k2 = 0.1;
      J: S => P; k*S;
      J2: P => ; k2*P;
    end
    """, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=61)
    sigma = {'[S]': 0.1, '[P]': 0.1}
    for kwargs in (dict(rel_step=1e-2), dict(method='forward')):
        ls = local_sensitivity(km, cond, ['k', 'k2'], outputs=['[S]', '[P]'],
                               **kwargs)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            ident = identifiability(ls, sigma)
        assert ident.numerical_rank == 2
        assert ident.fd_error_level > ident.rank_rtol
        assert ident.rank_test_doubtful is False
    decay = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    for kwargs in (dict(method='forward'), dict(tolerances=(1e-12, 1e-7))):
        ls = local_sensitivity(decay, nsk.Condition(t_end=10, n_steps=11),
                               ['k', '[S]'], outputs=['[S]', '[P]'],
                               **kwargs)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            ident = identifiability(ls, sigma)
        assert ident.rank_test_doubtful is False


def test_a_lowered_relative_floor_keeps_the_error_level_warning():
    # Under a relative_floor below the default, near-zero rows amplify the
    # finite-difference error beyond its estimate: the sum toy reads full
    # rank at s_min / s_max ~3e3 x the estimate, and the warning stays.
    km = nsk.KineticModel.from_antimony(SUM_TOY, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=60, n_steps=121),
                           ['k1', 'k2'], outputs=['[S]', '[P]'],
                           rel_step=2e-3)
    with pytest.warns(UserWarning, match='error level') as record:
        ident = identifiability(ls, ('relative', 0.1), relative_floor=1e-9)
    assert ident.numerical_rank == 2
    assert 'relative_floor=1e-09, below the default 0.001' in str(
        record[0].message)
    assert ident.rank_test_doubtful is True
    # At the defaults of local_sensitivity the unscaled level (~6e-7) is
    # below rank_rtol, yet precise data (rel = 0.01) at a floor of 1e-7
    # read both k1 and k2 identifiable at s_min / s_max ~8e-4: the level
    # scaled by the realized weight amplification (1e-3 / 1e-7, at the t0
    # row of P, which is 0) flags it.
    km = nsk.KineticModel.from_antimony(DEPLETED_TOY, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=30, n_steps=61),
                           ['k1', 'k2'], outputs=['[S]', '[P]'])
    with pytest.warns(UserWarning, match='error level') as record:
        precise = identifiability(ls, ('relative', 0.01), relative_floor=1e-7)
    assert len(record) == 1
    message = str(record[0].message)
    assert 'relative_floor=1e-07, below the default 0.001' in message
    assert 'Raise relative_floor to the default 0.001' in message
    level = _fd.error_level(ls.flags, ls.tolerances).total
    assert level < precise.rank_rtol < precise.fd_error_level
    assert precise.fd_error_level == pytest.approx(level * 1e4, rel=1e-12)
    assert precise.numerical_rank == 2
    assert precise.table['identifiable'].all()
    assert precise.rank_test_doubtful is True
    # No false alarm: at the default floor the same data raise no doubt,
    # and a well-posed problem stays clear at a floor of 1e-4.
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        default = identifiability(ls, ('relative', 0.01))
    assert default.numerical_rank == 1
    assert default.rank_test_doubtful is False
    km = nsk.KineticModel.from_antimony("""
    model d()
      S = 10; P = 0; k = 0.5; k2 = 0.1;
      J: S => P; k*S;
      J2: P => ; k2*P;
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=6, n_steps=61),
                           ['k', 'k2'], outputs=['[S]', '[P]'])
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        clear = identifiability(ls, ('relative', 0.1), relative_floor=1e-4)
    assert clear.numerical_rank == 2 and clear.rank_test_doubtful is False


def test_a_parameter_at_zero_makes_the_others_depend_on_its_unit():
    # P0 = 0 gets a 'linear' column per mol/L of P0: at S0 = 1e-4 mol/L it
    # dominates the matrix, k1 and kd lose their finite se and drop out of
    # best_subset. A warning names it; without it the others are judged
    # unit-invariantly, with the se of S0 = 10.
    cond = nsk.Condition(t_end=10, n_steps=41)
    se = {}
    for S0 in (10.0, 1e-4):
        km = nsk.KineticModel.from_antimony(f"""
        model s()
          S = {S0!r}; P = 0; k1 = 0.5; kd = 0.2;
          J: S => P; k1*S;
          J2: P => ; kd*P;
        end
        """, units={'time': 'h', 'conc': 'mol/L'})
        ls = local_sensitivity(km, cond, ['k1', '[P]', 'kd'],
                               outputs=['[S]', '[P]'],
                               abs_step={'[P]': 1e-3 * S0})
        with pytest.warns(UserWarning, match=r"\['\[P\]'\] are at p <= 0"
                          ) as record:
            full = identifiability(ls, ('relative', 0.1))
        message = str(record[0].message)
        assert len(record) == 1 and full.space[1] == 'linear'
        assert 'per model unit of p' in message and 'params=' in message
        assert 'relative_floor' in message           # the t0 weight
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            sub = identifiability(ls, ('relative', 0.1), params=['k1', 'kd'])
        assert sub.table['identifiable'].all()
        se[S0] = sub.se.to_numpy()
    assert full.numerical_rank < 3
    assert not full.table.loc[['k1', 'kd'], 'identifiable'].any()
    assert np.allclose(se[10.0], se[1e-4], rtol=1e-4)
    with pytest.warns(UserWarning, match='are at p <= 0') as record:
        identifiability(ls, {'[S]': 1e-5, '[P]': 1e-5})
    assert 'relative_floor' not in str(record[0].message)   # absolute


def test_a_linear_parameter_at_zero_passes_the_precision_check():
    # [P] at 0 has no relative scale (rel_se NaN): it is judged on its
    # finite se alone, in table['identifiable'] and in best_subset.
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=10, n_steps=11),
                           ['k', '[P]'], outputs=['[S]', '[P]'],
                           abs_step={'[P]': 1e-3}, tolerances=TIGHT)
    with pytest.warns(UserWarning, match='are at p <= 0'):
        ident = identifiability(ls, {'[S]': 0.5, '[P]': 0.5})
    table = ident.table
    assert table.loc['[P]', 'space'] == 'linear'
    assert np.isnan(table.loc['[P]', 'rel_se'])
    assert np.isfinite(table.loc['[P]', 'se'])
    assert table.loc['[P]', 'identifiable'] is np.True_
    assert '[P]' in ident.best_subset()
    assert repr(ident).startswith("IdentifiabilityResult(parameters=['k', "
                                  "'[P]'], numerical_rank=2/2")
    assert 'unreliable=0' in repr(ident)


def test_identifiability_is_independent_of_the_concentration_scale():
    # RoadRunner scales a scalar atol by each state's start value, so a
    # decay in mol/L from 1 mM is as accurate as one from 10: the relative
    # model once raised 'No measured rows', and an absolute sigma scaled by
    # S0 warned of a finite-difference error level of 1e-5 of s_max.
    se = {}
    for S0 in (10.0, 1e-3):
        km = nsk.KineticModel.from_antimony(f"""
        model d()
          S = {S0!r}; P = 0; k = 0.5;
          J: S => P; k*S;
        end
        """, units={'time': 'h', 'conc': 'mol/L'})
        ls = local_sensitivity(km, nsk.Condition(t_end=6, n_steps=61),
                               ['k', '[S]'], outputs=['[S]'])
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            relative = identifiability(ls, ('relative', 0.1))
            absolute = identifiability(ls, {'[S]': 0.1 * S0})
        assert relative.row_outputs == ['[S]'] * 61
        se[S0] = relative.se.to_numpy(), absolute.se.to_numpy()
    assert np.allclose(se[10.0][0], se[1e-3][0], rtol=1e-4)
    assert np.allclose(se[10.0][1], se[1e-3][1], rtol=1e-4)


def _record(fn):
    """``(result, messages)`` of ``fn()``, with every warning recorded."""
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter('always')
        out = fn()
    return out, [str(w.message) for w in record]


def test_rows_of_a_product_are_kept_at_any_concentration_scale():
    # [P] starts at 0. A noise-based row filter once dropped every [P] row
    # in mol/L, although they are accurate to ~1e-5: k2, which acts only
    # through P, came out non-identifiable. Under the bounded relative
    # model every row counts, and nothing is warned about, at any scale
    # and with or without the default atol scaling.
    sigmas = {'all': lambda S0: ('relative', 0.1),
              'P': lambda S0: {'[P]': ('relative', 0.1)},
              'absolute': lambda S0: {'[S]': 0.01 * S0, '[P]': 0.01 * S0}}
    found = {}
    for S0 in (10.0, 1e-3):
        km = nsk.KineticModel.from_antimony(f"""
        model t()
          S = {S0!r}; P = 0; k = 0.5; k2 = 0.2;
          J1: S => P; k*S;
          J2: P => ; k2*P;
        end
        """, units={'time': 'h', 'conc': 'mol/L'})
        for tolerances in (None, TIGHT):
            ls = local_sensitivity(km, nsk.Condition(t_end=10, n_steps=21),
                                   ['k', 'k2', '[S]'], outputs=['[S]', '[P]'],
                                   tolerances=tolerances)
            for key, sigma in sigmas.items():
                ident, messages = _record(
                    lambda: identifiability(ls, sigma(S0)))
                assert not messages
                assert ident.numerical_rank == 3
                assert ident.table['identifiable'].all()
                found[S0, tolerances, key] = (len(ident.row_outputs),
                                              ident.se.to_numpy())
    for tolerances in (None, TIGHT):
        for key in sigmas:
            n, se = found[10.0, tolerances, key]
            n_molar, se_molar = found[1e-3, tolerances, key]
            assert n == n_molar
            assert np.allclose(se, se_molar, rtol=1e-3)
    assert found[10.0, None, 'P'][0] == 21     # P = 0 at t0 counts too


CHAIN_TOY = """
model chain()
  A = {A0!r}; B = 0; C = 0; k1 = 1; k2 = 0.25;
  J1: A => B; k1*A;
  J2: B => C; k2*B;
end
"""


def test_outputs_starting_at_zero_are_unit_invariant():
    # B and C start at 0: RoadRunner integrates them to a fixed atol. The
    # default tolerances scale atol to their peak, so their log
    # sensitivities, and the identifiability under a relative sigma, are
    # the same at the micro- or nanomolar scale, without any warning.
    scaled, found = {}, {}
    for A0, conc in ((10.0, 'g/L'), (1e-3, 'mol/L'), (1e-8, 'mol/L')):
        km = nsk.KineticModel.from_antimony(CHAIN_TOY.format(A0=A0),
                                            units={'time': 'h', 'conc': conc})
        ls = local_sensitivity(km, nsk.Condition(t_end=10, n_steps=21),
                               ['k1', 'k2'], outputs=['[B]', '[C]'])
        peak = min(1.0, ls.nominal[:, 0].max(), ls.nominal[:, 1].max())
        assert ls.tolerances[0] == pytest.approx(1e-12 * peak, rel=1e-6)
        assert np.isfinite(ls.scaled[1:]).all()
        scaled[A0] = ls.scaled[1:]
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            ident = identifiability(ls, ('relative', 0.1))
        found[A0] = len(ident.row_outputs), ident.numerical_rank, ident.se
    for A0 in (1e-3, 1e-8):
        assert np.allclose(scaled[A0], scaled[10.0], rtol=0, atol=1e-6)
        assert found[A0][:2] == found[10.0][:2] == (42, 2)
        assert np.allclose(found[A0][2], found[10.0][2], rtol=1e-4)


def test_rows_of_a_trace_by_product_are_kept():
    # Q is a by-product at mg/L levels in a g/L model: a noise-based row
    # filter once dropped its 60 rows (below atol / rtol = 0.01), and
    # se[kq] came out 129 instead of ~0.013 (an absolute sigma gives
    # 0.0135). Every row counts.
    km = nsk.KineticModel.from_antimony("""
    model d()
      S = 10; P = 0; Q = 0; k = 0.5; k2 = 0.3; kq = 1e-4;
      J: S => P; k*S;
      J2: S => Q; kq*S;
      J3: P => ; k2*P;
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=6, n_steps=61),
                           ['k', 'k2', 'kq'], outputs=['[S]', '[P]', '[Q]'])
    assert ls.nominal[:, 2].max() < 2e-3
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ident = identifiability(ls, ('relative', 0.1))
    assert ident.row_outputs.count('[Q]') == 61
    assert ident.numerical_rank == 3
    assert 0 < ident.se['kq'] < 0.1


def test_identifiable_is_the_practical_verdict():
    # Michaelis-Menten in its linear regime (S << Km): only Vmax/Km is well
    # determined. Neither parameter is on an exact null direction, but each
    # has a relative SE of ~12, above REL_SE_MAX = 1: not identifiable in
    # practice.
    km = nsk.KineticModel.from_antimony("""
    model mm()
      S = 1; P = 0; Vmax = 1; Km = 100;
      J: S => P; Vmax*S/(Km + S);
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=50, n_steps=51),
                           ['Vmax', 'Km'], outputs=['S', 'P'])
    ident = identifiability(ls, {'S': 0.01, 'P': 0.01})
    assert ident.numerical_rank == 2
    table = ident.table
    assert np.all(np.isfinite(ident.se)) and np.all(ident.se > 1)
    assert np.allclose(table['rel_se'], ident.se)        # log space
    assert not table['null_direction'].any()
    assert not table['identifiable'].any()
    assert ident.best_subset() == ['Vmax']


def test_best_subset_excludes_a_parameter_at_the_noise_level():
    # Ki = 1e15: the inhibition term changes the rate by 1e-15, so the Ki
    # column is finite-difference noise (delta_msqr ~ 7e-8 against 39 for
    # k). Normalized, its noise direction is nearly orthogonal to k
    # (gamma(k, Ki) ~ 2.5), so the collinearity cap alone kept it and
    # dropped the well-determined Km (gamma(k, Km) ~ 13.8).
    km = nsk.KineticModel.from_antimony("""
    model mmi()
      S = 10; P = 0; I = 1; k = 0.5; Ki = 1e15; Km = 2;
      J1: S => P; k*S/(Km + S)/(1 + I/Ki);
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=20, n_steps=21),
                           ['k', 'Km', 'Ki'], outputs=['[S]', '[P]'])
    ident = identifiability(ls, {'[S]': 0.1, '[P]': 0.1})
    table = ident.table
    assert list(table['identifiable']) == [True, True, False]
    assert ident.gamma(['k', 'Ki']) < 10 < ident.gamma(['k', 'Km']) < 20
    for gamma_max in (10, 20, math.inf):
        assert 'Ki' not in ident.best_subset(gamma_max)
    # Km is left out at 10 by its collinearity with k, not because of Ki.
    assert ident.best_subset(10) == ['k']
    assert ident.best_subset(20) == ['k', 'Km']
    assert ident.best_subset(math.inf) == ['k', 'Km']


def test_times_and_params_select_rows_and_columns(decay_ls):
    ident = identifiability(decay_ls, sigma={'[S]': 0.5}, times=[6, 2, 4],
                            params=['[S]'])
    assert ident.parameters == ['[S]']
    assert np.array_equal(ident.row_times, [2.0, 4.0, 6.0])
    # d [S] / d ln S0 = S, divided by sigma.
    assert np.allclose(ident.sensitivity[:, 0],
                       10 * np.exp(-0.5 * np.array([2, 4, 6])) / 0.5,
                       rtol=1e-5)
    with pytest.raises(ValueError, match='not an output time'):
        identifiability(decay_ls, sigma={'[S]': 0.5}, times=[2.5])
    with pytest.raises(ValueError, match='more than once'):
        identifiability(decay_ls, sigma={'[S]': 0.5}, times=[2, 2])


def test_brun_measures_and_greedy_search_on_a_known_matrix():
    # Columns: a and b nearly parallel, c orthogonal, z without effect.
    S = np.array([[1.0, 1.0, 0.0, 0.0],
                  [0.0, 0.01, 0.0, 0.0],
                  [0.0, 0.0, 2.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0]])
    ident = IdentifiabilityResult(
        parameters=['a', 'b', 'c', 'z'], values=[1.0, 2.0, 3.0, 4.0],
        space=['log'] * 4, sensitivity=S, row_outputs=['y'] * 4,
        row_times=[0.0, 1.0, 2.0, 3.0], reliable=[True] * 4)
    assert np.allclose(ident.delta_msqr,
                       [0.5, math.sqrt(1.0001) / 2, 1.0, 0.0])
    cos = 1 / math.sqrt(1.0001)
    assert ident.gamma(['a', 'b']) == pytest.approx(1 / math.sqrt(1 - cos),
                                                    rel=1e-9)
    assert ident.gamma(['a', 'c']) == pytest.approx(1.0)
    assert math.isinf(ident.gamma(['c', 'z']))
    assert math.isinf(ident.gamma(['z']))
    assert ident.numerical_rank == 3 and math.isinf(ident.condition_number)
    assert math.isinf(ident.se['z']) and np.all(np.isfinite(ident.se[:3]))
    full = ident.collinearity(2)
    assert full.attrs['search'] == 'exhaustive' and len(full) == 6
    # rho of (a, c): det(diag(1, 4) / 4)^(1/4).
    rho = full.set_index('parameters')['rho']
    assert rho[('a', 'c')] == pytest.approx(0.25 ** 0.25)
    assert rho[('a', 'z')] == 0.0
    # gamma of the pairs: the three with z are infinite, then (a, b).
    assert np.isinf(full['gamma'].iloc[:3]).all()
    assert full['parameters'].iloc[3] == ('a', 'b')
    greedy = ident.collinearity(2, max_subsets=5)
    assert greedy.attrs['search'] == 'greedy'
    assert set(greedy['parameters']) <= set(full['parameters'])
    assert ('a', 'b') in set(greedy['parameters'])
    assert ident.collinearity(3)['gamma'].notna().all()
    # The subset skips z (no effect) and one of the near-parallel pair.
    assert ident.best_subset(gamma_max=10) == ['c', 'b']
    # Without a collinearity cap a and b still do not go together: within
    # the subset their relative standard errors are ~100 (> REL_SE_MAX).
    assert ident.best_subset(gamma_max=math.inf) == ['c', 'b']
    # More parameters than rows: rank-deficient, padded singular values.
    wide = IdentifiabilityResult(
        parameters=['a', 'b'], values=[1.0, 1.0], space=['log', 'linear'],
        sensitivity=[[1.0, 2.0]], row_outputs=['y'], row_times=[1.0],
        reliable=[True, False])
    assert wide.numerical_rank == 1 and wide.singular_values.shape == (2,)
    assert wide.singular_values[1] == 0 and np.isinf(wide.se).all()


def test_result_constructor_validates_its_inputs():
    good = dict(parameters=['a', 'b'], values=[1.0, 2.0],
                space=['log', 'log'], sensitivity=np.eye(2),
                row_outputs=['y', 'y'], row_times=[0.0, 1.0],
                reliable=[True, True])
    IdentifiabilityResult(**good)
    for change, match in ((dict(sensitivity=np.eye(3)), 'shapes'),
                          (dict(space=['log', 'ln']), 'space'),
                          (dict(parameters=['a', 'a']), 'Duplicate'),
                          (dict(rank_rtol=2.0), 'rank_rtol'),
                          (dict(relative_floor=0.0), 'relative_floor'),
                          (dict(relative_floor=1.5), 'relative_floor'),
                          (dict(fd_error_level=-1e-9), 'fd_error_level'),
                          (dict(fd_error_level=math.nan), 'fd_error_level'),
                          (dict(sensitivity=[[1.0, np.nan], [0.0, 1.0]]),
                           'finite'),
                          # ln p needs p > 0; a value must be finite
                          (dict(values=[-1.0, 1.0]), r'p > 0'),
                          (dict(values=[0.0, 1.0]), r'p > 0'),
                          (dict(values=[math.nan, 1.0]), 'finite'),
                          (dict(values=[math.inf, 1.0],
                                space=['linear', 'log']), 'finite')):
        with pytest.raises(ValueError, match=match):
            IdentifiabilityResult(**{**good, **change})
    # A 'linear' column takes any finite value; with an se of 1000 a
    # parameter is never identifiable (rel_se uses |p|).
    for value in (0.0, -1.0):
        result = IdentifiabilityResult(**{**good, 'values': [value, 1.0],
                                          'space': ['linear', 'log']})
        assert result.values[0] == value
    wide = IdentifiabilityResult(**{**good, 'values': [-1.0, 1.0],
                                    'space': ['linear', 'log'],
                                    'sensitivity': np.diag([1e-3, 1.0])})
    assert wide.se['a'] == pytest.approx(1e3)
    assert wide.table.loc['a', 'rel_se'] == pytest.approx(1e3)
    assert not wide.table.loc['a', 'identifiable']
    assert wide.best_subset() == ['b']
    with pytest.raises(TypeError, match='fd_error_level'):
        IdentifiabilityResult(**good, fd_error_level='1e-7')
    # A matrix built directly has no known error level.
    result = IdentifiabilityResult(**good)
    assert result.fd_error_level is None
    assert result.rank_test_doubtful is None
    assert 'rank_test' not in repr(result)
    # The result keeps copies: changing the caller's arrays afterwards
    # leaves it (and its derived attributes) consistent.
    arrays = dict(sensitivity=np.eye(2), values=np.array([1.0, 2.0]),
                  row_times=np.array([0.0, 1.0]),
                  reliable=np.array([True, True]))
    result = IdentifiabilityResult(**{**good, **arrays})
    arrays['sensitivity'][0, 0] = 5.0
    arrays['values'][0] = 7.0
    arrays['row_times'][0] = 3.0
    arrays['reliable'][0] = False
    assert result.sensitivity[0, 0] == 1.0 == result.fim.iloc[0, 0]
    assert result.values[0] == 1.0 and result.row_times[0] == 0.0
    assert result.reliable[0]


def test_svd_core_matches_numpy_on_a_full_rank_jacobian():
    rng = np.random.default_rng(0)
    J = rng.normal(size=(30, 4))
    svd = _svd_analysis(J)
    assert svd.numerical_rank == 4
    assert np.allclose(svd.covariance, np.linalg.inv(J.T @ J))
    assert np.allclose(svd.se, np.sqrt(np.diag(np.linalg.inv(J.T @ J))))
    assert svd.condition_number == pytest.approx(np.linalg.cond(J))
    assert np.allclose(svd.null_loading, 0)


def test_svd_core_keeps_estimable_parameters_off_a_noisy_null_direction():
    # Columns a and a + 1e-5 (e + 2 b) are collinear up to "noise" that
    # partly lies along b: a null direction that is not bit-exact and spills
    # ~1e-5 onto b. b is estimable and must keep the SE it has without the
    # collinear pair: pinv(M^T G M)_bb = (G^-1)_bb.
    rng = np.random.default_rng(1)
    a, e, b = rng.normal(size=(3, 50))
    J = np.column_stack([a, a + 1e-5 * (e + 2 * b), b])
    svd = _svd_analysis(J, rank_rtol=1e-4)
    assert svd.numerical_rank == 2
    assert np.isinf(svd.se[0]) and np.isinf(svd.se[1])
    assert np.isfinite(svd.se[2])
    reference = _svd_analysis(np.column_stack([a, b]))
    assert svd.se[2] == pytest.approx(reference.se[1], rel=1e-3)
    # the spill-over onto b is what a 1e-6 cutoff used to call "null"
    assert 1e-6 < svd.null_loading[2] < NULL_LOADING


SUM_TOY = """
model s()
  S = 10; P = 0; k1 = 0.2; k2 = 0.3;
  J: S => P; (k1 + k2)*S;
end
"""


def test_sum_rate_is_non_identifiable_at_default_settings(tmp_path):
    # Only k1 + k2 matters: in ln p the null direction is (k2, -k1) / |k|.
    # A sum gives no bit-exact cancellation, so the verdict rests on the FD
    # error level; the default tolerances put it far below rank_rtol.
    km = nsk.KineticModel.from_antimony(SUM_TOY, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=61)
    ls = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]', '[P]'])
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ident = identifiability(ls, sigma={'[S]': 0.1, '[P]': 0.1})
    assert ident.numerical_rank == 1
    assert np.isinf(ident.se).all() and not ident.table['identifiable'].any()
    v = ident.null_directions['v2'].to_numpy()
    expected = np.array([0.3, -0.2]) / math.hypot(0.2, 0.3)
    assert np.allclose(np.abs(v), np.abs(expected), atol=1e-6)
    assert v[0] * v[1] < 0
    level = _fd.error_level(ls.flags, ls.tolerances).total
    assert ident.fd_error_level == level < ident.rank_rtol
    assert ident.rank_test_doubtful is False
    assert 'rank_test' not in repr(ident)
    # At CVODE's own tolerances the noise hides it, and a warning says so,
    # with the rtol that brings noise plus truncation below rank_rtol.
    loose = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]', '[P]'],
                              tolerances=CVODE)
    with pytest.warns(UserWarning, match='error level') as record:
        missed = identifiability(loose, sigma={'[S]': 0.1, '[P]': 0.1})
    assert missed.numerical_rank == 2
    assert len(record) == 1
    message = str(record[0].message)
    assert 'tolerances=(1e-12, rtol), rtol below 2e-10' in message
    assert 'atol=' not in message                  # the atol is the default
    assert 'rel_step / 3' in message and 'falls about 9x' in message
    assert '(1 / condition_number at full rank)' in message
    # The noise signature of the step check is named as well.
    assert 'rises about 3x for one lifted by integration noise' in message
    assert 'singular_values[1] / singular_values[0] = ' in message
    assert 'gamma' not in message                  # that caveat is in Notes
    # The result keeps the doubt, through a CSV round trip too.
    assert missed.rank_test_doubtful is True
    assert repr(missed).endswith(', rank_test=doubtful)')
    missed.to_csv(tmp_path / 'missed.csv')
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        back = IdentifiabilityResult.from_csv(tmp_path / 'missed.csv')
    assert back.fd_error_level == missed.fd_error_level > missed.rank_rtol
    assert back.rank_test_doubtful is True


def test_the_rank_doubt_is_judged_on_the_retained_singular_values():
    # A singular value the rank test already counted as zero (se = inf) is
    # the conservative verdict, and finite-difference error cannot have
    # lifted a value it kept ~1e3 x above the error level: no doubt, no
    # warning. Here the product k1*k2 cancels to ~1e-16 at any step.
    km = nsk.KineticModel.from_antimony(PRODUCT_TOY, units=UNITS)
    cond = nsk.Condition(t_end=5, n_steps=21)
    sigma = {'[S]': 0.1, '[P]': 0.1}
    for kwargs in (dict(tolerances=CVODE), dict(method='forward')):
        ls = local_sensitivity(km, cond, ['k1', 'k2', 'k3'],
                               outputs=['[S]', '[P]'], **kwargs)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            ident = identifiability(ls, sigma)
        assert ident.numerical_rank == 2
        assert np.isinf(ident.se[['k1', 'k2']]).all()
        assert ident.fd_error_level > ident.rank_rtol
        assert ident.rank_test_doubtful is False
    # A retained value at the error level still raises the doubt, and the
    # warning quotes it rather than the exact null below it.
    mix = nsk.KineticModel.from_antimony("""
    model mix()
      S = 10; P = 0; k1 = 0.5; k2 = 2; k3 = 0.3; k4 = 0.2;
      J1: S => P; k1*k2*S;
      J2: P => ; (k3 + k4)*P;
    end
    """, units=UNITS)
    ls = local_sensitivity(mix, cond, ['k1', 'k2', 'k3', 'k4'],
                           outputs=['[S]', '[P]'], method='forward')
    with pytest.warns(UserWarning, match='error level') as record:
        ident = identifiability(ls, sigma)
    s = ident.singular_values
    assert len(record) == 1 and ident.numerical_rank == 3
    assert ident.rank_test_doubtful is True
    assert s[-1] / s[0] < 1e-14 < 1e-6 < s[2] / s[0] < 1e-4
    message = str(record[0].message)
    assert (f'singular_values[2] / singular_values[0] = {s[2] / s[0]:.1g}'
            in message)
    # An exactly inert parameter (no inhibitor: Ki has no effect) gives an
    # exact zero, never a doubt; a lifted null beside one still does.
    mmi = nsk.KineticModel.from_antimony("""
    model mmi()
      S = 10; P = 0; I = 0; Vmax = 1; Km = 2; Ki = 0.5;
      J1: S => P; Vmax*S/(Km + S)/(1 + I/Ki);
    end
    """, units=UNITS)
    cond = nsk.Condition(t_end=20, n_steps=41)
    for lkw, ikw in ((dict(method='forward'), {}),
                     ({}, dict(relative_floor=1e-4))):
        ls = local_sensitivity(mmi, cond, ['Vmax', 'Km', 'Ki'],
                               outputs=['[S]', '[P]'], **lkw)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            ident = identifiability(ls, ('relative', 0.1), **ikw)
        assert ident.numerical_rank == 2 and ident.singular_values[2] == 0
        assert np.isinf(ident.se['Ki'])
        assert ident.fd_error_level > ident.rank_rtol
        assert ident.rank_test_doubtful is False
    inert = nsk.KineticModel.from_antimony(
        SUM_TOY.replace('k2 = 0.3;', 'k2 = 0.3; z = 1;'), units=UNITS)
    ls = local_sensitivity(inert, nsk.Condition(t_end=6, n_steps=61),
                           ['k1', 'k2', 'z'], outputs=['[S]', '[P]'],
                           method='forward')
    with pytest.warns(UserWarning, match='error level') as record:
        ident = identifiability(ls, sigma)
    s = ident.singular_values
    assert len(record) == 1 and ident.rank_test_doubtful is True
    assert ident.numerical_rank == 2 and s[2] == 0 < s[1]
    assert (f'singular_values[1] / singular_values[0] = {s[1] / s[0]:.1g}'
            in str(record[0].message))
    # The gate itself, and a matrix with fewer rows than parameters (its
    # zero padding is no lifted null).
    assert _rank_test_doubtful(1e-3, 1e-6, [1.0, 1e-5, 0.0], 2) is True
    assert _rank_test_doubtful(1e-3, 1e-6, [1.0, 0.5, 0.0], 2) is False
    assert _rank_test_doubtful(1e-3, 1e-6, [1.0, 0.0], 1) is False
    assert _rank_test_doubtful(1e-3, 1e-6, [0.0, 0.0], 0) is False
    assert _rank_test_doubtful(1e-7, 1e-6, [1.0, 1e-7], 2) is False
    wide = IdentifiabilityResult(
        parameters=['a', 'b'], values=[1.0, 1.0], space=['log', 'log'],
        sensitivity=[[1.0, 2.0]], row_outputs=['y'], row_times=[1.0],
        reliable=[True, True], fd_error_level=1e-3)
    assert wide.numerical_rank == 1 and wide.rank_test_doubtful is False


def test_a_lowered_floor_that_never_binds_raises_no_doubt():
    # With Km far above S, S stays above 0.8 of its peak: no row reaches
    # a floor of 1e-3, so lower floors give the same matrix, the same
    # error level and the same verdict.
    km = nsk.KineticModel.from_antimony("""
    model mm()
      S = 1; P = 0; Vmax = 1; Km = 50;
      J: S => P; Vmax*S/(Km + S);
    end
    """, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=10, n_steps=21),
                           ['Vmax', 'Km'], outputs=['[S]'])
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        default = identifiability(ls, ('relative', 0.1))
        for floor in (1e-4, 1e-5):
            low = identifiability(ls, ('relative', 0.1),
                                  relative_floor=floor)
            assert np.array_equal(low.sensitivity, default.sensitivity)
            assert low.fd_error_level == default.fd_error_level
            assert low.rank_test_doubtful is False
    assert default.rank_test_doubtful is False


def test_an_error_level_that_overflows_is_unknown():
    # A subnormal abs_step overflows the noise estimate 3 rtol / h: the
    # level is recorded as unknown (with a warning), not passed on as inf.
    km = nsk.KineticModel.from_antimony("""
    model d()
      S = 10; P = 0; k = 0.5; k2 = 0.2;
      J: S => P; k*S;
      J2: P => ; k2*P;
    end
    """, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=13)
    with pytest.warns(UserWarning, match='use a step') as record:
        ls = local_sensitivity(km, cond, ['k', 'k2'], outputs=['[S]', '[P]'],
                               abs_step={'k': 1e-320})
    assert all(issubclass(w.category, UserWarning) for w in record)
    with pytest.warns(UserWarning, match='could not be estimated') as record:
        ident = identifiability(ls, {'[S]': 0.1, '[P]': 0.1})
    assert len(record) == 1
    assert ident.fd_error_level is None and ident.rank_test_doubtful is None


def test_truncation_error_of_the_differences_is_reported():
    # A sum k1 + k2 does not cancel the truncation error of a difference:
    # forward differences (first order) lift the null singular value to
    # ~3e-5 of the largest at h = 1e-3, central ones (second order) to
    # ~8e-5 at h = 0.1. Both look full rank, and both are warned about,
    # with or without known tolerances.
    km = nsk.KineticModel.from_antimony(SUM_TOY, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=61)
    sigma = {'[S]': 0.1, '[P]': 0.1}
    forward = local_sensitivity(km, cond, ['k1', 'k2'],
                                outputs=['[S]', '[P]'], method='forward')
    assert list(forward.flags['scheme']) == ['forward', 'forward']
    one_sided, central, names = _fd.truncation_error(forward.flags)
    assert one_sided == pytest.approx(0.5e-3) and central == 0.0
    assert names == ['k1', 'k2']
    for sens in (forward, dataclasses.replace(forward, tolerances=None)):
        with pytest.warns(UserWarning, match='truncation error') as record:
            ident = identifiability(sens, sigma)
        assert len(record) == 1
        message = str(record[0].message)
        assert ("one-sided differences of ['k1', 'k2']" in message
                and "method='central'" in message)
        assert ident.numerical_rank == 2
        assert np.isfinite(ident.se).all()
    s = ident.singular_values
    assert 1e-5 < s[1] / s[0] < 1e-4
    # The same on the round-3 sum toy under a relative sigma.
    km2 = nsk.KineticModel.from_antimony("""
    model s()
      S = 10; P = 0; k1 = 0.5; k2 = 1;
      J: S => P; (k1 + k2)*S;
    end
    """, units=UNITS)
    ls2 = local_sensitivity(km2, nsk.Condition(t_end=10, n_steps=61),
                            ['k1', 'k2'], outputs=['[S]', '[P]'],
                            method='forward')
    with pytest.warns(UserWarning, match='first-order accurate'):
        identifiability(ls2, ('relative', 0.05))
    # Central differences: 0.3 h^2 is 3e-7 at the default step (silent,
    # see test_sum_rate_is_non_identifiable_at_default_settings), 2.7e-6
    # at 3e-3 and 3e-3 at 0.1. At 3e-3 the toy above reaches full rank
    # (2.4e-6) under a relative sigma and is warned about; this one keeps
    # rank 1 there, the conservative verdict, with no warning.
    ls = local_sensitivity(km2, nsk.Condition(t_end=10, n_steps=61),
                           ['k1', 'k2'], outputs=['[S]', '[P]'],
                           rel_step=3e-3)
    with pytest.warns(UserWarning, match='Central differences at a relative '
                                         'step of 0.003') as record:
        ident = identifiability(ls, ('relative', 0.05))
    assert 'rel_step at most 1e-3' in str(record[0].message)
    assert ident.numerical_rank == 2
    ls = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]', '[P]'],
                           rel_step=3e-3)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ident = identifiability(ls, sigma)
    assert ident.numerical_rank == 1 and ident.rank_test_doubtful is False
    assert ident.fd_error_level > ident.rank_rtol
    ls = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]', '[P]'],
                           rel_step=0.1)
    with pytest.warns(UserWarning, match='Central differences at a relative '
                                         'step of 0.1') as record:
        ident = identifiability(ls, sigma)
    assert 'rel_step at most 1e-3' in str(record[0].message)
    assert ident.numerical_rank == 2
    assert ls.flags['scheme'].eq('central').all()
    one_sided, central, names = _fd.truncation_error(ls.flags)
    assert one_sided == 0.0 and not names
    assert central == pytest.approx(0.3 * 0.1 ** 2)
    level = _fd.error_level(ls.flags, ls.tolerances)
    assert level.central == central and level.one_sided == 0.0
    assert level.h_min == pytest.approx(0.1)
    assert level.total == pytest.approx(level.noise + central)
    # A large step set by abs_step (0.05 at k1 = 0.2: a relative step of
    # 0.25) is named, with abs_step as the remedy, not rel_step (already
    # 1e-3).
    ls = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]', '[P]'],
                           abs_step={'k1': 0.05})
    with pytest.warns(UserWarning, match='error level') as record:
        identifiability(ls, sigma)
    message = str(record[0].message)
    assert len(record) == 1
    assert "central differences of ['k1'] use an abs_step" in message
    assert 'a smaller abs_step' in message
    assert 'rel_step at most' not in message
    # Below the default rank_rtol, central differences at the default step
    # are too coarse: no rel_step advice (it is already 1e-3).
    ls = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]', '[P]'])
    with pytest.warns(UserWarning, match='error level') as record:
        identifiability(ls, sigma, rank_rtol=1e-8)
    message = str(record[0].message)
    assert ("central differences of ['k1', 'k2'], at relative steps of at "
            'most 0.001') in message
    assert 'too large for rank_rtol=1e-08' in message
    assert 'rel_step at most' not in message


def test_a_one_sided_fallback_is_reported_with_its_own_remedy():
    # The minus side of a = 1.0005 makes sqrt(a - 1) undefined, so 'a' falls
    # back to a reliable forward difference (first order) under
    # method='central': the warning names it, and leaving it out (not
    # method='central', already in use) is the remedy.
    km = nsk.KineticModel.from_antimony("""
    model f()
      S = 4; P = 0; a = 1.0005; b = 0.3;
      J: S => P; (sqrt(a - 1) + b)*S;
    end
    """, units=UNITS)
    with pytest.warns(UserWarning, match=r"central difference of \['a'\] "
                                         'failed'):
        ls = local_sensitivity(km, nsk.Condition(t_end=4, n_steps=41),
                               ['a', 'b'], outputs=['[S]', '[P]'])
    assert ls.flags.loc['a', 'scheme'] == 'forward'
    assert ls.flags.loc['a', 'reliable']
    with pytest.warns(UserWarning, match='fell back to one side') as record:
        identifiability(ls, {'[S]': 0.1, '[P]': 0.1})
    message = str(record[0].message)
    assert len(record) == 1
    assert "one-sided differences of ['a']" in message
    assert 'leave them out with params=' in message
    assert "method='central'" not in message


# Monod growth: only m1 + m2 matters, so m1 and m2 are exactly
# non-identifiable; S is depleted in a steep front around t = 20.
MONOD_TOY = """
model monod()
  X = 0.01; S = 20; m1 = 0.1; m2 = 0.3; Ks = 0.5;
  J1: S => X; (m1 + m2)*S/(Ks + S)*X;
end
"""


def test_a_depletion_front_needs_the_step_convergence_check():
    # At a depletion front a step in ln p shifts the front in time, so the
    # truncation error is far above the 0.3 h^2 estimate: at the defaults
    # the m1/m2 null reads full rank, identifiable, with gamma below 1e3
    # and no warning (the documented limitation). Recomputing at
    # rel_step / 3 exposes it: s_min / s_max falls about 9x (as h^2),
    # while that of the well-posed pair (m2, Ks) barely changes.
    km = nsk.KineticModel.from_antimony(MONOD_TOY, units=UNITS)
    cond = nsk.Condition(t_end=30, n_steps=301)
    ratio = {}
    for rel_step in (1e-3, 1e-3 / 3):
        ls = local_sensitivity(km, cond, ['m1', 'm2', 'Ks'],
                               outputs=['[X]', '[S]'], rel_step=rel_step)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            for params in (['m1', 'm2', 'Ks'], ['m2', 'Ks']):
                ident = identifiability(ls, ('relative', 0.1),
                                        params=params)
                s = ident.singular_values
                ratio[rel_step, len(params)] = s[-1] / s[0]
                if rel_step == 1e-3 and len(params) == 3:
                    assert ident.numerical_rank == 3
                    assert ident.table['identifiable'].all()
                    assert ident.gamma(['m1', 'm2']) < 1e3
                    # The estimate does not flag the front (the documented
                    # limitation): False here is no proof, and this
                    # full-rank verdict is wrong.
                    assert ident.rank_test_doubtful is False
    assert 6 < ratio[1e-3, 3] / ratio[1e-3 / 3, 3] < 12
    assert 0.9 < ratio[1e-3, 2] / ratio[1e-3 / 3, 2] < 1.1


def test_a_loose_atol_is_warned_about():
    # The error level models rtol only. An explicit atol of 1e-6 lets the
    # depleted rows of S lift the (k1 + k2) null to full rank with a finite
    # se, with no rank warning: a separate warning says so.
    km = nsk.KineticModel.from_antimony("""
    model s()
      S = 10; P = 0; k1 = 0.5; k2 = 1;
      J: S => P; (k1 + k2)*S;
    end
    """, units=UNITS)
    cond = nsk.Condition(t_end=10, n_steps=61)
    ls = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]', '[P]'],
                           tolerances=(1e-6, 1e-10))
    with pytest.warns(UserWarning, match='atol') as record:
        loose = identifiability(ls, ('relative', 0.1))
    message = str(record[0].message)
    assert len(record) == 1 and 'atol=1e-06' in message
    assert 'atol <= 1e-08' in message
    assert 'tolerances=(1e-12, rtol)' in message
    # The estimate does not model this atol: the result records the level,
    # and so the doubt, as unknown, not as "no doubt".
    assert loose.numerical_rank == 2
    assert loose.fd_error_level is None and loose.rank_test_doubtful is None
    # One parameter has no combination to hide: no warning, and never
    # doubtful, though the level is still unknown.
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        one = identifiability(ls, ('relative', 0.1), params=['k1'])
    assert one.fd_error_level is None and one.rank_test_doubtful is False
    # At a loose rtol too, the rank warning's advice does not echo the
    # loose atol.
    ls = local_sensitivity(km, cond, ['k1', 'k2'], outputs=['[S]', '[P]'],
                           tolerances=(1e-6, 1e-6))
    with pytest.warns(UserWarning) as record:
        both = identifiability(ls, ('relative', 0.1))
    messages = [str(w.message) for w in record]
    rank = [m for m in messages if m.startswith('The estimated')]
    assert len(messages) == 2 and len(rank) == 1
    assert both.rank_test_doubtful is None
    assert 'tolerances=(1e-12, rtol), rtol below' in rank[0]
    assert 'atol=1e-06 must come down too' in rank[0]


def test_symmetric_chain_keeps_the_estimable_parameters():
    # C(t) of A -> B -> C is symmetric in (k1, k2) at k1 = k2: (1, -1)/sqrt2
    # is a null direction that finite differences do not cancel exactly.
    # k3 and [A] stay estimable, with the SEs of the analysis without k2.
    km = nsk.KineticModel.from_antimony("""
    model chain()
      A = 10; B = 0; C = 0; k1 = 0.5; k2 = 0.5; k3 = 0.2; n = 0;
      J1: A => B; k1*A;
      J2: B => C; k2*B;
      J3: C => ; k3*C;
      E1: at time > 2: A = A + 5, n = n + 1;
      E2: at time > 4: A = A + 5, n = n + 1;
    end
    """, units=UNITS)
    cond = nsk.Condition(t_end=10, n_steps=101)
    ls = local_sensitivity(km, cond, ['k1', 'k2', 'k3', '[A]'],
                           outputs=['C'])
    ident = identifiability(ls, sigma={'C': 0.1}, rank_rtol=1e-4)
    assert ident.numerical_rank == 3
    table = ident.table
    assert list(table['identifiable']) == [False, False, True, True]
    assert np.isinf(ident.se[['k1', 'k2']]).all()
    ref = identifiability(local_sensitivity(km, cond, ['k1', 'k3', '[A]'],
                                            outputs=['C']),
                          sigma={'C': 0.1})
    assert np.allclose(ident.se[['k3', '[A]']], ref.se[['k3', '[A]']],
                       rtol=1e-2)
    assert ident.se['k3'] == pytest.approx(0.0306, rel=1e-2)
    assert ident.se['[A]'] == pytest.approx(0.0572, rel=1e-2)
    # At CVODE's own tolerances FD noise spills ~1e-5 onto k3 and [A]; they
    # must still count as estimable (and the error level is reported).
    loose = local_sensitivity(km, cond, ['k1', 'k2', 'k3', '[A]'],
                              outputs=['C'], tolerances=CVODE)
    with pytest.warns(UserWarning, match='error level'):
        noisy = identifiability(loose, sigma={'C': 0.1}, rank_rtol=1e-4)
    assert noisy.numerical_rank == 3
    assert (noisy.table['null_loading'][['k3', '[A]']] > 1e-6).all()
    assert list(noisy.table['identifiable']) == [False, False, True, True]
    assert np.allclose(noisy.se[['k3', '[A]']], ref.se[['k3', '[A]']],
                       rtol=1e-2)


def test_unsorted_times_keep_their_per_time_sigma(decay_ls, tmp_path):
    # Rows are sorted by time; a per-time sigma array follows its times.
    a = identifiability(decay_ls, sigma={'[S]': [0.01, 1.0]}, times=[6, 2])
    b = identifiability(decay_ls, sigma={'[S]': [1.0, 0.01]}, times=[2, 6])
    assert list(a.row_times) == [2.0, 6.0]
    assert np.array_equal(a.sensitivity, b.sensitivity)
    pd.testing.assert_series_equal(a.se, b.se, check_exact=True)
    assert a.sigma == b.sigma == {'[S]': [1.0, 0.01]}
    a.to_csv(tmp_path / 'a.csv')
    _assert_same(a, IdentifiabilityResult.from_csv(tmp_path / 'a.csv'))


def test_unreliable_columns_warn_and_nonfinite_columns_raise(decay_ls):
    # At k = 4 ln2 / 6 the 4th top-up lands at t = 6: k is flagged.
    km = nsk.KineticModel.from_antimony(EVENT_TOY, units=UNITS)
    cond = nsk.Condition(t_end=6, n_steps=13,
                         parameters={'k': 4 * math.log(2) / 6})
    ls = local_sensitivity(km, cond, ['k', 'z'], outputs=['[S]'],
                           tolerances=TIGHT)
    assert list(ls.flags['reliable']) == [False, True]
    with pytest.warns(UserWarning, match=r"\['k'\] are unreliable"):
        ident = identifiability(ls, sigma={'[S]': 0.2})
    assert list(ident.reliable) == [False, True]
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        inert = identifiability(ls, sigma={'[S]': 0.2}, params=['z'])
    assert inert.numerical_rank == 0 and np.isinf(inert.se['z'])
    raw = decay_ls.raw.copy()
    raw[:, :, 0] = np.nan                   # a failed perturbation of k
    broken = dataclasses.replace(decay_ls, raw=raw)
    with pytest.raises(ValueError,
                       match=r"Non-finite sensitivities for \['k'\]"):
        identifiability(broken, sigma={'[S]': 0.5})
    assert identifiability(broken, sigma={'[S]': 0.5},
                           params=['[S]']).numerical_rank == 1


def test_a_result_without_output_time_courses_is_refused():
    # identifiability uses output time courses only, not metrics
    from nskinetics.analysis import Final
    km = nsk.KineticModel.from_antimony(DECAY_TOY, units=UNITS)
    ls = local_sensitivity(km, nsk.Condition(t_end=2, n_steps=5), ['k'],
                           outputs=[], metrics={'P': Final('[P]')},
                           tolerances=TIGHT)
    for sigma in (('relative', 0.1), {'P': 0.5}):
        with pytest.raises(ValueError, match='needs output time courses'):
            identifiability(ls, sigma=sigma)


@pytest.mark.parametrize('kwargs, error, match', [
    (dict(sigma={'[X]': 1.0}), ValueError, 'not an output'),
    (dict(sigma={'[S]': 0.0}), ValueError, 'positive'),
    (dict(sigma={'[S]': -1.0}), ValueError, 'positive'),
    (dict(sigma={'[S]': [1.0, 2.0]}), ValueError, 'one value per'),
    (dict(sigma={'[S]': 'a'}), TypeError, 'positive number'),
    (dict(sigma={}), ValueError, 'at least one'),
    (dict(sigma=0.5), TypeError, 'mapping'),
    (dict(sigma=('absolute', 0.5)), TypeError, 'mapping'),
    (dict(sigma=('relative', 0)), ValueError, 'positive'),
    # a scalar rel names the argument, not the first output
    (dict(sigma=('relative', -1)), ValueError,
     r"rel in sigma=\('relative', rel\) must be positive"),
    (dict(sigma=('relative', 'x')), TypeError, 'rel in sigma'),
    # per-output relative sigmas are written as a mapping
    (dict(sigma=('relative', {'[S]': 0.1})), TypeError,
     r"written as a mapping, e\.g\. \{'\[S\]': \('relative', 0\.1\)"),
    (dict(sigma={'[X]': ('relative', 0.1)}), ValueError, 'not an output'),
    (dict(sigma={'[S]': 1.0}, params=['q']), ValueError, 'not parameters'),
    (dict(sigma={'[S]': 1.0}, params=['k', 'k']), ValueError,
     'more than once'),
    (dict(sigma={'[S]': 1.0}, rank_rtol=0), ValueError, 'rank_rtol'),
    # relative_floor is checked whatever the noise model, before any work
    (dict(sigma={'[S]': 1.0}, relative_floor=0), ValueError,
     r'relative_floor must be in \(0, 1\]'),
    (dict(sigma=('relative', 0.1), relative_floor=-1e-3), ValueError,
     'relative_floor'),
    (dict(sigma=('relative', 0.1), relative_floor=1.0001), ValueError,
     'relative_floor'),
    (dict(sigma=('relative', 0.1), relative_floor=math.nan), ValueError,
     'relative_floor must be finite'),
    (dict(sigma=('relative', 0.1), relative_floor=math.inf), ValueError,
     'relative_floor must be finite'),
    (dict(sigma=('relative', 0.1), relative_floor='1e-3'), TypeError,
     'relative_floor must be a real number'),
    (dict(sigma=('relative', 0.1), relative_floor=True), TypeError,
     'relative_floor must be a real number'),
    (dict(sigma=('relative', 0.1), relative_floor=0.0, times=[2.5]),
     ValueError, 'relative_floor'),
    (dict(sigma={'[S]': 1.0}, times='a'), TypeError, 'times'),
])
def test_invalid_arguments(decay_ls, kwargs, error, match):
    with pytest.raises(error, match=match):
        identifiability(decay_ls, **kwargs)


def test_invalid_method_arguments(decay_ls):
    with pytest.raises(TypeError, match='LocalSensitivityResult'):
        identifiability('ls', sigma={'[S]': 1.0})
    ident = identifiability(decay_ls, sigma={'[S]': 0.5})
    for k in (0, 3):
        with pytest.raises(ValueError, match='subset_size'):
            ident.collinearity(k)
    with pytest.raises(TypeError, match='subset_size'):
        ident.collinearity(2.0)
    with pytest.raises(ValueError, match='max_subsets'):
        ident.collinearity(2, max_subsets=0)
    with pytest.raises(ValueError, match='gamma_max'):
        ident.best_subset(gamma_max=0.5)
    with pytest.raises(ValueError, match='Unknown'):
        ident.gamma(['q'])
    with pytest.raises(ValueError, match='kind'):
        ident.plot(kind='pie')


# --- (j) CSV round trip ------------------------------------------------------
def _assert_same(a, b):
    assert type(a) is type(b)
    for name in ('parameters', 'space', 'row_outputs', 'sigma', 'rank_rtol',
                 'relative_floor', 'fd_error_level', 'rank_test_doubtful',
                 'units', 'numerical_rank', 'condition_number'):
        assert getattr(a, name) == getattr(b, name), name
    for name in ('values', 'sensitivity', 'row_times', 'reliable',
                 'singular_values'):
        assert np.array_equal(getattr(a, name), getattr(b, name),
                              equal_nan=True), name
    for name in ('fim', 'directions', 'null_directions', 'covariance',
                 'table', 'correlation'):
        pd.testing.assert_frame_equal(getattr(a, name), getattr(b, name),
                                      check_exact=True)
    for name in ('se', 'delta_msqr'):
        pd.testing.assert_series_equal(getattr(a, name), getattr(b, name),
                                       check_exact=True)


def test_csv_round_trip_is_exact(product_ls, decay_ls, tmp_path):
    ident = identifiability(product_ls, sigma={'[S]': 0.5, '[P]': 0.25})
    path = tmp_path / 'ident.csv'
    ident.to_csv(path)
    back = IdentifiabilityResult.from_csv(path)
    _assert_same(ident, back)
    pd.testing.assert_frame_equal(back.collinearity(2), ident.collinearity(2),
                                  check_exact=True)
    assert back.best_subset() == ident.best_subset()
    back.to_csv(tmp_path / 'again.csv')
    assert (tmp_path / 'again.csv').read_bytes() == path.read_bytes()
    assert back.relative_floor is None
    assert back.fd_error_level == ident.fd_error_level > 0
    assert repr(back) == repr(ident)
    rel = identifiability(decay_ls, sigma=('relative', 0.05), times=[0, 4, 8],
                          relative_floor=0.02)
    rel.to_csv(tmp_path / 'rel.csv')
    back = IdentifiabilityResult.from_csv(tmp_path / 'rel.csv')
    _assert_same(rel, back)
    assert back.relative_floor == rel.relative_floor == 0.02
    # The floor is part of the noise model: it changes the matrix.
    assert not np.array_equal(
        back.sensitivity, identifiability(decay_ls, ('relative', 0.05),
                                          times=[0, 4, 8]).sensitivity)


def test_csv_round_trip_of_directly_built_results(tmp_path):
    # Calibration builds results from a weighted residual Jacobian: no
    # sigma, a Fortran-ordered matrix, zero rows, or a given error level.
    cases = {
        'plain': dict(parameters=['a', 'b'], values=[1.0, -2.0],
                      space=['log', 'linear'],
                      sensitivity=np.asfortranarray([[1.0, 0.3], [0.2, 2.0],
                                                     [0.1, 0.0]]),
                      row_outputs=['y', 'y', 'z'], row_times=[0, 1, 1],
                      reliable=[True, False]),
        'empty': dict(parameters=['a'], values=[1.0], space=['log'],
                      sensitivity=np.zeros((0, 1)), row_outputs=[],
                      row_times=[], reliable=[True]),
        'level': dict(parameters=['a', 'b'], values=[1.0, 2.0],
                      space=['log', 'log'],
                      sensitivity=[[1.0, 1.0], [1.0, 1.0 + 1e-4]],
                      row_outputs=['y', 'y'], row_times=[0.0, 1.0],
                      reliable=[True, True], sigma=('relative', 0.1),
                      relative_floor=1e-4, fd_error_level=3e-6),
    }
    found = {}
    for name, kwargs in cases.items():
        result = IdentifiabilityResult(**kwargs)
        path = tmp_path / f'{name}.csv'
        result.to_csv(path)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            back = found[name] = IdentifiabilityResult.from_csv(path)
        _assert_same(result, back)
        back.to_csv(tmp_path / f'{name}_again.csv')
        assert (tmp_path / f'{name}_again.csv').read_bytes() == \
            path.read_bytes()
    plain = found['plain']
    assert plain.sigma is None and plain.fd_error_level is None
    assert plain.rank_test_doubtful is None
    assert found['empty'].numerical_rank == 0
    assert found['empty'].rank_test_doubtful is False      # one parameter
    level = found['level']
    assert level.numerical_rank == 2 and level.rank_test_doubtful is True
    assert level.relative_floor == 1e-4 and level.fd_error_level == 3e-6


def test_from_csv_rejects_other_files(decay_ls, tmp_path):
    decay_ls.to_csv(tmp_path / 'ls.csv')
    with pytest.raises(ValueError, match='IdentifiabilityResult'):
        IdentifiabilityResult.from_csv(tmp_path / 'ls.csv')
    (tmp_path / 'x.csv').write_text('a,b\n1,2\n')
    with pytest.raises(ValueError, match='IdentifiabilityResult'):
        IdentifiabilityResult.from_csv(tmp_path / 'x.csv')


def test_from_csv_rejects_incomplete_files(decay_ls, tmp_path):
    # A truncated or edited file is refused, never filled with
    # uninitialized memory.
    ident = identifiability(decay_ls, sigma={'[S]': 0.5, '[P]': 0.25})
    path = tmp_path / 'ident.csv'
    ident.to_csv(path)
    lines = path.read_text().splitlines(keepends=True)

    def without(predicate, name):
        kept = [line for line in lines if not predicate(line)]
        assert len(kept) == len(lines) - 1
        (tmp_path / name).write_text(''.join(kept))
        return tmp_path / name

    for quantity in ('sensitivity', 'parameter_value', 'reliable', 'row'):
        first = next(line for line in lines
                     if line.startswith(quantity + ','))
        with pytest.raises(ValueError, match=f'incomplete.*{quantity}'):
            IdentifiabilityResult.from_csv(
                without(lambda line, first=first: line == first,
                        f'{quantity}.csv'))
    with pytest.raises(ValueError, match=r"incomplete.*meta keys "
                                         r"\['relative_floor'\]"):
        IdentifiabilityResult.from_csv(without(
            lambda line: line.startswith('meta,relative_floor,'),
            'meta.csv'))
    # A duplicated entry in place of another is refused too.
    first = next(i for i, line in enumerate(lines)
                 if line.startswith('sensitivity,'))
    edited = lines[:first] + [lines[first + 1]] + lines[first + 1:]
    (tmp_path / 'dup.csv').write_text(''.join(edited))
    with pytest.raises(ValueError, match='incomplete.*duplicate'):
        IdentifiabilityResult.from_csv(tmp_path / 'dup.csv')


# --- (k) plots ---------------------------------------------------------------
def test_plot_smoke(product_ls):
    import matplotlib.pyplot as plt
    ident = identifiability(product_ls, sigma={'[S]': 0.5, '[P]': 0.5})
    fig, ax = ident.plot()
    assert len(ax.patches) == 3
    assert [t.get_text() for t in ax.get_yticklabels()][-1] in ('k1', 'k2')
    assert 'rank 2/3' in ax.get_title()
    fig2, ax2 = ident.plot(kind='singular_values')
    assert len(ax2.lines) == 2 and ax2.get_yscale() == 'log'
    plt.close(fig)
    plt.close(fig2)
    for kind in ('delta_msqr', 'singular_values'):
        fig0, ax0 = plt.subplots()
        assert ident.plot(kind=kind, ax=ax0) == (fig0, ax0)
    plt.close('all')
