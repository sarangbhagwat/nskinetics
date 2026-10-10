# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of ``sweep_utils.py``, the plot helpers shared by the 2-D grid sweeps
(``analyses/kinetic_sweeps/evaluate_EtOH_*.py`` and
``analyses/feeding_strategies/evaluate_feeding_strategies.py``).

The helpers were copied verbatim from those scripts, quirks included, so the
sweep figures do not change; these tests pin that behaviour (tick and level
values, including :func:`get_contour_info_from_metric_data`'s base-5 rounding
of its ``median + 0.2 stdev`` tick). The module is loaded by file path; its
import is probed in a fresh interpreter.
"""
import json
import os

import numpy as np
import pytest

from ._support import PKG_DIR, assert_ran_ok, load_by_path, run_isolated, skip_unless_modules

pytestmark = pytest.mark.isobutanol

SWEEP_UTILS = os.path.join(PKG_DIR, 'sweep_utils.py')
_SENTINEL = 'NSK_MODULES '


@pytest.fixture(scope='module')
def su():
    return load_by_path('sweep_utils', SWEEP_UTILS)


def test_import_is_light():
    """Loading sweep_utils by path imports numpy and matplotlib only: no
    nskinetics / biosteam / thermosteam / biorefinery module (the colormaps
    import biosteam's colours only when called). Probed in a fresh
    interpreter."""
    probe = (
        "import importlib.util, json, sys\n"
        f"spec = importlib.util.spec_from_file_location('sweep_utils', {SWEEP_UTILS!r})\n"
        "m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)\n"
        f"print({_SENTINEL!r} + json.dumps(sorted(sys.modules)))\n")
    proc = run_isolated(['-c', probe], timeout=120, simulation=False)
    out = assert_ran_ok(proc, _SENTINEL)
    line = [ln for ln in out.splitlines() if ln.startswith(_SENTINEL)][-1]
    modules = json.loads(line[len(_SENTINEL):])
    heavy = sorted(m for m in modules if m.split('.')[0] in (
        'nskinetics', 'biosteam', 'thermosteam', 'biorefineries', 'hensmith',
        'tellurium', 'roadrunner'))
    assert heavy == []
    assert 'numpy' in modules and 'matplotlib.colors' in modules


def test_exports(su):
    assert set(su.__all__) == {'JBEI_UCB_colormap', 'CABBI_green_colormap',
                               'tickmarks_from_data', 'tickmarks', 'get_OOM',
                               'my_round', 'get_contour_info_from_metric_data'}
    assert all(callable(getattr(su, name)) for name in su.__all__)


def test_tickmarks(su):
    """Evenly spaced ticks between the bounds rounded outwards to a multiple
    of `accuracy`."""
    assert su.tickmarks(12, 487) == [0.0, 125.0, 250.0, 375.0, 500.0]
    assert np.allclose(su.tickmarks(0.12, 0.87, accuracy=0.25, N_points=4),
                       [0.0, 1/3, 2/3, 1.0], rtol=0, atol=1e-15)
    assert su.tickmarks_from_data(np.array([[3.0, 140.0], [77.0, 251.0]])) == [
        0.0, 75.0, 150.0, 225.0, 300.0]


def test_get_OOM_and_my_round(su):
    """get_OOM = floor(log10(x + 1)); my_round rounds to a multiple of base
    (default 5)."""
    # (values away from the decade boundaries, where log(x, 10) round-off
    # is platform-dependent)
    assert [su.get_OOM(x) for x in (0.0, 0.5, 8.0, 50.0, 97.0, 500.0)] == [0, 0, 0, 1, 1, 2]
    assert [su.my_round(x) for x in (7, 8, 12.4, -3)] == [5, 10, 10, -5]
    assert abs(su.my_round(0.26, 0.05) - 0.25) < 1e-15
    assert su.my_round(13, 2) == 12


def test_contour_info_integer_scale(su):
    """median 30, stdev sqrt(200): levels 20..40 by 1, ticks every 5."""
    w_levels, w_ticks, cbar_ticks = su.get_contour_info_from_metric_data(
        np.array([10., 20., 30., 40., 50.]))
    assert np.array_equal(w_levels, np.arange(20.0, 41.0, 1.0))
    assert w_ticks == [20.0, 40.0]
    assert cbar_ticks.tolist() == [20.0, 25.0, 30.0, 35.0, 40.0]


def test_contour_info_ignores_nan_and_keeps_the_base5_quirk(su):
    """NaNs are dropped (median 3, stdev sqrt(2)); the median + 0.2 stdev
    tick is rounded to base 5 (-> 5.0), unlike its neighbours, and ticks
    within 80 % of their successor are dropped."""
    w_levels, w_ticks, cbar_ticks = su.get_contour_info_from_metric_data(
        np.array([[1.0, 2.0, np.nan], [3.0, 4.0, 5.0]]))
    assert len(w_levels) == 21
    assert np.allclose(w_levels, np.linspace(2.0, 4.0, 21), rtol=0, atol=1e-12)
    assert [float(t) for t in w_ticks] == [2.0, 4.0, 5.0]
    assert cbar_ticks.tolist() == [2.0, 2.5, 3.0, 3.5, 4.0]


def test_contour_info_hundreds(su):
    """A list input (converted to an array); levels every 10, ticks every 50."""
    w_levels, w_ticks, cbar_ticks = su.get_contour_info_from_metric_data(
        [100., 200., 300., 400., 500., 600.])
    assert np.array_equal(w_levels, np.arange(200.0, 501.0, 10.0))
    assert [float(t) for t in w_ticks] == [200.0, 300.0, 500.0]
    assert cbar_ticks.tolist() == [200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0]


@skip_unless_modules('biosteam')
def test_colormaps(su):
    """The JBEI/UCB map runs UCB yellow -> ... -> dark grey (reversible); the
    CABBI green map runs CABBI orange -> ... -> dark grey; N_levels sets the
    number of colours."""
    from biosteam.utils import colors
    jbei = su.JBEI_UCB_colormap()
    assert jbei.N == 90 and jbei.name == 'CABBI'
    assert np.allclose(jbei(0.0)[:3], (253/255, 181/255, 21/255))
    assert np.allclose(jbei(1.0)[:3], colors.grey_dark.RGBn)
    rev = su.JBEI_UCB_colormap(N_levels=10, reverse=True)
    assert rev.N == 10
    assert np.allclose(rev(0.0)[:3], colors.grey_dark.RGBn)
    assert np.allclose(rev(1.0)[:3], (253/255, 181/255, 21/255))
    green = su.CABBI_green_colormap()
    assert green.N == 90
    assert np.allclose(green(0.0)[:3], colors.CABBI_orange.RGBn)
    assert np.allclose(green(1.0)[:3], colors.grey_dark.RGBn)
