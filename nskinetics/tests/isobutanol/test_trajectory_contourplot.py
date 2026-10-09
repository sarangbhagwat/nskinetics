# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of ``plots/trajectory_contourplot.py``.

The module is loaded by file path (inside a fixture), so the biorefinery
package is never imported (no import-time simulation, no numba cache). It
needs the optional ``contourplots`` package; the tests skip without it. Uses
the Agg backend.

Moved from ``plots/test_trajectory_contourplot.py`` (the standalone
``__main__`` runner was dropped; the 14 tests are unchanged).
"""
import os

import numpy as np
import pytest

from ._support import load_by_path, load_paths, skip_unless_modules

pytestmark = [pytest.mark.isobutanol, skip_unless_modules('contourplots')]


@pytest.fixture(scope='module')
def tc():
    import matplotlib
    matplotlib.use('Agg')
    plots_dir = load_paths().PLOTS_DIR
    return load_by_path('trajectory_contourplot',
                        os.path.join(plots_dir, 'trajectory_contourplot.py'))


def test_squeeze_grid_accepts_2d_and_3d(tc):
    a2 = np.arange(6.0).reshape(2, 3)
    assert tc._squeeze_grid(a2).shape == (2, 3)
    a3 = np.arange(6.0).reshape(1, 2, 3)
    assert tc._squeeze_grid(a3).shape == (2, 3)
    assert np.array_equal(tc._squeeze_grid(a3), a2)


def test_snap_index_picks_nearest(tc):
    spec = np.array([0.0, 10.0, 20.0, 30.0, 40.0])
    assert tc._snap_index(spec, 0.0) == 0
    assert tc._snap_index(spec, 12.0) == 1
    assert tc._snap_index(spec, 26.0) == 3


def test_greedy_climb_walks_diagonally_to_min(tc):
    # bowl with min at (4,4); 8-neighbor climb from (0,0) goes diagonally.
    iy, ix = np.mgrid[0:5, 0:5]
    grid = (iy - 4.0) ** 2 + (ix - 4.0) ** 2
    path = tc._greedy_climb(grid, (0, 0), 'min')
    assert path[0] == (0, 0)
    assert path[-1] == (4, 4)
    assert path == [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]
    vals = [grid[c] for c in path]
    assert all(b < a for a, b in zip(vals, vals[1:]))  # strictly improving


def test_greedy_climb_max_sense(tc):
    iy, ix = np.mgrid[0:5, 0:5]
    grid = -((iy - 4.0) ** 2 + (ix - 4.0) ** 2)  # max at (4,4)
    path = tc._greedy_climb(grid, (0, 0), 'max')
    assert path[-1] == (4, 4)


def test_greedy_climb_skips_nan_neighbors(tc):
    iy, ix = np.mgrid[0:5, 0:5]
    grid = (iy - 4.0) ** 2 + (ix - 4.0) ** 2
    grid[1, 1] = np.nan  # block the diagonal; climb must route around it
    path = tc._greedy_climb(grid, (0, 0), 'min')
    assert (1, 1) not in path
    assert path[-1] == (4, 4)


def test_greedy_climb_start_at_optimum_is_single_point(tc):
    iy, ix = np.mgrid[0:5, 0:5]
    grid = (iy - 4.0) ** 2 + (ix - 4.0) ** 2
    path = tc._greedy_climb(grid, (4, 4), 'min')
    assert path == [(4, 4)]


def test_greedy_climb_nan_start_raises(tc):
    grid = np.zeros((3, 3))
    grid[0, 0] = np.nan
    with pytest.raises(ValueError) as e:
        tc._greedy_climb(grid, (0, 0), 'min')
    assert 'NaN' in str(e.value)


def test_greedy_climb_invalid_sense_raises(tc):
    with pytest.raises(ValueError):
        tc._greedy_climb(np.zeros((2, 2)), (0, 0), 'minimize')


def test_greedy_climb_min_rel_improvement_ignores_round_off_plateau(tc):
    # a flat column (differences ~1e-8 relative, like a converged simulation's
    # noise) next to a genuinely better cell two columns away that the climb
    # cannot reach without first taking a noise-level step.
    grid = np.full((3, 3), 100.0)
    grid[:, 0] = [100.0, 100.0 + 1e-6, 100.0 + 2e-6]   # plateau, rising "up"
    grid[2, 1] = 100.0 + 3e-6
    grid[2, 2] = 150.0
    # default (any strict improvement): walks the plateau and reaches 150
    path = tc._greedy_climb(grid, (0, 0), 'max')
    assert path[-1] == (2, 2)
    # with a threshold above the noise the start is already a local optimum
    # ((1, 1) is exactly equal, the plateau cells are within tolerance)
    assert tc._greedy_climb(grid, (0, 0), 'max', min_rel_improvement=1e-4) == [(0, 0)]
    # a real improvement still passes the threshold
    assert tc._greedy_climb(grid, (2, 1), 'max', min_rel_improvement=1e-4) == [(2, 1), (2, 2)]
    # zero current value: any strict improvement counts
    grid0 = np.zeros((1, 2)); grid0[0, 1] = 1e-9
    assert tc._greedy_climb(grid0, (0, 0), 'max', min_rel_improvement=1e-4) == [(0, 0), (0, 1)]


def _synthetic_results():
    # 6x7 grid (ny=6 over spec_2, nx=7 over spec_1). MPSP-like bowl (min at a
    # corner), IRR-like plane (max at the opposite corner), plus a constant.
    spec_1 = np.linspace(0.0, 30.0, 7)   # x
    spec_2 = np.linspace(0.0, 0.5, 6)    # y
    iy, ix = np.mgrid[0:6, 0:7]
    mpsp = (iy - 0.0) ** 2 + (ix - 0.0) ** 2 + 1.0   # min at (0,0)
    irr = iy.astype(float) + ix.astype(float)         # max at (5,6)
    const = np.full((6, 7), 3.0)
    results = {
        'MPSP': mpsp[None, :, :],
        'IRR': irr[None, :, :],
        'Const': const[None, :, :],
    }
    return results, spec_1, spec_2


def test_derive_levels_ranges_and_errors(tc):
    lv, ct = tc._derive_levels(np.array([1.0, np.nan, 3.0]))
    assert lv[0] == 1.0 and lv[-1] == 3.0
    assert ct[0] == 1.0 and ct[-1] == 3.0
    for bad in (np.array([np.nan, np.nan]), np.array([2.0, 2.0])):
        with pytest.raises(ValueError):
            tc._derive_levels(bad)


def test_plot_returns_fig_ax_and_trajectory_data(tc):
    results, s1, s2 = _synthetic_results()
    fig, ax, traj = tc.plot_metric_with_trajectories(
        results, s1, s2, color_metric='MPSP', baseline_point=(30.0, 0.5),
        trajectory_metrics=['MPSP', 'IRR'], senses={'MPSP': 'min', 'IRR': 'max'})
    assert set(traj) == {'MPSP', 'IRR'}
    # baseline snapped to the far corner (ix=6, iy=5)
    assert traj['MPSP']['path_ij'][0] == (5, 6)
    assert traj['MPSP']['optimum_xy'] == (0.0, 0.0)      # min bowl corner
    assert traj['IRR']['n_steps'] == 0                    # baseline already max
    # each path strictly improves per its sense
    mpsp = tc._squeeze_grid(results['MPSP'])
    vals = [mpsp[c] for c in traj['MPSP']['path_ij']]
    assert all(b < a for a, b in zip(vals, vals[1:]))
    # at least one trajectory polyline was drawn, plus a legend
    assert len(ax.lines) >= 1
    assert ax.get_legend() is not None
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_missing_sense_raises(tc):
    results, s1, s2 = _synthetic_results()
    with pytest.raises(KeyError):
        tc.plot_metric_with_trajectories(
            results, s1, s2, color_metric='MPSP', baseline_point=(0.0, 0.0),
            trajectory_metrics=['IRR'], senses={})


def test_plot_constant_color_metric_raises(tc):
    results, s1, s2 = _synthetic_results()
    with pytest.raises(ValueError):
        tc.plot_metric_with_trajectories(
            results, s1, s2, color_metric='Const', baseline_point=(0.0, 0.0),
            trajectory_metrics=['MPSP'], senses={'MPSP': 'min'})


def test_plot_trajectory_styles_and_legend_kwargs(tc):
    results, s1, s2 = _synthetic_results()
    fig, ax, traj = tc.plot_metric_with_trajectories(
        results, s1, s2, color_metric='MPSP', baseline_point=(30.0, 0.5),
        trajectory_metrics=['MPSP', 'IRR'], senses={'MPSP': 'min', 'IRR': 'max'},
        trajectory_linestyles={'MPSP': '--'}, trajectory_linewidths={'MPSP': 3.0},
        trajectory_markers={'MPSP': 'o'}, trajectory_marker_sizes={'MPSP': 7.0},
        legend_kwargs={'loc': 'lower left'})
    # the MPSP polyline (multi-point line) carries the requested style/width
    polylines = [l for l in ax.lines if len(l.get_xdata()) > 1]
    assert len(polylines) == 1                       # IRR has n_steps == 0
    assert polylines[0].get_linestyle() == '--'
    assert polylines[0].get_linewidth() == 3.0
    # legend entry mirrors the style; unspecified metric keeps the defaults
    by_label = {h.get_label(): h for h in ax.get_legend().legend_handles}
    assert by_label['MPSP'].get_linestyle() == '--'
    assert by_label['MPSP'].get_linewidth() == 3.0
    assert by_label['IRR'].get_linestyle() == '-'
    assert by_label['IRR'].get_linewidth() == 1.6
    assert by_label['MPSP'].get_marker() == 'o'
    assert by_label['IRR'].get_marker() == '*'      # default optimum_marker
    # the optimum marker drawn on the axes uses the per-metric marker/size too
    circles = [l for l in ax.lines if l.get_marker() == 'o']
    assert circles and all(l.get_markersize() == 7.0 for l in circles)
    assert by_label['MPSP'].get_markersize() == 7.0
    assert by_label['IRR'].get_markersize() == 12   # default optimum_marker_size
    assert ax.get_legend()._loc == 3                 # matplotlib code for 'lower left'
    import matplotlib.pyplot as plt
    plt.close(fig)
