# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Temperature-time heat maps of MA "selectivity" and "yield" (Fig. 11 of Kim
et al., Appl. Catal. B: Environ. Energy 386 (2026) 126403), computed with
the NSKinetics models of this package, and the measures used to compare a
computed map with the contour lines digitized from the published figure.

The two mapped quantities are those of the authors' heat-map script:
"selectivity" is the water-free liquid fraction of MA, ``X_MA``, and "yield"
is ``(1 - X_HP3)*X_MA``. The grid is the script's: 30 to 100 C in steps of
1 C and 300 plotted times from 10 to 1680 min. In the published figure the
plotted time is the reaction time plus 10 min (the script puts the pure
3-HP feed at the first plotted time); see ``NOTES.md``.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .pd_c_3hp_oxidation import (
    load_o2_model, load_o2_simplified_reversible_model, simulate_at_times)

__all__ = ('FIG11_T_C', 'FIG11_TIMES', 'FIG11_TIME_OFFSET', 'FIG11_LEVELS',
           'FIG11_MIN_PER_PX', 'FIG11_C_PER_PX',
           'simulate_heat_maps', 'fig11_heat_maps', 'table3_o2_heat_maps',
           'heat_map_maximum', 'heat_map_ridge',
           'contour_offsets', 'contour_shifts',)

#%% The grid of the published figure
#: Temperatures of the heat-map rows in C (the script's 303.15:1:373.15 K).
FIG11_T_C = 30. + np.arange(71.)

#: Plotted times of the heat-map columns in min (the script's
#: ``linspace(10, 1680, 300)``; 1680 min is the last 50 C sampling time).
FIG11_TIMES = np.linspace(10., 1680., 300)

#: Plotted time minus reaction time in the published Fig. 11, in min.
FIG11_TIME_OFFSET = 10.

#: Contour levels drawn in Fig. 11.
FIG11_LEVELS = (0.1, 0.3, 0.5)

#: Scale of the main panels of the published raster (1604 x 596 px), from
#: the tick calibration of the digitization: min and C per pixel.
FIG11_MIN_PER_PX = 3.0723
FIG11_C_PER_PX = 0.1476


#%% Heat maps
def simulate_heat_maps(te_r, T_C=None, times=None, time_offset=0.):
    """Compute the two Fig. 11 quantities of a model over temperature and time.

    Parameters
    ----------
    te_r : nskinetics.KineticModel
        A model with a temperature parameter ``T`` (K) and the water-free
        liquid fractions ``X_HP3`` and ``X_MA``: one from
        :func:`.load_o2_simplified_reversible_model` or
        :func:`.load_o2_model`. Its ``T`` is restored on return.
    T_C : array-like, optional
        Temperatures in C, one per row. Defaults to :data:`FIG11_T_C`.
    times : array-like, optional
        Increasing plotted times in min, one per column. Defaults to
        :data:`FIG11_TIMES`.
    time_offset : float, optional
        Plotted time minus reaction time, in min. Every row is one
        integration that starts from the feed at reaction time 0, i.e. at
        plotted time ``time_offset``. Default 0 (plotted time is reaction
        time); the published Fig. 11 has :data:`FIG11_TIME_OFFSET`.

    Returns
    -------
    T_C : numpy.ndarray
        Row temperatures in C.
    times : numpy.ndarray
        Column plotted times in min.
    selectivity : numpy.ndarray
        ``X_MA``, shape ``(len(T_C), len(times))``.
    MA_yield : numpy.ndarray
        ``(1 - X_HP3)*X_MA``, same shape.
    """
    T_C = FIG11_T_C if T_C is None else np.asarray(T_C, dtype=float)
    times = FIG11_TIMES if times is None else np.asarray(times, dtype=float)
    selectivity = np.empty((len(T_C), len(times)))
    MA_yield = np.empty_like(selectivity)
    T_before = te_r.get_value('T')
    try:
        for i, T in enumerate(T_C):
            te_r.set_value('T', T + 273.15)
            sim = simulate_at_times(te_r, times - time_offset,
                                    ['X_HP3', 'X_MA'])
            selectivity[i] = sim['X_MA']
            MA_yield[i] = (1. - sim['X_HP3'])*sim['X_MA']
    finally:
        te_r.set_value('T', T_before)
        te_r.reset()
    return T_C, times, selectivity, MA_yield


def fig11_heat_maps():
    """Compute the heat maps of the published Fig. 11.

    The simplified reversible model is run as in the authors' heat-map
    script (unit mole fractions at plotted time 10 min). Each panel uses
    its own constant set: ``'fig11a'`` for the selectivity map and
    ``'fig11b'`` for the yield map (see ``NOTES.md``).

    Returns
    -------
    T_C, times, selectivity, MA_yield
        As :func:`simulate_heat_maps`, on the 71 x 300 grid of the figure.
    """
    maps = [simulate_heat_maps(
                load_o2_simplified_reversible_model(parameter_set=name),
                time_offset=FIG11_TIME_OFFSET)
            for name in ('fig11a', 'fig11b')]
    T_C, times, selectivity, _ = maps[0]
    return T_C, times, selectivity, maps[1][3]


def table3_o2_heat_maps(time_offset=0.):
    """Compute the Fig. 11 quantities with the Table 3 O2 model.

    This is what the text of the paper describes for Fig. 11 (the O2 model
    of Fig. 9b extended with the Arrhenius equation), computed with
    :func:`.load_o2_model` on the grid of the figure. It is not how the
    published maps were made, and it differs from them (see ``NOTES.md``).

    Parameters
    ----------
    time_offset : float, optional
        Plotted time minus reaction time, in min. Default 0: the time axis
        is the reaction time.

    Returns
    -------
    T_C, times, selectivity, MA_yield
        As :func:`simulate_heat_maps`, on the 71 x 300 grid of the figure.
    """
    return simulate_heat_maps(load_o2_model(), time_offset=time_offset)


def heat_map_maximum(T_C, times, matrix):
    """Return the largest cell of a heat map as (value, temperature, time)."""
    i, j = np.unravel_index(np.argmax(matrix), matrix.shape)
    return matrix[i, j], T_C[i], times[j]


def heat_map_ridge(times, matrix):
    """Return, per row, the time of the row maximum (the red dashed line)."""
    return np.asarray(times)[np.argmax(matrix, axis=1)]


#%% Comparison with digitized contour lines
def contour_offsets(T_C, times, matrix, level, t, T,
                    min_per_px=FIG11_MIN_PER_PX, C_per_px=FIG11_C_PER_PX):
    """Offsets of digitized contour samples from a map's contour, in pixels.

    For each sample the offset is ``(value - level)/|gradient|``, with the
    map interpolated bilinearly (as a contour routine does) and the
    gradient taken per pixel of the published raster. For a sample close
    to the map's contour this is its distance from that contour along the
    normal; it is positive where the map is above ``level`` at the sample.
    It is a linearization, so it is not a distance for samples many pixels
    away.

    Parameters
    ----------
    T_C, times, matrix
        A heat map, as returned by :func:`simulate_heat_maps`.
    level : float
        Level of the digitized contour line.
    t, T : array-like
        Plotted times (min) and temperatures (C) of the digitized samples.
    min_per_px, C_per_px : float, optional
        Scale of the raster the samples were read from.

    Returns
    -------
    numpy.ndarray
        One offset per sample that lies inside the grid.
    """
    t, T = np.asarray(t, dtype=float), np.asarray(T, dtype=float)
    inside = ((t >= times[0]) & (t <= times[-1])
              & (T >= T_C[0]) & (T <= T_C[-1]))
    t, T = t[inside], T[inside]
    value = RegularGridInterpolator((T_C, times), matrix)

    def at(T, t):
        return value(np.column_stack((np.clip(T, T_C[0], T_C[-1]),
                                      np.clip(t, times[0], times[-1]))))

    # Central differences: over one pixel in time, and over one row (1 C)
    # in temperature, where the bilinear interpolant is piecewise linear.
    h = 0.5*min_per_px
    per_px_t = (at(T, t + h) - at(T, t - h))
    per_px_T = (at(T + 0.5, t) - at(T - 0.5, t))*C_per_px
    return (at(T, t) - level)/np.hypot(per_px_t, per_px_T)


def contour_shifts(T_C, times, matrix, level, branch, t, T, along='time'):
    """Shifts of a map's contour from digitized contour samples.

    Unlike :func:`contour_offsets` this locates the contour of the map
    itself, so it also holds for a map that is far from the samples.

    Parameters
    ----------
    T_C, times, matrix
        A heat map, as returned by :func:`simulate_heat_maps`.
    level : float
        Level of the digitized contour line.
    branch : {'early', 'late'}
        Side of the ridge the samples are on: ``'early'`` is the rising
        flank of the map in time (left of the red dashed line), ``'late'``
        the falling flank.
    t, T : array-like
        Plotted times (min) and temperatures (C) of the digitized samples.
    along : {'time', 'T'}, optional
        ``'time'`` (default): shift in min at the temperature of each
        sample. ``'T'``: shift in C at the time of each sample.

    Returns
    -------
    numpy.ndarray
        Map contour position minus sample position, one value per sample;
        ``nan`` where the map has no such crossing of ``level`` at the
        sample's temperature (or time). Where it has several, the one
        nearest to the sample is used.
    """
    t, T = np.asarray(t, dtype=float), np.asarray(T, dtype=float)
    if along == 'time':
        # Rows of the map at the sample temperatures.
        axis, position = times, t
        lines, inside = _interpolate(T_C, matrix, T)
    else:
        # Columns of the map at the sample times.
        axis, position = T_C, T
        lines, inside = _interpolate(times, matrix.T, t)
    # The early branch is where the map rises through the level with time
    # and, at a fixed time, with temperature; the late branch where it falls.
    sign = 1. if branch == 'early' else -1.
    shifts = np.full(len(position), np.nan)
    for n in np.flatnonzero(inside):
        d = sign*(lines[n] - level)
        cross = np.flatnonzero((d[:-1] < 0.) & (d[1:] >= 0.))
        if len(cross):
            where = (axis[cross] + (axis[cross + 1] - axis[cross])
                     * d[cross]/(d[cross] - d[cross + 1]))
            nearest = np.argmin(np.abs(where - position[n]))
            shifts[n] = where[nearest] - position[n]
    return shifts


def _interpolate(grid, matrix, points):
    """Rows of ``matrix`` interpolated linearly along its first axis, and
    which of the points lie inside the grid."""
    inside = (points >= grid[0]) & (points <= grid[-1])
    i = np.clip(np.searchsorted(grid, points) - 1, 0, len(grid) - 2)
    w = ((points - grid[i])/(grid[i + 1] - grid[i]))[:, None]
    return (1. - w)*matrix[i] + w*matrix[i + 1], inside
