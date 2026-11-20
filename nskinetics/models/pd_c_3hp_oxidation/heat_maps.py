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
           'simulate_heat_maps', 'fig11_panel_heat_map', 'fig11_heat_maps',
           'table3_o2_heat_maps', 'heat_map_maximum', 'heat_map_ridge',
           'contour_offsets', 'contour_shifts',)

#%% The grid of the published figure
#: Temperatures of the heat-map rows in degrees Celsius, 30 to 100 in steps
#: of 1 (the script's 303.15 to 373.15 K in steps of 1 K).
FIG11_T_C = 30. + np.arange(71.)

#: Plotted times of the heat-map columns in min (the script's
#: ``linspace(10, 1680, 300)``; 1680 min is the last 50 C sampling time).
FIG11_TIMES = np.linspace(10., 1680., 300)

#: Plotted time minus reaction time in the published Fig. 11, in min.
FIG11_TIME_OFFSET = 10.

#: Contour levels drawn in Fig. 11.
FIG11_LEVELS = (0.1, 0.3, 0.5)

#: Minutes per pixel in the main panels of the published raster of Fig. 11
#: (1604 x 596 px), from the positions of the tick marks of panel a (the
#: axis calibration is given in ``reference_data/README.md``).
FIG11_MIN_PER_PX = 3.0723

#: Degrees Celsius per pixel in the main panels of the published raster of
#: Fig. 11, from the positions of the tick marks of panel a.
FIG11_C_PER_PX = 0.1476

# Constant set of the simplified reversible model that reproduces each panel
# of Fig. 11, and the index of the panel's quantity in the tuple returned by
# simulate_heat_maps.
_FIG11_PANELS = {'a': ('fig11a', 2), 'b': ('fig11b', 3)}

# Range of temperatures in degrees Celsius that simulate_heat_maps accepts.
# It only guards against temperatures given in kelvin.
_T_C_RANGE = (-50., 200.)


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
        Temperatures in degrees Celsius (not K), one per row. Defaults to
        :data:`FIG11_T_C`.
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
        Row temperatures in degrees Celsius.
    times : numpy.ndarray
        Column plotted times in min.
    selectivity : numpy.ndarray
        ``X_MA``, shape ``(len(T_C), len(times))``.
    MA_yield : numpy.ndarray
        ``(1 - X_HP3)*X_MA``, same shape.

    Raises
    ------
    ValueError
        If the model has no temperature parameter ``T`` (the H2O2 model and
        the reversible O2 variant have none). Also if a temperature is
        outside -50 to 200 C, which is what temperatures given in kelvin
        look like; that range is a guard against that mistake only, not a
        range of validity.
    """
    if 'T' not in te_r._te.model.getGlobalParameterIds():
        raise ValueError(
            'simulate_heat_maps needs a model with a temperature parameter '
            'T: one from load_o2_model or load_o2_simplified_reversible_model '
            '(the H2O2 model and the reversible O2 variant have none)')
    T_C = FIG11_T_C if T_C is None else np.asarray(T_C, dtype=float)
    if not ((T_C >= _T_C_RANGE[0]) & (T_C <= _T_C_RANGE[1])).all():
        raise ValueError(
            f'T_C must lie within {_T_C_RANGE[0]:g} to {_T_C_RANGE[1]:g}: '
            'T_C is in degrees Celsius, not kelvin (for 348.15 K pass 75).')
    times = FIG11_TIMES if times is None else np.asarray(times, dtype=float)
    selectivity = np.empty((len(T_C), len(times)))
    MA_yield = np.empty_like(selectivity)
    T_before = te_r.get_value('T')
    try:
        for i, row_T_C in enumerate(T_C):
            te_r.set_value('T', row_T_C + 273.15)
            sim = simulate_at_times(te_r, times - time_offset,
                                    ['X_HP3', 'X_MA'])
            selectivity[i] = sim['X_MA']
            MA_yield[i] = (1. - sim['X_HP3'])*sim['X_MA']
    finally:
        te_r.set_value('T', T_before)
        te_r.reset()
    return T_C, times, selectivity, MA_yield


def fig11_panel_heat_map(panel):
    """Compute the heat map of one panel of the published Fig. 11.

    The simplified reversible model is run as in the authors' heat-map
    script (unit mole fractions at plotted time 10 min), with the constant
    set that reproduces the panel: ``'fig11a'`` for the selectivity map and
    ``'fig11b'`` for the yield map (see ``NOTES.md``).

    Parameters
    ----------
    panel : {'a', 'b'}
        ``'a'`` for the MA "selectivity" map, ``'b'`` for the MA "yield"
        map.

    Returns
    -------
    T_C : numpy.ndarray
        Row temperatures in degrees Celsius (:data:`FIG11_T_C`).
    times : numpy.ndarray
        Column plotted times in min (:data:`FIG11_TIMES`), which are
        reaction times plus :data:`FIG11_TIME_OFFSET`.
    matrix : numpy.ndarray
        The mapped quantity, shape ``(71, 300)``.

    Raises
    ------
    ValueError
        If ``panel`` is not ``'a'`` or ``'b'``.
    """
    if panel not in _FIG11_PANELS:
        raise ValueError(f"panel must be 'a' or 'b', not {panel!r}")
    parameter_set, index = _FIG11_PANELS[panel]
    maps = simulate_heat_maps(
        load_o2_simplified_reversible_model(parameter_set=parameter_set),
        time_offset=FIG11_TIME_OFFSET)
    return maps[0], maps[1], maps[index]


def fig11_heat_maps():
    """Compute the two heat maps of the published Fig. 11.

    Each panel uses its own constant set; see
    :func:`fig11_panel_heat_map`. This takes 2 x 71 integrations: one per
    temperature, as in the authors' script, for each of the two sets.

    Returns
    -------
    T_C, times, selectivity, MA_yield
        As :func:`simulate_heat_maps`, on the 71 x 300 grid of the figure.
        ``times`` are plotted times, i.e. reaction times plus
        :data:`FIG11_TIME_OFFSET`.
    """
    T_C, times, selectivity = fig11_panel_heat_map('a')
    return T_C, times, selectivity, fig11_panel_heat_map('b')[2]


def table3_o2_heat_maps(time_offset=0.):
    """Compute the Fig. 11 quantities with the Table 3 O2 model.

    This is what the text of the paper describes for Fig. 11 (the O2 model
    extended with the Arrhenius relation), computed with
    :func:`.load_o2_model` on the grid of the figure. It is not what the
    authors' heat-map script computes, and it differs from the published
    maps (see ``NOTES.md``).

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
    """Return the largest cell of a heat map.

    Parameters
    ----------
    T_C, times, matrix
        A heat map: row temperatures in degrees Celsius, column times in
        min and one of the two matrices returned by
        :func:`simulate_heat_maps`.

    Returns
    -------
    value : float
        The largest value of ``matrix``.
    T_C_max : float
        Temperature of its row, in degrees Celsius.
    time_max : float
        Time of its column in min, on the time axis of ``times``.
    """
    i, j = np.unravel_index(np.argmax(matrix), matrix.shape)
    return matrix[i, j], T_C[i], times[j]


def heat_map_ridge(times, matrix):
    """Return, per row, the time of the row maximum (the red dashed line).

    Parameters
    ----------
    times, matrix
        Column times in min and one of the two matrices returned by
        :func:`simulate_heat_maps`.

    Returns
    -------
    numpy.ndarray
        One time per row of ``matrix``, on the time axis of ``times``.
    """
    return np.asarray(times)[np.argmax(matrix, axis=1)]


#%% Comparison with digitized contour lines
def contour_offsets(T_C, times, matrix, level, sample_times, sample_T_C,
                    min_per_px=FIG11_MIN_PER_PX, C_per_px=FIG11_C_PER_PX,
                    max_offset=None):
    """Offsets of digitized contour samples from a map's contour, in pixels.

    For each sample the offset is ``(value - level)/|gradient|``, with the
    map interpolated bilinearly (as a contour routine does) and the
    gradient taken per pixel of the published raster. For a sample close
    to the map's contour this is its distance from that contour along the
    normal; it is positive where the map is above ``level`` at the sample.
    It is a linearization: for a sample more than about 1 px away it is
    an indication, not a distance.

    Parameters
    ----------
    T_C, times, matrix
        A heat map, as returned by :func:`simulate_heat_maps`.
    level : float
        Level of the digitized contour line.
    sample_times, sample_T_C : array-like
        Plotted times (min) and temperatures (degrees Celsius) of the
        digitized samples.
    min_per_px, C_per_px : float, optional
        Scale of the raster the samples were read from.
    max_offset : float, optional
        If given, offsets are clipped to ``-max_offset .. max_offset`` px,
        so that samples where the map is nearly flat do not dominate an
        rms. The rms offsets quoted in ``NOTES.md`` and on the
        documentation page use 30; for maps within about 5 px of the
        samples the clip changes nothing.

    Returns
    -------
    numpy.ndarray
        One offset per sample that lies inside the grid.
    """
    t = np.asarray(sample_times, dtype=float)
    T = np.asarray(sample_T_C, dtype=float)
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
    offsets = (at(T, t) - level)/np.hypot(per_px_t, per_px_T)
    if max_offset is not None:
        offsets = np.clip(offsets, -max_offset, max_offset)
    return offsets


def contour_shifts(T_C, times, matrix, level, branch, sample_times,
                   sample_T_C, along='time'):
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
    sample_times, sample_T_C : array-like
        Plotted times (min) and temperatures (degrees Celsius) of the
        digitized samples.
    along : {'time', 'T_C'}, optional
        ``'time'`` (default): shift in min at the temperature of each
        sample. ``'T_C'``: shift in degrees Celsius at the time of each
        sample.

    Returns
    -------
    numpy.ndarray
        Map contour position minus sample position, one value per sample;
        ``nan`` where the map has no such crossing of ``level`` at the
        sample's temperature (or time). Where it has several, the one
        nearest to the sample is used.
    """
    t = np.asarray(sample_times, dtype=float)
    T = np.asarray(sample_T_C, dtype=float)
    if along == 'time':
        # Rows of the map at the sample temperatures.
        axis, position = times, t
        lines, inside = _interpolate(T_C, matrix, T)
    elif along == 'T_C':
        # Columns of the map at the sample times.
        axis, position = T_C, T
        lines, inside = _interpolate(times, matrix.T, t)
    else:
        raise ValueError(f"along must be 'time' or 'T_C', not {along!r}")
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
