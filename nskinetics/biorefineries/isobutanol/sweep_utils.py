# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Plot helpers shared by the 2-D grid sweeps.

The contour-plot colormaps, tick helpers and the colorbar-level heuristic
that ``analyses/kinetic_sweeps/evaluate_EtOH_*.py`` and
``analyses/feeding_strategies/evaluate_feeding_strategies.py`` duplicated.
Only :func:`JBEI_UCB_colormap` is called live, by the kinetic sweeps
(``evaluate_feeding_strategies.py`` keeps its own variant of it); the rest
back those scripts' commented-out contour-level and tick cells. The function
bodies are copied verbatim from those scripts, quirks included, so the
figures do not change; in particular
:func:`get_contour_info_from_metric_data` rounds its ``median + 0.2 stdev``
tick with :func:`my_round`'s default base (5), not ``w_ticks_round_to_base``
like its neighbours.

Importing this module imports only numpy, matplotlib and :mod:`math`; the
colormaps import ``biosteam.utils.colors`` when called.
"""
from math import floor, ceil, log

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

__all__ = ('JBEI_UCB_colormap', 'CABBI_green_colormap',
           'tickmarks_from_data', 'tickmarks',
           'get_OOM', 'my_round', 'get_contour_info_from_metric_data')

get_median = np.median
np_round = np.round


#%% Colormaps

def JBEI_UCB_colormap(N_levels=90, reverse=False):
    """Return the JBEI/UCB contour colormap (UCB yellow, JBEI orange, UCB
    blue, dark grey).

    Parameters
    ----------
    N_levels : int, optional
        Number of discrete colors. Default 90.
    reverse : bool, optional
        If True, reverse the color order. Default False.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
    """
    from biosteam.utils import colors
    JBEI_orange = (233/255, 83/255, 39/255)
    UCB_blue = (0/255, 38/255, 118/255)
    UCB_yellow = (253/255, 181/255, 21/255)
    cmap_colors = [
                    UCB_yellow,
                    JBEI_orange,
                    UCB_blue,
                    # colors.CABBI_teal_green.shade(50).RGBn,
                    colors.grey_dark.RGBn]
    if reverse: cmap_colors.reverse()
    return LinearSegmentedColormap.from_list('CABBI', cmap_colors, N_levels)

def CABBI_green_colormap(N_levels=90):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.

    Parameters
    ----------
    N_levels : int, optional
        Number of discrete colors. Default 90.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
    """
    from biosteam.utils import colors
    CABBI_colors = (colors.CABBI_orange.RGBn,
                    colors.CABBI_yellow.RGBn,

                    colors.CABBI_green.RGBn,
                    # colors.CABBI_teal_green.shade(50).RGBn,
                    colors.grey_dark.RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)


#%% Tickmark utils

def tickmarks_from_data(data, accuracy=50, N_points=5):
    """Evenly spaced tick values spanning ``data`` (see :func:`tickmarks`).

    Parameters
    ----------
    data : numpy.ndarray
        Values to span (``data.min()`` to ``data.max()``).
    accuracy : float, optional
        The bounds are rounded outwards to a multiple of this. Default 50.
    N_points : int, optional
        Number of ticks. Default 5.

    Returns
    -------
    list[float]
    """
    dmin = data.min()
    dmax = data.max()
    return tickmarks(dmin, dmax, accuracy, N_points)

def tickmarks(dmin, dmax, accuracy=50, N_points=5):
    """``N_points`` evenly spaced ticks from ``dmin`` rounded down to
    ``dmax`` rounded up, both to a multiple of ``accuracy``.

    Parameters
    ----------
    dmin, dmax : float
        Range to span.
    accuracy : float, optional
        Rounding multiple of the bounds. Default 50.
    N_points : int, optional
        Number of ticks. Default 5.

    Returns
    -------
    list[float]
    """
    dmin = floor(dmin/accuracy) * accuracy
    dmax = ceil(dmax/accuracy) * accuracy
    step = (dmax - dmin) / (N_points - 1)
    return [dmin + step * i for i in range(N_points)]


#%% Contour levels and colorbar ticks

def get_OOM(number):
    """Order of magnitude of ``number``: ``floor(log10(number + 1))``.

    Parameters
    ----------
    number : float

    Returns
    -------
    int
    """
    return floor(log(number+1., 10))

def my_round(x, base=5):
    """Round ``x`` to the nearest multiple of ``base``.

    Parameters
    ----------
    x : float
    base : float, optional
        Default 5.

    Returns
    -------
    float
    """
    return base * round(x/base)

def get_contour_info_from_metric_data(
                                        metric_data, # numpy array
                                        round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM=1,
                                        n_levels_between_cbar_ticks=5,
                                        n_stdevs_for_bounds=1,
                                        lb=None,
                                        ub=None,
                                        multiply_step_size_by=0.5,
                                        w_ticks_round_to_base_divisor=2,
                                        remove_w_ticks_if_greater_than_this_fraction_of_successor=0.8,
                                      ):
    """Contour levels, labelled contour ticks and colorbar ticks for a
    metric's grid, centred on its median +/- ``n_stdevs_for_bounds`` standard
    deviations (NaNs ignored).

    Parameters
    ----------
    metric_data : array_like
        The metric's values over the grid.
    round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM : int, optional
        Tick rounding relative to the data's order of magnitude. Default 1.
    n_levels_between_cbar_ticks : int, optional
        Contour levels per colorbar-tick step. Default 5.
    n_stdevs_for_bounds : float, optional
        Half-width of the level range in standard deviations. Default 1.
    lb, ub : float, optional
        Explicit level bounds (default: from the median and standard
        deviation).
    multiply_step_size_by : float, optional
        Scale of the colorbar-tick step. Default 0.5.
    w_ticks_round_to_base_divisor : float, optional
        Unused (kept for call compatibility).
    remove_w_ticks_if_greater_than_this_fraction_of_successor : float, optional
        Drop a labelled tick whose ratio to the next one exceeds this.
        Default 0.8.

    Returns
    -------
    w_levels : numpy.ndarray
        Contour levels.
    w_ticks : list[float]
        Labelled contour values.
    cbar_ticks : numpy.ndarray
        Colorbar ticks.
    """
    if not type(metric_data) == np.ndarray: metric_data = np.array(metric_data)
    metric_data = metric_data[~np.isnan(metric_data)]
    median, stdev = get_median(metric_data), metric_data.std()
    bound_diff = n_stdevs_for_bounds*stdev
    ub_temp = max(abs(median-bound_diff), abs(median+bound_diff))
    # breakpoint()
    OOM = get_OOM(ub_temp)
    log_ub_temp = log(ub_temp, 10)

    round_to_decimal_place = None

    if log_ub_temp < 1.:
        round_to_decimal_place = int(round(ub_temp/(10.*round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM), 0))
    else:
        round_to_decimal_place = -OOM -1 + round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM

    if lb==None: lb = np_round(median - bound_diff, round_to_decimal_place)
    if ub==None: ub = np_round(median + bound_diff, round_to_decimal_place)

    cbar_ticks_step_size = (10**OOM)/(round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM) * multiply_step_size_by
    w_levels_step_size = cbar_ticks_step_size/n_levels_between_cbar_ticks
    cbar_ticks = np.arange(lb, ub+cbar_ticks_step_size, cbar_ticks_step_size)
    w_levels = np.arange(lb, ub+w_levels_step_size, w_levels_step_size)

    w_ticks_round_to_base = 10**-round_to_decimal_place / 2
    w_ticks = [*set([lb,
               my_round(median-0.6*stdev, w_ticks_round_to_base),
               my_round(median-0.4*stdev, w_ticks_round_to_base),
               my_round(median-0.2*stdev, w_ticks_round_to_base),
               np_round(median, round_to_decimal_place),
               my_round(median+0.2*stdev, ),
               my_round(median+0.4*stdev, w_ticks_round_to_base),
               my_round(median+0.6*stdev, w_ticks_round_to_base),
               ub])]

    # w_ticks.sort()
    # w_ticks_to_remove = []
    # for i in range(len(w_ticks)-1):
    #     if w_ticks[i]/w_ticks[i+1] > remove_w_ticks_if_greater_than_this_fraction_of_successor:
    #         w_ticks_to_remove.append(w_ticks[i])
    # for j in w_ticks_to_remove:
    #     w_ticks.remove(j)
    #     w_ticks.append(j-w_levels_step_size)

    w_ticks.sort()
    w_ticks_to_remove = []
    for i in range(len(w_ticks)-1):
        if w_ticks[i]/w_ticks[i+1] > remove_w_ticks_if_greater_than_this_fraction_of_successor:
            w_ticks_to_remove.append(w_ticks[i])
    for j in w_ticks_to_remove:
        w_ticks.remove(j)

    w_ticks.sort()
    return w_levels, w_ticks, cbar_ticks
