# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
The O2 experiments at 25, 75 and 100 C and the Arrhenius plot (Fig. 10 of
Kim et al., Appl. Catal. B: Environ. Energy 386 (2026) 126403): the measures
used to compare a model with the experimental points of Fig. 10a-c, and the
straight lines of Fig. 10d.

The paper does not say of which model the points of Fig. 10d are rate
constants. They cannot be recomputed and are read from the deposited table
(``read_fig10d_arrhenius``, which also explains their names); the lines
drawn through them are recomputed here.

The points are taken to be rate constants of the authors' simplified
reversible model (``load_o2_simplified_reversible_model``), and the 50 C
point supplies three of the default constants of that model. Taking them so
means that the plotted constants follow the rate law of the authors'
heat-map script, which the reasons given in ``NOTES.md`` support but do not
establish. The model curves of Fig. 10a-c are not reproduced by this
package.
"""
import numpy as np
import pandas as pd

from .pd_c_3hp_oxidation import NORMALIZED
from .reference import read_fig10d_arrhenius

__all__ = ('FIG10_TEMPERATURES_C', 'FIG10D_REGRESSION_TEMPERATURES_C',
           'TEXT_ACTIVATION_ENERGIES',
           'pooled_r2', 'rms_difference', 'fig10d_points',
           'fig10d_regression_lines', 'fig10d_text_lines',)

#: Temperatures in C of the experiments of Fig. 10a, b and c.
FIG10_TEMPERATURES_C = (25, 75, 100)

#: Temperatures in C of the points each line of Fig. 10d is fitted to. Every
#: line leaves out one of the four plotted points: 50 C for ``k_1``, 25 C
#: for ``k_2`` and ``k_5``.
FIG10D_REGRESSION_TEMPERATURES_C = {'k_1': (25, 75, 100),
                                    'k_2': (50, 75, 100),
                                    'k_5': (50, 75, 100)}

#: Activation energies in J/mol that the text of the paper quotes for the
#: three constants of Fig. 10d (59.8, 91.6 and 59.5 kJ/mol), with the digits
#: of the authors' heat-map script. Two of them (``k_2``, ``k_5``) are not
#: the slopes of the lines drawn in Fig. 10d (68.96 and 52.09 kJ/mol); the
#: ``k_1`` value agrees with its line (59.93 kJ/mol) to 0.15 kJ/mol. The
#: keys are those of ``read_fig10d_arrhenius`` (3-HP to MA, MA to AAc, MA to
#: D); the simplified reversible model has these values as ``Ea_1``,
#: ``Ea_2`` and ``Ea_5``.
TEXT_ACTIVATION_ENERGIES = {'k_1': 59781.2, 'k_2': 91594.6, 'k_5': 59549.8}

_R_GAS = 8.314          # J/(mol K), as in the authors' heat-map script
_T_REF = 323.15         # K


#%% Agreement of a model with a table of normalized concentrations
def pooled_r2(observed, modeled):
    """Coefficient of determination over all species and rows together.

    Parameters
    ----------
    observed, modeled : pandas.DataFrame
        Tables with the columns :data:`.NORMALIZED` and the same rows (for
        example an ``'experiment'`` table of ``reference`` and the output
        of ``simulate_at_times`` at its times).

    Returns
    -------
    float
        ``1 - SS_res/SS_tot`` with both sums taken over the five species
        and all rows, and ``SS_tot`` about the mean of all observed values
        (unweighted, as in the authors' script).
    """
    observed = observed[list(NORMALIZED)].to_numpy()
    modeled = modeled[list(NORMALIZED)].to_numpy()
    ss_res = ((observed - modeled)**2).sum()
    ss_tot = ((observed - observed.mean())**2).sum()
    return 1. - ss_res/ss_tot


def rms_difference(observed, modeled):
    """Root-mean-square difference over all species and rows together.

    Parameters
    ----------
    observed, modeled : pandas.DataFrame
        As for :func:`pooled_r2`.

    Returns
    -------
    float
    """
    difference = (observed[list(NORMALIZED)].to_numpy()
                  - modeled[list(NORMALIZED)].to_numpy())
    return np.sqrt((difference**2).mean())


#%% Arrhenius plot (Fig. 10d)
def fig10d_points():
    """The points of Fig. 10d with their temperatures.

    Returns
    -------
    pandas.DataFrame
        One row per temperature: ``T_C`` (25, 50, 75, 100), ``inv_T``
        (``1/(T_C + 273.15)`` in 1/K, unrounded) and ``ln_k_1``,
        ``ln_k_2``, ``ln_k_5`` as deposited (natural logarithms of three
        rate constants, taken to be in 1/min and named as by
        :func:`.read_fig10d_arrhenius`: 3-HP to MA, MA to AAc and MA to D).
        They are taken to be ``k_1``, ``k_2`` and ``k_5`` of the simplified
        reversible model.

    Notes
    -----
    The deposited table prints 1/T to five decimals only; the temperatures
    are recovered from it and the abscissa recomputed.
    """
    points = read_fig10d_arrhenius()['points']
    T_C = 25.*np.round((1./points['inv_T'] - 273.15)/25.)
    table = pd.DataFrame({'T_C': T_C, 'inv_T': 1./(T_C + 273.15)})
    for name in TEXT_ACTIVATION_ENERGIES:
        table['ln_' + name] = points['ln_' + name]
    return table.sort_values('T_C').reset_index(drop=True)


def _line(slope, intercept, temperatures_C):
    return {'slope': slope, 'intercept': intercept,
            'Ea': -slope*_R_GAS, 'temperatures_C': tuple(temperatures_C)}


def fig10d_regression_lines():
    """The straight lines drawn in Fig. 10d, recomputed from the points.

    Each line is the ordinary least-squares fit of ``ln k`` on
    ``1/(T_C + 273.15)`` through three of the four points
    (:data:`FIG10D_REGRESSION_TEMPERATURES_C`).

    Returns
    -------
    dict of dict
        For ``'k_1'``, ``'k_2'`` and ``'k_5'``: ``'slope'`` (K),
        ``'intercept'``, ``'Ea'`` (the activation energy ``-slope*R`` in
        J/mol, with R = 8.314 J/(mol K)) and ``'temperatures_C'`` (the
        temperatures of the fitted points).
    """
    points = fig10d_points()
    lines = {}
    for name, temperatures_C in FIG10D_REGRESSION_TEMPERATURES_C.items():
        fitted = points[points['T_C'].isin(temperatures_C)]
        slope, intercept = np.polyfit(fitted['inv_T'], fitted['ln_' + name],
                                      1)
        lines[name] = _line(slope, intercept, temperatures_C)
    return lines


def fig10d_text_lines():
    """The lines that the activation energies quoted in the text imply.

    Each line has the slope ``-Ea/R`` of :data:`TEXT_ACTIVATION_ENERGIES`
    and passes through the 50 C point of Fig. 10d. These are the constants
    ``k_1(T)``, ``k_2(T)`` and ``k_5(T)`` of the simplified reversible
    model with its default constants
    (``load_o2_simplified_reversible_model``): its 50 C values are the 50 C
    points, and its activation energies, those of the authors' heat-map
    script, are the ones the text quotes. The lines are therefore how that
    model, like the script, extrapolates the three constants from 50 C.
    They miss the points at 25, 75 and 100 C (``NOTES.md``).

    Returns
    -------
    dict of dict
        As returned by :func:`fig10d_regression_lines`; ``'temperatures_C'``
        is ``(50,)``, the point each line passes through.
    """
    points = fig10d_points()
    at_50C = points[points['T_C'] == 50.].iloc[0]
    lines = {}
    for name, Ea in TEXT_ACTIVATION_ENERGIES.items():
        slope = -Ea/_R_GAS
        lines[name] = _line(slope, at_50C['ln_' + name] - slope/_T_REF,
                            (50,))
    return lines
