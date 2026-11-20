# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Least-squares fit of constants of the simplified reversible O2 model to the
50 C data, set up as in the heat-map script that the authors of Kim et al.,
Appl. Catal. B: Environ. Energy 386 (2026) 126403, deposited.

The script fits the ten 50 C constants of its model to the table of 16
samples it carries (``read_heatmap_script_data``). The fitted constants are
not in the deposit, and the script's own procedure (a particle swarm of
three iterations followed by five least-squares iterations, integrated with
``ode23s`` at default tolerances) is random and truncated, so its result
cannot be recomputed. :func:`fit_simplified_reversible_constants` is the fit
this package makes instead: the same data, the same initial condition and
the same residuals, carried to convergence with an accurate solver, for the
constants the deposit gives no value for. The default constants ``k_3`` and
``k_4`` of the model are its result (``NOTES.md``, "Constants at 50 C").

Nothing is loaded or fitted on import.
"""
import numpy as np
from scipy.optimize import least_squares

from .pd_c_3hp_oxidation import (
    NORMALIZED, SIMPLIFIED_REVERSIBLE_PARAMETER_SETS,
    load_o2_simplified_reversible_model, simulate_at_times)
from .reference import read_heatmap_script_data

__all__ = ('fit_simplified_reversible_constants',)

# Starting value of a constant that the script's vector k0 prints as zero
# (k_3, printed 0.0000, i.e. below 5e-5): a fit cannot be scaled by 0.
_START_IF_ZERO = 1e-5


def fit_simplified_reversible_constants(free=('k_3', 'k_4'), values=None,
                                        start=None, data=None):
    """Fit constants of the simplified reversible model to the 50 C data.

    The model (:func:`.load_o2_simplified_reversible_model`) is run as the
    authors' heat-map script runs it when it fits: at 50 C, from an initial
    3-HP "mole fraction" of 0.000339 (the ``'fitting'`` initial condition).
    The residuals are the differences between the five water-free liquid
    fractions ``X_HP3, X_MA, X_AAc, X_OXA, X_D`` of the model and of the
    data at the sample times, all of them, unweighted (for the script's
    table 16 x 5 = 80 values, the t = 0 row included, where they are zero).
    The constants in ``free`` are varied to minimize the sum of their
    squares; the others are held.

    Parameters
    ----------
    free : sequence of str, optional
        Names of the 50 C constants to fit, out of ``'k_1'`` to ``'k_5'``
        (rate constants, 1/min) and ``'K_6'`` to ``'K_10'`` (equilibrium
        constants). The names of the model parameters, ``'k_1_ref'`` and so
        on, are accepted as well. Default ``('k_3', 'k_4')``: the two
        constants for which the deposit gives no value.
    values : dict, optional
        The constants that are held, and any others, keyed like
        :data:`.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS` (``'k_1_ref'`` ...).
        Constants not given take the value of the ``'deposit'`` set.
    start : sequence of float, optional
        Starting values of the free constants, in the order of ``free``;
        all positive. By default the vector ``k0`` printed in the script
        (``'k0_printed'``), with 1e-5 where it prints zero (``k_3``).
    data : pandas.DataFrame, optional
        Samples to fit: ``time`` (min, starting at 0) and the five
        fractions. By default the table of the script
        (:func:`.read_heatmap_script_data`).

    Returns
    -------
    dict
        ``'fitted'``
            The fitted constants, ``{name: value}`` with the names as in
            ``values`` (``'k_3_ref'`` ...).
        ``'standard_errors'``
            Their standard errors, same keys: the square roots of the
            diagonal of ``s2*inv(J.T @ J)``, with ``J`` the Jacobian of the
            residuals at the optimum and ``s2`` the sum of squares over
            the number of residuals minus the number of fitted constants.
            A linearized measure of how well the data determine each
            constant with the others held; ``nan`` if the Jacobian is
            singular.
        ``'sum_of_squares'``
            The sum of the squared residuals at the optimum.
        ``'n_residuals'``
            The number of residuals.
        ``'values'``
            The complete constant set, held and fitted, keyed like
            :data:`.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS`.

    Raises
    ------
    ValueError
        If a name in ``free`` or ``values`` is not a constant of the model,
        a name is given twice, or ``start`` does not have one positive
        value per free constant.

    Notes
    -----
    This is not the authors' procedure, which cannot be repeated (see the
    module docstring). It differs from it in three ways: the least-squares
    iteration is carried to convergence; the integration is accurate (CVODE
    at a relative tolerance of 1e-12, where the script uses ``ode23s`` at
    its default 1e-3); and each trajectory is one integration from t = 0
    through all sample times, where the script restarts its solver at every
    sample. The constants are fitted as multiples of their starting
    values, bounded below by zero as in the script, with a central-
    difference Jacobian. The function is deterministic: the same call
    returns the same numbers. From other starting points the optimum is
    found again to about 1e-5 of ``k_3``, the less well determined of the
    two default constants. With the defaults the fit takes less than a
    second and returns ``k_3`` = 1.4200e-5 1/min and ``k_4`` = 2.3204e-3
    1/min, the values of the ``'deposit'`` set, with standard errors of
    31 % and 3.7 %.
    """
    names = tuple(SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['deposit'])
    free = [name if name.endswith('_ref') else name + '_ref'
            for name in free]
    constants = dict(SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['deposit'])
    constants.update(values or {})
    unknown = [name for name in (*free, *constants) if name not in names]
    if unknown:
        raise ValueError(
            f'unknown constants {unknown} of the simplified reversible '
            f'model: choose from {list(names)}')
    if len(set(free)) != len(free) or not free:
        raise ValueError(f'free must name each constant once: {free}')
    if start is None:
        printed = SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['k0_printed']
        start = [printed[name] or _START_IF_ZERO for name in free]
    start = np.asarray(start, dtype=float)
    if start.shape != (len(free),) or not (start > 0.).all():
        raise ValueError(
            'start must hold one positive value per free constant, '
            f'not {start!r}')
    if data is None:
        data = read_heatmap_script_data()
    times = data['time'].to_numpy(dtype=float)
    observed = data[list(NORMALIZED)].to_numpy(dtype=float)

    km = load_o2_simplified_reversible_model(initial_condition='fitting')
    integrator = km._te.integrator
    integrator.relative_tolerance = 1e-12
    integrator.absolute_tolerance = 1e-16
    for name, value in constants.items():
        km.set_value(name, float(value))

    def residuals(multiples):
        for name, value in zip(free, multiples*start):
            km.set_value(name, float(value))
        modeled = simulate_at_times(km, times)[list(NORMALIZED)].to_numpy()
        return (modeled - observed).ravel()

    fit = least_squares(residuals, np.ones(len(free)), bounds=(0., np.inf),
                        jac='3-point', diff_step=1e-4, xtol=1e-12,
                        ftol=1e-14, gtol=1e-14)
    fitted = fit.x*start
    sum_of_squares = float(2.*fit.cost)
    n = fit.fun.size
    jacobian = fit.jac/start        # with respect to the constants
    try:
        covariance = (np.linalg.inv(jacobian.T @ jacobian)
                      * sum_of_squares/(n - len(free)))
        errors = np.sqrt(np.diag(covariance))
    except np.linalg.LinAlgError:
        errors = np.full(len(free), np.nan)
    constants.update(zip(free, map(float, fitted)))
    return {'fitted': dict(zip(free, map(float, fitted))),
            'standard_errors': dict(zip(free, map(float, errors))),
            'sum_of_squares': sum_of_squares,
            'n_residuals': int(n),
            'values': {name: float(constants[name]) for name in names}}
