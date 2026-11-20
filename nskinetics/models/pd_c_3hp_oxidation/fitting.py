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
samples it carries (``read_heatmap_script_data``). The constants its run
fitted are not in the deposit, and the script's own procedure (a particle
swarm of three iterations followed by five least-squares iterations,
integrated with ``ode23s`` at default tolerances) is random and truncated,
so its result cannot be recomputed. The package uses every constant the
deposit gives as given, and where the deposit gives one twice, the
full-precision value: ``k_1``, ``k_2`` and ``k_5`` from the deposited
Arrhenius table, not the four-decimal values of the script's vector ``k0``
(``NOTES.md``, "Constants at 50 C"). The one constant the deposit does not
give a value for, ``k_3`` (``k0`` writes it as 0.0000), is fitted here:
:func:`fit_simplified_reversible_constants` takes the script's data, its
initial condition and its residuals and carries the fit to convergence with
an accurate solver. The default ``k_3`` of the model is its result.

Nothing is loaded or fitted on import.
"""
import numpy as np
from scipy.optimize import least_squares

from .pd_c_3hp_oxidation import (
    NORMALIZED, SIMPLIFIED_REVERSIBLE_PARAMETER_SETS,
    load_o2_simplified_reversible_model, simulate_at_times)
from .reference import read_heatmap_script_data

__all__ = ('fit_simplified_reversible_constants',)

# Starting value of a constant that the script's vector k0 writes as zero
# (k_3, 0.0000, i.e. below 5e-5): a fit cannot be scaled by 0.
_START_IF_ZERO = 1e-5


def _model_parameter_name(name):
    """``'k_3'`` -> ``'k_3_ref'``; a name ending in ``_ref`` is kept."""
    return name if name.endswith('_ref') else name + '_ref'


def fit_simplified_reversible_constants(free=('k_3',), values=None,
                                        start=None, data=None):
    """Fit constants of the simplified reversible model to the 50 C data.

    The model (:func:`.load_o2_simplified_reversible_model`) is run as the
    authors' heat-map script runs it when it fits: at its reference
    temperature, 50 C, from an initial 3-HP "mole fraction" of 0.000339
    (the ``'fitting'`` initial condition). The residuals are the
    differences between the five water-free liquid fractions ``X_HP3,
    X_MA, X_AAc, X_OXA, X_D`` of the model and of the data at the sample
    times, all of them, unweighted (for the script's table 16 x 5 = 80
    values, the t = 0 row included, where they are zero). The constants in
    ``free`` are varied to minimize the sum of their squares; the others
    are held.

    Parameters
    ----------
    free : str or sequence of str, optional
        Names of the 50 C constants to fit, out of ``'k_1'`` to ``'k_5'``
        (rate constants, 1/min) and ``'K_6'`` to ``'K_10'`` (equilibrium
        constants); the names of the model parameters, ``'k_1_ref'`` and so
        on, are accepted as well. A single name may be given as a string.
        Default ``('k_3',)``: the one constant for which the deposit gives
        no value (the script writes it as 0.0000).
    values : dict, optional
        Values of constants, keyed by either form of their names. Constants
        not in ``free`` are held at these values, and a constant in
        ``free`` is fitted starting from its value here (unless ``start``
        is given). Constants not given take the value of the ``'deposit'``
        set of :data:`.SIMPLIFIED_REVERSIBLE_PARAMETER_SETS`.
    start : sequence of float, optional
        Starting values of the free constants, in the order of ``free``;
        all finite and positive. By default the value in ``values`` if
        there is one, else that of the vector ``k0`` written in the script
        (``'k0_printed'``), with 1e-5 where it is zero (``k_3``).
    data : pandas.DataFrame, optional
        50 C samples to fit: ``time`` (min, starting at 0) and the five
        fractions. By default the table of the script
        (:func:`.read_heatmap_script_data`). The model is run at 50 C
        whatever the data: samples taken at another temperature give the
        constants at that temperature, which belong on a model loaded at
        the default ``T`` (setting them on a model loaded at that
        temperature would apply the Arrhenius factors twice).

    Returns
    -------
    dict
        ``'fitted'``
            The fitted constants, ``{name: value}`` keyed by the model
            parameter names (``'k_3_ref'`` ...).
        ``'standard_errors'``
            Their standard errors, same keys: the square roots of the
            diagonal of ``s2*inv(J.T @ J)``, with ``J`` the Jacobian of the
            residuals at the optimum and ``s2`` the sum of squares over
            the number of residuals minus the number of fitted constants
            (the script divides by twice the number of residuals minus
            ten). A linearized measure of how well the data determine each
            constant with the others held; ``nan`` where ``J.T @ J`` is
            singular or nearly so.
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
        If ``free`` is empty or names a constant twice, if a name in
        ``free`` or ``values`` is not a constant of the model (the message
        gives it as it was passed), or if ``start`` does not have one
        finite, positive value per free constant.

    Notes
    -----
    This is not the authors' procedure, which cannot be repeated (see the
    module docstring). It differs from it in four ways: the least-squares
    iteration is carried to convergence; the integration is accurate
    (CVODE at a relative tolerance of 1e-12, where the script uses
    ``ode23s`` at its default 1e-3) and runs through all sample times at
    once, where the script restarts its solver at every sample; by default
    only ``k_3`` is free, where the script frees all ten constants; and it
    is deterministic: the same call returns the same numbers. The constants
    are fitted as multiples of their starting values, bounded below by zero
    as in the script, with a central-difference (three-point) Jacobian;
    with scipy's default two-point Jacobian the iteration stops early
    (``k_3`` 2.8e-4 of its value lower from the default start, up to
    3.6e-3 away from others). From other starting points the optimum is
    found again to 1e-5 of ``k_3``. With the defaults the fit takes less
    than a second and returns ``k_3`` = 1.4215e-5 1/min, the value of the
    ``'deposit'`` set, with a standard error of 31 %.
    """
    names = tuple(SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['deposit'])
    if isinstance(free, str):
        free = (free,)
    given = list(free) + list(values or {})
    unknown = [name for name in given
               if _model_parameter_name(name) not in names]
    if unknown:
        raise ValueError(
            f'unknown constants {unknown} of the simplified reversible '
            f'model: choose from {list(names)} (or the same names without '
            f'"_ref")')
    if not free:
        raise ValueError('free must name at least one constant')
    if len({_model_parameter_name(name) for name in free}) != len(free):
        raise ValueError(f'free must name each constant once: {list(free)}')
    if values and len({_model_parameter_name(name)
                       for name in values}) != len(values):
        raise ValueError(
            f'values must give each constant once: {list(values)}')
    free = [_model_parameter_name(name) for name in free]
    values = {_model_parameter_name(name): value
              for name, value in (values or {}).items()}
    constants = dict(SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['deposit'])
    constants.update(values)
    if start is None:
        printed = SIMPLIFIED_REVERSIBLE_PARAMETER_SETS['k0_printed']
        start = [values[name] if name in values
                 else printed[name] or _START_IF_ZERO for name in free]
    start = np.asarray(start, dtype=float)
    if (start.shape != (len(free),) or not np.isfinite(start).all()
            or not (start > 0.).all()):
        raise ValueError(
            'start must hold one finite, positive value per free constant, '
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

    # The three-point Jacobian is needed for the iteration to converge to
    # the fifth digit of k_3 (see the Notes above).
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
