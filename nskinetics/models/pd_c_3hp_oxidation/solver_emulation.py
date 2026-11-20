# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Emulation of the solver run behind the model curves that the authors of
Kim et al., Appl. Catal. B: Environ. Energy 386 (2026) 126403, deposited
for Fig. 9a, Fig. 9b and SI Fig. S12a (H2O2).

The deposited model curves are not converged solutions: they have 19 or
23 unevenly spaced solver steps (20 or 24 nodes with the one at t = 0) that
lie up to 2.0e-3 from the solution, and they are reproduced step for step
by the method of MATLAB ``ode23s`` at its default tolerances, the solver
the paper names. An accurate solver can therefore reproduce them only to
that error. This module repeats the authors' run instead: it integrates
the state vector of their posted script with the same method and step-size
control, taking every reaction rate from a shipped model. With the constant
sets named in the functions below it returns the deposited nodes (times and
values) to the digits printed in the deposited files. It is used for
verification only (the package tests and ``NOTES.md``); to simulate a
model, use the model.

:func:`ode23s` is a port of the modified Rosenbrock (2,3) pair of
L. F. Shampine and M. W. Reichelt, "The MATLAB ODE Suite", SIAM J. Sci.
Comput. 18 (1997) 1-22 (https://doi.org/10.1137/S1064827594276424), written
from that paper and from the documented behaviour of the MATLAB function.
It covers what these runs need (an autonomous system, default options, a
finite-difference Jacobian, output at every step) and is not a
general-purpose solver. MATLAB itself was not run.

Nothing is loaded or simulated on import.
"""
import numpy as np
import pandas as pd

from .pd_c_3hp_oxidation import (
    NORMALIZED, SPECIES, load_h2o2_model, load_o2_model)
from .reference import read_fig9a_h2o2, read_fig9b_o2

__all__ = ('ode23s', 'emulate_authors_run', 'emulate_fig9a', 'emulate_fig9b',
           'node_deviations',)


#%% Solver
def _jacobian(fun, y, f0, threshold):
    """Forward-difference Jacobian of ``fun`` at ``y``, where it is ``f0``.

    The increment of state ``j`` is ``sqrt(eps)*max(|y_j|, threshold)``,
    the starting increment of MATLAB's ``numjac``. (``numjac`` goes on to
    adapt its increments; that is not repeated here. The right-hand sides
    of these runs are linear or bilinear in the states, so that a forward
    difference is exact up to rounding.)
    """
    jac = np.empty((y.size, y.size))
    for j in range(y.size):
        step = np.sqrt(np.finfo(float).eps)*max(abs(y[j]), threshold)
        shifted = y.copy()
        shifted[j] += step
        jac[:, j] = (fun(shifted) - f0)/step
    return jac


def ode23s(fun, t_end, y0, rtol=1e-3, atol=1e-6, max_step=None):
    """Integrate ``y' = fun(y)`` from 0 to ``t_end`` as MATLAB ``ode23s`` does.

    A port of the modified Rosenbrock (2,3) pair of Shampine and Reichelt
    (SIAM J. Sci. Comput. 18 (1997) 1-22) with the step-size control of
    the MATLAB function, for the emulation of the authors' runs; see the
    module docstring. It is not a general-purpose solver: the system must
    be autonomous, the Jacobian is always taken by finite differences, and
    every step is returned (as MATLAB does for ``tspan = [0 t_end]``).

    Parameters
    ----------
    fun : callable
        Right-hand side ``fun(y)``, returning an array like ``y``.
    t_end : float
        End of the integration, which starts at 0.
    y0 : array-like
        Initial state.
    rtol : float, optional
        Relative tolerance, by default 1e-3 (MATLAB ``RelTol``).
    atol : float, optional
        Absolute tolerance, by default 1e-6 (MATLAB ``AbsTol``).
    max_step : float, optional
        Largest step, by default ``0.1*t_end`` (MATLAB ``MaxStep``).

    Returns
    -------
    t : numpy.ndarray
        The times of the steps taken, from 0 to ``t_end``.
    y : numpy.ndarray
        The states at those times, one row per time.
    info : dict
        ``'first_step'``: the first step size tried; ``'failed'``: one
        tuple (time, step size, error over ``rtol``) per rejected step;
        ``'error'``: the error estimate over ``rtol`` of each step taken
        (a step is taken when it is at most 1).

    Notes
    -----
    One step from ``y`` with ``W = I - h*d*J``, ``d = 1/(2 + sqrt(2))``,
    ``e32 = 6 + sqrt(2)`` and ``J`` the Jacobian at ``y``::

        k1 = W \\ f(y)
        k2 = W \\ (f(y + h*k1/2) - k1) + k1
        y_new = y + h*k2
        k3 = W \\ (f(y_new) - e32*(k2 - f(y + h*k1/2)) - 2*(k1 - f(y)))
        error = h/6*(k1 - 2*k2 + k3)

    The step is taken if the largest component of ``error``, each divided
    by ``max(|y|, |y_new|, atol/rtol)``, is at most ``rtol``. (The
    deposited runs do not test the ``|y_new|`` term of this weight: without
    it all three take the same steps. The package tests check it on
    ``y' = y``.) Step-size rules:

    * first step: ``1/h = 1.25*sqrt(0.5*max(|J f|/w))/rtol**(1/3)`` with
      ``w = max(|y0|, atol/rtol)``, at most ``max_step``;
    * after a step without a rejection the next step is the last one
      divided by ``1.25*(error/rtol)**(1/3)``, at most five times it (in
      the three deposited runs the error is small enough for this limit
      only in the last step and in steps that are already of the largest
      size, ``max_step``, a tenth of the run; the package tests check it on
      ``y' = -100 y``);
    * a rejected step is retried with the step multiplied by
      ``max(0.1, 0.8*(rtol/error)**(1/3))``, a second rejection onwards
      with half the step, and the step size is then kept for the next
      step (no step of the three deposited runs is rejected; the package
      tests check these rules on ``y' = cos(10 t)``);
    * a step that would end within 10 % of its own length before
      ``t_end`` is stretched to end there (the last step of each of the
      three deposited runs is cut short at ``t_end`` instead; the package
      tests check the stretching on ``y' = -y``).
    """
    y = np.array(y0, dtype=float)
    identity = np.eye(y.size)
    threshold = atol/rtol
    if max_step is None:
        max_step = 0.1*t_end
    d = 1./(2. + np.sqrt(2.))
    e32 = 6. + np.sqrt(2.)

    t = 0.
    f0 = fun(y)
    jac = _jacobian(fun, y, f0, threshold)
    # First step, from an estimate of the second derivative, J f.
    weight = np.maximum(np.abs(y), threshold)
    rh = 1.25*np.sqrt(0.5*np.max(np.abs(jac @ f0)/weight))/rtol**(1./3.)
    absh = min(max_step, t_end)
    if absh*rh > 1.:
        absh = 1./rh
    info = {'first_step': absh, 'failed': [], 'error': []}

    times, states = [t], [y]
    done = False
    while not done:
        hmin = 16.*np.spacing(t)
        absh = min(max_step, max(hmin, absh))
        h = absh
        if 1.1*absh >= t_end - t:
            h = absh = t_end - t
            done = True
        rejected = False
        while True:
            W = identity - h*d*jac
            k1 = np.linalg.solve(W, f0)
            f1 = fun(y + 0.5*h*k1)
            k2 = np.linalg.solve(W, f1 - k1) + k1
            y_new = y + h*k2
            f2 = fun(y_new)
            k3 = np.linalg.solve(W, f2 - e32*(k2 - f1) - 2.*(k1 - f0))
            scale = np.maximum(np.maximum(np.abs(y), np.abs(y_new)),
                               threshold)
            error = absh/6.*np.max(np.abs(k1 - 2.*k2 + k3)/scale)
            if error <= rtol:
                break
            info['failed'].append((t, absh, error/rtol))
            if absh <= hmin:
                raise RuntimeError(
                    f'ode23s: step size below {hmin:g} at t = {t:g}')
            if rejected:
                absh = max(hmin, 0.5*absh)
            else:
                absh = max(hmin, absh*max(0.1, 0.8*(rtol/error)**(1./3.)))
                rejected = True
            h = absh
            done = False
        info['error'].append(error/rtol)
        if not rejected:
            factor = 1.25*(error/rtol)**(1./3.)
            absh = absh/factor if factor > 0.2 else 5.*absh
        t = t_end if done else t + h
        # The last function value of a step is the first of the next one.
        y, f0 = y_new, f2
        times.append(t)
        states.append(y)
        if not done:
            jac = _jacobian(fun, y, f0, threshold)
    return np.array(times), np.array(states), info


#%% The authors' system, with the rates of a shipped model
# The authors' posted O2 script integrates ten states: the liquid mole
# fractions of 3-HP, OPA, MA, AAc, OXA and D, then H2O, H2, O2 and CO2. The
# O2 state is ``yO2`` of the O2 model. H2O, H2 and CO2, which that model
# does not have, are integrated here as the script does: the rate of change
# is the row of the script's stoichiometric matrix (reaction: coefficient)
# applied to the model's reaction rates. These states do not feed back into
# any rate; they are kept because the error estimate of the solver, which
# sets the step sizes, runs over all states.
_SCRIPT_ROWS = {'H2O': {'R1': 1.},
                'H2': {'R4': 1.},
                'CO2': {'R3': 1., 'R4': 1.}}


def _state_names(te_r):
    """Names of the states of the authors' run for a model, in their order.

    A name is a species of the model (integrated as its mole fraction on
    the authors' scale, ``x_per_M`` times the concentration), ``'yO2'``
    (the O2 state of the O2 model) or a key of ``_SCRIPT_ROWS``.
    """
    species = tuple(te_r._te.model.getFloatingSpeciesIds())
    variables = te_r._te.model.getGlobalParameterIds()
    if species == SPECIES + ('H2O2',):
        # The authors posted the O2 version of their script only. The H2O2
        # run is taken as the six liquid species and H2O2.
        return species
    if species == SPECIES and 'yO2' in variables:
        return species + ('H2O', 'H2', 'yO2', 'CO2')
    raise ValueError(
        'no deposited curve to emulate for this model: pass a model from '
        'load_h2o2_model or load_o2_model')


def emulate_authors_run(te_r, t_end, rtol=1e-3, atol=1e-6,
                        script_states=True):
    """Repeat the authors' solver run with the rates of a shipped model.

    The state vector of the authors' script (liquid mole fractions on their
    scale, 0.000339 for the 60 mM feed, and the gas states of the script)
    is integrated from the initial state of the model with :func:`ode23s`.
    Every rate of change is evaluated by the model itself at the state the
    solver asks for: the rate laws, the constants, the stoichiometry and,
    for the O2 model, the rate rule of ``yO2`` are those of ``te_r``. The
    plotted quantities are read from the model at each step.

    Parameters
    ----------
    te_r : nskinetics.KineticModel
        A model from :func:`~.pd_c_3hp_oxidation.load_h2o2_model` or
        :func:`~.pd_c_3hp_oxidation.load_o2_model`, with any values set on
        it. It is reset before the run and again after it;
        values set with ``set_value`` are kept.
    t_end : float
        End of the run in min. The authors' script integrates to 1.1
        times the last sampling time of the experiment.
    rtol, atol : float, optional
        Tolerances of :func:`ode23s`, by default those of MATLAB.
    script_states : bool, optional
        Whether the states of the posted script that the model does not
        have are integrated along (default). They enter no rate, and
        leaving them out changes no step of the three deposited runs (see
        the Notes).

    Returns
    -------
    pandas.DataFrame
        One row per step of the solver, the initial state included:
        ``time`` (min), the plotted fractions ``X_HP3, X_MA, X_AAc, X_OXA,
        X_D`` and ``x_OPA``, for the H2O2 model also ``eq_H2O2``, and the
        integrated states in columns named ``state_<name>``. The
        information returned by :func:`ode23s` is in ``.attrs['ode23s']``.

    Raises
    ------
    ValueError
        If ``te_r`` is not one of the two models named above (there is
        no deposited curve of the simplified reversible model).

    Notes
    -----
    The integrated states are, in this order: for the O2 model the ten
    states of the posted script, 3-HP, OPA, MA, AAc, OXA, D, H2O, H2, O2 and
    CO2, with ``yO2`` of the model as the O2 state; for the H2O2 model
    the six liquid species and H2O2. States that the model does not have
    (H2O, H2, CO2) follow the stoichiometric matrix of the script and
    do not enter any rate. They could matter only through the error
    estimate of the solver, which runs over all states; in the deposited
    runs they never set it, and with ``script_states=False`` the steps are
    the same to 1e-12.
    """
    r = te_r._te
    names = _state_names(te_r)
    if not script_states:
        names = tuple(name for name in names if name not in _SCRIPT_ROWS)
    species = r.model.getFloatingSpeciesIds()
    reactions = r.model.getReactionIds()
    x_per_M = te_r.get_value('x_per_M')
    # Model quantity behind each state (item access), its value per unit of
    # the state, and the script's row for the states the model lacks.
    selections = [('[%s]' % name, 1./x_per_M) if name in species
                  else (name, 1.) if name == 'yO2' else None
                  for name in names]
    rows = [None if selection is not None
            else np.array([_SCRIPT_ROWS[name].get(reaction, 0.)
                           for reaction in reactions])
            for name, selection in zip(names, selections)]

    def set_state(y):
        for selection, value in zip(selections, y):
            if selection is not None:
                r[selection[0]] = value*selection[1]

    def fun(y):
        set_state(y)
        # Reaction rates in mole fraction per min (the compartment is 1).
        rates = x_per_M*np.array([r[reaction] for reaction in reactions])
        return np.array([
            row @ rates if selection is None
            else r[selection[0].strip('[]') + "'"]/selection[1]
            for selection, row in zip(selections, rows)])

    te_r.reset()
    y0 = np.array([0. if selection is None else r[selection[0]]/selection[1]
                   for selection in selections])
    if 'H2O' in names:
        # Script: the liquid mole fractions add up to 1.
        y0[names.index('H2O')] = 1. - y0[:len(SPECIES)].sum()
    times, states, info = ode23s(fun, t_end, y0, rtol=rtol, atol=atol)

    plotted = list(NORMALIZED) + ['x_OPA']
    if 'H2O2' in names:
        plotted.append('eq_H2O2')
    values = []
    for y in states:
        set_state(y)
        values.append([r[name] for name in plotted])
    te_r.reset()
    result = pd.DataFrame(np.column_stack([times, values, states]),
                          columns=['time'] + plotted
                          + ['state_' + name for name in names])
    result.attrs['ode23s'] = info
    return result


#%% The deposited curves
def _t_end(reference):
    """End of the authors' run: 1.1 times the last sampling time."""
    return 1.1*reference['experiment']['time'].iloc[-1]


def emulate_fig9a(parameter_set='fig9a_nodes'):
    """Repeat the solver run behind the deposited Fig. 9a curve (H2O2).

    Parameters
    ----------
    parameter_set : str, optional
        Constant set of :func:`~.pd_c_3hp_oxidation.load_h2o2_model`.
        ``'fig9a_nodes'`` (default) returns the deposited nodes of the
        five species to their printed digits. ``'figS12a_h2o2_nodes'``
        returns the deposited H2O2 trace of SI Fig. S12a, which is
        reproduced only as a run of its own: its column ``eq_H2O2`` is
        then the trace, and the deposited file holds all of its steps but
        the last. ``'table3'`` shows how close the constants of Table 3 as
        printed come.

    Returns
    -------
    pandas.DataFrame
        The steps of the run from 0 to 1584 min, as returned by
        :func:`emulate_authors_run`.
    """
    return emulate_authors_run(load_h2o2_model(parameter_set),
                               _t_end(read_fig9a_h2o2()))


def emulate_fig9b(parameter_set='fig9b_nodes'):
    """Repeat the solver run behind the deposited Fig. 9b curve (O2, 50 C).

    Parameters
    ----------
    parameter_set : str, optional
        Constant set of :func:`~.pd_c_3hp_oxidation.load_o2_model`.
        ``'fig9b_nodes'`` (default) returns the deposited nodes to their
        printed digits; ``'table3'`` shows how close the constants of
        Table 3 as printed come.

    Returns
    -------
    pandas.DataFrame
        The steps of the run from 0 to 1848 min, as returned by
        :func:`emulate_authors_run`.
    """
    return emulate_authors_run(load_o2_model(parameter_set),
                               _t_end(read_fig9b_o2()))


def node_deviations(emulated, deposited, columns=NORMALIZED):
    """Compare the steps of an emulated run with deposited nodes, row by row.

    Parameters
    ----------
    emulated : pandas.DataFrame
        Steps returned by one of the functions of this module.
    deposited : pandas.DataFrame
        Deposited nodes, a table of ``reference`` (``'model'`` or
        ``'h2o2_trace'``).
    columns : sequence of str, optional
        Quantities to compare, by default the five plotted fractions.

    Returns
    -------
    dict
        ``'steps'`` and ``'deposited_steps'``: the number of steps of the
        two; over the rows both have, ``'time'``: the largest difference
        between the times of corresponding rows in min, and ``'rms'`` and
        ``'max'``: the root mean square and the largest absolute value of
        the differences in ``columns``.
    """
    columns = list(columns)
    n = min(len(emulated), len(deposited))
    difference = (emulated[columns].to_numpy()[:n]
                  - deposited[columns].to_numpy()[:n])
    return {'steps': len(emulated) - 1,
            'deposited_steps': len(deposited) - 1,
            'time': float(np.abs(emulated['time'].to_numpy()[:n]
                                 - deposited['time'].to_numpy()[:n]).max()),
            'rms': float(np.sqrt((difference**2).mean())),
            'max': float(np.abs(difference).max())}
