# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Kinetic models for the selective oxidation of 3-hydroxypropionic acid (3-HP)
to malonic acid (MA) over 10 % Pd/C, from S. Y. Kim, J. E. Chavarrio,
J. Install and G. W. Huber, Appl. Catal. B: Environ. Energy 386 (2026)
126403 (https://doi.org/10.1016/j.apcatb.2026.126403).

Two models are provided, each as an Antimony file next to this module:

* the H2O2 model (Table 3, H2O2 column), which reproduces Fig. 9a. Its
  predicted H2O2 is compared with SI Fig. S12a, whose H2O2 curve lies up
  to 0.38 equivalents lower and is matched only with other constants;
* the O2 model (Table 3, O2 column), which reproduces Fig. 9b and SI
  Fig. S12b, with an Arrhenius temperature dependence that reduces to the
  published 50 C model at the default ``T = 323.15`` K. The temperature
  dependence is an NSKinetics addition: away from 50 C it does not
  reproduce the experiments at 25, 75 and 100 C.

The rate laws of the two Table 3 models have the form of SI Fig. S2 (rate
expressions), ``r_j = k_jf*prod(a_i^n_ij)``, with every order ``n`` equal
to 1, unit activity coefficients and the liquid activities on the authors'
mole-fraction scale. The O2 model departs from Fig. S2 in one place: the gas
activity ``phi*y*P/f0`` of O2 is replaced by the inferred constant
``x_O2_ref`` times ``yO2``. ``NOTES.md`` (section "SI Fig. S2 and Eq. S1")
maps both models to Fig. S2 term by term.

Time is in min, liquid concentrations are in M and temperatures are in K
(arguments and parameters named ``T``) unless a name ends in ``_C``
(degrees Celsius). Nothing is loaded or simulated on import: call
:func:`load_h2o2_model` or :func:`load_o2_model`.
``NOTES.md`` in this folder records what is printed in the paper, what was
inferred from the authors' data deposit, and where the two disagree.

The model curves the authors deposited are solver steps and lie up to
2.0e-3 from the solutions of these models. ``solver_emulation``
repeats that solver run with a model loaded here; the constant sets
:data:`H2O2_NODE_PARAMETER_SETS` and :data:`O2_NODE_PARAMETER_SETS`, which
the loaders of the two Table 3 models accept by name, exist for that
comparison only.
"""
import os

import numpy as np
import pandas as pd
import tellurium as te

import nskinetics as nsk

__all__ = ('load_h2o2_model', 'load_o2_model', 'simulate_at_times',
           'h2o2_antimony_filename', 'o2_antimony_filename',
           'SPECIES', 'NORMALIZED', 'TABLE_S4',
           'H2O2_NODE_PARAMETER_SETS', 'O2_NODE_PARAMETER_SETS',)

#%% Model files
h2o2_antimony_filename = 'pd_c_3hp_oxidation_h2o2_antimony.txt'
o2_antimony_filename = 'pd_c_3hp_oxidation_o2_antimony.txt'

#: Liquid species of the two Table 3 models (concentrations in M), in the
#: order 3-HP, 3-oxopropanoic acid, malonic acid, acetic acid, oxalic acid
#: and the lumped other by-products.
SPECIES = ('HP3', 'OPA', 'MA', 'AAc', 'OXA', 'D')

#: Assignment rules holding the "Normalized Concentration" plotted in Fig. 9:
#: each measured species over the sum of the five measured species (OPA is
#: not in the sum).
NORMALIZED = ('X_HP3', 'X_MA', 'X_AAc', 'X_OXA', 'X_D')

#: SI Table S4 ("Calculated reaction orders and equilibrium constant"), kept
#: for documentation only. Neither the orders ``n`` nor the equilibrium
#: constants ``Keq`` enter the rate laws. The paper itself treats all
#: reactions as irreversible because the ``Keq`` are large (Section 3.4).
#: That the orders are exactly 1 is inferred: the model curves the authors
#: deposited are reproduced by first-order steps and not with these orders,
#: and the authors' posted script has no orders either (see ``NOTES.md``).
#: ``None`` marks a dash in the table. The order of H2O2 reaction 7 is
#: cited from the literature there.
TABLE_S4 = {
    'H2O2': {
        'reaction': ('3-HP + H2O2 -> OPA + 2 H2O',
                     'OPA + H2O2 -> MA + H2O',
                     'MA -> AAc + CO2',
                     'MA + H2O2 -> OXA + 2 H2 + CO2',
                     '3-HP -> D',
                     'MA -> D',
                     'H2O2 -> H2O + 0.5 O2'),
        'n': (0.99, 0.99, 1.49, 1.14, 1.02, 1.01, 1.0),
        'Keq': (3.44e47, 1.53e59, 5.75e4, 5.76e27, None, None, 2.52e17),
        },
    'O2': {
        'reaction': ('3-HP + 0.5 O2 -> OPA + H2O',
                     'OPA + 0.5 O2 -> MA',
                     'MA -> AAc + CO2',
                     'MA + 1.5 O2 -> OXA + H2O + CO2',
                     '3-HP -> D',
                     'MA -> D'),
        'n': (0.99, 0.95, 0.78, 0.87, 0.98, 1.01),
        'Keq': (1.37e30, 6.09e41, 5.75e4, 2.44e92, None, None),
        },
    }

#: Constant sets of the H2O2 model that return the deposited model nodes
#: when the authors' solver run is repeated (``solver_emulation``). They are
#: **not** the published constants: those are Table 3 as printed, the default
#: of the Antimony file and of :func:`load_h2o2_model`. The sets were
#: inferred from the deposited nodes by fits that are not in the repository
#: and are given to eight significant digits, which the deposited time
#: grids (printed to five decimals) need. ``k_1`` to ``k_6`` in
#: 1/(g_cat min), ``k_7`` in 1/min; ``k_6`` is not identifiable and is kept
#: as printed.
#:
#: ``'fig9a_nodes'``
#:     Reproduces the 24 nodes of the deposited Fig. 9a model curve. Over
#:     Table 3 as printed: 1.00242, 0.99399, 1.00522, 0.99756, 0.99827, 1
#:     and 1.00129 for ``k_1`` to ``k_7``, more than rounding to three
#:     digits allows for ``k_1``, ``k_2``, ``k_3``, ``k_5`` and ``k_7``.
#: ``'figS12a_h2o2_nodes'``
#:     Reproduces the 24 nodes of the deposited H2O2 trace of SI Fig. S12a,
#:     which is not the H2O2 of the Fig. 9a curve. Over Table 3 as printed:
#:     1.04354, 1.05228, 0.91569, 0.96122, 0.78103, 1 and 1.07045.
H2O2_NODE_PARAMETER_SETS = {
    name: dict(zip(('k_1', 'k_2', 'k_3', 'k_4', 'k_5', 'k_6', 'k_7'), values))
    for name, values in (
        ('fig9a_nodes', (5.3228596, 30.913113, 7.3582251e-4, 1.1471899,
                         0.0030247652, 4.10e-9, 0.0064883315)),
        ('figS12a_h2o2_nodes', (5.5412148, 32.725866, 6.7028262e-4,
                                1.1053986, 0.0023665310, 4.10e-9,
                                0.0069365260)),
        )}

#: Constant set of the O2 model that returns the 20 nodes of the deposited
#: Fig. 9b model curve when the authors' solver run is repeated
#: (``solver_emulation``): ``k_1_ref`` to ``k_6_ref`` in 1/(g_cat min), for
#: ``m_cat`` = 0.2 g and ``x_O2_ref`` = 0.005085 (only the products
#: ``k_i*m_cat*x_O2_ref`` of reactions 1, 2 and 4 and ``k_i*m_cat`` of the
#: others are determined). It is **not** the published set: that is Table 3
#: as printed, the default of the Antimony file and of
#: :func:`load_o2_model`. The set was inferred from the deposited nodes by a
#: fit that is not in the repository and is given to eight significant
#: digits. Over Table 3 as printed: 1.00007, 0.99770, 1.00028, 1.00013,
#: 0.99743 and 1.00023, all within the rounding of Table 3 to three digits.
O2_NODE_PARAMETER_SETS = {
    'fig9b_nodes': dict(zip(
        ('k_1_ref', 'k_2_ref', 'k_3_ref', 'k_4_ref', 'k_5_ref', 'k_6_ref'),
        (7.1405267, 138.68010, 3.1608819e-4, 0.058607778, 0.011270915,
         0.0055212537))),
    }

# Range of temperatures in K that the loader accepts. It only guards against
# a temperature given in degrees Celsius; it is not a range of validity.
_T_RANGE = (250., 450.)


#%% Loaders
def _load(antimony_filename, values=None):
    """Build a model from an Antimony file, optionally with other values.

    ``values`` maps parameter names to the values a loader selects. They
    are written into the model text as further assignments, so that they
    are the initial values of the model built from it and not changes made
    to a loaded model: ``compile_events()``, which regenerates a model from
    its own definition, therefore keeps them.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        antimony_filename)
    with open(path, encoding='utf-8') as f:
        text = f.read()
    if values:
        assignments = ''.join(f'  {name} = {float(value)!r};\n'
                              for name, value in values.items())
        # A later assignment replaces an earlier one; these go last in the
        # model block, before the line that closes it. (Antimony keeps 15
        # significant digits of a number.)
        body, closing = text.rsplit('\nend\n', 1)
        text = body + '\n' + assignments + 'end\n' + closing
    r = te.loada(text)
    te_r = nsk.KineticModel(r, units={'time': 'min', 'conc': 'M'})
    integrator = r.integrator
    integrator.absolute_tolerance = 1e-10
    integrator.relative_tolerance = 1e-9
    return te_r


def _check_T(T):
    """Raise ValueError for a temperature that cannot be meant in kelvin."""
    if not _T_RANGE[0] <= T <= _T_RANGE[1]:
        raise ValueError(
            f'T = {T!r} is outside {_T_RANGE[0]:g} to {_T_RANGE[1]:g} K. '
            'T is in kelvin: for 75 C pass T=348.15 (75 + 273.15).')


def _table3_or_node_set(parameter_set, node_sets, model):
    """Values of a constant set of a Table 3 model: none for ``'table3'``
    (the Antimony file as it is), else the named set of ``node_sets``."""
    if parameter_set == 'table3':
        return {}
    try:
        return dict(node_sets[parameter_set])
    except (KeyError, TypeError):
        names = ', '.join(repr(name) for name in ('table3', *node_sets))
        raise ValueError(
            f'unknown parameter set {parameter_set!r} of the {model} model: '
            f'choose from {names}') from None


def load_h2o2_model(parameter_set='table3'):
    """Load the H2O2 model (Table 3, H2O2 column; Fig. 9a) at 50 C.

    Reproduces Fig. 9a. Its predicted H2O2 (``eq_H2O2``) is compared with
    SI Fig. S12a, whose H2O2 curve lies up to 0.38 equivalents below this
    model and is matched only with other constants (see ``NOTES.md``).
    Source: the paper (constants of Table 3 as printed; rate form of SI
    Fig. S2 with all orders equal to 1). How H2O2 enters the rates was
    inferred from the deposited model curve.

    Parameters
    ----------
    parameter_set : str, optional
        ``'table3'`` (default): the constants of Table 3 as printed, the
        published model. ``'fig9a_nodes'`` or ``'figS12a_h2o2_nodes'``: a
        key of :data:`H2O2_NODE_PARAMETER_SETS`, constants inferred from
        the deposited nodes with which the repeated solver run of the
        authors (``solver_emulation``) returns the deposited Fig. 9a curve
        or the deposited H2O2 trace to their printed digits. They are for
        that verification only.

    Returns
    -------
    nskinetics.KineticModel
        A new model with 60 mM 3-HP, 15 equivalents (0.9 M) of H2O2 and
        0.2 g of catalyst. Species: ``HP3, OPA, MA, AAc, OXA, D, H2O2``.
        Plotted quantities (select them by name): ``X_HP3, X_MA, X_AAc,
        X_OXA, X_D`` (Fig. 9a), ``eq_H2O2`` (SI Fig. S12a inset) and
        ``x_OPA``.

    Raises
    ------
    ValueError
        If ``parameter_set`` is not one of the three names; the message
        lists them.

    Notes
    -----
    Parameters meant to be changed with
    :meth:`~nskinetics.KineticModel.set_value` (all positive):

    ``m_cat``
        Catalyst mass **of the paper's batch** in g of 10 % Pd/C (default
        0.2, the value of the Fig. 9 caption). It is not a mass per litre
        of the model's compartment, and the model has no volume term: the
        rates of reactions 1 to 6 are proportional to ``m_cat``, so
        ``m_cat = 0.4`` stands for twice the catalyst loading of the paper,
        whatever the liquid volume. The paper's batch is 0.2 g in 10 mL of
        liquid (20 g/L); it adds that substrate, NaOH and water were
        tripled for time-resolved runs and does not say whether the
        catalyst was, so the loading of the kinetic runs is 20 g/L or
        6.7 g/L (Section 2.3; see ``NOTES.md``).
    ``C0``
        Initial 3-HP concentration in M (default 0.06).
    ``eq_H2O2_0``
        Initial H2O2 in equivalents of ``C0`` (default 15, i.e. 0.9 M).
        The initial H2O2 concentration is ``eq_H2O2_0*C0``, so changing
        ``C0`` alone also changes it.
    ``k_1`` to ``k_7``
        Rate constants: ``k_1`` to ``k_6`` in 1/(g_cat min), ``k_7`` in
        1/min.

    ``C0`` and ``eq_H2O2_0`` define the initial state and take effect at
    the next ``reset()`` (:func:`simulate_at_times` resets). The constants
    of Table 3 were fitted to one experiment (60 mM 3-HP, 0.2 g of
    catalyst, 15 equivalents of H2O2, 50 C, samples up to 1440 min); any
    other value of these parameters is an extrapolation. For ``eq_H2O2_0``
    the paper's own data speak against it: in its 8 h screening experiments
    (Fig. 2a) the conversion falls above 15 equivalents, whereas in this
    model it rises with every equivalent added (see ``NOTES.md``).
    ``set_value`` does not check its input: zero or negative values of
    ``m_cat``, ``C0`` or ``eq_H2O2_0`` give meaningless results without an
    error.

    Each call loads a new, independent model. Read and set model values
    through :meth:`~nskinetics.KineticModel.get_value` and
    :meth:`~nskinetics.KineticModel.set_value` (item access), not through
    attributes of the underlying roadrunner object: attribute access is
    unreliable once more than one model has been loaded in a process.
    ``reset()`` restores the species to the initial state and keeps
    parameter values set with ``set_value``.
    :meth:`~nskinetics.KineticModel.compile_events` regenerates the model
    and discards every value set with ``set_value`` before it: set values
    after compiling events. A constant set selected with ``parameter_set``
    is part of the model that is built and is kept by both.
    """
    return _load(h2o2_antimony_filename, _table3_or_node_set(
        parameter_set, H2O2_NODE_PARAMETER_SETS, 'H2O2'))


def load_o2_model(T=None, parameter_set='table3'):
    """Load the O2 model (Table 3, O2 column; Fig. 9b).

    Reproduces Fig. 9b and SI Fig. S12b at the default 50 C. Source: the
    paper (constants of Table 3 as printed; rate form of SI Fig. S2 with
    all orders equal to 1), with the factor ``x_O2_ref``, which stands in
    for the O2 activity of Fig. S2, inferred from the deposited model
    curve; the temperature dependence is an NSKinetics addition.

    Parameters
    ----------
    T : float, optional
        Reaction temperature in K (for 75 C pass ``T=348.15``). Defaults
        to the model's reference temperature, 323.15 K (50 C), at which
        the model is the published one. Other temperatures scale the six
        rate constants by Arrhenius factors (the Arrhenius relation of the
        paper, the second of its two equations numbered 6, with the
        activation energies of the authors' heat-map script); see the
        warning below.
    parameter_set : str, optional
        ``'table3'`` (default): the constants of Table 3 as printed, the
        published model. ``'fig9b_nodes'``, the key of
        :data:`O2_NODE_PARAMETER_SETS`: the 50 C constants inferred from
        the deposited nodes, with which the repeated solver run of the
        authors (``solver_emulation``) returns the deposited Fig. 9b curve
        to its printed digits. It is for that verification only.

    Returns
    -------
    nskinetics.KineticModel
        A new model with 60 mM 3-HP and 0.2 g of catalyst. Species:
        ``HP3, OPA, MA, AAc, OXA, D``; ``yO2`` is the dimensionless O2
        activity. Plotted quantities (select them by name): ``X_HP3, X_MA,
        X_AAc, X_OXA, X_D`` (Fig. 9b) and ``x_OPA`` (SI Fig. S12b inset).

    Raises
    ------
    ValueError
        If ``T`` is outside 250 to 450 K, which is what a temperature
        given in degrees Celsius looks like. The range is a guard against
        that mistake only, not a range of validity. Also if
        ``parameter_set`` is not one of the two names; the message lists
        them.

    Warnings
    --------
    Only the 50 C model is a published one. Away from 50 C the model is an
    extrapolation: it is too fast for the experiments at 25, 75 and 100 C
    (pooled R2 0.579, 0.644 and 0.746). See ``NOTES.md`` and the
    documentation page of this package.

    Notes
    -----
    Parameters meant to be changed with
    :meth:`~nskinetics.KineticModel.set_value`:

    ``m_cat``
        Catalyst mass **of the paper's batch** in g of 10 % Pd/C (default
        0.2); all six rates are proportional to it. It is not a mass per
        litre; see :func:`load_h2o2_model`.
    ``C0``
        Initial 3-HP concentration in M (default 0.06); takes effect at
        the next ``reset()``. The plotted fractions ``X_i`` hardly depend
        on it, because every rate is first order in its liquid species:
        halving or doubling ``C0`` changes them by less than 1e-4, through
        the O2 depletion.
    ``T``
        Temperature in K (default 323.15). ``set_value`` does not check
        it: the guard of this loader applies to its argument only.
    ``k_1_ref`` to ``k_6_ref``, ``Ea_1`` to ``Ea_6``
        Rate constants at 50 C in 1/(g_cat min) and activation energies in
        J/mol. ``k_1`` to ``k_6`` are computed from them and cannot be
        set.
    ``x_O2_ref``
        The inferred oxidant factor (default 0.005085, for 3 bar of O2).
        The model has no O2 pressure parameter, and this factor is not
        one: scaling it makes reactions 1, 2 and 4 proportionally faster,
        whereas in the paper's screening experiments conversion and MA
        selectivity fall above 3 bar (Fig. 3a; see ``NOTES.md``).

    The constants of Table 3 were fitted to one experiment (60 mM 3-HP,
    0.2 g of catalyst, 3 bar O2, 50 C, samples up to 1680 min); any other
    value of these parameters is an extrapolation. For ``T`` the
    comparison with the experiments of the paper is poor (see the warning
    above; ``NOTES.md`` also has the paper's 8 h screening experiments,
    which the models were not built on). ``set_value`` does not check its
    input: zero or negative values of ``m_cat`` or ``C0`` give meaningless
    results without an error.

    The temperature and the constant set given here are part of the model
    that is built, so ``reset()`` and
    :meth:`~nskinetics.KineticModel.compile_events` keep them. A value set
    later with ``set_value`` is kept by ``reset()``, which restores only
    the species and ``yO2``, and is discarded by ``compile_events()``,
    which regenerates the model: set values after compiling events. See
    :func:`load_h2o2_model` on item access.
    """
    values = _table3_or_node_set(parameter_set, O2_NODE_PARAMETER_SETS, 'O2')
    if T is not None:
        _check_T(T)
        values['T'] = T
    return _load(o2_antimony_filename, values)


#%% Simulation at given times
def simulate_at_times(te_r, times, selections=None):
    """Reset a model and return its trajectory at exactly the given times.

    Parameters
    ----------
    te_r : nskinetics.KineticModel
        Model from one of the loaders of this module.
    times : array-like
        Increasing times in min. The integration starts at t = 0 whether
        or not 0 is listed.
    selections : list of str, optional
        Quantities to record. Defaults to the Fig. 9 normalized
        concentrations (:data:`NORMALIZED`).

    Returns
    -------
    pandas.DataFrame
        Column ``time`` followed by one column per selection, one row per
        requested time.
    """
    times = np.asarray(times, dtype=float)
    if selections is None:
        selections = list(NORMALIZED)
    start_at_zero = times[0] == 0.
    grid = times if start_at_zero else np.concatenate(([0.], times))
    te_r.reset()
    result = te_r._te.simulate(times=list(grid),
                               selections=['time'] + list(selections))
    df = pd.DataFrame(np.array(result), columns=['time'] + list(selections))
    if not start_at_zero:
        df = df.iloc[1:].reset_index(drop=True)
    return df
