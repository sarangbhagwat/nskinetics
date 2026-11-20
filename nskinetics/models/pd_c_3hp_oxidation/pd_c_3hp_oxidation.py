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

Three models are provided, each as an Antimony file next to this module:

* the H2O2 model (Table 3, H2O2 column), which reproduces Fig. 9a. Its
  predicted H2O2 is compared with SI Fig. S12a, whose H2O2 curve lies up
  to 0.38 equivalents lower and is matched only with other constants;
* the O2 model (Table 3, O2 column), which reproduces Fig. 9b and SI
  Fig. S12b. It is the published model at 50 C and has no temperature
  parameter;
* the five-reaction O2 model (no OPA) of the authors' deposited heat-map
  script (Fig. 11), called here the simplified reversible model. It is not
  in the paper. Every constant the deposit gives is used as given, and
  where the deposit gives one twice, the full-precision value: three rate
  constants from the 50 C row of the deposited Arrhenius table of
  Fig. 10d (not the four-decimal values of ``k0``), and ``k_4`` and the
  five equilibrium constants from the vector ``k0`` written in the
  script. The one constant without a value in the deposit, ``k_3``
  (``k0`` writes 0.0000), is fitted to the 50 C data of the script
  (``fitting``). The temperature dependence is the script's.

The rate laws of the two Table 3 models have the form of SI Fig. S2 (rate
expressions), ``r_j = k_jf*prod(a_i^n_ij)``, with every order ``n`` equal
to 1, unit activity coefficients and the liquid activities on the authors'
mole-fraction scale. The O2 model departs from Fig. S2 in one place: the gas
activity ``phi*y*P/f0`` of O2 is replaced by the inferred constant
``x_O2_ref`` times ``yO2``. ``NOTES.md`` (section "SI Fig. S2 and Eq. S1")
maps each of the three models to Fig. S2 term by term.

Time is in min, liquid concentrations are in M and temperatures are in K
(arguments and parameters named ``T``) unless a name ends in ``_C``
(degrees Celsius). Nothing is loaded or simulated on import: call
:func:`load_h2o2_model`, :func:`load_o2_model` or
:func:`load_o2_simplified_reversible_model`.
``NOTES.md`` in this folder records what is printed in the paper, what was
inferred from the authors' data deposit, and where the two disagree.

The model curves the authors deposited for Fig. 9 are solver steps and lie
up to 2.0e-3 from the solutions of these models. ``solver_emulation``
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

__all__ = ('load_h2o2_model', 'load_o2_model',
           'load_o2_simplified_reversible_model', 'simulate_at_times',
           'h2o2_antimony_filename', 'o2_antimony_filename',
           'o2_simplified_reversible_antimony_filename',
           'SPECIES', 'NORMALIZED', 'TABLE_S4',
           'H2O2_NODE_PARAMETER_SETS', 'O2_NODE_PARAMETER_SETS',
           'SIMPLIFIED_REVERSIBLE_PARAMETER_SETS',
           'SIMPLIFIED_REVERSIBLE_INITIAL_CONDITIONS',)

#%% Model files
h2o2_antimony_filename = 'pd_c_3hp_oxidation_h2o2_antimony.txt'
o2_antimony_filename = 'pd_c_3hp_oxidation_o2_antimony.txt'
o2_simplified_reversible_antimony_filename = (
    'pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt')

#: Liquid species of the two Table 3 models (concentrations in M), in the
#: order 3-HP, 3-oxopropanoic acid, malonic acid, acetic acid, oxalic acid
#: and the lumped other by-products. The simplified reversible model has
#: the same without OPA.
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
#: (``solver_emulation``): ``k_1`` to ``k_6`` in 1/(g_cat min), for
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
        ('k_1', 'k_2', 'k_3', 'k_4', 'k_5', 'k_6'),
        (7.1405267, 138.68010, 3.1608819e-4, 0.058607778, 0.011270915,
         0.0055212537))),
    }

_SIMPLIFIED_REVERSIBLE_NAMES = (
    'k_1_ref', 'k_2_ref', 'k_3_ref', 'k_4_ref', 'k_5_ref',
    'K_6_ref', 'K_7_ref', 'K_8_ref', 'K_9_ref', 'K_10_ref')

#: Constant sets of the simplified reversible model at 50 C, in the order of
#: the authors' script: five rate constants (1/min) and five equilibrium
#: constants. The constants the script's run fitted were not deposited (see
#: ``NOTES.md``).
#:
#: ``'deposit'``
#:     The default of the Antimony file and of the loader: every constant
#:     the deposit gives, used as given (where it gives one twice, the
#:     full-precision value), and the one it does not give fitted. ``k_1``,
#:     ``k_2`` and ``k_5`` are ``exp()`` of the 50 C row of the deposited
#:     ``11-4_Arrhenius.csv``, the full-precision source of the three (the
#:     vector ``k0`` written in the script has them to four
#:     decimals, 0.0021, 0.0001 and 0.0010; the last disagrees with the
#:     table's 1.0625e-3, which would round to 0.0011). ``k_4`` = 0.0023 (two
#:     significant digits) and ``K_6`` to ``K_10`` are those of ``k0``,
#:     their only source. ``k_3``, which ``k0`` writes as 0.0000 (zero as
#:     written; below 5e-5 if ``k0`` is a four-decimal print of an earlier
#:     value, which the deposit does not state), is the least-squares fit
#:     to the 50 C data table of the script with the other nine held
#:     (``fit_simplified_reversible_constants``), to five significant
#:     digits.
#: ``'k0_printed'``
#:     The vector ``k0`` written in the script, to its four decimals
#:     (``k_3`` reads 0). The script uses it only for the bounds of its
#:     particle swarm (0 to 1.5 times ``k0``); whether it holds the result
#:     of an earlier run is not stated. Its five equilibrium constants and
#:     its ``k_4`` are those of the ``'deposit'`` set. Run as a constant set,
#:     it is far from the published Fig. 11.
SIMPLIFIED_REVERSIBLE_PARAMETER_SETS = {
    name: dict(zip(_SIMPLIFIED_REVERSIBLE_NAMES, values))
    for name, values in (
        ('deposit', (2.139999e-3, 5.199984e-5, 1.4215e-5, 0.0023,
                     1.062501e-3, 4.7840, 0.5918, 3.5116, 81.6870, 52.0472)),
        ('k0_printed', (0.0021, 0.0001, 0.0000, 0.0023, 0.0010,
                        4.7840, 0.5918, 3.5116, 81.6870, 52.0472)),
        )}

#: Initial 3-HP "mole fraction" of the two initial conditions the authors'
#: script uses for the simplified reversible model: ``'heatmap'`` (unit mole
#: fractions, used for Fig. 11) and ``'fitting'`` (0.000339, used to fit the
#: constants to the 50 C data). The model parameter ``x_per_M`` is this value
#: divided by the 60 mM feed concentration.
SIMPLIFIED_REVERSIBLE_INITIAL_CONDITIONS = {'heatmap': 1., 'fitting': 0.000339}

# Names that say what the two initial conditions are, not what the authors'
# script uses them for; accepted by the loader as aliases.
_INITIAL_CONDITION_ALIASES = {'unscaled': 'heatmap', 'dilute': 'fitting'}

# Initial 3-HP concentration held by the three model files, in M.
_C0 = 0.06

# Range of temperatures in K that the loader of the simplified reversible
# model accepts. It only guards against a temperature given in degrees
# Celsius; it is not a range of validity.
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
    """Raise ValueError for a temperature that cannot be meant in kelvin.

    A positive value that is in the range once 273.15 is added (up to
    176.85) is taken to be degrees Celsius, and the message gives it in
    kelvin.
    """
    if not _T_RANGE[0] <= T <= _T_RANGE[1]:
        if 0. < T and _T_RANGE[0] <= T + 273.15 <= _T_RANGE[1]:
            hint = f': for {T:g} C pass T={T + 273.15:g} ({T:g} + 273.15)'
        else:
            hint = ' (degrees Celsius + 273.15)'
        raise ValueError(
            f'T = {T!r} is outside {_T_RANGE[0]:g} to {_T_RANGE[1]:g} K. '
            f'T is in kelvin{hint}.')


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


def load_o2_model(parameter_set='table3'):
    """Load the O2 model (Table 3, O2 column; Fig. 9b) at 50 C.

    Reproduces Fig. 9b and SI Fig. S12b. Source: the paper (constants of
    Table 3 as printed; rate form of SI Fig. S2 with all orders equal to
    1), with the factor ``x_O2_ref``, which stands in for the O2 activity
    of Fig. S2, inferred from the deposited model curve.

    Parameters
    ----------
    parameter_set : str, optional
        ``'table3'`` (default): the constants of Table 3 as printed, the
        published model. ``'fig9b_nodes'``, the key of
        :data:`O2_NODE_PARAMETER_SETS`: the constants inferred from the
        deposited nodes, with which the repeated solver run of the authors
        (``solver_emulation``) returns the deposited Fig. 9b curve to its
        printed digits. It is for that verification only.

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
        If ``parameter_set`` is not one of the two names; the message
        lists them.

    Notes
    -----
    The model is the published one, at 50 C, and has no temperature
    parameter: the paper fitted the constants of Table 3 at that
    temperature only. For other temperatures the package has the authors'
    simplified reversible model
    (:func:`load_o2_simplified_reversible_model`), which is a different
    model with its own limits.

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
    ``k_1`` to ``k_6``
        Rate constants at 50 C in 1/(g_cat min).
    ``x_O2_ref``
        The inferred oxidant factor (default 0.005085, for 3 bar of O2).
        The model has no O2 pressure parameter, and this factor is not
        one: scaling it makes reactions 1, 2 and 4 proportionally faster,
        whereas in the paper's screening experiments conversion and MA
        selectivity fall above 3 bar (Fig. 3a; see ``NOTES.md``).

    The constants of Table 3 were fitted to one experiment (60 mM 3-HP,
    0.2 g of catalyst, 3 bar O2, 50 C, samples up to 1680 min); any other
    value of these parameters is an extrapolation (``NOTES.md`` has the
    paper's 8 h screening experiments, which the models were not built
    on). ``set_value`` does not check its input: zero or negative values
    of ``m_cat`` or ``C0`` give meaningless results without an error.

    The constant set given here is part of the model that is built, so
    ``reset()`` and :meth:`~nskinetics.KineticModel.compile_events` keep
    it. A value set later with ``set_value`` is kept by ``reset()``, which
    restores only the species and ``yO2``, and is discarded by
    ``compile_events()``, which regenerates the model: set values after
    compiling events. See :func:`load_h2o2_model` on item access.
    """
    return _load(o2_antimony_filename, _table3_or_node_set(
        parameter_set, O2_NODE_PARAMETER_SETS, 'O2'))


def load_o2_simplified_reversible_model(T=None, parameter_set='deposit',
                                        initial_condition='heatmap'):
    """Load the simplified reversible O2 model (the model of Fig. 11).

    Source: the authors' deposited heat-map script, not the paper; the name
    "simplified reversible model" is this package's. Every constant the
    deposit gives is used as given, and where the deposit gives one twice,
    the full-precision value: ``k_1``, ``k_2`` and ``k_5`` from the 50 C
    row of the deposited Arrhenius table of Fig. 10d (not the four-decimal
    values of ``k0``), ``k_4`` = 0.0023 and ``K_6`` to ``K_10`` from the
    vector ``k0`` written in the script, their only source.
    The one constant without a value in the deposit, ``k_3`` (``k0`` writes
    0.0000), is fitted to the 50 C data of the script
    (``fit_simplified_reversible_constants``); see
    :data:`SIMPLIFIED_REVERSIBLE_PARAMETER_SETS`.

    This is the five-reaction model of the authors' heat-map script, with
    no OPA intermediate and a reverse term in every step. It is not the
    Table 3 / Fig. 9b model (:func:`load_o2_model`) and is not written out
    in the paper. Run as the script runs it, with the default constants,
    its contour lines lie 0.30 px (Fig. 11b) and 1.05 px (Fig. 11a) rms
    from those of the published heat maps, whose reading uncertainty is
    about 0.3 px (see ``NOTES.md``).

    Parameters
    ----------
    T : float, optional
        Reaction temperature in K (for 75 C pass ``T=348.15``). Defaults
        to the reference temperature, 323.15 K (50 C). Other temperatures
        scale all ten constants (rate and equilibrium constants) by the
        Arrhenius factors of the script.
    parameter_set : str, optional
        Key of :data:`SIMPLIFIED_REVERSIBLE_PARAMETER_SETS`: ``'deposit'``
        (default) or ``'k0_printed'``.
    initial_condition : str, optional
        Key of :data:`SIMPLIFIED_REVERSIBLE_INITIAL_CONDITIONS`.
        ``'heatmap'`` (default) starts from unit mole fractions, as the
        script does for Fig. 11: water is then absent, the first step is
        irreversible and the O2 fraction falls as by-products form.
        ``'fitting'`` starts from a 3-HP mole fraction of 0.000339, as the
        script does when it fits the 50 C data. The two give different
        curves for the same constants. ``'unscaled'`` and ``'dilute'`` are
        accepted as other names of ``'heatmap'`` and ``'fitting'``.

    Returns
    -------
    nskinetics.KineticModel
        A new model with 60 mM 3-HP. Species: ``HP3, MA, AAc, OXA, D``;
        ``yH2`` and ``yCO2`` are the script's gas fractions. Plotted
        quantities (select them by name): ``X_HP3, X_MA, X_AAc, X_OXA,
        X_D`` (water-free liquid fractions) and the two quantities mapped
        in Fig. 11, ``fig11a_selectivity`` (= ``X_MA``) and ``fig11b_yield``
        (= ``(1 - X_HP3)*X_MA``).

    Raises
    ------
    ValueError
        If ``T`` is outside 250 to 450 K (a guard against a temperature
        given in degrees Celsius, not a range of validity), or if
        ``parameter_set`` or ``initial_condition`` is not one of the names
        above; the message lists the valid names.

    Notes
    -----
    Model time is reaction time; the time axis of Fig. 11 is reaction time
    plus 10 min (see ``heat_maps``).

    The rate constants ``k_1_ref`` to ``k_5_ref`` are in 1/min, lumped for
    the paper's batch (0.2 g of catalyst); ``K_6_ref`` to ``K_10_ref`` are
    dimensionless. There is no catalyst-mass parameter. ``T`` (K) and the
    pressure of the script, ``P`` (bar, default 3.5), can be changed with
    ``set_value``, which does not check them. Do not use ``P`` to explore
    the effect of pressure: the forward rates of steps 1 and 3 are
    proportional to it, so the conversion only rises with it, whereas in
    the paper's 8 h screening experiments (Fig. 3a; 1 to 20 bar of O2)
    conversion and MA selectivity are highest at 3 bar and fall above it
    (see ``NOTES.md``).

    ``x_per_M`` holds the initial condition as the script's initial 3-HP
    mole fraction divided by the 60 mM feed; if ``C0`` is changed, set
    ``x_per_M`` to that mole fraction over the new ``C0`` as well, or the
    model no longer is the script's. The script fits its constants to the
    50 C experiment only, and the Fig. 11 grid runs to 1680 min; at other
    temperatures the model is the script's extrapolation. That
    extrapolation does not describe the paper's experiments at 25, 75 and
    100 C: the pooled R2 against the deposited points of Fig. 10a-c is
    0.653, 0.641 and 0.757 with the ``'fitting'`` initial condition, where
    it is 0.9939 at 50 C (see ``NOTES.md``, "Limits of the temperature
    dependence").

    The temperature, constant set and initial condition given here are
    part of the model that is built, so ``reset()`` and
    :meth:`~nskinetics.KineticModel.compile_events` keep them. A value set
    later with ``set_value`` is kept by ``reset()``, which restores only
    the species and gas fractions, and is discarded by
    ``compile_events()``: set values after compiling events. See
    :func:`load_h2o2_model` on item access.
    """
    try:
        values = dict(SIMPLIFIED_REVERSIBLE_PARAMETER_SETS[parameter_set])
    except (KeyError, TypeError):
        names = ', '.join(repr(name)
                          for name in SIMPLIFIED_REVERSIBLE_PARAMETER_SETS)
        raise ValueError(
            f'unknown parameter set {parameter_set!r} of the simplified '
            f'reversible model: choose from {names}') from None
    try:
        x_3HP = SIMPLIFIED_REVERSIBLE_INITIAL_CONDITIONS[
            _INITIAL_CONDITION_ALIASES.get(initial_condition,
                                           initial_condition)]
    except (KeyError, TypeError):
        names = ' and '.join(
            f'{name!r} (or {alias!r})'
            for alias, name in _INITIAL_CONDITION_ALIASES.items())
        raise ValueError(
            f'unknown initial condition {initial_condition!r} of the '
            f'simplified reversible model: choose from {names}') from None
    values['x_per_M'] = x_3HP/_C0
    if T is not None:
        _check_T(T)
        values['T'] = T
    return _load(o2_simplified_reversible_antimony_filename, values)


#%% Simulation at given times
def simulate_at_times(te_r, times, selections=None):
    """Reset a model and return its trajectory at exactly the given times.

    Parameters
    ----------
    te_r : nskinetics.KineticModel
        Model from one of the loaders of this module.
    times : array-like
        Increasing times in min, at least one of them after t = 0. The
        integration starts at t = 0 whether or not 0 is listed.
    selections : list of str, optional
        Quantities to record. Defaults to the Fig. 9 normalized
        concentrations (:data:`NORMALIZED`).

    Returns
    -------
    pandas.DataFrame
        Column ``time`` followed by one column per selection, one row per
        requested time.

    Raises
    ------
    ValueError
        If ``times`` is not a one-dimensional sequence (a single number
        included), is empty or holds only t = 0 (the solver needs an
        interval to integrate over).
    """
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or not (times > 0.).any():
        raise ValueError('times must be a one-dimensional sequence with at '
                         'least one time after t = 0, not '
                         f'{times.tolist()!r}')
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
