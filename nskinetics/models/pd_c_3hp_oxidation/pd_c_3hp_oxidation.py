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

Four models are provided, each as an Antimony file next to this module:

* the H2O2 model (Table 3, H2O2 column), which reproduces Fig. 9a. Its
  predicted H2O2 is compared with SI Fig. S12a, whose H2O2 curve lies up
  to 0.38 equivalents lower and is matched only with other constants;
* the O2 model (Table 3, O2 column), which reproduces Fig. 9b and SI
  Fig. S12b, with an Arrhenius temperature dependence that reduces to the
  published 50 C model at the default ``T = 323.15`` K. The temperature
  dependence is an NSKinetics addition: away from 50 C it reproduces
  neither the published heat maps nor the experiments at 25, 75 and 100 C;
* the reversible six-reaction O2 variant, which reproduces the model curves
  of Fig. 10a-c (25, 75 and 100 C) with one constant set per curve. It is
  not in the paper or the authors' scripts: structure and constants were
  inferred from the deposited curves;
* the authors' simplified reversible five-reaction O2 model (no OPA), which
  is not in the paper but is the model of their deposited heat-map script
  and reproduces the heat maps (Fig. 11); the 50 C point of the Arrhenius
  plot (Fig. 10d) supplies three of its rate constants.

The rate laws of the two Table 3 models have the form of SI Fig. S2 (rate
expressions), ``r_j = k_jf*prod(a_i^n_ij)``, with every order ``n`` equal
to 1, unit activity coefficients and the liquid activities on the authors'
mole-fraction scale. The O2 model departs from Fig. S2 in one place: the gas
activity ``phi*y*P/f0`` of O2 is replaced by the inferred constant
``x_O2_ref`` times ``yO2``. ``NOTES.md`` (section "SI Fig. S2 and Eq. S1")
maps each of the four models to Fig. S2 term by term.

Time is in min, liquid concentrations are in M and temperatures are in K
(arguments and parameters named ``T``) unless a name ends in ``_C``
(degrees Celsius). Nothing is loaded or simulated on import: call
:func:`load_h2o2_model`, :func:`load_o2_model`,
:func:`load_o2_reversible_model` or
:func:`load_o2_simplified_reversible_model`.
``NOTES.md`` in this folder records what is printed in the paper, what was
inferred from the authors' data deposit, and where the two disagree.

The model curves the authors deposited are solver steps and lie up to
2.7e-3 from the solutions of these models. ``solver_emulation``
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

__all__ = ('load_h2o2_model', 'load_o2_model', 'load_o2_reversible_model',
           'load_o2_simplified_reversible_model', 'simulate_at_times',
           'h2o2_antimony_filename', 'o2_antimony_filename',
           'o2_reversible_antimony_filename',
           'o2_simplified_reversible_antimony_filename',
           'SPECIES', 'NORMALIZED', 'TABLE_S4',
           'H2O2_NODE_PARAMETER_SETS', 'O2_NODE_PARAMETER_SETS',
           'O2_REVERSIBLE_PARAMETER_SETS',
           'SIMPLIFIED_REVERSIBLE_PARAMETER_SETS',
           'SIMPLIFIED_REVERSIBLE_INITIAL_CONDITIONS',)

#%% Model files
h2o2_antimony_filename = 'pd_c_3hp_oxidation_h2o2_antimony.txt'
o2_antimony_filename = 'pd_c_3hp_oxidation_o2_antimony.txt'
o2_reversible_antimony_filename = (
    'pd_c_3hp_oxidation_o2_reversible_antimony.txt')
o2_simplified_reversible_antimony_filename = (
    'pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt')

#: Liquid species of the two Table 3 models (concentrations in M), in the
#: order 3-HP, 3-oxopropanoic acid, malonic acid, acetic acid, oxalic acid
#: and the lumped other by-products. The reversible O2 variant has the same
#: plus ``W`` (the water formed by reaction); the simplified reversible
#: model has the same without OPA.
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

_O2_REVERSIBLE_NAMES = ('k_1', 'k_1r', 'k_2', 'k_2r', 'k_3', 'k_4',
                        'k_5', 'k_5r', 'k_6', 'k_6r')

#: Constant sets of the reversible O2 variant, keyed by a label that names
#: the experiment whose deposited model curve the set reproduces: ``'25C'``
#: (Fig. 10a), ``'75C'`` (Fig. 10b) and ``'100C'`` (Fig. 10c). Forward
#: constants ``k_i`` and reverse constants ``k_ir`` in 1/min for the 0.2 g
#: batch; ``k_1r`` is per unit of ``[W]/C0``. All values were inferred from
#: the deposited curves (see ``NOTES.md``); the three sets are separate fits
#: and follow no temperature law, so the labels are names, not temperatures
#: that could be interpolated. There is no set for 50 C.
O2_REVERSIBLE_PARAMETER_SETS = {
    label: dict(zip(_O2_REVERSIBLE_NAMES, values))
    for label, values in (
        ('25C', (0.0465036, 51.0644, 0.0088933, 0.00136441, 2.43424e-4,
                 6.17632e-5, 1.01563e-6, 1.47217e-7, 0.00284329,
                 3.12354e-4)),
        ('75C', (0.00918893, 1.88278e-5, 0.0490036, 0.00193586, 4.9301e-4,
                 7.97777e-6, 0.00160383, 2.86632e-5, 0.00613222,
                 2.56106e-5)),
        ('100C', (0.0414885, 0.00474984, 0.122878, 0.0103261, 0.00169028,
                  7.86978e-6, 9.32644e-7, 0., 0.0157177, 1.05066e-4)),
        )}

_SIMPLIFIED_REVERSIBLE_NAMES = (
    'k_1_ref', 'k_2_ref', 'k_3_ref', 'k_4_ref', 'k_5_ref',
    'K_6_ref', 'K_7_ref', 'K_8_ref', 'K_9_ref', 'K_10_ref')

#: Constant sets of the simplified reversible model at 50 C, in the order of
#: the authors' script: five rate constants (1/min) and five equilibrium
#: constants. The fitted constants were not deposited, and the two panels of
#: Fig. 11 need different sets (see ``NOTES.md``).
#:
#: ``'fig11b'``
#:     Reproduces Fig. 11b and is the default of the Antimony file. ``k_1``,
#:     ``k_2`` and ``k_5`` are ``exp()`` of the 50 C row of the deposited
#:     ``11-4_Arrhenius.csv``; ``k_3`` and ``k_4`` were inferred from the
#:     contour lines digitized from Fig. 11b; ``K_6`` to ``K_10`` are the
#:     values printed in the script.
#: ``'fig11a'``
#:     Reproduces Fig. 11a. As ``'fig11b'`` except for ``k_3``, ``k_4``,
#:     ``k_5`` and ``K_10``, inferred from the contour lines of Fig. 11a.
#: ``'k0_printed'``
#:     The starting vector printed in the script, to its four decimals
#:     (``k_3`` reads 0). For reference only: it does not reproduce Fig. 11.
SIMPLIFIED_REVERSIBLE_PARAMETER_SETS = {
    name: dict(zip(_SIMPLIFIED_REVERSIBLE_NAMES, values))
    for name, values in (
        ('fig11b', (2.139999e-3, 5.199984e-5, 1.435e-5, 2.3040e-3,
                    1.062501e-3, 4.7840, 0.5918, 3.5116, 81.6870, 52.0472)),
        ('fig11a', (2.139999e-3, 5.199984e-5, 1.4912e-5, 2.2397e-3,
                    1.1110e-3, 4.7840, 0.5918, 3.5116, 81.6870, 23.28)),
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

# Initial 3-HP concentration held by the four model files, in M.
_C0 = 0.06

# Range of temperatures in K that the loaders accept. It only guards against
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
    other value of these parameters is an extrapolation that was never
    compared with data. ``set_value`` does not check its input: zero or
    negative values of ``m_cat``, ``C0`` or ``eq_H2O2_0`` give meaningless
    results without an error.

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
    extrapolation that reproduces no published item: its heat maps miss the
    contour lines of Fig. 11 (by about 7.6 px rms, against 0.17 and 0.27 px
    for the simplified reversible model that reproduces them), and it is
    too fast for the experiments at 25, 75 and 100 C (pooled R2 0.579,
    0.644 and 0.746, against 0.9996, 0.990 and 0.993 for the
    per-temperature constant sets of :func:`load_o2_reversible_model`).
    See ``NOTES.md`` and the documentation page of this package.

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
        The model has no O2 pressure parameter.

    The constants of Table 3 were fitted to one experiment (60 mM 3-HP,
    0.2 g of catalyst, 3 bar O2, 50 C, samples up to 1680 min); any other
    value of these parameters is an extrapolation that was never compared
    with data, except for ``T``, for which the comparison is poor (see the
    warning above). ``set_value`` does not check its input: zero or
    negative values of ``m_cat`` or ``C0`` give meaningless results
    without an error.

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


def load_o2_reversible_model(parameter_set='75C'):
    """Load the reversible O2 variant that reproduces the curves of Fig. 10a-c.

    Reproduces the model curves of Fig. 10a-c (25, 75 and 100 C). Source:
    inferred from the deposited curves; it is in neither the paper nor the
    authors' scripts.

    This is the six-reaction network of the Table 3 O2 model with
    reactions 1, 2, 5 and 6 reversible. The text of the paper treats all
    reactions as irreversible; the structure and the constants of this
    variant were inferred from the model curves deposited for the
    experiments at 25, 75 and 100 C. See ``NOTES.md``.

    Parameters
    ----------
    parameter_set : {'25C', '75C', '100C'}, optional
        Key of :data:`O2_REVERSIBLE_PARAMETER_SETS`: the constant set that
        reproduces the curve of the experiment at 25 C (Fig. 10a), 75 C
        (Fig. 10b, the default and the set held by the Antimony file) or
        100 C (Fig. 10c). The model has no temperature parameter and no
        temperature law: the three sets are separate fits, and a label
        selects one of them.

    Returns
    -------
    nskinetics.KineticModel
        A new model with 60 mM 3-HP. Species: ``HP3, OPA, MA, AAc, OXA,
        D`` and ``W``, the water formed by reactions 1 and 4, which starts
        at zero and enters the reverse of reaction 1. Plotted quantities
        (select them by name): ``X_HP3, X_MA, X_AAc, X_OXA, X_D``
        (Fig. 10a-c) and ``x_OPA``.

    Raises
    ------
    ValueError
        If ``parameter_set`` is not one of the three labels; the message
        lists them. In particular there is no set for 50 C: the 50 C curve
        (Fig. 9b) is the irreversible Table 3 model, :func:`load_o2_model`.

    Notes
    -----
    The constants ``k_1, k_1r, k_2, k_2r, k_3, k_4, k_5, k_5r, k_6, k_6r``
    are effective first-order constants in 1/min, lumped for the paper's
    batch (0.2 g of catalyst, 3 bar O2), and are not in the units of
    Table 3. There is no catalyst-mass, pressure or temperature parameter.
    ``C0`` (M, default 0.06) can be changed and takes effect at the next
    ``reset()``; the plotted fractions do not depend on it. Each set
    belongs to one experiment only (samples up to 1890, 1440 and 1560 min
    at 25, 75 and 100 C).

    The deposited curves are solver steps, not converged solutions (they
    are reproduced by the method of MATLAB ``ode23s`` at default
    tolerances); this model, solved accurately, differs from them by up to
    1.8e-4 (25 C), 2.5e-3 (75 C) and 2.7e-3 (100 C). With that solver run
    repeated (``solver_emulation.emulate_fig10abc``) it returns the
    deposited nodes within 1e-5, the last decimal the files print.

    The selected constant set is part of the model that is built, so
    ``reset()`` and :meth:`~nskinetics.KineticModel.compile_events` keep
    it. A value set later with ``set_value`` is kept by ``reset()`` and
    discarded by ``compile_events()``: set values after compiling events.
    See :func:`load_h2o2_model` on item access.
    """
    try:
        values = O2_REVERSIBLE_PARAMETER_SETS[parameter_set]
    except (KeyError, TypeError):
        labels = ', '.join(repr(label)
                           for label in O2_REVERSIBLE_PARAMETER_SETS)
        raise ValueError(
            f'unknown parameter set {parameter_set!r} of the reversible O2 '
            f'variant: choose from {labels} (separate fits to the '
            'experiments at 25, 75 and 100 C, with no temperature law). '
            'For 50 C (Fig. 9b) use load_o2_model, the irreversible '
            'Table 3 model.') from None
    return _load(o2_reversible_antimony_filename, values)


def load_o2_simplified_reversible_model(T=None, parameter_set='fig11b',
                                        initial_condition='heatmap'):
    """Load the authors' simplified reversible O2 model (the model of Fig. 11).

    Reproduces Fig. 11; the 50 C point of Fig. 10d supplies three of its
    rate constants. Source: the authors' deposited heat-map script, not the
    paper; the constants that are not in the deposit (``k_3`` and ``k_4``;
    for Fig. 11a also ``k_5`` and ``K_10``) were inferred from the published
    contour lines.

    This is the five-reaction model of the authors' heat-map script, with
    no OPA intermediate and a reverse term in every step. It is not the
    Table 3 / Fig. 9b model (:func:`load_o2_model`) and is not written out
    in the paper; run as the script runs it, it reproduces the published
    heat maps (Fig. 11).

    Parameters
    ----------
    T : float, optional
        Reaction temperature in K (for 75 C pass ``T=348.15``). Defaults
        to the reference temperature, 323.15 K (50 C). Other temperatures
        scale all ten constants (rate and equilibrium constants) by the
        Arrhenius factors of the script.
    parameter_set : str, optional
        Key of :data:`SIMPLIFIED_REVERSIBLE_PARAMETER_SETS`: ``'fig11b'``
        (default), ``'fig11a'`` or ``'k0_printed'``.
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
    ``set_value``, which does not check them. ``x_per_M`` holds the
    initial condition as the script's initial 3-HP mole fraction divided
    by the 60 mM feed; if ``C0`` is changed, set ``x_per_M`` to that mole
    fraction over the new ``C0`` as well, or the model no longer is the
    script's. The script fits its constants to the 50 C experiment only,
    and the Fig. 11 grid runs to 1680 min; at other temperatures the model
    is the script's extrapolation.

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
