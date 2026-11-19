# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Kinetic models for the selective oxidation of 3-hydroxypropionic acid (3-HP)
to malonic acid (MA) over 10 % Pd/C, from Kim, Chavarrio, Install and Huber,
Appl. Catal. B: Environ. Energy 386 (2026) 126403
(https://doi.org/10.1016/j.apcatb.2026.126403).

Four models are provided, each as an Antimony file next to this module:

* the H2O2 model (Table 3, H2O2 column), which reproduces Fig. 9a and SI
  Fig. S12a;
* the O2 model (Table 3, O2 column), which reproduces Fig. 9b and SI
  Fig. S12b, with an Arrhenius temperature dependence that reduces to the
  published 50 C model at the default ``T = 323.15`` K. The temperature
  dependence is an NSKinetics addition: away from 50 C it reproduces
  neither the published heat maps nor the experiments at 25, 75 and 100 C;
* the reversible six-reaction O2 variant, which reproduces the model curves
  of Fig. 10a-c (25, 75 and 100 C). It is not in the paper or the authors'
  scripts: structure and constants were inferred from the deposited curves;
* the authors' simplified reversible five-reaction O2 model (no OPA), which
  is not in the paper but generated its heat maps (Fig. 11); three of its
  rate constants are the ones plotted in the Arrhenius plot (Fig. 10d).

Time is in min and liquid concentrations are in M. Nothing is loaded or
simulated on import: call :func:`load_h2o2_model`, :func:`load_o2_model`,
:func:`load_o2_reversible_model` or
:func:`load_o2_simplified_reversible_model`.
``NOTES.md`` in this folder records what is printed in the paper, what was
inferred from the authors' data deposit, and where the two disagree.
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

#: Liquid species of the two Table 3 models, in M: 3-HP, 3-oxopropanoic
#: acid, malonic acid, acetic acid, oxalic acid and the lumped other
#: by-products. The reversible O2 variant has the same plus ``W`` (the
#: water formed by reaction); the simplified reversible model has the same
#: without OPA.
SPECIES = ('HP3', 'OPA', 'MA', 'AAc', 'OXA', 'D')

#: Assignment rules holding the "Normalized Concentration" plotted in Fig. 9:
#: each measured species over the sum of the five measured species (OPA is
#: not in the sum).
NORMALIZED = ('X_HP3', 'X_MA', 'X_AAc', 'X_OXA', 'X_D')

#: SI Table S4 ("Calculated reaction orders and equilibrium constant"), kept
#: for documentation only. Neither the orders ``n`` nor the equilibrium
#: constants ``Keq`` are used by the rate laws: the model curves the authors
#: deposited are reproduced only by first-order, irreversible steps, and the
#: authors' posted script does not use them either. ``None`` marks a dash in
#: the table. The H2O2 reaction 7 order is cited from the literature there.
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

_O2_REVERSIBLE_NAMES = ('k_1', 'k_1r', 'k_2', 'k_2r', 'k_3', 'k_4',
                        'k_5', 'k_5r', 'k_6', 'k_6r')

#: Constant sets of the reversible O2 variant, keyed by the temperature in C
#: of the experiment whose deposited model curve they reproduce (Fig. 10a,
#: b and c). Forward constants ``k_i`` and reverse constants ``k_ir`` in
#: 1/min for the 0.2 g batch; ``k_1r`` is per unit of ``[W]/C0``. All values
#: were inferred from the deposited curves (see ``NOTES.md``); the three
#: sets are separate fits and follow no temperature law. There is no set
#: for 50 C.
O2_REVERSIBLE_PARAMETER_SETS = {
    temperature_C: dict(zip(_O2_REVERSIBLE_NAMES, values))
    for temperature_C, values in (
        (25, (0.0465036, 51.0644, 0.0088933, 0.00136441, 2.43424e-4,
              6.17632e-5, 1.01563e-6, 1.47217e-7, 0.00284329, 3.12354e-4)),
        (75, (0.00918893, 1.88278e-5, 0.0490036, 0.00193586, 4.9301e-4,
              7.97777e-6, 0.00160383, 2.86632e-5, 0.00613222, 2.56106e-5)),
        (100, (0.0414885, 0.00474984, 0.122878, 0.0103261, 0.00169028,
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


#%% Loaders
def _load(antimony_filename):
    r = te.loadAntimonyModel(
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     antimony_filename))
    te_r = nsk.KineticModel(r, units={'time': 'min', 'conc': 'M'})
    integrator = r.integrator
    integrator.absolute_tolerance = 1e-10
    integrator.relative_tolerance = 1e-9
    return te_r


def load_h2o2_model():
    """Load the H2O2 model (Table 3, H2O2 column; Fig. 9a) at 50 C.

    Reproduces Fig. 9a and SI Fig. S12a. Source: the paper (constants of
    Table 3 as printed); the rate laws were reconstructed from the deposited
    model curve.

    Returns
    -------
    nskinetics.KineticModel
        A new model with 60 mM 3-HP, 15 equivalents (0.9 M) of H2O2 and
        0.2 g of catalyst. Species: ``HP3, OPA, MA, AAc, OXA, D, H2O2``.
        Plotted quantities (select them by name): ``X_HP3, X_MA, X_AAc,
        X_OXA, X_D`` (Fig. 9a), ``eq_H2O2`` (SI Fig. S12a inset) and
        ``x_OPA``.

    Notes
    -----
    Each call loads a new, independent model. Read and set model values
    through :meth:`~nskinetics.KineticModel.get_value` and
    :meth:`~nskinetics.KineticModel.set_value` (item access), not through
    attributes of the underlying roadrunner object: attribute access is
    unreliable once more than one model has been loaded in a process.
    """
    return _load(h2o2_antimony_filename)


def load_o2_model(T=None):
    """Load the O2 model (Table 3, O2 column; Fig. 9b).

    Reproduces Fig. 9b and SI Fig. S12b at the default 50 C. Source: the
    paper (constants of Table 3 as printed), with the factor ``x_O2_ref``
    inferred from the deposited model curve; the temperature dependence is
    an NSKinetics addition.

    Parameters
    ----------
    T : float, optional
        Reaction temperature in K. Defaults to the model's reference
        temperature, 323.15 K (50 C), at which the model is the published
        one. Other temperatures scale the six rate constants by Arrhenius
        factors (the form of the paper's Eq. 6 with the activation energies
        of the authors' heat-map script); see the warning below.

    Returns
    -------
    nskinetics.KineticModel
        A new model with 60 mM 3-HP and 0.2 g of catalyst. Species:
        ``HP3, OPA, MA, AAc, OXA, D``; ``yO2`` is the dimensionless O2
        activity. Plotted quantities (select them by name): ``X_HP3, X_MA,
        X_AAc, X_OXA, X_D`` (Fig. 9b) and ``x_OPA`` (SI Fig. S12b inset).

    Warnings
    --------
    Only the 50 C model is a published one. Away from 50 C the model is an
    extrapolation that reproduces no published item: its heat maps miss the
    contour lines of Fig. 11 (by about 7.6 px rms, against 0.17 and 0.27 px
    for the simplified reversible model that drew them), and it is too fast
    for the experiments at 25, 75 and 100 C (pooled R2 0.579, 0.644 and
    0.746, against 0.9996, 0.990 and 0.993 for the per-temperature fits of
    :func:`load_o2_reversible_model`). See ``NOTES.md`` and the
    documentation page of this package.

    Notes
    -----
    See :func:`load_h2o2_model` on item access. ``reset()`` restores the
    species and ``yO2`` but keeps changed parameters such as ``T``.
    """
    te_r = _load(o2_antimony_filename)
    if T is not None:
        te_r.set_value('T', T)
    return te_r


def load_o2_reversible_model(temperature_C=75):
    """Load the reversible O2 variant behind the curves of Fig. 10a-c.

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
    temperature_C : {25, 75, 100}, optional
        Temperature in C of the experiment; selects one of the three
        constant sets of :data:`O2_REVERSIBLE_PARAMETER_SETS`. Default 75
        (Fig. 10b), the set held by the Antimony file. The sets were fitted
        separately and there is no temperature law, so no other value is
        accepted.

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
        If no constant set exists for ``temperature_C``. In particular
        there is none for 50 C: the 50 C curve (Fig. 9b) is the
        irreversible Table 3 model, :func:`load_o2_model`.

    Notes
    -----
    The constants are effective first-order constants in 1/min for the
    0.2 g batch and are not in the units of Table 3. The deposited curves
    are MATLAB ``ode23s`` output at default tolerances; this model, solved
    accurately, differs from them by up to 1.8e-4 (25 C), 2.5e-3 (75 C)
    and 2.7e-3 (100 C). See :func:`load_h2o2_model` on item access.
    """
    if temperature_C not in O2_REVERSIBLE_PARAMETER_SETS:
        raise ValueError(
            'no parameter set of the reversible O2 variant exists for '
            f'{temperature_C!r} C: the sets are separate fits at '
            f'{", ".join(str(T) for T in O2_REVERSIBLE_PARAMETER_SETS)} C '
            'with no temperature law. For 50 C (Fig. 9b) use load_o2_model, '
            'the irreversible Table 3 model.')
    te_r = _load(o2_reversible_antimony_filename)
    for name, value in O2_REVERSIBLE_PARAMETER_SETS[temperature_C].items():
        te_r.set_value(name, value)
    return te_r


def load_o2_simplified_reversible_model(T=None, parameter_set='fig11b',
                                        initial_condition='heatmap'):
    """Load the authors' simplified reversible O2 model (behind Fig. 11).

    Reproduces Fig. 11; three of its rate constants are the points of
    Fig. 10d. Source: the authors' deposited heat-map script, not the paper;
    the constants that are not in the deposit (``k_3`` and ``k_4``; for
    Fig. 11a also ``k_5`` and ``K_10``) were inferred from the published
    contour lines.

    This is the five-reaction model of the authors' heat-map script, with
    no OPA intermediate and a reverse term in every step. It is not the
    Table 3 / Fig. 9b model (:func:`load_o2_model`) and is not written out
    in the paper, but it generated the published heat maps (Fig. 11).

    Parameters
    ----------
    T : float, optional
        Reaction temperature in K. Defaults to the reference temperature,
        323.15 K (50 C). Other temperatures scale all ten constants (rate
        and equilibrium constants) by the Arrhenius factors of the script.
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

    Notes
    -----
    Model time is reaction time; the time axis of Fig. 11 is reaction time
    plus 10 min (see ``heat_maps``). See :func:`load_h2o2_model` on item
    access. ``reset()`` restores the species and gas fractions but keeps
    changed parameters such as ``T``.
    """
    te_r = _load(o2_simplified_reversible_antimony_filename)
    for name, value in SIMPLIFIED_REVERSIBLE_PARAMETER_SETS[
            parameter_set].items():
        te_r.set_value(name, value)
    initial_condition = _INITIAL_CONDITION_ALIASES.get(initial_condition,
                                                       initial_condition)
    te_r.set_value(
        'x_per_M',
        SIMPLIFIED_REVERSIBLE_INITIAL_CONDITIONS[initial_condition]
        / te_r.get_value('C0'))
    if T is not None:
        te_r.set_value('T', T)
    return te_r


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
