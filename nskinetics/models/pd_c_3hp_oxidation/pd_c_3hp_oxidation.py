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

Two models are provided, each as an Antimony file next to this module:

* the H2O2 model (Table 3, H2O2 column; Fig. 9a; SI Fig. S12a), and
* the O2 model (Table 3, O2 column; Fig. 9b; SI Fig. S12b), with an Arrhenius
  temperature dependence that reduces to the published 50 C model at the
  default ``T = 323.15`` K.

Time is in min and liquid concentrations are in M. Nothing is loaded or
simulated on import: call :func:`load_h2o2_model` or :func:`load_o2_model`.
``NOTES.md`` in this folder records what is printed in the paper, what was
inferred from the authors' data deposit, and where the two disagree.
"""
import os

import numpy as np
import pandas as pd
import tellurium as te

import nskinetics as nsk

__all__ = ('load_h2o2_model', 'load_o2_model', 'simulate_at_times',
           'h2o2_antimony_filename', 'o2_antimony_filename',
           'SPECIES', 'NORMALIZED', 'TABLE_S4',)

#%% Model files
h2o2_antimony_filename = 'pd_c_3hp_oxidation_h2o2_antimony.txt'
o2_antimony_filename = 'pd_c_3hp_oxidation_o2_antimony.txt'

#: Liquid species shared by the models, in M: 3-HP, 3-oxopropanoic acid,
#: malonic acid, acetic acid, oxalic acid and the lumped other by-products.
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

    Parameters
    ----------
    T : float, optional
        Reaction temperature in K. Defaults to the model's reference
        temperature, 323.15 K (50 C), at which the model is the published
        one. Other temperatures scale the six rate constants by Arrhenius
        factors; this is an extrapolation the paper does not show for this
        model (see ``NOTES.md``).

    Returns
    -------
    nskinetics.KineticModel
        A new model with 60 mM 3-HP and 0.2 g of catalyst. Species:
        ``HP3, OPA, MA, AAc, OXA, D``; ``yO2`` is the dimensionless O2
        activity. Plotted quantities (select them by name): ``X_HP3, X_MA,
        X_AAc, X_OXA, X_D`` (Fig. 9b) and ``x_OPA`` (SI Fig. S12b inset).

    Notes
    -----
    See :func:`load_h2o2_model` on item access. ``reset()`` restores the
    species and ``yO2`` but keeps changed parameters such as ``T``.
    """
    te_r = _load(o2_antimony_filename)
    if T is not None:
        te_r.set_value('T', T)
    return te_r


#%% Simulation at given times
def simulate_at_times(te_r, times, selections=None):
    """Reset a model and return its trajectory at exactly the given times.

    Parameters
    ----------
    te_r : nskinetics.KineticModel
        Model from :func:`load_h2o2_model` or :func:`load_o2_model`.
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
