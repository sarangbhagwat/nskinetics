# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Readers for the reference data shipped in ``reference_data/``: the authors'
deposited Fig. 9 tables (Dryad, https://doi.org/10.5061/dryad.jsxksn0q4, CC0,
copied verbatim) and two curves digitized from SI Fig. S12 of Kim et al.,
Appl. Catal. B: Environ. Energy 386 (2026) 126403.

The readers return tidy tables whose column names are the model's selection
names (``X_HP3``, ``X_MA``, ``X_AAc``, ``X_OXA``, ``X_D``, ``eq_H2O2``,
``x_OPA``) and undo the known quirks of the deposited files; the CSV files
themselves are left untouched. ``reference_data/README.md`` gives the
provenance of every file.
"""
import csv
import os

import pandas as pd

__all__ = ('reference_data_dir',
           'read_fig9a_h2o2', 'read_fig9b_o2',
           'read_figS12a_h2o2_inset_digitized',
           'read_figS12b_opa_inset_digitized',)

reference_data_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'reference_data')

_COLUMNS = ['time', 'X_HP3', 'X_MA', 'X_AAc', 'X_OXA', 'X_D']


def _rows(filename):
    # utf-8-sig: the deposited files start with a byte-order mark.
    with open(os.path.join(reference_data_dir, filename),
              encoding='utf-8-sig', newline='') as f:
        rows = list(csv.reader(f))
    # Data start two rows below the 'Time,3HP,...' header (a label row follows).
    start = next(i for i, row in enumerate(rows) if row[0] == 'Time') + 2
    return rows[start:]


def _block(rows, first, columns):
    """Columns ``first .. first+len(columns)`` of the rows where they are filled."""
    last = first + len(columns)
    data = [[float(cell) for cell in row[first:last]]
            for row in rows if len(row) >= last and row[first] != '']
    return pd.DataFrame(data, columns=columns)


def read_fig9a_h2o2():
    """Read ``10-1_kinetic_model_H2O2.csv`` (Fig. 9a and SI Fig. S12a, H2O2).

    Returns
    -------
    dict of pandas.DataFrame
        ``'experiment'``
            13 sampling times (0 to 1440 min) with the measured normalized
            concentrations ``X_HP3, X_MA, X_AAc, X_OXA, X_D``.
        ``'model'``
            The authors' model curve of Fig. 9a: 24 solver nodes
            (0 to 1584 min), same columns.
        ``'h2o2_trace'``
            The authors' predicted H2O2 of SI Fig. S12a in equivalents of
            the initial 3-HP (``time``, ``eq_H2O2``): 24 nodes
            (0 to 1537 min) on a time grid of its own.

    Notes
    -----
    Three quirks of the deposited file are undone here:

    * the block labelled ``exp`` (columns 1-6) holds the model nodes and
      the block labelled ``model`` (columns 9-14) holds the experimental
      points;
    * the experimental t = 0 row is printed as zeros and is returned as
      the feed, (1, 0, 0, 0, 0);
    * the t = 0 cell of the H2O2 trace (columns 7-8) reads 1 and is
      returned as 15, the initial charge.

    The H2O2 trace was generated with a different parameter set than the
    model block and Table 3 (see ``NOTES.md``).
    """
    rows = _rows('10-1_kinetic_model_H2O2.csv')
    model = _block(rows, 0, _COLUMNS)
    h2o2_trace = _block(rows, 6, ['time', 'eq_H2O2'])
    experiment = _block(rows, 8, _COLUMNS)
    experiment.loc[0, _COLUMNS[1:]] = [1., 0., 0., 0., 0.]
    h2o2_trace.loc[0, 'eq_H2O2'] = 15.
    return {'experiment': experiment, 'model': model,
            'h2o2_trace': h2o2_trace}


def read_fig9b_o2():
    """Read ``10-2_kinetic_model_O2.csv`` (Fig. 9b, O2).

    Returns
    -------
    dict of pandas.DataFrame
        ``'experiment'``
            16 sampling times (0 to 1680 min) with the measured normalized
            concentrations ``X_HP3, X_MA, X_AAc, X_OXA, X_D``.
        ``'model'``
            The authors' model curve of Fig. 9b: 20 solver nodes
            (0 to 1848 min), same columns.

    Notes
    -----
    This file is labelled correctly (columns 1-6 experiment, 7-12 model).
    """
    rows = _rows('10-2_kinetic_model_O2.csv')
    return {'experiment': _block(rows, 0, _COLUMNS),
            'model': _block(rows, 6, _COLUMNS)}


def read_figS12a_h2o2_inset_digitized():
    """Read the predicted H2O2 curve digitized from SI Fig. S12a.

    Returns
    -------
    pandas.DataFrame
        ``time`` (min) and ``eq_H2O2`` (equivalents of the initial 3-HP),
        with the digitization columns ``source_axes`` (``'inset'``, 0-15
        scale, or ``'main'``, where the curve is drawn below 1), ``kind``
        and ``flag``.

    Notes
    -----
    Digitized from the published figure, not deposited data. Reading
    uncertainty: about 0.02 equivalents (or 1.5 min) on the inset and
    0.3 % of the value on the main axes.
    """
    df = pd.read_csv(os.path.join(reference_data_dir,
                                  'figS12a_h2o2_inset_digitized.csv'),
                     keep_default_na=False)
    return df.rename(columns={'t_min': 'time', 'h2o2_norm': 'eq_H2O2'})


def read_figS12b_opa_inset_digitized():
    """Read the predicted OPA curve digitized from the SI Fig. S12b inset.

    Returns
    -------
    pandas.DataFrame
        The 13 vertices of the plotted polyline: ``time`` (min),
        ``time_unc`` (its reading uncertainty, min), ``x_OPA`` (OPA liquid
        mole fraction; the inset axis is in units of 1e-5), ``X_HP3`` (the
        3-HP model line of the main axes at the same times) and
        ``ratio_OPA_over_3HP``.

    Notes
    -----
    Digitized from the published figure, not deposited data. Reading
    uncertainty of ``x_OPA``: 0.006e-5 on the decay and 0.015e-5 around
    the peak (6 to 26 min).
    """
    df = pd.read_csv(os.path.join(reference_data_dir,
                                  'figS12b_opa_inset_digitized.csv'))
    return df.rename(columns={'t_min': 'time', 't_unc': 'time_unc',
                              'OPA': 'x_OPA', 'HP3': 'X_HP3'})
