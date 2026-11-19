# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Readers for the reference data shipped in ``reference_data/``: the authors'
deposited Fig. 9 and Fig. 10 tables (Dryad,
https://doi.org/10.5061/dryad.jsxksn0q4, CC0, copied verbatim), two curves
digitized from SI Fig. S12 and the contour lines digitized from Fig. 11 of
Kim et al., Appl. Catal. B: Environ. Energy 386 (2026) 126403.

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
           'read_figS12b_opa_inset_digitized',
           'read_fig10_o2', 'read_fig10d_arrhenius',
           'read_fig11a_selectivity_digitized',
           'read_fig11b_yield_digitized',)

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
            The authors' model curve of Fig. 9a: 24 solver nodes (23
            steps, 0 to 1584 min; here and below the node at t = 0 is
            counted), same columns.
        ``'h2o2_trace'``
            The authors' predicted H2O2 of SI Fig. S12a in equivalents of
            the initial 3-HP (``time``, ``eq_H2O2``): 24 nodes (23 steps,
            0 to 1537 min) on a time grid of its own.

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

    The H2O2 trace is not the H2O2 of the Table 3 model, which lies up to
    0.38 equivalents above it (see ``NOTES.md``).
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
            The authors' model curve of Fig. 9b: 20 solver nodes (19
            steps, 0 to 1848 min), same columns.

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


_FIG10_FILES = {25: '11-1_25C_Kinetic.csv', 75: '11-2_75C_Kinetic.csv',
                100: '11-3_100C_Kinetic.csv'}
# The labels of the constant sets of the reversible O2 variant ('25C', ...)
# name the same three experiments.
_FIG10_LABELS = {'%dC' % T_C: T_C for T_C in _FIG10_FILES}


def _fig10_temperature(temperature_C):
    """Return 25, 75 or 100 for the temperature or label of a Fig. 10 panel."""
    try:
        T_C = _FIG10_LABELS.get(temperature_C, temperature_C)
        known = T_C in _FIG10_FILES
    except TypeError:               # not hashable, e.g. a list
        known = False
    if not known:
        raise ValueError(
            f'no Fig. 10 data for {temperature_C!r}: the panels are at 25, '
            "75 and 100 C (pass the number or the label '25C', '75C' or "
            "'100C'; for 50 C, Fig. 9b, use read_fig9b_o2)")
    return int(T_C)


def read_fig10_o2(temperature_C):
    """Read the deposited table of one panel of Fig. 10a-c (O2).

    Parameters
    ----------
    temperature_C : {25, 75, 100} or {'25C', '75C', '100C'}
        Temperature of the experiment in degrees Celsius: Fig. 10a
        (``11-1_25C_Kinetic.csv``), Fig. 10b (``11-2_75C_Kinetic.csv``) or
        Fig. 10c (``11-3_100C_Kinetic.csv``). The label of the constant set
        of that experiment, as ``load_o2_reversible_model`` takes it, is
        accepted as well.

    Returns
    -------
    dict of pandas.DataFrame
        ``'experiment'``
            The sampling times with the measured normalized concentrations
            ``X_HP3, X_MA, X_AAc, X_OXA, X_D``: 14 rows (0 to 1890 min) at
            25 C, 12 (0 to 1440 min) at 75 C, 11 (0 to 1560 min) at 100 C.
        ``'model'``
            The authors' model curve, same columns: the steps of a
            solver run, 19 nodes (18 steps, 0 to 2079 min) at 25 C, 20
            nodes (19 steps, 0 to 1584 min) at 75 C and 25 nodes (24
            steps, 0 to 1716 min) at 100 C, the node at t = 0 included.
            The last time is 1.1 times the last sampling time.

    Raises
    ------
    ValueError
        For any other temperature or label.

    Notes
    -----
    The files are labelled correctly (columns 1-6 experiment, 7-12 model)
    and need no correction. The model curves are coarse polylines of 19
    to 25 nodes (18 to 24 solver steps). The reversible O2 variant
    (``load_o2_reversible_model``), solved accurately, lies within 2.7e-3
    of the nodes; with the authors' solver run repeated
    (``solver_emulation``) it returns them within 1e-5 (see ``NOTES.md``).
    """
    rows = _rows(_FIG10_FILES[_fig10_temperature(temperature_C)])
    return {'experiment': _block(rows, 0, _COLUMNS),
            'model': _block(rows, 6, _COLUMNS)}


def read_fig10d_arrhenius():
    """Read ``11-4_Arrhenius.csv`` (Arrhenius plot of Fig. 10d).

    Returns
    -------
    dict of pandas.DataFrame
        ``'points'``
            The four fitted points (25, 50, 75 and 100 C): ``inv_T`` (1/K,
            as printed, to five decimals) and the natural logarithms
            ``ln_k_1``, ``ln_k_2`` and ``ln_k_5`` of three rate constants
            in 1/min, named as described in the Notes.
        ``'line_k_1'``, ``'line_k_2'``, ``'line_k_5'``
            The straight line of the published figure for each constant:
            ``inv_T`` and the logarithm.

    Notes
    -----
    The file labels the three constants k1, k2 and k3, and the paper does
    not say of which model they are. They are taken to be ``k_1`` (3-HP to
    MA), ``k_2`` (MA to AAc) and ``k_5`` (MA to D) of the simplified
    reversible model (see ``NOTES.md``), and the column labelled k3 is
    renamed accordingly. The 50 C row supplies ``k_1``, ``k_2`` and ``k_5``
    of the ``'fig11b'`` constant set of that model.
    """
    with open(os.path.join(reference_data_dir, '11-4_Arrhenius.csv'),
              encoding='utf-8-sig', newline='') as f:
        rows = list(csv.reader(f))[2:]      # two header rows
    tables = {'points': _block(rows, 0, ['inv_T', 'ln_k_1', 'ln_k_2',
                                         'ln_k_5'])}
    for first, name in ((4, 'k_1'), (6, 'k_2'), (8, 'k_5')):
        tables['line_' + name] = _block(rows, first, ['inv_T', 'ln_' + name])
    return tables


def _read_fig11_digitized(filename):
    df = pd.read_csv(os.path.join(reference_data_dir, filename))
    df = df.rename(columns={'t_min': 'time'})
    contours = df[df['series'] == 'contour']
    stars = df[df['series'].isin(['red_star', 'purple_star'])]
    return {
        'contours': contours[['level', 'branch', 'time', 'T_C']].reset_index(
            drop=True),
        'ridge': df.loc[df['series'] == 'ridge', ['time', 'T_C']].reset_index(
            drop=True),
        'stars': stars.set_index('series')[['time', 'T_C']].rename_axis(None),
        }


def read_fig11a_selectivity_digitized():
    """Read the lines digitized from the published Fig. 11a (selectivity map).

    Returns
    -------
    dict of pandas.DataFrame
        ``'contours'``
            Samples along the 0.1, 0.3 and 0.5 contour lines of the main
            panel: ``level``, ``branch`` (``'early'``: left of the red
            dashed line, ``'late'``: right of it), ``time`` (plotted time,
            min) and ``T_C`` (C).
        ``'ridge'``
            Samples along the red dashed line: ``time`` and ``T_C``.
        ``'stars'``
            Centres of the red and purple stars, indexed ``'red_star'``
            and ``'purple_star'``: ``time`` and ``T_C``.

    Notes
    -----
    Digitized from the published figure, not deposited data; every second
    sample of the digitization is kept, and samples the digitization marked
    as uncertain (crowded lines by the left frame and around the red star)
    are left out. Reading uncertainty: about 0.3 px normal to a line, i.e.
    1 min or 0.05 C. The plotted time is the reaction time plus 10 min.
    """
    return _read_fig11_digitized('fig11a_selectivity_digitized.csv')


def read_fig11b_yield_digitized():
    """Read the lines digitized from the published Fig. 11b (yield map).

    Returns
    -------
    dict of pandas.DataFrame
        ``'contours'``, ``'ridge'`` and ``'stars'``, as returned by
        :func:`read_fig11a_selectivity_digitized`.

    Notes
    -----
    See :func:`read_fig11a_selectivity_digitized`.
    """
    return _read_fig11_digitized('fig11b_yield_digitized.csv')
