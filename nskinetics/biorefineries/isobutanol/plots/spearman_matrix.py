# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Ellipse correlation matrix of the Spearman's rank correlations in an
uncertainty-analysis workbook (``*_1_full_evaluation.xlsx``, sheets
'Spearman' and 'Spearman p-values').

Run as a script to draw the newest scenario-A workbook in the results root::

    python plots/spearman_matrix.py [WORKBOOK] [--out-dir DIR]

Created on Fri Feb 20 01:51:07 2026

@author: saran
"""
import argparse
import glob
import importlib.util
import os

import pandas as pd
import matplotlib.pyplot as plt
import contourplots

__all__ = ('plot_spearman_matrix',)


def _strip_units(names):
    """'k_1 [1/h]' -> 'k_1 ' (the text before the first '['), per name."""
    return [n[:n.index('[')] if '[' in n else n for n in names]


#%%
def plot_spearman_matrix(data, data_p, remove_units=True):
    """Draw the ellipse correlation matrix of a Spearman sheet.

    Parameters
    ----------
    data : str or pandas.DataFrame
        Workbook path (its 'Spearman' sheet is read) or the sheet's frame.
    data_p : str or pandas.DataFrame
        Workbook path (its 'Spearman p-values' sheet is read) or the frame.
    remove_units : bool, optional
        Strip the ' [unit]' suffix from the parameter names. Frames passed in
        are never modified.
    """
    df, df_p = data, data_p
    if not isinstance(data, pd.DataFrame):
        df = pd.read_excel(data, sheet_name='Spearman')
    if not isinstance(data_p, pd.DataFrame):
        df_p = pd.read_excel(data_p, sheet_name='Spearman p-values')
    if remove_units:
        # whole-column assignment on copies (a chained df['col'][i] = ...
        # never reaches df under pandas 3 copy-on-write)
        df = df.copy()
        df_p = df_p.copy()
        df['Parameter'] = _strip_units(df['Parameter'])
        df_p['Parameter'] = _strip_units(df_p['Parameter'])
    contourplots.ellipse_correlation_matrix_plot(df=df, param_col='Parameter',
                                                 category_col='Element', corr_prefix='Correlation with ',
                                                 small_corr_threshold=0.05,
                                                 exclude_parameters_all_below_threshold=True,
                                                 legend_max_cols=2,
                                                 format_kinetic_param_names=True,
                                                 )


#%%
def _load_paths():
    """The package's stdlib-only paths.py, loaded by file path (sim-safe)."""
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'paths.py')
    spec = importlib.util.spec_from_file_location('_ibo_paths', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _newest_workbook(paths, scenario='A'):
    pattern = os.path.join(glob.escape(paths.results_dir()), '*_1_full_evaluation.xlsx')
    files = [f for f in glob.glob(pattern) if f"['{scenario}']" in os.path.basename(f)]
    if not files:
        paths.require_input(
            paths.results_dir(f"*['{scenario}']*_1_full_evaluation.xlsx"),
            f'the scenario-{scenario} uncertainty analysis '
            '(analyses/uncertainty/uncertainties_IBO_EtOH.py)')
    return max(files, key=os.path.getmtime)


def main(argv=None):
    paths = _load_paths()
    ap = argparse.ArgumentParser(description='Spearman ellipse correlation matrix.')
    ap.add_argument('workbook', nargs='?', default=None,
                    help='*_1_full_evaluation.xlsx (default: the newest '
                         'scenario-A workbook in the results root)')
    ap.add_argument('--out-dir', default=paths.publication_dir('Uncertainty'))
    args = ap.parse_args(argv)
    data = args.workbook or _newest_workbook(paths)
    plot_spearman_matrix(data, data)
    os.makedirs(args.out_dir, exist_ok=True)
    out = os.path.join(args.out_dir, 'spearman_matrix_A_0.05.png')
    plt.savefig(
        fname=out,
        transparent=False,
        facecolor='white',
        bbox_inches='tight',
        dpi=600,)
    print(f'Saved {out}')
    return out


if __name__ == '__main__':
    main()
