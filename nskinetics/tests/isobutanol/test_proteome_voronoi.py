# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Offline pure-logic tests of ``plots/plot_proteome_voronoi.py`` (no
``load()``, no simulation, no Node render).

Loads the Stage-1 plotter by file path (it loads ``ps`` and ``eb`` by path)
and checks the baseline tile's allocation against the ``eb`` constants and
the spec's verification numbers, the derated / infeasible tiles, the
document shape and the JSON writer. The baseline record is read from the
scenario-A workbook, so every test needs ``openpyxl``.

Converted from the script ``analyses/test_proteome_voronoi_offline.py``; its
35 checks are kept, grouped by tile / document.
"""
import json
import os

import pytest

from ._support import load_by_path, load_paths, skip_unless_modules

pytestmark = [pytest.mark.isobutanol, skip_unless_modules('pandas', 'openpyxl')]


@pytest.fixture(scope='module')
def pv():
    plots_dir = load_paths().PLOTS_DIR
    return load_by_path('plot_proteome_voronoi',
                        os.path.join(plots_dir, 'plot_proteome_voronoi.py'))


@pytest.fixture(scope='module')
def eb(pv):
    return pv.eb


@pytest.fixture(scope='module')
def rec(pv):
    """The scenario-A baseline record."""
    return pv.ps.baseline_set()


def _synthetic(eb, rec, phi_m, phi_t_demand, label):
    """A record shaped like a set row with all metabolic pool in Glycolysis
    (r1) so category-sum Phi_M == phi_m, and burden_factor = the model's
    derating d = clip((F_flex - Phi_M)/phi_T, 0, 1)."""
    s = dict(rec)
    for st in eb.STEP_ORDER:
        s[f'pool_{st}'] = 0.0
    s['pool_r1'] = phi_m
    s['Phi_M'] = phi_m
    s['phi_T'] = phi_t_demand
    s['F_flex'] = eb.F_FLEX
    d = (eb.F_FLEX - phi_m) / phi_t_demand
    s['burden_factor'] = max(0.0, min(1.0, d))
    s['label'] = label
    s['is_baseline'] = False
    return s


def test_baseline_tile(pv, eb, rec):
    """The baseline tile: four top-level pieces; housekeeping 0.5*P = 0.245;
    translation = record phi_T ~= PHI_T_WT; four metabolic categories summing
    to the record's Phi_M; the isobutanol category = the native Adh6 pool r17
    (Ehrlich r13-r16 off at A); slack = F_flex - Phi_M - phi_T >= 0; cells
    sum to PROTEIN_CONTENT (0.49); no warning."""
    tile = pv.tile_from_record(rec)
    kids = {c['piece']: c for c in tile['children']}
    PC = eb.PROTEIN_CONTENT

    assert len(tile['children']) == 4
    assert abs(kids['housekeeping']['value'] - PC * eb.HOUSEKEEPING_FRACTION) < 1e-9
    assert abs(kids['housekeeping']['value'] - 0.245) < 1e-6
    assert abs(kids['translation']['value'] - float(rec['phi_T'])) < 1e-9
    assert abs(kids['translation']['value'] - eb.PHI_T_WT) < 1e-3   # ~0.110 (g=1)

    metabolic = kids['metabolic']
    mkids = {c['piece']: c for c in metabolic['children']}
    assert len(metabolic['children']) == 4
    assert (abs(metabolic['value'] - sum(c['value'] for c in metabolic['children']))
            < 1e-12)
    assert abs(metabolic['value'] - float(rec['Phi_M'])) < 1e-4
    # The Ehrlich isobutanol enzymes (r13-r16) are off in scenario A, but since
    # the 2026-09-15 r16/r17 split the isobutanol category also carries the NATIVE
    # Adh6 step r17 (constitutively expressed, so nonzero at baseline A): the
    # category equals the native r17 pool, not zero.
    assert abs(mkids['cat_isobutanol']['value'] - float(rec['pool_r17'])) < 1e-9
    assert float(rec['pool_r17']) > 1e-9                 # Adh6 constitutive

    slack = kids['slack']['value']
    assert abs(slack - (eb.F_FLEX - metabolic['value']
                        - kids['translation']['value'])) < 1e-9
    assert slack > -1e-9

    total = (kids['housekeeping']['value'] + metabolic['value']
             + kids['translation']['value'] + slack)
    assert abs(total - PC) < 1e-6
    assert 'warning' not in tile


def test_derated_growth_tile(pv, eb, rec):
    """A derated-growth trial (feasible: Phi_M < F_flex, but Phi_M + demand >
    F_flex so d < 1): the translation cell is the derated allocation
    (F_flex - Phi_M < demand), slack 0, no warning, cells still sum to 0.49."""
    PC = eb.PROTEIN_CONTENT
    dr = pv.tile_from_record(_synthetic(eb, rec, 0.20, 0.10, 'Derated'))
    dk = {c['piece']: c for c in dr['children']}
    assert abs(dk['translation']['value'] - (eb.F_FLEX - 0.20)) < 1e-9
    assert dk['translation']['value'] < 0.10 - 1e-9
    assert abs(dk['slack']['value']) < 1e-9
    assert 'warning' not in dr
    assert abs(sum(c['value'] for c in dr['children']) - PC) < 1e-6


def test_infeasible_tile(pv, eb, rec):
    """A truly infeasible trial (Phi_M >= F_flex, d = 0): no room for any
    translation; slack clamped to 0 and a warning."""
    inf = pv.tile_from_record(_synthetic(eb, rec, 0.30, 0.10, 'Infeasible'))
    ik = {c['piece']: c for c in inf['children']}
    assert abs(ik['translation']['value']) < 1e-12
    assert abs(ik['slack']['value']) < 1e-12
    assert 'warning' in inf


def test_document_shape(pv, eb, rec):
    """meta.protein_content 0.49, meta.F_flex == eb.F_FLEX, one tile,
    7-entry piece_order and legend: the four metabolic categories first (in
    ps.BURDEN_CATEGORIES order, each with its steps), then unallocated /
    translation / housekeeping; the legend steps partition eb.STEP_ORDER."""
    doc = pv.build_document([rec], band_campaign=None)
    assert abs(doc['meta']['protein_content'] - 0.49) < 1e-9
    assert abs(doc['meta']['F_flex'] - eb.F_FLEX) < 1e-9
    assert len(doc['tiles']) == 1
    assert len(doc['meta']['piece_order']) == 7
    leg = doc['meta']['legend']
    assert len(leg) == 7
    assert ([e['piece'] for e in leg[:4]]
            == [pv.CATEGORY_PIECES[n] for n, _ in pv.ps.BURDEN_CATEGORIES]
            and all(e['steps'] for e in leg[:4]))
    assert ([e['piece'] for e in leg[4:]] == ['slack', 'translation', 'housekeeping']
            and all(e['steps'] is None for e in leg[4:]))
    assert sorted(st for e in leg[:4] for st in e['steps']) == sorted(eb.STEP_ORDER)


def test_tile_base_colors(pv, eb, rec):
    """Per-tile base colour matches the parameter-sets figure: the baseline
    tile is grey, the first campaign tile takes ps.HUE_COLORS[0]."""
    doc = pv.build_document([rec], band_campaign=None)
    assert doc['tiles'][0]['base_color'] == pv.ps.BASELINE_COLOR
    two = pv.build_document([pv.ps.baseline_set(),
                             _synthetic(eb, rec, 0.20, 0.10, 'Campaign')],
                            band_campaign=None)
    assert two['tiles'][0]['base_color'] == pv.ps.BASELINE_COLOR
    assert two['tiles'][1]['base_color'] == pv.ps.HUE_COLORS[0]


def test_write_document_round_trips_json(pv, tmp_path):
    """write_document writes a valid JSON document (the CLI's --no-render
    output) whose tile cells sum to 0.49."""
    # a baseline-only run needs --set-free defaults dropped; force one tile via
    # a direct document write instead of the study-CSV default path:
    td = str(tmp_path)
    doc = pv.build_document([pv.ps.baseline_set()], band_campaign=None)
    path = pv.write_document(doc, td, 'unit', 'stamp')
    with open(path, encoding='utf-8') as fh:
        loaded = json.load(fh)
    assert loaded['tiles'][0]['label'] == pv.ps.baseline_set()['label']
    assert abs(sum((c['value'] for c in loaded['tiles'][0]['children']))
               - 0.49) < 1e-6
