# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Offline pure-logic + render tests of ``plots/_pathway_lever_map.py`` (no
``load()``, no simulation).

Loads the parent plotter (``ps`` = ``plots/plot_kin_opt_parameter_sets.py``,
which itself loads ``ko``/``eb`` by path), ``eb`` and the companion (``plm``)
by file path, builds a real scenario-A baseline record (from the A workbook,
needs ``openpyxl``) plus synthetic campaign records under the split_12d
layout, and checks the lever / outcome / proteome helpers, the
objective-keyed colour map, and that the figure renders to a temp PNG/PDF.

Converted from the script ``analyses/test_pathway_lever_map_offline.py``; its
24 checks are kept, grouped by the helper they exercise.
"""
import os

import numpy as np
import pytest

from ._support import load_by_path, load_paths, skip_unless_modules

pytestmark = pytest.mark.isobutanol

needs_workbooks = skip_unless_modules('pandas', 'openpyxl')


@pytest.fixture(scope='module')
def ps():
    """The parameter-sets plotter, loaded by path, in the split_12d layout."""
    plots_dir = load_paths().PLOTS_DIR
    ps = load_by_path('ps', os.path.join(plots_dir, 'plot_kin_opt_parameter_sets.py'))
    ps.use_split_12d_layout()   # populate the split_12d panel-B globals
    return ps


@pytest.fixture(scope='module')
def plm():
    plots_dir = load_paths().PLOTS_DIR
    return load_by_path('_pathway_lever_map',
                        os.path.join(plots_dir, '_pathway_lever_map.py'))


@pytest.fixture(scope='module')
def base(ps):
    """Real scenario-A record (no simulation)."""
    return ps.baseline_set()


def _campaign(base, label, **over):
    """A synthetic campaign record: copy of the baseline with selected levers
    and outcomes overridden. is_baseline False so it takes a hue."""
    s = dict(base)
    s['label'] = label
    s['campaign'] = 'kin_opt_ethanol_isobutanol_metabolic_split_12d_ibo_titer_gp'
    s['is_baseline'] = False
    s['objective'] = 'IBO titer'
    s['trial_number'] = 1
    s.update(over)
    return s


@pytest.fixture(scope='module')
def sets(base):
    """The set list used by the colour-map checks and the render."""
    return [base,
            _campaign(base, 'IBO titer', objective='IBO titer',
                      k_13=3.2, k_14=1.4, k_15=1.2, k_16_rel=1.5, glycolysis=1.6,
                      inhib_ethanol=0.8, **{'IBO titer': 35.0, 'EtOH titer': 40.0,
                                            'IRR': 0.10}),
            _campaign(base, 'EtOH titer', objective='EtOH titer',
                      k_6_rel=2.1, glycolysis=1.8,
                      **{'IBO titer': 1.0, 'EtOH titer': 190.0, 'IRR': 0.05}),
            _campaign(base, 'Financial', objective='PI',
                      k_13=1.0, k_16_rel=1.2, glycolysis=1.3,
                      **{'IBO titer': 33.0, 'EtOH titer': 63.0, 'IRR': 0.27})]


@needs_workbooks
def test_lever_value(plm, base):
    """Baseline fold levers are 1.0; the absolute Ehrlich rate k_13 is 0.0; a
    missing key gives nan."""
    assert abs(plm.lever_value(base, 'glycolysis') - 1.0) < 1e-9
    assert abs(plm.lever_value(base, 'k_3_rel') - 1.0) < 1e-9
    assert abs(plm.lever_value(base, 'k_13')) < 1e-12
    assert np.isnan(plm.lever_value(base, 'nonexistent'))


@needs_workbooks
def test_is_off_in_wildtype(plm, base):
    """k_13 / k_14 are off in the wild type; an installed k_13 is not."""
    assert plm.is_off_in_wildtype(base, 'k_13') is True
    assert plm.is_off_in_wildtype(base, 'k_14') is True
    inst = _campaign(base, 'IBO', k_13=3.0, k_14=1.2, k_15=1.1)
    assert plm.is_off_in_wildtype(base, 'k_13') and not plm.is_off_in_wildtype(inst, 'k_13')


def test_clamp_outcome(ps, plm):
    """IRR treatment: a finite loss and -inf clamp to the floor, nan is
    omitted, a positive IRR passes through; a non-clamped titer passes
    through unchanged."""
    cn = ps.CLAMP_NEG_TO_ZERO
    assert plm.clamp_outcome('IRR', -0.2, 0.0, cn) == 0.0
    assert plm.clamp_outcome('IRR', float('-inf'), 0.0, cn) == 0.0
    assert np.isnan(plm.clamp_outcome('IRR', np.nan, 0.0, cn))
    assert abs(plm.clamp_outcome('IRR', 0.27, 0.0, cn) - 0.27) < 1e-12
    assert plm.clamp_outcome('IBO titer', -1.0, 0.0, cn) == -1.0


@needs_workbooks
def test_proteome_segments(ps, plm, base):
    """Segments sum to PROTEIN_CONTENT (0.49); housekeeping == 0.5*PC; the
    baseline isobutanol segment is the native Adh6 pool r17 (Ehrlich off);
    isobutanol + rest == Phi_M."""
    eb = ps.eb
    seg = plm.proteome_segments(base, eb, ps.BURDEN_CATEGORIES)
    PC = float(eb.PROTEIN_CONTENT)
    assert abs(sum(seg.values()) - PC) < 1e-6
    assert abs(seg['housekeeping'] - PC * 0.5) < 1e-9
    assert abs(seg['isobutanol'] - float(base['pool_r17'])) < 1e-9
    assert abs(seg['isobutanol'] + seg['rest'] - float(base['Phi_M'])) < 1e-9


@needs_workbooks
def test_objective_keyed_colour_map(ps, base, sets):
    """The baseline is grey; a financial record (objective 'PI') maps to the
    financial hue via the OUTCOME_OWN_OBJECTIVES alias, not by set order; the
    IBO-titer and EtOH-titer campaigns take their own hues."""
    cmap = ps.pathway_colors(sets)
    assert cmap[id(base)] == ps.BASELINE_COLOR
    assert cmap[id(sets[3])] == ps.HUE_COLORS[0]     # cyan
    assert cmap[id(sets[1])] == ps.HUE_COLORS[2]     # green
    assert cmap[id(sets[2])] == ps.HUE_COLORS[5]     # yellow


def test_pathway_objectives(ps):
    """PATHWAY_OBJECTIVES has the three story objectives."""
    assert ([o for _l, o in ps.PATHWAY_OBJECTIVES]
            == ['PI (log-tail)', 'IBO titer', 'EtOH titer'])


@needs_workbooks
def test_renders_png_and_pdf(ps, plm, sets, tmp_path):
    """The lever map renders to a temp PNG + PDF without error and returns
    out_stem."""
    colors = ps.pathway_colors(sets)
    helpers = ps.pathway_helpers()
    stem = os.path.join(str(tmp_path), 'lever_map')
    out = plm.plot_pathway_lever_map(sets, colors, stem, ko=ps.ko, eb=ps.eb,
                                     helpers=helpers, dpi=90)
    assert out == stem
    assert os.path.isfile(stem + '.png')
    assert os.path.isfile(stem + '.pdf')
