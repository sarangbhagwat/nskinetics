# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests for the r14-r17 rate laws, affinity anchoring and Ehrlich
stoichiometry bookkeeping of the shipped *S. cerevisiae* ethanol/isobutanol
model (see
``docs/superpowers/specs/2026-09-15-r16-split-aro10-adh6-design.md`` and its
2026-09-13 predecessor).

The Ehrlich branch's terminal step is split in two: r16 is the Aro10
2-keto-acid decarboxylase (irreversible Michaelis-Menten in KIV, releasing
isobutyraldehyde + CO2, no cofactor and no cross-product terms, mirroring
r3/Pdc), and r17 is the Adh6 NADPH reductase (a full structural mirror of
r6/Adh1: Haldane reverse ``k_17r``, competitive isobutanol ``K_17e``, and
acetate/ethanol cross-inhibition ``k_17ia``/``k_17ie``). ``K_17``/``K_17e``
are anchored to the fitted ``K_6``/``K_6e`` by the in-vitro Adh6:Adh1 Km
ratio. ``k_16r``/``K_16i`` and the lumped step's cross-product coefficients
``k_16ia``/``k_16ie`` remain declared at 0 and inert on the clean
decarboxylase, so the isobutanol repo's workbooks, which ``exec``
``K_16i = x`` / ``k_16ia = x`` / ``k_16ie = x`` onto the model, keep
resolving.

Marked ``slow`` as a whole and deliberately not registered in
``nskinetics/tests/__init__.py``: the shipped subpackage builds ``te_r`` on
import, and ``test_processes_contract.test_import_is_lightweight`` asserts
that nothing named for the shipped model is in ``sys.modules`` during the
``-m "not slow"`` run. Every import of it (and of libsbml) therefore lives
inside a test body. Never load a second Antimony model in this process: doing
so perturbs the shipped model. Every test that applies scenario B or writes a
parameter restores scenario A and the default in ``finally`` --
``reset_nsk_kinetic_model`` resets state, not parameters.
"""

import contextlib
import os

import pytest

pytestmark = pytest.mark.slow

# --- expected coefficients ---------------------------------------------------
# The split is taken at the molar node, so at steady flux it is mass-identical
# to the lumped r16 it replaces: 0.621 * 1.028 = 0.638 g IBO and
# 0.621 * 0.222 = 0.138 g $Red per g KIV, the lumped law's own coefficients.

#: g isobutyraldehyde per g KIV on r16: 72.11/116.12 (the rest leaves as CO2).
R16_ALD = 0.621

#: g Red per g isobutyraldehyde on r17: one NADPH per mole aldehyde,
#: 16/72.11 (Lei's $Red convention: 16 g O-equivalents per mole).
R17_RED = 0.222

#: ``{reaction: (reactants, products)}`` as ``{species_id: stoichiometry}``.
#: CO2: r13 44.01/(2*88.06) per g pyruvate, r16 44.01/116.12 per g KIV;
#: Red on r14: one NADPH per mole acetolactate, 16/132.11.
STOICH = {
    'r13': ({'s_pyr': 1.0}, {'s_AL': 0.750, 'CO2': 0.250}),
    'r14': ({'s_AL': 1.0, 'Red': 0.121}, {'s_DHI': 1.015}),
    'r15': ({'s_DHI': 1.0}, {'s_KIV': 0.866}),
    'r16': ({'s_KIV': 1.0}, {'s_isobutald': R16_ALD, 'CO2': 0.379}),
    'r17': ({'s_isobutald': 1.0, 'Red': R17_RED}, {'s_IBO': 1.028}),
}

#: Signed coefficients of ``qO2*x*env*32/1000`` (g O2-equivalents/h): the
#: rule subtracts the $Red consumers r6, r14 and r17. The NADPH credit moved
#: from r16 to r17 with the split, since r17 is the step that consumes it.
QO2_TERMS = {'r1': 0.178, 'r2': 0.908, 'r4': 0.363, 'r5': 1.066,
             'r6': -0.363, 'r17': -R17_RED, 'r14': -0.121,
             'r7': 0.063, 'r8': 0.214}

#: Coefficients of ``qCO2*x*env*44.01/1000`` (g CO2/h). The pre-existing
#: r3/r5 mismatches against the reactions (0.33 vs 0.5; 1.446 vs 1.466) are
#: Lei's and are reproduced, not corrected. r16 keeps the decarboxylation.
QCO2_TERMS = {'r2': 1.499, 'r3': 0.5, 'r5': 1.466, 'r7': 0.127, 'r8': 0.325,
              'r13': 0.250, 'r16': 0.379}

RATE_IDS = ('r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8',
            'r13', 'r14', 'r15', 'r16', 'r17')


# --- helpers ---------------------------------------------------------------

def _model():
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        te_r, reset_nsk_kinetic_model, apply_scenario_A, apply_scenario_B)
    return te_r, reset_nsk_kinetic_model, apply_scenario_A, apply_scenario_B


@contextlib.contextmanager
def _scenario_B_state(**species):
    """Reset the shipped model, apply scenario B, write ``species`` (g/L;
    ``env`` is 1 after a reset, so amount == concentration) and yield the raw
    RoadRunner object; restores scenario A and resets afterwards."""
    te_r, reset, apply_A, apply_B = _model()
    r = te_r._te
    try:
        reset(te_r)
        apply_B(te_r)
        r.is_aerobic = 1
        for name, value in species.items():
            r[name] = value
        yield r
    finally:
        apply_A(te_r)
        reset(te_r)


def _sbml_model(sbml_string):
    import libsbml
    doc = libsbml.readSBMLFromString(sbml_string)
    assert doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) == 0
    return doc.getModel()


def _formula(math_container):
    import libsbml
    return libsbml.formulaToL3String(math_container.getMath())


def _reaction_sides(m, rid):
    rx = m.getReaction(rid)
    assert rx is not None, rid
    reactants = {s.getSpecies(): s.getStoichiometry()
                 for s in rx.getListOfReactants()}
    products = {s.getSpecies(): s.getStoichiometry()
                for s in rx.getListOfProducts()}
    return rx, reactants, products


# --- 1. defaults -----------------------------------------------------------

def test_defaults():
    te_r, *_ = _model()
    r = te_r._te
    assert r.K_13 == pytest.approx(0.1)
    assert r.K_14 == pytest.approx(0.017)
    assert r.K_15 == pytest.approx(0.080)
    assert r.K_16 == pytest.approx(0.27)
    assert r.K_16i == 0.0
    assert r.k_16r == 0.0
    # still declared: the isobutanol workbooks exec `K_16i = x` onto the model
    ids = set(r.getGlobalParameterIds())
    assert {'K_16i', 'k_16r'} <= ids
    # Adh6 (r17): K_17 anchors to the fitted K_6 by the Adh6:Adh1 aldehyde-Km
    # ratio with the molar-mass correction, K_17e to K_6e through the same
    # 0.70 in-vitro -> in-vivo transfer factor r6 itself implies.
    assert r.K_17 == pytest.approx(0.0086)
    assert r.K_17e == pytest.approx(0.020)
    assert r.k_17r == pytest.approx(0.00025)
    assert r.k_17ia == pytest.approx(0.06)
    assert r.k_17ie == pytest.approx(0.02)
    # Aro10 (k_16) and Adh6 (k_17) are both native constitutive enzymes, so
    # both vmaxes default nonzero; the branch is gated by the engineered
    # upstream block k_13-k_15, not by these terminal steps. Each vmax anchors
    # to its fitted ethanol-branch analog -- k_16 to the Pdc vmax k_3, k_17 to
    # the Adh1 vmax k_6 -- by the kcat, abundance and substrate-MW ratios (a
    # vmax is kcat x [E] x MW_sub); see the anchor tests below.
    assert r.k_16 == pytest.approx(0.02115)
    assert r.k_17 == pytest.approx(0.1077)
    assert {'k_17', 'K_17', 'k_17r', 'K_17e', 'k_17ia', 'k_17ie'} <= ids
    # the lumped law's own cross-product coefficients stay declared, at 0 and
    # inert, for the same workbook reason as K_16i/k_16r (rows 74-75 of the
    # scenario-B and opt_* workbooks exec them)
    assert {'k_16ia', 'k_16ie'} <= ids
    assert r.k_16ia == 0.0 and r.k_16ie == 0.0


def test_K_17_is_anchored_to_the_fitted_K_6():
    # The anchoring rule, recomputed from its inputs rather than restated:
    # K_17 = K_6 * (Km_Adh6,ald / Km_Adh1,ald) * (MW_isobutald / MW_acetald),
    # with Adh6 0.17 mM (Larroy 2002, 3-methylbutanal surrogate) and Adh1
    # 1.1 mM (Ganzhorn 1987).
    te_r, *_ = _model()
    r = te_r._te
    expected = r.K_6 * (0.17 / 1.1) * (72.11 / 44.05)
    assert r.K_17 == pytest.approx(expected, rel=0.03)
    # K_17e = K_17 / Ki_IBO,model, Ki_IBO,model = 0.70 * 50 * 0.17 mM in g/L;
    # the 0.70 is r6's own transfer factor, K_6/K_6e = 0.60 g/L vs 18 mM.
    ki_ibo = 0.70 * 50 * 0.17e-3 * 74.12          # mol/L -> g/L isobutanol
    assert r.K_17e == pytest.approx(r.K_17 / ki_ibo, rel=0.05)


def test_k_17_is_anchored_to_the_fitted_k_6_by_kcat_abundance_and_mw():
    # k_17 (Adh6 vmax) mirrors r6 the way K_17/k_17r do, but a vmax is
    # kcat x [E] x MW_substrate, so the transfer from the fitted Adh1 vmax k_6
    # carries all three ratios -- not the kcat ratio alone, and not the old
    # (296/19) x scenario-B k_16 basis:
    #   kcat_Adh6/kcat_Adh1 = 296/1800      (Larroy 2002 / Ganzhorn 1987)
    #   [Adh6]/[Adh1]       = 14717/103727  (SGD proteomics medians, ~7:1)
    #   MW_ibald/MW_acetald = 72.11/44.05   (the same MW term K_17 carries)
    te_r, *_ = _model()
    r = te_r._te
    expected = r.k_6 * (296.0/1800.0) * (14717.0/103727.0) * (72.11/44.05)
    assert r.k_17 == pytest.approx(expected, rel=0.01)


def test_k_16_is_anchored_to_the_fitted_k_3_by_kcat_abundance_and_mw():
    # k_16 (Aro10 vmax) is set exactly as k_17 was: a vmax is
    # kcat x [E] x MW_substrate, so the transfer from the fitted Pdc vmax k_3
    # to the native constitutive Aro10 carries all three ratios, staying
    # consistent with K_16's own anchoring to K_3:
    #   kcat_Aro10/kcat_Pdc1 = 19/60        (Kneen 2011 / Balakrishnan 2012)
    #   [Aro10]/[Pdc1]       = 5068/581219  (SGD proteomics medians, ~1:115)
    #   MW_KIV/MW_pyr        = 116.12/88.06 (the same MW term K_16 carries)
    te_r, *_ = _model()
    r = te_r._te
    expected = r.k_3 * (19.0/60.0) * (5068.0/581219.0) * (116.12/88.06)
    assert r.k_16 == pytest.approx(expected, rel=0.01)


def test_r16_ignores_its_retired_cross_product_coefficients():
    # k_16ia/k_16ie are declared only for the isobutanol workbooks that exec
    # them; the decarboxylase law has no carrier for them, so a workbook
    # writing 0.06/0.02 (its baseline rows) must not touch r16 -- and r17,
    # which carries the live k_17ia/k_17ie, must not read them either.
    with _scenario_B_state(s_KIV=0.1, s_isobutald=0.01, s_acetate=1.0,
                           s_EtOH=20.0, s_IBO=26.0, x=10.0) as r:
        try:
            r16_default, r17_default = r['r16'], r['r17']
            r.k_16ia = 0.06
            r.k_16ie = 0.02
            r16_armed, r17_armed = r['r16'], r['r17']
        finally:
            r.k_16ia = 0.0
            r.k_16ie = 0.0
    assert r16_default > 0.0
    assert r16_armed == pytest.approx(r16_default, rel=1e-12)
    assert r17_armed == pytest.approx(r17_default, rel=1e-12)


# --- 2. r16 is insensitive to the retired terms -----------------------------

def test_r16_ignores_K_16i_and_k_16r():
    # s_KIV = 0.1 < 0.0125*s_IBO = 0.325: under the old law this state parked
    # r16 at a negative rate (the k_16r reverse term) and the K_16i term
    # suppressed the forward rate. Under the new law neither parameter reaches
    # the rate, and the rate is strictly positive.
    with _scenario_B_state(s_KIV=0.1, s_IBO=26.0, x=10.0) as r:
        try:
            r.K_16i = 0.0
            r.k_16r = 0.0
            rate_default = r['r16']
            r.K_16i = 0.057
            r.k_16r = 0.0125
            rate_armed = r['r16']
        finally:
            r.K_16i = 0.0
            r.k_16r = 0.0
    assert rate_default > 0.0
    assert rate_armed == pytest.approx(rate_default, rel=1e-12)


# --- 2b. the split: a clean decarboxylase and a full ADH mirror --------------

def test_r16_is_a_clean_decarboxylase_law():
    # Aro10 mirrors r3/Pdc: no cofactor, no reverse, no product or
    # cross-product terms. Everything ADH-specific lives on r17 now.
    te_r, *_ = _model()
    m = _sbml_model(te_r._te.getSBML())
    law = _formula(m.getReaction('r16').getKineticLaw())
    for absent in ('k_16r', 'K_16i', 'k_16ia', 'k_16ie', 'k_17ia', 'k_17ie',
                   's_IBO', 's_EtOH', 's_acetate', 'Red'):
        assert absent not in law, absent
    assert 'k_16' in law and 'K_16' in law and 's_KIV' in law


def test_r17_mirrors_the_r6_ADH_law():
    # Structural mirror of r6 with isobutanol as the self-product: the same
    # Haldane numerator, competitive denominator and two cross-product
    # exponentials, on the aldehyde r16 now makes.
    te_r, *_ = _model()
    m = _sbml_model(te_r._te.getSBML())
    law = _formula(m.getReaction('r17').getKineticLaw())
    for present in ('k_17', 'k_17r', 'K_17', 'K_17e', 'k_17ia', 'k_17ie',
                    's_isobutald', 's_IBO', 's_acetate', 's_EtOH'):
        assert present in law, present
    # r17 must not inhibit itself on the product it does not make
    assert 'k_17ii' not in law
    assert m.getReaction('r17').getReversible() is True


def test_r17_reverse_and_product_terms_are_live():
    # Unlike the retired k_16r/K_16i, these reach the rate: the ADH step is a
    # real carrier for them. Raising K_17e must slow r17 at a high titer.
    with _scenario_B_state(s_isobutald=0.01, s_IBO=26.0, x=10.0) as r:
        try:
            base = r['r17']
            r.K_17e = 0.20
            suppressed = r['r17']
            r.K_17e = 0.020
            r.k_17r = 0.05
            reversed_ = r['r17']
        finally:
            r.K_17e = 0.020
            r.k_17r = 0.00025
    assert base > 0.0
    assert suppressed < base
    # k_17r * s_IBO = 1.3 > s_isobutald = 0.01, so the net rate goes negative
    assert reversed_ < 0.0


# --- 3. stoichiometry -------------------------------------------------------

def test_stoichiometry_and_reversibility():
    te_r, *_ = _model()
    m = _sbml_model(te_r._te.getSBML())
    for rid, (reactants, products) in STOICH.items():
        rx, got_reactants, got_products = _reaction_sides(m, rid)
        assert got_reactants == pytest.approx(reactants), rid
        assert got_products == pytest.approx(products), rid
        assert rx.getReversible() is (rid == 'r17'), rid
    # r6 and r17 are the model's two reversible steps (Haldane ADH forms)
    assert m.getReaction('r6').getReversible() is True
    # the retired terms are gone from the r16 law itself, not just zeroed
    law = _formula(m.getReaction('r16').getKineticLaw())
    assert 'k_16r' not in law and 'K_16i' not in law


def test_the_split_conserves_the_lumped_coefficients():
    # At steady flux through the node, r16 + r17 must move exactly what the
    # single lumped r16 moved per g KIV: 0.638 g isobutanol and 0.138 g $Red.
    te_r, *_ = _model()
    m = _sbml_model(te_r._te.getSBML())
    _, r16_in, r16_out = _reaction_sides(m, 'r16')
    _, r17_in, r17_out = _reaction_sides(m, 'r17')
    ald_per_kiv = r16_out['s_isobutald']
    assert ald_per_kiv * r17_out['s_IBO'] == pytest.approx(0.638, abs=5e-4)
    assert ald_per_kiv * r17_in['Red'] == pytest.approx(0.138, abs=5e-4)
    # and the decarboxylation CO2 is unchanged by the split
    assert r16_out['CO2'] == pytest.approx(0.379)
    # carbon closes on r16: isobutyraldehyde + CO2 = 1 g per g KIV
    assert ald_per_kiv + r16_out['CO2'] == pytest.approx(1.0, abs=1e-3)


def test_isobutaldehyde_is_an_internal_intermediate():
    # Like s_AL / s_DHI / s_KIV: made and consumed inside the branch, with no
    # dilution outflow of its own (D = 0 in fed-batch use, but the asymmetry
    # would still be wrong) and no biosteam chemical mapping.
    te_r, *_ = _model()
    r = te_r._te
    assert 's_isobutald' in set(r.getFloatingSpeciesIds())
    assert 's_isobutald_out' not in set(r.getReactionIds())
    m = _sbml_model(r.getSBML())
    consumers = [rx.getId() for rx in m.getListOfReactions()
                 if any(s.getSpecies() == 's_isobutald'
                        for s in rx.getListOfReactants())]
    producers = [rx.getId() for rx in m.getListOfReactions()
                 if any(s.getSpecies() == 's_isobutald'
                        for s in rx.getListOfProducts())]
    assert producers == ['r16'] and consumers == ['r17']


# --- 4. rule coefficients, behaviourally ------------------------------------

def test_qO2_and_qCO2_rules_carry_the_expected_coefficients():
    # Every reaction rate nonzero (r6 may be negative here: acetaldehyde is
    # below k_6r*s_EtOH, which is fine -- the rule is linear in the rates).
    with _scenario_B_state(s_glu=50.0, s_pyr=0.05, s_acetate=1.0,
                           s_acetald=0.05, s_EtOH=20.0, s_AL=0.01,
                           s_DHI=0.02, s_KIV=0.1, s_IBO=26.0, x=10.0) as r:
        rates = {rid: float(r[rid]) for rid in RATE_IDS}
        x, env = float(r['[x]']), float(r.env)
        qO2, qCO2 = float(r.qO2), float(r.qCO2)
    for rid, v in rates.items():
        assert v != 0.0, rid
    expected_O2 = sum(c * rates[rid] for rid, c in QO2_TERMS.items())
    expected_CO2 = sum(c * rates[rid] for rid, c in QCO2_TERMS.items())
    assert qO2 * x * env * 32 / 1000 == pytest.approx(expected_O2, rel=1e-9)
    assert qCO2 * x * env * 44.01 / 1000 == pytest.approx(expected_CO2, rel=1e-9)


# --- 5. scenario dicts ------------------------------------------------------

def test_scenario_dicts_carry_only_the_four_capacities():
    from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
        SCENARIO_A_EHRLICH, SCENARIO_B_EHRLICH)
    keys = {'k_13', 'k_14', 'k_15', 'k_16'}
    assert set(SCENARIO_A_EHRLICH) == keys
    assert set(SCENARIO_B_EHRLICH) == keys
    assert all(v == 0.0 for v in SCENARIO_A_EHRLICH.values())
    assert SCENARIO_B_EHRLICH == {'k_13': 5.81, 'k_14': 4.8, 'k_15': 4.8,
                                  'k_16': 2.82}


# --- 6. shipped SBML in sync ------------------------------------------------

def test_shipped_sbml_matches_the_live_model():
    import libsbml
    import nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo as pkg
    te_r, *_ = _model()
    xml = os.path.join(os.path.dirname(os.path.abspath(pkg.__file__)),
                       's_cerevisiae_ferm_fb_inhib_mod_ibo_sbml.xml')
    doc = libsbml.readSBML(xml)
    assert doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) == 0
    shipped = doc.getModel()
    live = _sbml_model(te_r._te.getSBML())      # as loaded, scenario-invariant
    for pid in ('K_13', 'K_14', 'K_15', 'K_16', 'K_16i', 'k_16r', 'k_16',
                'k_17', 'K_17', 'k_17r', 'K_17e', 'k_17ia', 'k_17ie',
                'anaerobic_growth_mult'):
        assert shipped.getParameter(pid).getValue() == pytest.approx(
            live.getParameter(pid).getValue()), pid
    for rid in ('r13', 'r14', 'r15', 'r16', 'r17'):
        s_rx, l_rx = shipped.getReaction(rid), live.getReaction(rid)
        assert s_rx.getReversible() == l_rx.getReversible(), rid
        assert _formula(s_rx.getKineticLaw()) == _formula(l_rx.getKineticLaw()), rid
        _, sr, sp = _reaction_sides(shipped, rid)
        _, lr, lp = _reaction_sides(live, rid)
        assert sr == pytest.approx(lr) and sp == pytest.approx(lp), rid
    for var in ('qO2', 'qCO2'):
        assert _formula(shipped.getRule(var)) == _formula(live.getRule(var)), var
