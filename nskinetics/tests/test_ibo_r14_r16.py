# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests for the r14-r16 rate laws, affinity anchoring and Ehrlich
stoichiometry bookkeeping of the shipped *S. cerevisiae* ethanol/isobutanol
model (see
``docs/superpowers/specs/2026-09-13-r14-r16-rate-laws-affinity-anchoring-design.md``).

r16 (lumped 2-keto-acid decarboxylase + ADH) is irreversible Michaelis-Menten
in KIV: the ``k_16r`` reverse term and the ``K_16i`` isobutanol product term
copied from r6 are gone from the law but stay declared at 0 so that the
isobutanol repo's workbooks, which ``exec`` ``K_16i = x`` onto the model, keep
resolving. ``K_14``/``K_15``/``K_16`` are anchored to the fitted core, and the
Ehrlich branch's NAD(P)H and CO2 are bookkept in the reactions and in the
``qO2``/``qCO2`` assignment rules.

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
# Commit 1 of the spec lands the r16 law change with r16's pre-existing
# 0.363 $Red left in place and the rules untouched; commit 2 moves these four
# constants to the corrected bookkeeping. Nothing else in this module changes
# between the two commits.

#: g Red per g KIV on r16 (commit 1: r6's 0.363, uncorrected).
R16_RED = 0.363

#: ``{reaction: (reactants, products)}`` as ``{species_id: stoichiometry}``.
STOICH = {
    'r13': ({'s_pyr': 1.0}, {'s_AL': 0.750}),
    'r14': ({'s_AL': 1.0}, {'s_DHI': 1.015}),
    'r15': ({'s_DHI': 1.0}, {'s_KIV': 0.866}),
    'r16': ({'s_KIV': 1.0, 'Red': R16_RED}, {'s_IBO': 0.638}),
}

#: Signed coefficients of ``qO2*x*env*32/1000`` (g O2-equivalents/h).
QO2_TERMS = {'r1': 0.178, 'r2': 0.908, 'r4': 0.363, 'r5': 1.066,
             'r6': -0.363, 'r16': -R16_RED, 'r7': 0.063, 'r8': 0.214}

#: Coefficients of ``qCO2*x*env*44.01/1000`` (g CO2/h).
QCO2_TERMS = {'r2': 1.499, 'r3': 0.5, 'r5': 1.466, 'r7': 0.127, 'r8': 0.325}

RATE_IDS = ('r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8',
            'r13', 'r14', 'r15', 'r16')


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


# --- 3. stoichiometry -------------------------------------------------------

def test_stoichiometry_and_reversibility():
    te_r, *_ = _model()
    m = _sbml_model(te_r._te.getSBML())
    for rid, (reactants, products) in STOICH.items():
        rx, got_reactants, got_products = _reaction_sides(m, rid)
        assert got_reactants == pytest.approx(reactants), rid
        assert got_products == pytest.approx(products), rid
        assert rx.getReversible() is False, rid
    # r6 is the model's one reversible step (Haldane form on Adh1)
    assert m.getReaction('r6').getReversible() is True
    # the retired terms are gone from the r16 law itself, not just zeroed
    law = _formula(m.getReaction('r16').getKineticLaw())
    assert 'k_16r' not in law and 'K_16i' not in law
    assert 'k_16ia' in law and 'k_16ie' in law


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
    for pid in ('K_13', 'K_14', 'K_15', 'K_16', 'K_16i', 'k_16r',
                'anaerobic_growth_mult'):
        assert shipped.getParameter(pid).getValue() == pytest.approx(
            live.getParameter(pid).getValue()), pid
    for rid in ('r13', 'r14', 'r15', 'r16'):
        s, l = shipped.getReaction(rid), live.getReaction(rid)
        assert s.getReversible() == l.getReversible(), rid
        assert _formula(s.getKineticLaw()) == _formula(l.getKineticLaw()), rid
        _, sr, sp = _reaction_sides(shipped, rid)
        _, lr, lp = _reaction_sides(live, rid)
        assert sr == pytest.approx(lr) and sp == pytest.approx(lp), rid
    for var in ('qO2', 'qCO2'):
        assert _formula(shipped.getRule(var)) == _formula(live.getRule(var)), var
