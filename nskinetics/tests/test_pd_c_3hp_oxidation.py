# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of the Pd/C 3-HP oxidation models (Kim et al., Appl. Catal. B:
Environ. Energy 386 (2026) 126403): ``nskinetics.models.pd_c_3hp_oxidation``.

The H2O2 and O2 models are checked against the reference data shipped with
the package: the model curves and experimental points the authors deposited
for Fig. 9, and the inset curves digitized from SI Fig. S12. Simulations are
evaluated at exactly the reference times; the deposited curves are coarse
solver polylines and must not be interpolated.

The tolerances are those of the reconstruction of these models. The
deposited curves carry the authors' solver error (MATLAB ode23s at default
tolerances) and Table 3 is rounded to three digits, so a converged solution
with the printed constants is expected to differ from them by up to 2e-3
(O2) and 1e-3 (H2O2), not to match them exactly. The limits are still tight
enough to fail a 2 % error in most rate constants, a missing catalyst-mass
or oxidant factor, the SI Table S4 reaction orders, or a missing reaction.

The two models are small and load in milliseconds, so this module is not
marked slow. It imports pytest at top level and is therefore not registered
in ``nskinetics/tests/__init__.py``.
"""
import hashlib
import os

import numpy as np
import pytest

import nskinetics as nsk
import nskinetics.models.pd_c_3hp_oxidation as pkg
from nskinetics.models.pd_c_3hp_oxidation import (
    NORMALIZED, SPECIES, load_h2o2_model, load_o2_model, simulate_at_times)

_PKG_DIR = os.path.dirname(os.path.abspath(pkg.__file__))
_NORMALIZED = list(NORMALIZED)
_CONCENTRATIONS = ['[%s]' % s for s in SPECIES]

# Paper Table 3 as printed. k_1 to k_6 in 1/(g_cat min), k_7 in 1/min.
TABLE_3_H2O2 = {'k_1': 5.31, 'k_2': 3.11e1, 'k_3': 7.32e-4, 'k_4': 1.15e0,
                'k_5': 3.03e-3, 'k_6': 4.10e-9, 'k_7': 6.48e-3}
TABLE_3_O2 = {'k_1_ref': 7.14, 'k_2_ref': 1.39e2, 'k_3_ref': 3.16e-4,
              'k_4_ref': 5.86e-2, 'k_5_ref': 1.13e-2, 'k_6_ref': 5.52e-3}
ACTIVATION_ENERGIES = {'Ea_1': 59781.2, 'Ea_2': 59781.2, 'Ea_3': 91594.6,
                       'Ea_4': 10000., 'Ea_5': 51788.73, 'Ea_6': 59549.8}
M_CAT = 0.2                     # g
X_3HP_FEED = 0.000339           # authors' mole fraction of the 60 mM feed
C0 = 0.06                       # M


@pytest.fixture(scope='module')
def h2o2_model():
    return load_h2o2_model()


@pytest.fixture(scope='module')
def o2_model():
    return load_o2_model()


@pytest.fixture(scope='module')
def fig9a():
    return pkg.read_fig9a_h2o2()


@pytest.fixture(scope='module')
def fig9b():
    return pkg.read_fig9b_o2()


def _rms_and_max(simulated, reference):
    diff = simulated[_NORMALIZED].to_numpy() - reference[_NORMALIZED].to_numpy()
    return np.sqrt((diff**2).mean(axis=0)), np.abs(diff).max(axis=0)


def _pooled_r2(observed, modeled):
    """Unweighted R2 over all species and rows, as in the authors' script."""
    observed = observed[_NORMALIZED].to_numpy()
    modeled = modeled[_NORMALIZED].to_numpy()
    ss_res = ((observed - modeled)**2).sum()
    ss_tot = ((observed - observed.mean())**2).sum()
    return 1. - ss_res/ss_tot


# --- 1. loaders -------------------------------------------------------------

def test_loaders_return_configured_kinetic_models():
    for load, species in ((load_h2o2_model, SPECIES + ('H2O2',)),
                          (load_o2_model, SPECIES)):
        km = load()
        assert isinstance(km, nsk.KineticModel)
        assert km._units == {'time': 'min', 'conc': 'M'}
        km.validate_units()
        assert tuple(km._te.getFloatingSpeciesIds()) == species
        assert km._te.integrator.absolute_tolerance == 1e-10
        assert km._te.integrator.relative_tolerance == 1e-9
        assert load() is not km         # every call builds a new model


def test_initial_state_is_the_fig9_feed(h2o2_model, o2_model):
    for km in (h2o2_model, o2_model):
        km.reset()
        assert km.get_value('[HP3]') == pytest.approx(C0)
        for product in ('OPA', 'MA', 'AAc', 'OXA', 'D'):
            assert km.get_value('[%s]' % product) == 0.
        assert km.get_value('m_cat') == M_CAT
        assert km.get_value('X_HP3') == pytest.approx(1.)
    assert h2o2_model.get_value('[H2O2]') == pytest.approx(15*C0)
    assert h2o2_model.get_value('eq_H2O2') == pytest.approx(15.)
    assert o2_model.get_value('yO2') == 1.
    assert o2_model.get_value('T') == o2_model.get_value('T_ref') == 323.15


def test_kinetic_model_simulate_runs(o2_model):
    o2_model.simulate(0., 1440., 1441, reset=True)
    df = o2_model.results_df
    assert len(df) == 1441
    assert df['[MA]'].max() == pytest.approx(0.0337, abs=2e-4)   # 56 % of C0
    o2_model.reset()


# --- 2. constants and rate laws ---------------------------------------------

def test_rate_constants_are_table_3_as_printed(h2o2_model, o2_model):
    for name, value in TABLE_3_H2O2.items():
        assert h2o2_model.get_value(name) == value, name
    for name, value in TABLE_3_O2.items():
        assert o2_model.get_value(name) == value, name
    for name, value in ACTIVATION_ENERGIES.items():
        assert o2_model.get_value(name) == value, name
    for km in (h2o2_model, o2_model):
        assert km.get_value('x_per_M') == pytest.approx(X_3HP_FEED/C0,
                                                        rel=1e-12)
    # the inferred oxidant reference factor of the O2 model
    assert o2_model.get_value('x_O2_ref') == pytest.approx(15*X_3HP_FEED,
                                                           rel=1e-12)


def test_o2_rate_laws():
    km = load_o2_model()
    hp3, opa, ma, y = 0.05, 0.004, 0.02, 0.9
    km.set_value('[HP3]', hp3)
    km.set_value('[OPA]', opa)
    km.set_value('[MA]', ma)
    km.set_value('yO2', y)
    k = [TABLE_3_O2['k_%d_ref' % i] for i in range(1, 7)]
    ox = 15*X_3HP_FEED*y            # oxidant factor on R1, R2 and R4 only
    expected = {'R1': k[0]*M_CAT*ox*hp3, 'R2': k[1]*M_CAT*ox*opa,
                'R3': k[2]*M_CAT*ma, 'R4': k[3]*M_CAT*ox*ma,
                'R5': k[4]*M_CAT*hp3, 'R6': k[5]*M_CAT*ma}
    for reaction, rate in expected.items():
        assert km.get_value(reaction) == pytest.approx(rate, rel=1e-12), reaction
    depletion = -(X_3HP_FEED/C0)*(expected['R1'] + expected['R2']
                                  + expected['R4'])
    assert km.get_value("yO2'") == pytest.approx(depletion, rel=1e-12)


def test_h2o2_rate_laws():
    km = load_h2o2_model()
    hp3, opa, ma, h2o2 = 0.05, 0.004, 0.02, 0.5
    km.set_value('[HP3]', hp3)
    km.set_value('[OPA]', opa)
    km.set_value('[MA]', ma)
    km.set_value('[H2O2]', h2o2)
    k = TABLE_3_H2O2
    x_h2o2 = (X_3HP_FEED/C0)*h2o2   # H2O2 mole fraction, on R1, R2 and R4
    expected = {'R1': k['k_1']*M_CAT*x_h2o2*hp3,
                'R2': k['k_2']*M_CAT*x_h2o2*opa,
                'R3': k['k_3']*M_CAT*ma,
                'R4': k['k_4']*M_CAT*x_h2o2*ma,
                'R5': k['k_5']*M_CAT*hp3,
                'R6': k['k_6']*M_CAT*ma,
                'R7': k['k_7']*h2o2}        # not multiplied by m_cat
    for reaction, rate in expected.items():
        assert km.get_value(reaction) == pytest.approx(rate, rel=1e-12), reaction
    # H2O2 is consumed 1:1 by R1, R2, R4 and by its decomposition R7
    consumption = -(expected['R1'] + expected['R2'] + expected['R4']
                    + expected['R7'])
    assert km.get_value("[H2O2]'") == pytest.approx(consumption, rel=1e-12)


@pytest.mark.parametrize('load, t_end', [(load_h2o2_model, 1584.),
                                         (load_o2_model, 1848.)])
def test_liquid_species_are_conserved(load, t_end):
    # Every tracked step converts one liquid species into one other, so the
    # six concentrations always add up to the initial 3-HP.
    times = np.linspace(0., t_end, 50)
    sim = simulate_at_times(load(), times, _CONCENTRATIONS)
    total = sim[_CONCENTRATIONS].sum(axis=1)
    assert np.allclose(total, C0, rtol=1e-6, atol=0.)
    assert (sim[_CONCENTRATIONS].to_numpy() > -1e-9).all()


# --- 3. O2 model against the deposit (Fig. 9b, SI Fig. S12b) ----------------

def test_o2_model_matches_the_deposited_fig9b_curve(o2_model, fig9b):
    nodes = fig9b['model']
    nodes = nodes[nodes['time'] > 0.]
    assert len(nodes) == 19
    rms, max_abs = _rms_and_max(simulate_at_times(o2_model, nodes['time']),
                                nodes)
    #                        3-HP    MA      AAc   OXA   D
    assert (rms <= np.array([1.5e-3, 1.5e-3, 5e-5, 5e-5, 3.5e-4])).all(), rms
    assert (max_abs <= np.array([3e-3, 3e-3, 1e-4, 1e-4, 6e-4])).all(), max_abs


def test_o2_model_opa_matches_the_published_inset(o2_model):
    published = pkg.read_figS12b_opa_inset_digitized()
    sim = simulate_at_times(o2_model, published['time'], ['x_OPA'])
    assert np.abs(sim['x_OPA'] - published['x_OPA']).max() <= 5e-7
    # peak of the smooth curve: 1.43e-5 at 20.5 min (4.2 % of the 3-HP feed)
    early = simulate_at_times(o2_model, np.arange(0., 60.01, 0.5), ['x_OPA'])
    peak = early['x_OPA'].idxmax()
    assert 18. <= early['time'][peak] <= 23.
    assert 1.40e-5 <= early['x_OPA'][peak] <= 1.47e-5


def test_o2_model_parity_with_experiment(o2_model, fig9b):
    experiment = fig9b['experiment']
    assert len(experiment) == 16
    sim = simulate_at_times(o2_model, experiment['time'])
    # 0.9945, not the 0.997 printed in SI Fig. S11a (see NOTES.md)
    assert _pooled_r2(experiment, sim) == pytest.approx(0.9945, abs=1e-3)


# --- 4. H2O2 model against the deposit (Fig. 9a, SI Fig. S12a) --------------

def test_h2o2_model_matches_the_deposited_fig9a_curve(h2o2_model, fig9a):
    nodes = fig9a['model']
    assert len(nodes) == 24
    rms, max_abs = _rms_and_max(simulate_at_times(h2o2_model, nodes['time']),
                                nodes)
    assert (rms <= 5e-4).all(), rms
    assert max_abs.max() <= 1.5e-3, max_abs


def test_h2o2_trace(h2o2_model, fig9a):
    trace = fig9a['h2o2_trace']
    sim = simulate_at_times(h2o2_model, trace['time'], ['eq_H2O2'])
    assert sim['eq_H2O2'][0] == pytest.approx(15.)
    # Known deposit inconsistency: the deposited trace comes from another
    # parameter set and lies up to 0.38 equivalents below the Table 3 model,
    # so only a loose bound holds. The t = 0 cell of the deposit is skipped.
    diff = (sim['eq_H2O2'] - trace['eq_H2O2'])[1:]
    assert np.abs(diff).max() <= 0.45
    assert diff.min() > 0.
    late = simulate_at_times(h2o2_model, [1000.], ['eq_H2O2'])
    assert 0. < late['eq_H2O2'][0] < 0.03


def test_h2o2_model_parity_with_experiment(h2o2_model, fig9a):
    experiment = fig9a['experiment']
    sim = simulate_at_times(h2o2_model, experiment['time'])
    # 0.945 is reached on the 11 rows that leave out t = 20 and 1205 min
    keep = ~experiment['time'].isin([20., 1205.])
    assert keep.sum() == 11
    assert _pooled_r2(experiment[keep], sim[keep]) == pytest.approx(0.945,
                                                                    abs=1e-3)


# --- 5. Arrhenius temperature dependence of the O2 model --------------------

def test_o2_model_at_reference_temperature_uses_the_reference_constants():
    km = load_o2_model()
    for i in range(1, 7):
        assert km.get_value('k_%d' % i) == km.get_value('k_%d_ref' % i)
    reference = simulate_at_times(km, [20., 240., 1440.])
    # With no temperature sensitivity left, T has no other way into the model.
    km = load_o2_model(T=348.15)
    for name in ACTIVATION_ENERGIES:
        km.set_value(name, 0.)
    assert simulate_at_times(km, [20., 240., 1440.]).equals(reference)


@pytest.mark.parametrize('T', [298.15, 348.15, 373.15])
def test_o2_rate_constants_follow_arrhenius(T):
    km = load_o2_model(T=T)
    assert km.get_value('T') == T
    for i in range(1, 7):
        factor = np.exp(-ACTIVATION_ENERGIES['Ea_%d' % i]/8.314
                        * (1./T - 1./323.15))
        ratio = km.get_value('k_%d' % i)/km.get_value('k_%d_ref' % i)
        assert ratio == pytest.approx(factor, rel=1e-12), i
    # reset() restores the state but keeps the temperature
    km.reset()
    assert km.get_value('T') == T


def test_o2_model_is_faster_at_higher_temperature():
    half_conversion = []
    for T in (298.15, 323.15, 348.15):
        sim = simulate_at_times(load_o2_model(T=T),
                                np.arange(0., 3000.1, 1.), ['X_HP3'])
        half_conversion.append(sim['time'][(sim['X_HP3'] < 0.5).idxmax()])
    assert half_conversion[0] > half_conversion[1] > half_conversion[2]


# --- 6. reference data ------------------------------------------------------

def test_deposited_files_are_verbatim_copies():
    # SHA-256 of the files in the Dryad deposit (doi:10.5061/dryad.jsxksn0q4)
    for filename, sha256 in (
            ('10-1_kinetic_model_H2O2.csv',
             'eab233b4e2e28fc00743f5637d002a8d87f9c640daa8939779bfe80145b81cdc'),
            ('10-2_kinetic_model_O2.csv',
             '3aff103e744a09925551c10116db6d987a6c8883d0c5114b9e23380c4ed7036e')):
        with open(os.path.join(pkg.reference_data_dir, filename), 'rb') as f:
            assert hashlib.sha256(f.read()).hexdigest() == sha256, filename


def test_readers_undo_the_deposit_quirks(fig9a, fig9b):
    assert {k: len(v) for k, v in fig9a.items()} == {
        'experiment': 13, 'model': 24, 'h2o2_trace': 24}
    assert {k: len(v) for k, v in fig9b.items()} == {
        'experiment': 16, 'model': 20}
    # the swapped 'exp'/'model' labels of 10-1: sampling times vs solver nodes
    assert list(fig9a['experiment']['time']) == [
        0, 20, 40, 60, 90, 120, 180, 240, 300, 420, 1205, 1370, 1440]
    assert fig9a['model']['time'].iloc[1] == 12.547
    assert fig9a['model']['time'].iloc[-1] == 1584.
    # t = 0 rows: feed instead of zeros, 15 equivalents instead of 1
    assert list(fig9a['experiment'].iloc[0]) == [0., 1., 0., 0., 0., 0.]
    assert list(fig9a['h2o2_trace'].iloc[0]) == [0., 15.]
    # every row is a set of fractions of the five measured species
    for table in (fig9a['experiment'], fig9a['model'],
                  fig9b['experiment'], fig9b['model']):
        assert list(table.columns) == ['time'] + _NORMALIZED
        assert np.allclose(table[_NORMALIZED].sum(axis=1), 1., atol=5e-5)


def test_digitized_inset_readers():
    h2o2 = pkg.read_figS12a_h2o2_inset_digitized()
    assert {'time', 'eq_H2O2', 'source_axes'} <= set(h2o2.columns)
    assert set(h2o2['source_axes']) == {'inset', 'main'}
    assert h2o2['eq_H2O2'].between(0., 15.1).all()
    opa = pkg.read_figS12b_opa_inset_digitized()
    assert len(opa) == 13
    assert opa['x_OPA'].max() == pytest.approx(1.42e-5, abs=2e-7)


# --- 7. shipped SBML in sync with the Antimony ------------------------------

_MODEL_FILES = [(pkg.h2o2_antimony_filename, load_h2o2_model),
                (pkg.o2_antimony_filename, load_o2_model)]


def _sbml_path(antimony_filename):
    return os.path.join(
        _PKG_DIR, antimony_filename.replace('_antimony.txt', '_sbml.xml'))


@pytest.mark.parametrize('antimony_filename, load', _MODEL_FILES)
def test_shipped_sbml_reproduces_the_antimony_model(antimony_filename, load):
    from_sbml = nsk.KineticModel.from_sbml(_sbml_path(antimony_filename))
    from_sbml._te.integrator.absolute_tolerance = 1e-10
    from_sbml._te.integrator.relative_tolerance = 1e-9
    times = [0., 20., 240., 1440.]
    selections = _NORMALIZED + _CONCENTRATIONS
    assert np.allclose(simulate_at_times(from_sbml, times, selections),
                       simulate_at_times(load(), times, selections),
                       rtol=1e-9, atol=1e-12)


@pytest.mark.parametrize('antimony_filename, load', _MODEL_FILES)
def test_shipped_sbml_matches_the_live_model(antimony_filename, load):
    libsbml = pytest.importorskip('libsbml')

    def formula(math_container):
        return libsbml.formulaToL3String(math_container.getMath())

    doc = libsbml.readSBML(_sbml_path(antimony_filename))
    assert doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) == 0
    shipped = doc.getModel()
    live_doc = libsbml.readSBMLFromString(load()._te.getSBML())
    live = live_doc.getModel()

    assert ([p.getId() for p in shipped.getListOfParameters()]
            == [p.getId() for p in live.getListOfParameters()])
    for parameter in live.getListOfParameters():
        if parameter.isSetValue():
            pid = parameter.getId()
            assert (shipped.getParameter(pid).getValue()
                    == parameter.getValue()), pid
    assert ([r.getId() for r in shipped.getListOfReactions()]
            == [r.getId() for r in live.getListOfReactions()])
    for l_rx in live.getListOfReactions():
        s_rx = shipped.getReaction(l_rx.getId())
        assert s_rx.getReversible() == l_rx.getReversible()
        assert formula(s_rx.getKineticLaw()) == formula(l_rx.getKineticLaw())
        for side in ('getListOfReactants', 'getListOfProducts'):
            assert ({s.getSpecies(): s.getStoichiometry()
                     for s in getattr(s_rx, side)()}
                    == {s.getSpecies(): s.getStoichiometry()
                        for s in getattr(l_rx, side)()})
    assert ({r.getVariable(): formula(r) for r in shipped.getListOfRules()}
            == {r.getVariable(): formula(r) for r in live.getListOfRules()})
