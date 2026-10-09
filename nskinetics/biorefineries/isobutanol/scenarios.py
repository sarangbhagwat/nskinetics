# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Scenario registry of the isobutanol biorefinery -- the single source of
truth for every scenario's baseline kinetics + uncertainty distributions
(from its parameter-distributions workbook), baseline feeding strategy and
smoke-test pins.

A "scenario" = (workbook: baseline kinetics + distributions)
             + (feeding strategy: max_n_spikes, threshold_conc, target_conc)
             + (expected pins for the smoke tests)
             + (enzyme-burden policy).

:func:`load_scenario` applies a scenario to the ALREADY-BUILT model (call
``nskinetics.biorefineries.isobutanol.load(...)`` first, choosing the
separation configuration); it never calls ``load()`` itself. The
``assert_*`` helpers are the smoke-test gates (see
:mod:`nskinetics.biorefineries.isobutanol.smoke`).

Importing this module is build-free: model/system objects,
``kinetic_optimization`` and ``enzyme_burden`` are imported lazily, inside
the functions that need them (valid only after ``load()``).
"""
import os
import math
from dataclasses import dataclass

from . import paths as _paths

__all__ = ('ScenarioSpec', 'SCENARIOS', 'load_scenario',
           'assert_MPSPs_pinned', 'assert_MPSPs_stable',
           'assert_objective_reproduced', 'assert_spike_feed_residual',
           'SPIKE_FEED_RESIDUAL_TOL')

_WORKBOOK_DIR = _paths.workbook_dir()


@dataclass(frozen=True)
class ScenarioSpec:
    """One registered scenario (an entry of :data:`SCENARIOS`).

    Parameters
    ----------
    name : str
        Scenario key, e.g. ``'A'``, ``'B'``, ``'opt_PI_TRY_informed'``.
    workbook : str
        Bare filename of the parameter-distribution workbook under
        ``data/parameter_distributions/`` (see :attr:`workbook_path`).
    max_n_spikes : int
        Glucose-spike cap of the baseline feeding strategy (0 = batch).
    threshold_conc, target_conc : float
        Baseline feeding strategy [g/L glucose]: spike when glucose falls to
        ``threshold_conc``, back up to ``target_conc``.
    expected : dict
        Smoke-test MPSP pins, ``{'ethanol': float, 'isobutanol': float}``
        (purity-adjusted, $/kg at 15 % IRR; nan for an empty product).
    objective_name : str, optional
        For opt_* scenarios, the ``kinetic_optimization.OBJECTIVE_REGISTRY``
        key of the campaign objective; None otherwise.
    objective_value : float, optional
        The recorded best-trial value of that objective.
    objective_tol : float, optional
        Relative tolerance of :func:`assert_objective_reproduced`. Default 0.02.
    spike_conc, stage_1_max_x : float, optional
        None for every current scenario (the ``load()`` defaults, 600 g/L
        and 5.0 g/L, apply); a non-None value is not yet supported (see
        :func:`load_scenario`).
    burden_default : bool, optional
        Enzyme-burden policy :func:`load_scenario` applies when its own
        ``burden`` argument is None. Default True (False for scenario B).
    """
    name: str
    workbook: str
    max_n_spikes: int
    threshold_conc: float
    target_conc: float
    expected: dict
    objective_name: str = None
    objective_value: float = None
    objective_tol: float = 0.02
    spike_conc: float = None
    stage_1_max_x: float = None
    burden_default: bool = True

    @property
    def workbook_path(self):
        """Absolute path of the scenario's parameter-distribution workbook."""
        return os.path.join(_WORKBOOK_DIR, self.workbook)


#: Every registered scenario by name (see :class:`ScenarioSpec`). The
#: ``expected`` pins gate smoke cases 1-4 and 9-16 (smoke.SMOKE_CASES).
SCENARIOS = {
    'A': ScenarioSpec(
        name='A',
        workbook='parameter-distributions_corn_IBO_EtOH_A.xlsx',
        max_n_spikes=16, threshold_conc=217.125, target_conc=221.25,
        expected={'ethanol': 0.86021, 'isobutanol': math.nan}),
    'B': ScenarioSpec(
        name='B',
        workbook='parameter-distributions_corn_IBO_EtOH_B.xlsx',
        max_n_spikes=0, threshold_conc=34.25, target_conc=140.0,
        expected={'ethanol': 0.59530, 'isobutanol': 1.2800},
        burden_default=False),
    # opt_* scenarios (since 2026-10-08) reproduce the best trial -- by each
    # campaign's OWN objective -- of the eight 2026-09-23/24 scenario-A-
    # anchored ethanol_isobutanol x metabolic_split_12d GP campaigns: the
    # seven seed-350 campaigns (_rs350; PI (log-tail), IBO / EtOH titer,
    # yield, productivity) and the PI (log-tail) relay over the six
    # process-level _rs350 donors (_rl15c111dc). Trials (CSV trial_number):
    # opt_PI_uninformed #103 (_rs350 PI campaign; never left the ethanol
    # basin), opt_PI_TRY_informed #1912 (relay; co-production IBO 39.6 +
    # EtOH 40.6 g/L), opt_IBO_titer #84, opt_IBO_yield #841,
    # opt_IBO_productivity #1472, opt_EtOH_titer #1820, opt_EtOH_yield
    # #1923, opt_EtOH_productivity #365. The two PI scenarios pin PI (the
    # campaigns optimized PI (log-tail), whose argmax is PI's). Workbooks
    # and pins are generated by build_opt_split12d_workbooks.py (analyses)
    # (ask-first): the workbook = scenario A's live kinetics (INTENDED
    # k_7 / k_8, every K_*) + the trial's k_3 / k_6 / k_13 / k_17 and its
    # glycolysis / ehrlich_downstream / inhib_* group members; the burden
    # (A-referenced, burden_default=True) derates k_7 / k_8 at the
    # load_simulate choke point, as during the campaigns. The pins are the
    # LIVE-model reproduction (hensmith 95e9362): every fermentation metric
    # reproduces the recorded trial to < 1e-5, while the TEA metrics
    # (IRR / TCI / PI) differ from the recorded values by ~1-5 % because
    # the campaigns ran under the hensmith 2b5b27d export (PI 0.0493 ->
    # 0.0761 for opt_PI_uninformed, 0.8083 -> 0.8559 for
    # opt_PI_TRY_informed). The near-zero-coproduct MPSPs (ethanol of the
    # IBO scenarios, e.g. 5760.9 for opt_IBO_yield; isobutanol nan for the
    # ethanol-basin scenarios) are purity-adjusted MPSP artifacts, not bugs.
    # Superseded: opt_IRR (trial 1602 of the 2026-09-16 PI (log-tail) GP
    # study) and the 2026-09-07/08 metabolic_minimal_subset optima (k_17 =
    # 44 workbooks); their re-pin history is in the git history of this file.
    'opt_PI_uninformed': ScenarioSpec(
        name='opt_PI_uninformed',
        workbook='parameter-distributions_corn_IBO_EtOH_opt_PI_uninformed.xlsx',
        max_n_spikes=43, threshold_conc=230.76862454834765, target_conc=235.76862454834765,
        expected={'ethanol': 0.8209082015734926, 'isobutanol': math.nan},
        objective_name='PI', objective_value=0.07605700208182188),
    'opt_PI_TRY_informed': ScenarioSpec(
        name='opt_PI_TRY_informed',
        workbook='parameter-distributions_corn_IBO_EtOH_opt_PI_TRY_informed.xlsx',
        max_n_spikes=1, threshold_conc=189.19244688495075, target_conc=194.19244688495075,
        expected={'ethanol': 0.4606577978165082, 'isobutanol': 1.1123449034049862},
        objective_name='PI', objective_value=0.8559120965194471),
    'opt_IBO_titer': ScenarioSpec(
        name='opt_IBO_titer',
        workbook='parameter-distributions_corn_IBO_EtOH_opt_IBO_titer.xlsx',
        max_n_spikes=46, threshold_conc=295.0, target_conc=300.0,
        expected={'ethanol': 207.4614507511063, 'isobutanol': 2.4624841465596883},
        objective_name='IBO titer', objective_value=46.0374309889061),
    'opt_IBO_yield': ScenarioSpec(
        name='opt_IBO_yield',
        workbook='parameter-distributions_corn_IBO_EtOH_opt_IBO_yield.xlsx',
        max_n_spikes=21, threshold_conc=35.0569044943704, target_conc=40.0569044943704,
        expected={'ethanol': 5760.9136657577355, 'isobutanol': 1.8532747760425223},
        objective_name='IBO yield', objective_value=0.3814326759972127),
    'opt_IBO_productivity': ScenarioSpec(
        name='opt_IBO_productivity',
        workbook='parameter-distributions_corn_IBO_EtOH_opt_IBO_productivity.xlsx',
        max_n_spikes=21, threshold_conc=295.0, target_conc=300.0,
        expected={'ethanol': 491.158340624422, 'isobutanol': 2.269370275493917},
        objective_name='IBO productivity', objective_value=2.062857742480627),
    'opt_EtOH_titer': ScenarioSpec(
        name='opt_EtOH_titer',
        workbook='parameter-distributions_corn_IBO_EtOH_opt_EtOH_titer.xlsx',
        max_n_spikes=35, threshold_conc=295.0, target_conc=300.0,
        expected={'ethanol': 0.8567173739935816, 'isobutanol': math.nan},
        objective_name='EtOH titer', objective_value=168.45160120071247),
    'opt_EtOH_yield': ScenarioSpec(
        name='opt_EtOH_yield',
        workbook='parameter-distributions_corn_IBO_EtOH_opt_EtOH_yield.xlsx',
        max_n_spikes=8, threshold_conc=85.18036240246147, target_conc=107.54474841058254,
        expected={'ethanol': 0.868382857165542, 'isobutanol': math.nan},
        objective_name='EtOH yield', objective_value=0.48466606785335037),
    'opt_EtOH_productivity': ScenarioSpec(
        name='opt_EtOH_productivity',
        workbook='parameter-distributions_corn_IBO_EtOH_opt_EtOH_productivity.xlsx',
        max_n_spikes=15, threshold_conc=295.0, target_conc=300.0,
        expected={'ethanol': 0.908130907666657, 'isobutanol': math.nan},
        objective_name='EtOH productivity', objective_value=6.896774756187426),
}


def _handles():
    """Live model/system handles. Valid only after the package load()."""
    from nskinetics.biorefineries import isobutanol
    m = isobutanol.models
    model = m.models_EtOH_IBO_corn.model
    return dict(model=model,
                namespace_dict=m.namespace_dict,
                fbs_spec=m.fbs_spec,
                model_specification=model.specification,
                solve_TEA=isobutanol.system.solve_TEA,
                V406=model.system.flowsheet.V406,
                tea=model.system.TEA)


#: Cached scenario-A reference kinetics (the burden's wild-type calibration
#: point; native pools inert there). Populated once per kernel from the live
#: model after the A workbook is applied. A's kinetics are fixed within a
#: process, so caching is safe.
_A_REFERENCE = [None]


def _load_workbook(model, spec, h):
    """Apply a scenario workbook's Baseline column to the live model
    (kinetics + distributions), matching the existing load_scenario recipe."""
    model.parameters = ()
    model.load_parameter_distributions(spec.workbook_path, h['namespace_dict'])
    model.metrics_at_baseline()     # workbook Baseline column SETS kinetics


def _snapshot_kinetics(h):
    """The live kinetic parameters ({name: value}) off r_te, via the same
    reader the optimizer uses (kinetic_optimization.discover_kinetic_parameters)."""
    from nskinetics.biorefineries.isobutanol import kinetic_optimization as ko
    r_te = h['V406'].nsk_kinetic_model._te
    return ko.discover_kinetic_parameters(r_te)


def _scenario_A_reference(model, h):
    """Scenario A's full live kinetics, cached. On a cache miss, apply the A
    workbook and snapshot (this leaves the model on A; the caller re-applies
    the requested scenario afterwards)."""
    if _A_REFERENCE[0] is None:
        _load_workbook(model, SCENARIOS['A'], h)
        _A_REFERENCE[0] = _snapshot_kinetics(h)
    return _A_REFERENCE[0]


def load_scenario(scenario, *, apply=True, burden=None):
    """Apply a registered scenario to the built model.

    Loads the scenario's baseline kinetics + distributions from its workbook
    (the workbook's Baseline column is applied by one
    ``model.metrics_at_baseline()`` simulation), sets the baseline feeding
    strategy, installs the enzyme burden per policy and, when ``apply`` is
    True, runs one baseline ``model_specification``. Requires the model to
    be built already (``nskinetics.biorefineries.isobutanol.load(...)``).

    Parameters
    ----------
    scenario : str
        Key of :data:`SCENARIOS`.
    apply : bool, optional
        If True (default), finish with one ``model_specification`` at the
        scenario's feeding strategy.
    burden : bool or None, optional
        None (default) applies the scenario's ``burden_default``; True/False
        overrides it. When effectively on, the enzyme ``BurdenModel`` is
        built from SCENARIO A's kinetics (the wild-type calibration where
        native pools are inert) -- NOT from this scenario's own baseline --
        and installed via ``system.set_active_burden``, so the
        ``load_simulate`` choke point derates k_7/k_8 for every simulation
        of this scenario.

    Returns
    -------
    dict
        Bundle with keys ``spec`` (the :class:`ScenarioSpec`),
        ``feeding_kwargs`` (``threshold_conc``/``target_conc``), ``expected``,
        ``model``, ``model_specification``, ``solve_TEA``, ``fbs_spec``,
        ``namespace_dict``, ``V406``, ``tea``, ``burden_model`` and
        ``burden_on``.

    Raises
    ------
    ValueError
        If `scenario` is not registered, or if the burden is on and the
        scenario's own baseline point is over the proteome cap (e.g.
        scenario B's Ehrlich constants) -- pass ``burden=False`` then.
    """
    if scenario not in SCENARIOS:
        raise ValueError(f'scenario {scenario!r} not in SCENARIOS: '
                         f'{sorted(SCENARIOS)}')
    spec = SCENARIOS[scenario]
    burden_on = spec.burden_default if burden is None else bool(burden)
    h = _handles()
    model = h['model']

    # Snapshot the scenario-A burden reference BEFORE loading a non-A
    # workbook (cache miss loads A, then we overwrite with the request).
    A_reference = _scenario_A_reference(model, h) if burden_on else None

    _load_workbook(model, spec, h)

    h['fbs_spec'].max_n_spikes = spec.max_n_spikes
    if spec.spike_conc is not None or spec.stage_1_max_x is not None:
        raise NotImplementedError(
            'load_scenario does not yet override spike_conc / stage_1_max_x; '
            'every registered scenario pins them at the load() defaults '
            '(600 g/L / 5.0 g/L).')
    feeding_kwargs = dict(threshold_conc=spec.threshold_conc,
                          target_conc=spec.target_conc)

    from nskinetics.biorefineries.isobutanol import system as _system
    burden_model = None
    if burden_on:
        from nskinetics.biorefineries.isobutanol import enzyme_burden as eb
        burden_model = eb.BurdenModel.from_reference(A_reference)
        point = _snapshot_kinetics(h)          # this scenario's baseline point
        result = burden_model.evaluate(point)
        if not result.feasible:
            raise ValueError(
                f'scenario {scenario!r} baseline is over the enzyme-burden '
                f'proteome cap (Phi_M {result.Phi_M:.4f} >= F_flex '
                f'{result.F_flex:.4f} g/gDCW) against the scenario-A '
                'reference; it cannot be simulated with the burden on. Pass '
                'burden=False to load it burden-free.')
        _system.set_active_burden(burden_model)
    else:
        _system.set_active_burden(None)

    if apply:
        h['model_specification'](**feeding_kwargs)

    return dict(spec=spec, feeding_kwargs=feeding_kwargs,
                expected=spec.expected, model=model,
                model_specification=h['model_specification'],
                solve_TEA=h['solve_TEA'], fbs_spec=h['fbs_spec'],
                namespace_dict=h['namespace_dict'], V406=h['V406'],
                tea=h['tea'], burden_model=burden_model, burden_on=burden_on)


def assert_MPSPs_pinned(expected, current, sim_number, rel_tol=0.01):
    """Assert the current MPSPs match their pins.

    Parameters
    ----------
    expected : dict
        Pinned MPSPs by stream ID (e.g. ``ScenarioSpec.expected``); a nan
        pin requires a nan MPSP (empty product).
    current : dict
        A ``solve_TEA`` result (its ``'MPSPs'`` entry is read).
    sim_number : int
        Simulation index, for the assertion message.
    rel_tol : float, optional
        Maximum relative deviation from a finite pin. Default 0.01.

    Raises
    ------
    AssertionError
        On the first MPSP outside its pin.
    """
    for ID, ref in expected.items():
        cur = current['MPSPs'][ID]
        if math.isnan(ref):
            assert math.isnan(cur), \
                f'sim {sim_number}: {ID} MPSP {cur} expected nan (empty product)'
        else:
            assert abs(cur - ref)/abs(ref) < rel_tol, \
                (f'sim {sim_number}: {ID} MPSP {cur} not within rel tol '
                 f'{rel_tol} of pinned {ref}')


def assert_MPSPs_stable(reference, current, sim_number, rel_tol=5e-3):
    """Assert the current MPSPs match a reference simulation's.

    Parameters
    ----------
    reference, current : dict
        ``solve_TEA`` results (their ``'MPSPs'`` entries are compared); nan
        must match nan.
    sim_number : int
        Simulation index, for the assertion message.
    rel_tol : float, optional
        Maximum relative drift. Default 5e-3 (~3 significant figures).

    Raises
    ------
    AssertionError
        On the first MPSP that drifted.
    """
    for ID, ref in reference['MPSPs'].items():
        cur = current['MPSPs'][ID]
        if math.isnan(ref) or math.isnan(cur):
            assert math.isnan(ref) and math.isnan(cur), \
                f'sim {sim_number}: {ID} MPSP {cur} vs first-sim {ref} (nan mismatch)'
        else:
            assert abs(cur - ref)/abs(ref) < rel_tol, \
                (f'sim {sim_number}: {ID} MPSP {cur} drifted from first-sim '
                 f'value {ref} beyond rel tol {rel_tol}')


#: Pin on the fermentor's reactor-boundary spike-feed residual, |(implied -
#: delivered)/delivered| of the spiked species (nskinetics
#: NSKBatchReactor.spike_feed_residual): the reconciler's own closure
#: tolerance (spike_feed_reconciliation_tol = 1e-3). nskinetics report
#: docs/reports/fed-batch-spike-feed-reconciliation.md section 6.4.4: a
#: flowsheet-level yield ceiling (ethanol/glucose <= 0.511) lets a 20 %
#: phantom feed through at a 0.40 g/g yield; a 0.1 % residual pin does not.
#: Baseline values: A -1.6e-4 (7 spikes), opt_* within +-4e-5 (all fed
#: batch, 1-15 spikes), the batch B exactly 0.
SPIKE_FEED_RESIDUAL_TOL = 1e-3

def assert_spike_feed_residual(sim_number, abs_tol=SPIKE_FEED_RESIDUAL_TOL):
    """Assert the fermentor V406 closed its spike-feed balance.

    Checks the simulation just run: the reconciler hook is on (the residual
    is not None -- the factory installs ``fbs_spec`` as the reconciler) and
    ``|V406.spike_feed_residual| <= abs_tol``. Reads the built system, so
    valid only after ``load()``.

    Parameters
    ----------
    sim_number : int
        Simulation index, for the assertion message.
    abs_tol : float, optional
        Maximum absolute residual. Default :data:`SPIKE_FEED_RESIDUAL_TOL`.

    Raises
    ------
    AssertionError
        If the hook is off or the residual exceeds `abs_tol`.
    """
    from nskinetics.biorefineries.isobutanol import system as _system
    reactor = _system.f.unit.V406
    residual = getattr(reactor, 'spike_feed_residual', None)
    assert residual is not None, \
        (f'sim {sim_number}: V406.spike_feed_residual is None -- the '
         'spike-feed reconciler hook is off (V406.spike_feed_reconciler = '
         f'{getattr(reactor, "spike_feed_reconciler", None)!r})')
    assert abs(residual) <= abs_tol, \
        (f'sim {sim_number}: V406 spike-feed residual {residual:+.3e} '
         f'exceeds {abs_tol:g} (implied {reactor.spike_feed_implied:.6g} '
         f'vs delivered {reactor.spike_feed_delivered:.6g} kg/hr after '
         f'{reactor.n_spike_feed_reconciliation_passes} reconciliation '
         'pass(es))')


def assert_objective_reproduced(spec, results, bundle, sim_number, rel_tol=None):
    """Assert an opt_* scenario reproduces its campaign's best objective.

    Recomputes the study's objective from the current simulated state (via
    ``kinetic_optimization.OBJECTIVE_REGISTRY``, the exact study definition)
    and compares it with ``spec.objective_value``. A no-op for scenarios
    without an objective (A, B).

    Parameters
    ----------
    spec : ScenarioSpec
        The scenario (``bundle['spec']``).
    results : dict
        The ``solve_TEA`` result of the simulation just run.
    bundle : dict
        The :func:`load_scenario` bundle (``V406`` and ``tea`` are read).
    sim_number : int
        Simulation index, for the assertion message.
    rel_tol : float, optional
        Relative tolerance; None (default) uses ``spec.objective_tol``.

    Raises
    ------
    AssertionError
        If the objective is outside the tolerance.
    """
    if spec.objective_name is None:
        return
    tol = spec.objective_tol if rel_tol is None else rel_tol
    from nskinetics.biorefineries.isobutanol import kinetic_optimization as ko
    getter = ko.OBJECTIVE_REGISTRY[spec.objective_name]['getter']
    h = {'V406': bundle['V406'], 'tea': bundle['tea'],
         'latest_TEA_solution': results}
    value = getter(h)
    ref = spec.objective_value
    assert abs(value - ref)/abs(ref) < tol, \
        (f'sim {sim_number}: {spec.name} objective {spec.objective_name} '
         f'{value} not within rel tol {tol} of recorded best-trial {ref}')
