# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Smoke tests of the isobutanol biorefinery: 16 pinned baseline cases.

Each case builds the biorefinery in one separation configuration
(``load(separation_processes=...)``), optionally re-gates the separation
splitter S201, applies a registered scenario
(:func:`scenarios.load_scenario`) and then simulates it ``n_sims`` times
(``model_specification`` + ``solve_TEA``), checking EVERY simulation:

* the purity-adjusted MPSPs are within 1 % of their pins
  (:func:`scenarios.assert_MPSPs_pinned`; a nan pin requires nan);
* V406 closed its spike-feed balance
  (:func:`scenarios.assert_spike_feed_residual`);
* opt_* scenarios reproduce their campaign's best objective within 2 %
  (:func:`scenarios.assert_objective_reproduced`; a no-op for A and B);
* every MPSP stays within 0.5 % of the first simulation's
  (:func:`scenarios.assert_MPSPs_stable`).

A violation raises ``AssertionError``. Each case must run in a FRESH
process: the build happens once per process (a repeat ``load()`` is a
no-op) and scenarios mutate global fermentation state. From a shell::

    python -m nskinetics.biorefineries.isobutanol.smoke 1
    python -m nskinetics.biorefineries.isobutanol.smoke 9 --json --n-sims 3

prints one line per simulation, ``sim1: IRR=<repr> MPSPs={'ethanol': <repr>,
'isobutanol': <repr>}``, then ``elapsed: <s> s`` (and, with ``--json``, one
final JSON line with the case and every result), and exits non-zero on any
failed check. Run the cases strictly one at a time (the numba on-disk caches
of thermosteam/biosteam do not tolerate concurrent processes). Importing
this module is build-free.

Cases (:data:`SMOKE_CASES`)::

    #   scenario               separation_processes   S201        pins
    1   A                      both (default)         1.0         registry
    2   B                      both (default)         1.0         registry
    3   A                      ('IBO_EtOH',)          --          registry
    4   B                      ('IBO_EtOH',)          --          registry
    5   A                      ('ethanol',)           --          0.83971 / nan
    6   B                      ('ethanol',)           --          1.6271 / nan
    7   A                      both                   re-gated 0  0.83971 / nan
    8   B                      both                   re-gated 0  1.6271 / nan
    9   opt_PI_uninformed      ('IBO_EtOH',)          --          registry + objective
    10  opt_PI_TRY_informed    ('IBO_EtOH',)          --          registry + objective
    11  opt_IBO_titer          ('IBO_EtOH',)          --          registry + objective
    12  opt_IBO_yield          ('IBO_EtOH',)          --          registry + objective
    13  opt_IBO_productivity   ('IBO_EtOH',)          --          registry + objective
    14  opt_EtOH_titer         ('IBO_EtOH',)          --          registry + objective
    15  opt_EtOH_yield         ('IBO_EtOH',)          --          registry + objective
    16  opt_EtOH_productivity  ('IBO_EtOH',)          --          registry + objective

"registry" pins are ``scenarios.SCENARIOS[scenario].expected``. Cases 3/4
check that dropping the (pass-through) gating splitter reproduces the
both-trains baselines 1/2. Cases 5-8 carry their own pins
(:data:`ETHANOL_PRIMARY_PINS`): the ethanol-primary train alone (5/6), and
both trains with all broth re-gated to it (7/8: an idle zero-flow IBO/EtOH
branch must be economically equivalent to an absent one, so 7/8 track 5/6
to ~6 significant figures). These replace the former ``analyses/
smoke_test_1.py`` ... ``smoke_test_16.py`` (same configurations and pins;
their re-pin histories are in the Bioindustrial-Park git history of those
files).
"""
import math
from dataclasses import dataclass

__all__ = ('SmokeCase', 'SMOKE_CASES', 'ETHANOL_PRIMARY_PINS',
           'get_case', 'load_simulate_baseline', 'main')

#: Separation configurations of :func:`load`.
BOTH_TRAINS = ('IBO_EtOH', 'ethanol')
IBO_EtOH_ONLY = ('IBO_EtOH',)
ETHANOL_ONLY = ('ethanol',)

#: Pins of the ethanol-primary configurations (cases 5-8; purity-adjusted
#: $/kg at 15 % IRR): the ethanol-only build and the both-trains build
#: re-gated to S201.split = 0.0. Not in the scenario registry, whose A/B pins
#: are those of the IBO/EtOH-train builds.
ETHANOL_PRIMARY_PINS = {
    'A': {'ethanol': 0.83971, 'isobutanol': math.nan},
    'B': {'ethanol': 1.6271, 'isobutanol': math.nan},
}


@dataclass(frozen=True)
class SmokeCase:
    """One smoke-test configuration (an entry of :data:`SMOKE_CASES`).

    Parameters
    ----------
    number : int
        Case number (1-16).
    scenario : str
        Key of ``scenarios.SCENARIOS``.
    separation_processes : tuple[str, ...]
        Passed to ``load(separation_processes=...)``.
    S201_split : float, optional
        If not None, ``system.sep_udct['S201'].split`` is set to this after
        ``load()`` and before the scenario is applied (both-trains builds
        only). Default None (the build's gate, 1.0, or no S201 at all).
    expected : dict, optional
        MPSP pins ``{'ethanol': float, 'isobutanol': float}``; None
        (default) uses the scenario registry's ``expected``.
    description : str, optional
        One-line summary.
    """
    number: int
    scenario: str
    separation_processes: tuple
    S201_split: float = None
    expected: dict = None
    description: str = ''

    @property
    def pins_source(self):
        """'registry' (scenario pins) or 'case' (this case's own pins)."""
        return 'registry' if self.expected is None else 'case'

    def pins(self):
        """Return the MPSP pins this case asserts against.

        Returns
        -------
        dict
            ``{'ethanol': float, 'isobutanol': float}`` (nan = empty product).
        """
        if self.expected is not None:
            return self.expected
        from . import scenarios
        return scenarios.SCENARIOS[self.scenario].expected


def _opt_case(number, scenario):
    return SmokeCase(number, scenario, IBO_EtOH_ONLY,
                     description=f'{scenario} (opt_* campaign optimum; '
                                 'IBO/EtOH train only; objective reproduced)')


#: The 16 smoke cases by number (see the module docstring table).
SMOKE_CASES = {
    1: SmokeCase(1, 'A', BOTH_TRAINS,
                 description='scenario A, both trains (default build)'),
    2: SmokeCase(2, 'B', BOTH_TRAINS,
                 description='scenario B, both trains (default build)'),
    3: SmokeCase(3, 'A', IBO_EtOH_ONLY,
                 description='scenario A, IBO/EtOH train only (reproduces 1)'),
    4: SmokeCase(4, 'B', IBO_EtOH_ONLY,
                 description='scenario B, IBO/EtOH train only (reproduces 2)'),
    5: SmokeCase(5, 'A', ETHANOL_ONLY, expected=ETHANOL_PRIMARY_PINS['A'],
                 description='scenario A, ethanol-primary train only'),
    6: SmokeCase(6, 'B', ETHANOL_ONLY, expected=ETHANOL_PRIMARY_PINS['B'],
                 description='scenario B, ethanol-primary train only'),
    7: SmokeCase(7, 'A', BOTH_TRAINS, S201_split=0.0,
                 expected=ETHANOL_PRIMARY_PINS['A'],
                 description='scenario A, both trains, all broth re-gated to '
                             'the ethanol-primary train (tracks 5)'),
    8: SmokeCase(8, 'B', BOTH_TRAINS, S201_split=0.0,
                 expected=ETHANOL_PRIMARY_PINS['B'],
                 description='scenario B, both trains, all broth re-gated to '
                             'the ethanol-primary train (tracks 6)'),
    9: _opt_case(9, 'opt_PI_uninformed'),
    10: _opt_case(10, 'opt_PI_TRY_informed'),
    11: _opt_case(11, 'opt_IBO_titer'),
    12: _opt_case(12, 'opt_IBO_yield'),
    13: _opt_case(13, 'opt_IBO_productivity'),
    14: _opt_case(14, 'opt_EtOH_titer'),
    15: _opt_case(15, 'opt_EtOH_yield'),
    16: _opt_case(16, 'opt_EtOH_productivity'),
}


def get_case(case):
    """Return the :class:`SmokeCase` for a case number (or a SmokeCase).

    Parameters
    ----------
    case : int, str or SmokeCase
        Case number 1-16 (an int or a digit string), or a case itself.

    Returns
    -------
    SmokeCase

    Raises
    ------
    KeyError
        If the number is not a registered case.
    """
    if isinstance(case, SmokeCase):
        return case
    number = int(case)
    if number not in SMOKE_CASES:
        raise KeyError(f'smoke case {case!r} not in SMOKE_CASES: '
                       f'{sorted(SMOKE_CASES)}')
    return SMOKE_CASES[number]


def load_simulate_baseline(case, plot=False,
                           stream_IDs=('ethanol', 'isobutanol'),
                           IRR_for_MPSP=0.15,
                           n_sims=3):
    """Build, apply and simulate one smoke case, asserting every simulation.

    Calls the package ``load()`` itself (once per process), so run each case
    in a fresh process.

    Parameters
    ----------
    case : int, str or SmokeCase
        The case (see :data:`SMOKE_CASES`).
    plot : bool, optional
        If True, plot the last simulation's kinetics
        (``models.plot_kinetic_results(xlim=(0, 80), ylim=(0, 250))``).
    stream_IDs : tuple[str, ...], optional
        Products whose MPSPs are solved. Default ('ethanol', 'isobutanol').
    IRR_for_MPSP : float, optional
        Fixed IRR at which the MPSPs are solved. Default 0.15.
    n_sims : int, optional
        Simulations in this process (the in-process stability check).
        Default 3.

    Returns
    -------
    list[dict]
        One ``solve_TEA`` result per simulation:
        ``{'IRR': float, 'MPSPs': {stream_ID: float, ...}}``.

    Raises
    ------
    AssertionError
        On the first failed check (see the module docstring).
    RuntimeError
        If this process already built a different separation configuration.
    """
    case = get_case(case)
    from . import load, scenarios, system
    if system._loaded:
        active = system._published.get('separation_processes')
        if tuple(case.separation_processes) != active:
            raise RuntimeError(
                f'smoke case {case.number} needs separation_processes='
                f'{case.separation_processes!r}, but this process already '
                f'built {active!r}; run each smoke case in a fresh process.')
    load(separation_processes=case.separation_processes)
    if case.S201_split is not None:
        # Re-gate the both-trains build before any scenario simulation.
        system.sep_udct['S201'].split = case.S201_split

    bundle = scenarios.load_scenario(case.scenario)
    model_specification = bundle['model_specification']
    solve_TEA = bundle['solve_TEA']
    feeding_kwargs = bundle['feeding_kwargs']
    expected_MPSPs = (bundle['expected'] if case.expected is None
                      else case.expected)
    spec = bundle['spec']

    all_results = []
    for i in range(n_sims):
        model_specification(**feeding_kwargs)
        results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
        scenarios.assert_MPSPs_pinned(expected_MPSPs, results, i+1)
        scenarios.assert_spike_feed_residual(i+1)
        scenarios.assert_objective_reproduced(spec, results, bundle, i+1)
        if all_results:
            scenarios.assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)

    if plot:
        from . import models
        models.plot_kinetic_results(xlim=(0, 80), ylim=(0, 250))

    return all_results


#: Top-level name of the pre-migration copy of this package (Bioindustrial-
#: Park), which may still be importable. Joined at run time so that source
#: scans for stale imports of it do not match this guard.
_PRE_MIGRATION_PACKAGE = '.'.join(('biorefineries', 'isobutanol'))


def _stale_copy_modules():
    """Names of any loaded module of the pre-migration package copy."""
    import sys
    return sorted(m for m in sys.modules
                  if m == _PRE_MIGRATION_PACKAGE
                  or m.startswith(_PRE_MIGRATION_PACKAGE + '.'))


def main(argv=None):
    """Command-line entry point (``python -m ...isobutanol.smoke N``).

    Parameters
    ----------
    argv : list[str], optional
        Arguments (default ``sys.argv[1:]``): the case number, plus
        ``--n-sims K`` (default 3) and ``--json``.

    Returns
    -------
    int
        0 on success; a failed check raises instead.
    """
    import argparse
    import json
    import time
    parser = argparse.ArgumentParser(
        prog='python -m nskinetics.biorefineries.isobutanol.smoke',
        description='Run one isobutanol-biorefinery smoke case (build + '
                    'scenario + pinned simulations) in this process.')
    parser.add_argument('case', type=int, choices=sorted(SMOKE_CASES),
                        help='smoke case number')
    parser.add_argument('--n-sims', type=int, default=3,
                        help='simulations to run and check (default 3)')
    parser.add_argument('--json', action='store_true',
                        help='also print one final JSON line with the results')
    args = parser.parse_args(argv)
    case = get_case(args.case)

    t0 = time.time()
    results = load_simulate_baseline(case, n_sims=args.n_sims)
    elapsed = time.time() - t0

    stale = _stale_copy_modules()
    if stale:
        raise RuntimeError(
            'the pre-migration copy of this package (Bioindustrial-Park) was '
            f'imported in this process ({stale}); results are not trustworthy.')

    rows = []
    for i, r in enumerate(results):
        m = {k: float(v) for k, v in r['MPSPs'].items()}
        IRR = float(r['IRR'])
        print(f"sim{i+1}: IRR={IRR!r} MPSPs={m!r}")
        rows.append({'IRR': IRR, 'MPSPs': m})
    print(f"elapsed: {elapsed:.1f} s")
    if args.json:
        print(json.dumps({'case': case.number,
                          'scenario': case.scenario,
                          'separation_processes': list(case.separation_processes),
                          'S201_split': case.S201_split,
                          'pins_source': case.pins_source,
                          'n_sims': args.n_sims,
                          'elapsed_s': round(elapsed, 1),
                          'results': rows}))
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
