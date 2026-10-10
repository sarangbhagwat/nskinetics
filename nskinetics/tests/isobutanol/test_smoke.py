# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""The 16 isobutanol smoke cases, each in a fresh subprocess.

Each case runs ``python -m nskinetics.biorefineries.isobutanol.smoke N
--json``: one build, the case's scenario, then 3 simulate + ``solve_TEA``
cycles. The smoke module itself asserts every simulation against its pins
(1 %), the in-process stability (0.5 % of the first simulation), the
spike-feed residual and, for opt_* cases, the reproduced objective -- so a
zero exit certifies the case. These tests add the process-level checks:
exit 0, one ``simN:`` line and one result per simulation, a parseable JSON
payload for the right case.

The ``canonical`` gate is cases 1 and 2 run TWICE (fresh processes) with the
two runs' sim1 MPSPs agreeing to 3 significant figures (in practice they are
identical)::

    pytest nskinetics/tests/isobutanol/test_smoke.py -m canonical -v

Run strictly sequentially (no pytest-xdist): two simulating processes at once
corrupt the numba on-disk caches.
"""
import json
import math
import re

import pytest

from ._support import assert_ran_ok, requires_simulation, run_isolated

#: Module the smoke cases are run from (``python -m``).
SMOKE_MODULE = 'nskinetics.biorefineries.isobutanol.smoke'
#: Smoke case numbers (``smoke.SMOKE_CASES``).
CASES = tuple(range(1, 17))
#: Cases of the canonical gate (scenario A and B, default both-trains build).
CANONICAL_CASES = (1, 2)
#: Simulations per process (the smoke CLI default).
N_SIMS = 3
#: Seconds before a smoke process is killed. The slowest case (13) takes
#: ~160 s warm; a cold numba cache adds minutes.
SMOKE_TIMEOUT_S = 1500
#: Relative tolerance of the cross-run comparison: "agree to 3 significant
#: figures" (half a unit in the 4th significant figure, at most).
CANONICAL_REL_TOL = 5e-4

pytestmark = [pytest.mark.slow, pytest.mark.isobutanol, pytest.mark.simulation,
              pytest.mark.smoke, requires_simulation]

#: case -> payloads of the runs made in this pytest session, so the
#: canonical repeat reuses the first run of cases 1 and 2.
_RUNS = {}


def _payload(stdout):
    """The smoke CLI's final JSON line (it may contain NaN tokens)."""
    for line in reversed(stdout.splitlines()):
        line = line.strip()
        if line.startswith('{'):
            return json.loads(line)
    pytest.fail('the smoke CLI printed no JSON line\n' + stdout[-4000:])


def run_smoke_case(case):
    """Run one smoke case in a fresh process; return its JSON payload.

    Parameters
    ----------
    case : int
        Smoke case number.

    Returns
    -------
    dict
        ``{'case', 'scenario', 'separation_processes', 'S201_split',
        'pins_source', 'n_sims', 'elapsed_s', 'results': [{'IRR', 'MPSPs'}]}``.
    """
    proc = run_isolated(['-m', SMOKE_MODULE, str(case), '--json'],
                        timeout=SMOKE_TIMEOUT_S)
    assert_ran_ok(proc)
    payload = _payload(proc.stdout)
    assert payload['case'] == case, payload
    assert payload['n_sims'] == N_SIMS, payload
    results = payload['results']
    assert len(results) == N_SIMS, payload
    for i, r in enumerate(results, 1):
        assert re.search(rf'^sim{i}: IRR=', proc.stdout, re.M), (
            f'no sim{i} line printed\n{proc.stdout[-4000:]}')
        assert set(r['MPSPs']) == {'ethanol', 'isobutanol'}, r
        assert math.isfinite(r['MPSPs']['ethanol']), r
    runs = _RUNS.setdefault(case, [])
    runs.append(payload)
    # The child's result lines, verbatim (shown with pytest -s or -rP).
    print(f'--- smoke case {case}, run {len(runs)} ---')
    for line in proc.stdout.splitlines():
        if re.match(r'(sim\d+|elapsed): ', line):
            print(line)
    return payload


def _same_to_3_sig_figs(a, b):
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return math.isclose(a, b, rel_tol=CANONICAL_REL_TOL, abs_tol=0.0)


@pytest.mark.parametrize('case', [
    pytest.param(n, id=f'case{n}',
                 marks=[pytest.mark.canonical] if n in CANONICAL_CASES else [])
    for n in CASES])
def test_smoke_case(case):
    run_smoke_case(case)


@pytest.mark.canonical
@pytest.mark.parametrize('case', CANONICAL_CASES, ids=lambda n: f'case{n}')
def test_canonical_rerun_reproduces(case):
    """A second fresh run of case 1 / 2 reproduces the first run's sim1 MPSPs.

    Reuses the first run from :func:`test_smoke_case` when it ran in this
    session; otherwise runs the case twice.
    """
    first = (_RUNS.get(case) or [None])[0] or run_smoke_case(case)
    second = run_smoke_case(case)
    mpsp_1 = first['results'][0]['MPSPs']
    mpsp_2 = second['results'][0]['MPSPs']
    bad = {k: (mpsp_1[k], mpsp_2[k]) for k in mpsp_1
           if not _same_to_3_sig_figs(mpsp_1[k], mpsp_2[k])}
    assert not bad, (f'case {case}: sim1 MPSPs differ between two fresh runs beyond '
                     f'3 significant figures (run 1, run 2): {bad}')
