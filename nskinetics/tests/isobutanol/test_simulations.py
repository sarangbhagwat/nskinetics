# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Simulation regressions of the isobutanol biorefinery, each in a fresh subprocess.

Every check is a script in ``scripts/`` (named ``sim_*.py`` so pytest never
collects it): the eight flowsheet regressions build the biorefinery
(``isobutanol.load()``) and end with ``ALL <n> CHECKS PASSED``; the two
separation-train checks load ``separations.py`` by path, simulate the
standalone trains and end with a ``===== SUMMARY =====`` block of PASS/FAIL
lines. A script is run directly with e.g.
``python nskinetics/tests/isobutanol/scripts/sim_unit_groups.py``.

Run strictly sequentially (no pytest-xdist)::

    pytest nskinetics/tests/isobutanol/test_simulations.py -m "slow or not slow" -v
"""
import os
import re

import pytest

from ._support import SCRIPTS_DIR, assert_ran_ok, requires_simulation, run_isolated

pytestmark = [pytest.mark.slow, pytest.mark.isobutanol, pytest.mark.simulation,
              requires_simulation]

#: script -> (exact success line, timeout [s]). Timeouts are ~3x the warm
#: runtimes measured 2026-10-09 (child process, numba caches warm).
FLOWSHEET_CHECKS = {
    'sim_ddgs_dryer_acid_split.py':    ('ALL 4 CHECKS PASSED', 180),   # 62 s
    'sim_fed_batch_volume_guard.py':   ('ALL 7 CHECKS PASSED', 105),   # 34 s
    'sim_kin_opt_effluent_metrics.py': ('ALL 3 CHECKS PASSED', 210),   # 68 s
    'sim_reset_after_fail.py':         ('ALL 5 CHECKS PASSED', 80),    # 26 s
    'sim_scenario_burden.py':          ('ALL 4 CHECKS PASSED', 75),    # 24 s
    'sim_scenario_burden_sweep.py':    ('ALL 3 CHECKS PASSED', 75),    # 25 s
    'sim_unit_groups.py':              ('ALL 11 CHECKS PASSED', 160),  # 53 s
    'sim_v406_burden_guard.py':        ('ALL 6 CHECKS PASSED', 70),    # 23 s
}

#: script -> (number of variants in its SUMMARY block, timeout [s]).
SEPARATION_CHECKS = {
    'sim_separation_system.py':                 (20, 65),  # 21 s
    'sim_separation_system_ethanol_primary.py': (14, 60),  # 19 s
}


def _run_script(script, timeout):
    return run_isolated([os.path.join(SCRIPTS_DIR, script)], timeout=timeout)


@pytest.mark.parametrize('script', list(FLOWSHEET_CHECKS),
                         ids=[s[len('sim_'):-len('.py')] for s in FLOWSHEET_CHECKS])
def test_flowsheet_regression(script):
    sentinel, timeout = FLOWSHEET_CHECKS[script]
    proc = _run_script(script, timeout)
    out = assert_ran_ok(proc)
    lines = [line.strip() for line in proc.stdout.splitlines()]
    assert sentinel in lines, (f'{script}: success line {sentinel!r} not printed\n'
                               + '\n'.join(out.splitlines()[-40:]))


@pytest.mark.parametrize('script', list(SEPARATION_CHECKS),
                         ids=[s[len('sim_'):-len('.py')] for s in SEPARATION_CHECKS])
def test_separation_train(script):
    n_variants, timeout = SEPARATION_CHECKS[script]
    proc = _run_script(script, timeout)
    out = assert_ran_ok(proc)
    stdout = proc.stdout
    assert '===== SUMMARY =====' in stdout, out[-4000:]
    summary = stdout.split('===== SUMMARY =====', 1)[1].strip().splitlines()
    verdicts = [re.search(r'\s(PASS|FAIL)\s*$', line) for line in summary]
    assert all(verdicts), f'{script}: unparsed SUMMARY lines\n' + '\n'.join(summary)
    failed = [line.strip() for line, v in zip(summary, verdicts) if v.group(1) == 'FAIL']
    assert not failed, f'{script}: failed variants {failed}\n' + out[-4000:]
    assert len(summary) == n_variants, (
        f'{script}: {len(summary)} SUMMARY lines, expected {n_variants}\n'
        + '\n'.join(summary))
