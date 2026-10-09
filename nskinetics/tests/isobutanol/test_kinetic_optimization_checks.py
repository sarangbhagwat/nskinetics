# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Offline checks of ``kinetic_optimization`` (``scripts/ko_offline_checks.py``).

The script (~8000 lines, ~100 sequential checks: search spaces, the objective
registry, trajectory CSVs, the optuna engine on fake handles, the BO driver /
supervisor source and CLI, output-path budgets) imports
``nskinetics.biorefineries.isobutanol``, so it runs in a fresh subprocess --
never in the pytest process, whose import-guard tests must not see the
biorefinery. It never simulates the flowsheet, so the child keeps the
conftest's ``NUMBA_DISABLE_JIT=1``: the same 103 passing checks, a faster run
(108 s vs 126 s with JIT on, measured 2026-10-09) and no numba on-disk cache
writes, so it is not a simulation-class process.

With optuna, torch, openpyxl and the scenario A/B workbooks present the run
must report exactly ``ALL <EXPECTED_N_CHECKS> CHECKS PASSED`` with no
``SKIP`` line; without them the script skips the dependent parts (``SKIP n``
lines) and only a clean exit is required.
"""
import os
import re
import shutil
import tempfile

import pytest

from ._support import (SCRIPTS_DIR, assert_ran_ok, load_paths, missing_modules,
                       run_isolated, skip_unless_modules)

#: The script-style check.
KO_SCRIPT = os.path.join(SCRIPTS_DIR, 'ko_offline_checks.py')
#: Checks the script runs when every optional dependency is present.
EXPECTED_N_CHECKS = 103
#: Seconds before the child is killed (~110 s warm).
KO_TIMEOUT_S = 330
#: Optional dependencies whose absence turns checks into SKIP lines.
FULL_RUN_MODULES = ('optuna', 'torch', 'openpyxl')

pytestmark = [pytest.mark.slow, pytest.mark.isobutanol,
              # importing the package needs the biorefinery stack
              skip_unless_modules('biosteam', 'thermosteam', 'biorefineries', 'hensmith')]


def _full_run_blockers():
    """Why the script would SKIP checks here (empty list = it must not)."""
    blockers = missing_modules(*FULL_RUN_MODULES)
    paths = load_paths()
    blockers += [os.path.basename(p) for p in map(paths.workbook_path, ('A', 'B'))
                 if not os.path.isfile(p)]
    return blockers


def test_kinetic_optimization_offline_checks():
    # The script leaves ~60 mkdtemp folders behind: give it a private TEMP
    # and delete it afterwards.
    tmp = tempfile.mkdtemp(prefix='nsk_ko_')
    try:
        proc = run_isolated([KO_SCRIPT], timeout=KO_TIMEOUT_S, simulation=False,
                            extra_env={'TMP': tmp, 'TEMP': tmp, 'TMPDIR': tmp})
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    out = assert_ran_ok(proc)
    lines = [line.strip() for line in proc.stdout.splitlines()]
    tail = '\n'.join(out.splitlines()[-40:])
    summary = [line for line in lines if re.fullmatch(r'ALL \d+ CHECKS PASSED', line)]
    assert summary, 'no "ALL <n> CHECKS PASSED" line\n' + tail
    skips = [line for line in lines if line.startswith('SKIP')]
    blockers = _full_run_blockers()
    if blockers:
        return  # partial run (SKIP lines expected); a clean exit is the check
    assert not skips, (f'{len(skips)} checks skipped although {FULL_RUN_MODULES} and the '
                       'A/B workbooks are available:\n' + '\n'.join(skips))
    assert summary == [f'ALL {EXPECTED_N_CHECKS} CHECKS PASSED'], (
        f'expected ALL {EXPECTED_N_CHECKS} CHECKS PASSED, got {summary}\n' + tail)
