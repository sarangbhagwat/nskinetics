# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Shared helpers for the isobutanol biorefinery tests.

Imports nothing heavy: test modules may import this at module level without
pulling biosteam, the biorefinery or the upstream ``biorefineries`` package
into the pytest process (the import-guard tests depend on that).
"""
import importlib.util
import os
import subprocess
import sys
import time

import pytest

__all__ = (
    'TESTS_DIR', 'SCRIPTS_DIR', 'FIXTURES_DIR', 'REPO_ROOT', 'NSK_DIR', 'PKG_DIR',
    'load_by_path', 'load_paths', 'missing_modules', 'skip_unless_modules',
    'simulation_deps_missing', 'requires_simulation', 'run_isolated', 'assert_ran_ok',
    'SETTLE_S', 'TIMEOUT_SCALE_ENV', 'timeout_scale',
)

#: This folder.
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
#: Script-style checks run in fresh subprocesses (names never match test_*.py).
SCRIPTS_DIR = os.path.join(TESTS_DIR, 'scripts')
#: Small committed inputs used by the tests.
FIXTURES_DIR = os.path.join(TESTS_DIR, 'fixtures')
#: Repository root (holds setup.py).
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(TESTS_DIR)))
#: The ``nskinetics`` package folder.
NSK_DIR = os.path.join(REPO_ROOT, 'nskinetics')
#: The ``nskinetics.biorefineries.isobutanol`` package folder.
PKG_DIR = os.path.join(NSK_DIR, 'biorefineries', 'isobutanol')

#: Seconds to wait after each simulation subprocess so a finishing process's
#: numba-cache writes settle before the next one imports (never run two
#: simulation processes at once).
SETTLE_S = 5.0

#: Environment variable multiplying every :func:`run_isolated` timeout, e.g.
#: ``2`` after clearing the numba caches (the first process recompiles the
#: kernels) or on a slower machine. Default 1.
TIMEOUT_SCALE_ENV = 'NSKINETICS_TEST_TIMEOUT_SCALE'

_NUMBA_CACHE_SIGNATURES = (
    'underlying object has vanished',
    'access violation reading 0xFFFFFFFFFFFFFFFF',
)

# Env vars the root conftest exports that must NOT reach simulation children:
# thermosteam treats any non-empty NUMBA_DISABLE_JIT (even '0') as disabled,
# and the pinned results were produced with JIT on.
_SIM_ENV_STRIP = ('NUMBA_DISABLE_JIT', 'DISABLE_PREFERENCES', 'FILTER_WARNINGS')


def load_by_path(name, path, register=False):
    """Import a module from a file path without importing its parent package."""
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    if register:
        sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_paths():
    """Return the package's stdlib-only ``paths`` module, loaded by path."""
    return load_by_path('_ibo_paths', os.path.join(PKG_DIR, 'paths.py'))


def _has_top_level(name):
    # find_spec on a TOP-LEVEL name never imports anything.
    return importlib.util.find_spec(name) is not None


def missing_modules(*names):
    """Return the subset of top-level module ``names`` that are not importable.

    ``'biorefineries'`` additionally requires the upstream ``corn``,
    ``cellulosic`` and ``ethanol`` subpackages to exist on disk (checked
    without importing them).
    """
    missing = []
    for name in names:
        if name == 'biorefineries':
            spec = importlib.util.find_spec('biorefineries')
            locs = list(spec.submodule_search_locations or []) if spec else []
            ok = bool(locs) and all(
                any(os.path.isdir(os.path.join(loc, sub)) for loc in locs)
                for sub in ('corn', 'cellulosic', 'ethanol'))
            if not ok:
                missing.append(name)
        elif not _has_top_level(name):
            missing.append(name)
    return missing


def skip_unless_modules(*names):
    """``pytest.mark.skipif`` marker that skips when any module is missing."""
    missing = missing_modules(*names)
    return pytest.mark.skipif(bool(missing),
                              reason=f'optional dependencies missing: {missing}')


#: Optional dependencies a full biorefinery simulation needs.
SIMULATION_DEPS = ('biosteam', 'thermosteam', 'biorefineries', 'hensmith',
                   'contourplots', 'openpyxl')


def simulation_deps_missing():
    """Return the missing simulation dependencies (empty list = can simulate)."""
    return missing_modules(*SIMULATION_DEPS)


#: Marker for tests that build/simulate the biorefinery in a subprocess.
requires_simulation = pytest.mark.skipif(
    bool(simulation_deps_missing()),
    reason=f'biorefinery simulation dependencies missing: {simulation_deps_missing()}')


def timeout_scale():
    """Return the timeout multiplier from :data:`TIMEOUT_SCALE_ENV` (default 1)."""
    raw = os.environ.get(TIMEOUT_SCALE_ENV, '').strip()
    try:
        scale = float(raw) if raw else 1.0
    except ValueError:
        raise ValueError(f'{TIMEOUT_SCALE_ENV}={raw!r} is not a number') from None
    if not scale > 0:
        raise ValueError(f'{TIMEOUT_SCALE_ENV}={raw!r} must be positive')
    return scale


def _text(data):
    if data is None:
        return ''
    if isinstance(data, bytes):
        return data.decode('utf-8', errors='replace')
    return data


def run_isolated(args, timeout, simulation=True, cwd=None, extra_env=None):
    """Run ``python <args>`` in a fresh interpreter and return the CompletedProcess.

    Parameters
    ----------
    args : sequence of str
        Arguments after the interpreter, e.g. ``['-m', 'nskinetics...smoke', '1']``
        or ``[script_path]``.
    timeout : float
        Seconds before the child is killed and the test fails (with the
        child's output so far). Multiplied by :func:`timeout_scale`.
    simulation : bool, optional
        If True (default), strip the conftest's JIT/preference env vars so the
        child simulates exactly like the canonical runner, and sleep
        :data:`SETTLE_S` afterwards. Simulation children never run under
        pytest-xdist (the test is skipped): two simulating processes at once
        corrupt the numba on-disk caches.
    cwd : str, optional
        Working directory (default: the repository root).
    extra_env : dict, optional
        Extra environment variables for the child.
    """
    if simulation and os.environ.get('PYTEST_XDIST_WORKER'):
        pytest.skip('isobutanol simulation subprocesses must run serially; '
                    'run without pytest-xdist (-n)')
    timeout = timeout*timeout_scale()
    env = dict(os.environ)
    if simulation:
        for key in _SIM_ENV_STRIP:
            env.pop(key, None)
    env['MPLBACKEND'] = 'Agg'
    env['PYTHONIOENCODING'] = 'utf-8'
    env['PYTHONUNBUFFERED'] = '1'
    env['PYTHONPATH'] = os.pathsep.join(
        [REPO_ROOT] + [p for p in env.get('PYTHONPATH', '').split(os.pathsep) if p])
    if extra_env:
        env.update(extra_env)
    try:
        return subprocess.run([sys.executable, *args], cwd=cwd or REPO_ROOT, env=env,
                              capture_output=True, text=True, encoding='utf-8',
                              errors='replace', timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        out = _text(exc.stdout) + _text(exc.stderr)
        last = '\n'.join(out.splitlines()[-60:])
        pytest.fail(f'child {list(args)} killed after {timeout:.0f} s (cold numba '
                    f'cache or slow machine? set {TIMEOUT_SCALE_ENV}=2 to allow '
                    f'more time)\n{last}')
    finally:
        if simulation:
            time.sleep(SETTLE_S)


def assert_ran_ok(proc, sentinel=None, tail=60):
    """Fail with the child's output tail unless it exited 0 (and printed ``sentinel``)."""
    out = (proc.stdout or '') + (proc.stderr or '')
    last = '\n'.join(out.splitlines()[-tail:])
    if proc.returncode != 0:
        if any(sig in out for sig in _NUMBA_CACHE_SIGNATURES):
            pytest.fail('numba on-disk cache corruption (NOT a model regression): clear '
                        'the *.nbc/*.nbi/*.pyc caches (see CLAUDE.md) and re-run.\n' + last)
        pytest.fail(f'exit code {proc.returncode}\n{last}')
    if sentinel is not None and sentinel not in out:
        pytest.fail(f'sentinel {sentinel!r} not printed\n{last}')
    return out
