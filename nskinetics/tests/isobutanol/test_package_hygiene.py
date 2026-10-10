# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Static hygiene checks of the migrated isobutanol biorefinery
(``nskinetics.biorefineries.isobutanol``).

* No stale reference to the pre-migration package (the Bioindustrial-Park
  copy is still on the developer's PYTHONPATH, so a missed rewrite would
  silently import the OLD code).
* Nothing under ``nskinetics/biorefineries`` is collectable by pytest
  (``pytest .`` would otherwise import -- and build -- the biorefinery).
* Every scenario workbook the registry names ships with the package.
* No Excel lock file (``~$*``) is tracked.
* ``paths.py`` is stdlib-only, and the modules loaded BY FILE PATH (the
  sim-safe ones) have no top-level heavy or relative imports.

All of these read source files only (nothing is imported), except the one
``slow`` + ``simulation`` probe, which imports the package in a fresh
interpreter.
"""
import ast
import json
import os
import re
import shutil
import subprocess
import sys

import pytest

from ._support import (NSK_DIR, PKG_DIR, REPO_ROOT, assert_ran_ok, requires_simulation,
                       run_isolated)

pytestmark = pytest.mark.isobutanol

#: Directories never scanned (generated / third-party / local-only content).
_SKIP_DIRS = {'__pycache__', 'node_modules', '.git', '.ipynb_checkpoints'}

#: The upstream package name, spelled so this file never matches its own scan.
_UPSTREAM = 'biorefineries'
#: A reference to the pre-migration package: ``<upstream>.isobutanol`` not
#: preceded by ``nskinetics.`` (or any other dotted prefix), or
#: ``from <upstream> import isobutanol``.
_STALE = re.compile(r'(?<![\w.])' + _UPSTREAM + r'\.isobutanol\b'
                    r'|(?<![\w.])' + _UPSTREAM + r'\s+import\s+isobutanol\b')

#: Top-level names a by-path (sim-safe) module must never import at module level.
_HEAVY_ROOTS = frozenset({'nskinetics', _UPSTREAM, 'biosteam', 'thermosteam',
                          'hensmith', 'tellurium', 'roadrunner'})

#: Modules loaded BY FILE PATH (relative to the package folder), besides every
#: module under plots/ (MIGRATION_SPEC non-negotiable 8).
_BY_PATH_MODULES = (
    'paths.py',
    'enzyme_burden.py',
    'kinetic_optimization.py',
    'sensitivity_analysis.py',
    'sweep_utils.py',
    os.path.join('analyses', 'supervise_sweep.py'),
    os.path.join('analyses', 'kinetic_optimization', 'optimize_kinetics_BO_supervised.py'),
    os.path.join('analyses', 'opt_scenarios', 'build_opt_split12d_workbooks.py'),
)


def _walk_files(root, suffix='.py'):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if d not in _SKIP_DIRS)
        for name in sorted(filenames):
            if name.endswith(suffix):
                yield os.path.join(dirpath, name)


def _parse(path):
    with open(path, encoding='utf-8') as fh:
        return ast.parse(fh.read(), path)


def _module_level_imports(tree):
    """(root name, level, line) of every import executed at module level
    (including inside module-level if / try / with blocks, but not inside
    functions or classes)."""
    found = []

    def visit(body):
        for node in body:
            if isinstance(node, ast.Import):
                found.extend((a.name.split('.')[0], 0, node.lineno) for a in node.names)
            elif isinstance(node, ast.ImportFrom):
                found.append(((node.module or '').split('.')[0], node.level, node.lineno))
            elif isinstance(node, (ast.If, ast.For, ast.While, ast.With, ast.Try)):
                visit(node.body)
                visit(getattr(node, 'orelse', []))
                visit(getattr(node, 'finalbody', []))
                for handler in getattr(node, 'handlers', []):
                    visit(handler.body)
    visit(tree.body)
    return found


def _scenario_workbooks():
    """{scenario key: workbook filename} parsed statically from the
    ``SCENARIOS = {...}`` registry in scenarios.py (no import)."""
    tree = _parse(os.path.join(PKG_DIR, 'scenarios.py'))
    for node in tree.body:
        targets = (node.targets if isinstance(node, ast.Assign)
                   else [node.target] if isinstance(node, ast.AnnAssign) else [])
        if any(isinstance(t, ast.Name) and t.id == 'SCENARIOS' for t in targets):
            assert isinstance(node.value, ast.Dict), 'SCENARIOS is not a dict literal'
            workbooks = {}
            for key, value in zip(node.value.keys, node.value.values):
                assert isinstance(key, ast.Constant), ast.dump(key)
                assert isinstance(value, ast.Call), ast.dump(value)
                kw = {k.arg: k.value for k in value.keywords}
                assert 'workbook' in kw and isinstance(kw['workbook'], ast.Constant), key.value
                workbooks[key.value] = kw['workbook'].value
            return workbooks
    raise AssertionError('no SCENARIOS assignment found in scenarios.py')


def _git_ls_files(*args):
    git = shutil.which('git')
    if git is None:
        pytest.skip('git is not available')
    proc = subprocess.run([git, 'ls-files', '-z', *args], cwd=REPO_ROOT,
                          capture_output=True, text=True, encoding='utf-8')
    if proc.returncode != 0:
        pytest.skip(f'not a git checkout: {proc.stderr.strip()}')
    return [p for p in proc.stdout.split('\0') if p]


def test_stale_reference_pattern():
    """The scan pattern itself: it flags the old package and spares the
    migrated one."""
    old = _UPSTREAM + '.isobutanol'
    for bad in (f'from {old} import system', f'import {old}', f"'{old}.scenarios'",
                f'from {_UPSTREAM} import isobutanol'):
        assert _STALE.search(bad), bad
    for good in (f'from nskinetics.{old} import system',
                 f'import nskinetics.{old}.paths',
                 f'from nskinetics.{_UPSTREAM} import isobutanol',
                 f'{_UPSTREAM}/isobutanol/analyses', f'{_UPSTREAM}.corn'):
        assert not _STALE.search(good), good


def test_no_stale_package_references():
    """No .py file in nskinetics/ refers to the pre-migration package (only
    nskinetics.<...>.isobutanol is allowed), in code, strings or comments."""
    hits = []
    for path in _walk_files(NSK_DIR):
        with open(path, encoding='utf-8', errors='replace') as fh:
            for lineno, line in enumerate(fh, 1):
                if _STALE.search(line):
                    hits.append(f'{os.path.relpath(path, REPO_ROOT)}:{lineno}: {line.strip()}')
    assert not hits, 'stale references to the pre-migration package:\n' + '\n'.join(hits)


def test_nothing_collectable_inside_the_package():
    """No test_*.py / *_test.py / conftest.py under nskinetics/biorefineries:
    `pytest .` would import it (building the biorefinery, or worse) at
    collection time. Script-style checks live in nskinetics/tests/isobutanol."""
    root = os.path.join(NSK_DIR, 'biorefineries')
    bad = [os.path.relpath(p, REPO_ROOT) for p in _walk_files(root)
           if re.fullmatch(r'test_.*\.py|.*_test\.py|conftest\.py', os.path.basename(p))]
    assert not bad, bad


def test_scripts_folder_is_not_collectable():
    """The fresh-process scripts under tests/isobutanol/scripts never match
    pytest's collection patterns."""
    scripts = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts')
    if not os.path.isdir(scripts):
        pytest.skip('no scripts folder')
    bad = [os.path.relpath(p, REPO_ROOT) for p in _walk_files(scripts)
           if re.fullmatch(r'test_.*\.py|.*_test\.py|conftest\.py', os.path.basename(p))]
    assert not bad, bad


def test_every_scenario_workbook_ships():
    """Every workbook SCENARIOS names exists under data/parameter_distributions
    and follows paths.workbook_path's naming for its key."""
    workbooks = _scenario_workbooks()
    assert {'A', 'B'} <= set(workbooks), sorted(workbooks)
    folder = os.path.join(PKG_DIR, 'data', 'parameter_distributions')
    missing = [wb for wb in workbooks.values() if not os.path.isfile(os.path.join(folder, wb))]
    assert not missing, missing
    misnamed = {key: wb for key, wb in workbooks.items()
                if wb != f'parameter-distributions_corn_IBO_EtOH_{key}.xlsx'}
    assert not misnamed, misnamed


def test_scenario_workbooks_are_tracked():
    """The registry's workbooks are committed (``*.xlsx`` is gitignored
    except for these, so a new one is easy to forget)."""
    tracked = set(_git_ls_files('--', 'nskinetics/biorefineries/isobutanol/data'))
    prefix = 'nskinetics/biorefineries/isobutanol/data/parameter_distributions/'
    untracked = [wb for wb in _scenario_workbooks().values() if prefix + wb not in tracked]
    assert not untracked, untracked


def test_no_office_lock_files_tracked():
    """No Excel/Word lock file (``~$...``) is tracked anywhere."""
    locks = [p for p in _git_ls_files() if os.path.basename(p).startswith('~$')]
    assert not locks, locks


def test_paths_module_is_stdlib_only():
    """paths.py (loaded by file path everywhere) imports only the standard
    library, anywhere in the module, and has no relative imports."""
    tree = _parse(os.path.join(PKG_DIR, 'paths.py'))
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend((a.name.split('.')[0], 0) for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(((node.module or '').split('.')[0], node.level))
    assert imports, 'no imports found: parse error?'
    relative = [name for name, level in imports if level]
    assert not relative, relative
    non_stdlib = sorted({name for name, _ in imports if name not in sys.stdlib_module_names})
    assert not non_stdlib, non_stdlib


def _by_path_modules():
    files = [os.path.join(PKG_DIR, rel) for rel in _BY_PATH_MODULES]
    files += list(_walk_files(os.path.join(PKG_DIR, 'plots')))
    return files


@pytest.mark.parametrize('path', _by_path_modules(),
                         ids=lambda p: os.path.relpath(p, PKG_DIR).replace(os.sep, '/'))
def test_by_path_module_contract(path):
    """A module loaded BY FILE PATH has no module-level import of nskinetics,
    biosteam, thermosteam, tellurium, hensmith or the upstream biorefineries
    package, and no relative import anywhere (lazy absolute imports inside
    functions are fine)."""
    assert os.path.isfile(path), path
    tree = _parse(path)
    heavy = [(name, line) for name, level, line in _module_level_imports(tree)
             if name in _HEAVY_ROOTS]
    assert not heavy, heavy
    relative = [node.lineno for node in ast.walk(tree)
                if isinstance(node, ast.ImportFrom) and node.level]
    assert not relative, f'relative imports at lines {relative}'


_SENTINEL = 'NSK_PROBE '


@pytest.mark.slow
@pytest.mark.simulation
@requires_simulation
def test_package_import_never_loads_the_old_copy():
    """Importing nskinetics.biorefineries.isobutanol (and calling nothing)
    loads no module of the pre-migration copy, stays build-free (system not
    loaded) and does not import the plots subpackage or contourplots. Runs in
    a fresh interpreter: the package import pulls biosteam builds and the
    upstream corn/cellulosic biorefineries, which must never enter this
    pytest process."""
    code = (
        "import json, sys\n"
        "import nskinetics.biorefineries.isobutanol as ibo\n"
        f"print({_SENTINEL!r} + json.dumps({{\n"
        "    'modules': sorted(sys.modules),\n"
        "    'file': ibo.__file__,\n"
        "    'loaded': getattr(ibo.system, '_loaded', 'missing'),\n"
        "}))\n")
    proc = run_isolated(['-c', code], timeout=900, simulation=True)
    out = assert_ran_ok(proc, _SENTINEL)
    line = [ln for ln in out.splitlines() if ln.startswith(_SENTINEL)][-1]
    probe = json.loads(line[len(_SENTINEL):])
    old = _UPSTREAM + '.isobutanol'
    stale = [m for m in probe['modules'] if m == old or m.startswith(old + '.')]
    assert not stale, stale
    assert os.path.normcase(os.path.dirname(probe['file'])) == os.path.normcase(PKG_DIR)
    assert probe['loaded'] is False
    assert 'nskinetics.biorefineries.isobutanol.plots' not in probe['modules']
    assert 'contourplots' not in probe['modules']
