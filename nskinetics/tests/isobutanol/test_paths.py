# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests of the isobutanol package's stdlib-only ``paths.py`` (loaded by file
path, as the sim-safe scripts load it): the results root and its
``NSKINETICS_IBO_RESULTS_DIR`` override, the publication / analyses /
workbook folders, the workbook naming, ``require_input``'s error message and
``load_by_path``.
"""
import os
import sys

import pytest

from ._support import PKG_DIR, load_paths

pytestmark = pytest.mark.isobutanol

ENV = 'NSKINETICS_IBO_RESULTS_DIR'


@pytest.fixture(scope='module')
def paths():
    return load_paths()


def _same(a, b):
    return os.path.normcase(os.path.normpath(a)) == os.path.normcase(os.path.normpath(b))


def test_package_folders(paths):
    """The module's folders are the package's own, and they exist."""
    assert _same(paths.PACKAGE_DIR, PKG_DIR)
    assert _same(paths.ANALYSES_DIR, os.path.join(PKG_DIR, 'analyses'))
    assert _same(paths.DATA_DIR, os.path.join(PKG_DIR, 'data'))
    assert _same(paths.PLOTS_DIR, os.path.join(PKG_DIR, 'plots'))
    for folder in (paths.PACKAGE_DIR, paths.ANALYSES_DIR, paths.DATA_DIR, paths.PLOTS_DIR):
        assert os.path.isdir(folder), folder
    assert paths.RESULTS_DIR_ENV == ENV
    assert _same(paths.DEFAULT_RESULTS_DIR, os.path.join(PKG_DIR, 'analyses', 'results'))
    assert all(hasattr(paths, name) for name in paths.__all__)


def test_results_dir_default(paths, monkeypatch):
    """Without the environment variable (or with it empty) the results root
    is <package>/analyses/results."""
    monkeypatch.delenv(ENV, raising=False)
    assert paths.results_dir() == os.path.join(paths.DEFAULT_RESULTS_DIR)
    assert paths.results_dir('a', 'b.csv') == os.path.join(paths.DEFAULT_RESULTS_DIR, 'a', 'b.csv')
    monkeypatch.setenv(ENV, '')
    assert paths.results_dir() == os.path.join(paths.DEFAULT_RESULTS_DIR)


def test_results_dir_honours_env(paths, monkeypatch, tmp_path):
    """NSKINETICS_IBO_RESULTS_DIR moves the results root, read at call time;
    create=True makes the directory."""
    root = str(tmp_path / 'ibo_results')
    monkeypatch.setenv(ENV, root)
    assert paths.results_dir() == os.path.join(root)
    assert paths.results_dir('sub', 'x.csv') == os.path.join(root, 'sub', 'x.csv')
    assert not os.path.exists(root)                        # create defaults to False
    made = paths.results_dir('campaign', 'run1', create=True)
    assert made == os.path.join(root, 'campaign', 'run1') and os.path.isdir(made)
    assert paths.results_dir('campaign', 'run1', create=True) == made   # exist_ok
    monkeypatch.delenv(ENV)
    assert paths.results_dir() == paths.DEFAULT_RESULTS_DIR


def test_results_dir_relative_env_is_absolute(paths, monkeypatch, tmp_path):
    """A relative NSKINETICS_IBO_RESULTS_DIR is resolved against the working
    directory at call time, so a script that builds its paths and then
    chdir()s into the results root still finds them."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv(ENV, os.path.join('rel', 'ibo_results'))
    root = paths.results_dir()
    assert os.path.isabs(root)
    assert _same(root, os.path.join(str(tmp_path), 'rel', 'ibo_results'))
    made = paths.results_dir('sub', create=True)
    assert os.path.isabs(made) and os.path.isdir(made)
    assert _same(made, os.path.join(root, 'sub'))
    assert _same(paths.publication_dir('Uncertainty'),
                 os.path.join(root, 'publication', 'Uncertainty'))
    checkpoint = os.path.join(made, 'x_checkpoint.csv')
    with open(checkpoint, 'w') as fh:
        fh.write('x\n')
    monkeypatch.chdir(root)                 # what the sweeps do mid-run
    assert os.path.exists(checkpoint)


def test_publication_dir(paths, monkeypatch, tmp_path):
    """publication_dir is <results root>/publication and follows the
    override."""
    monkeypatch.delenv(ENV, raising=False)
    assert paths.publication_dir() == os.path.join(paths.DEFAULT_RESULTS_DIR, 'publication')
    assert paths.publication_dir('Uncertainty', 'f.png') == os.path.join(
        paths.DEFAULT_RESULTS_DIR, 'publication', 'Uncertainty', 'f.png')
    root = str(tmp_path)
    monkeypatch.setenv(ENV, root)
    assert paths.publication_dir('Uncertainty') == os.path.join(root, 'publication', 'Uncertainty')
    made = paths.publication_dir('Uncertainty', create=True)
    assert os.path.isdir(made)


def test_analyses_dir(paths):
    """analyses_dir joins under the analyses folder; the thematic subfolders
    exist."""
    assert paths.analyses_dir() == os.path.join(paths.ANALYSES_DIR)
    assert paths.analyses_dir('sensitivity', 'x.py') == os.path.join(
        paths.ANALYSES_DIR, 'sensitivity', 'x.py')
    for sub in ('feeding_strategies', 'kinetic_sweeps', 'objective_landscape',
                'kinetic_optimization', 'opt_scenarios', 'sensitivity',
                'tea_breakdowns', 'uncertainty'):
        assert os.path.isdir(paths.analyses_dir(sub)), sub


def test_workbook_naming(paths):
    """workbook_path follows the registry naming under
    data/parameter_distributions; the A and B workbooks ship."""
    assert paths.workbook_dir() == os.path.join(paths.DATA_DIR, 'parameter_distributions')
    for scenario in ('A', 'B', 'opt_PI_TRY_informed'):
        path = paths.workbook_path(scenario)
        assert os.path.dirname(path) == paths.workbook_dir()
        assert os.path.basename(path) == (
            f'parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx')
    assert os.path.isfile(paths.workbook_path('A'))
    assert os.path.isfile(paths.workbook_path('B'))


def test_require_input(paths, tmp_path):
    """require_input returns an existing path unchanged and otherwise raises
    FileNotFoundError naming the path, its producer and the override."""
    present = tmp_path / 'present.csv'
    present.write_text('x\n')
    assert paths.require_input(str(present), 'nothing') == str(present)
    assert paths.require_input(str(tmp_path), 'nothing') == str(tmp_path)  # directories too
    missing = str(tmp_path / 'missing.csv')
    producer = 'analyses/sensitivity/evaluate_sobol_split12d.py'
    with pytest.raises(FileNotFoundError) as e:
        paths.require_input(missing, producer)
    msg = str(e.value)
    assert repr(missing) in msg
    assert f'produced by {producer}' in msg
    assert 'Run that analysis first' in msg and ENV in msg


def test_load_by_path(paths, tmp_path):
    """load_by_path imports a file under the given name without importing a
    parent package; register=True adds it to sys.modules before executing
    it."""
    src = tmp_path / 'some_module_file.py'
    src.write_text('import sys\n'
                   'VALUE = 42\n'
                   'REGISTERED = __name__ in sys.modules\n')
    name = '_nsk_test_paths_probe_module'
    assert name not in sys.modules
    try:
        mod = paths.load_by_path(name, str(src))
        assert mod.__name__ == name and mod.VALUE == 42
        assert mod.__file__ == str(src)
        assert mod.REGISTERED is False and name not in sys.modules
        mod2 = paths.load_by_path(name, str(src), register=True)
        assert mod2.REGISTERED is True and sys.modules[name] is mod2
        assert mod2 is not mod
    finally:
        sys.modules.pop(name, None)
