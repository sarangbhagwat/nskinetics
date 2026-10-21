# -*- coding: utf-8 -*-
"""
Configuration for pytest to run with and without numba JIT compiled functions.
"""
import pytest
import os
import sys

@pytest.fixture(autouse=True)
def _restore_roadrunner_properties():
    """Undo libroadrunner's class-level property churn after every test.

    libroadrunner installs one Python property per model id (``r.k``,
    ``r.s_glu``...) on the ``roadrunner.roadrunner.RoadRunner`` class and
    deletes the previous model's on every load. A test that loads a toy model
    would otherwise make later attribute-style writes on a long-lived model
    (e.g. the shipped ``te_r``, ``r.env = 1.``) silently create instance
    attributes, so results would depend on test order. Properties removed
    during a test are put back afterwards (ones the test's own models define
    are left alone). roadrunner is looked up lazily so this file never
    imports it.
    """
    rrm = sys.modules.get('roadrunner.roadrunner')
    rr = getattr(rrm, 'RoadRunner', None)
    saved = ({n: rr.__dict__[n] for n in list(getattr(rr, '_properties', ()))
              if n in rr.__dict__} if rr is not None else {})
    yield
    for name, prop in saved.items():
        if name not in rr.__dict__:
            setattr(rr, name, prop)
            rr._properties.append(name)

def pytest_ignore_collect(collection_path):
    # pytest >= 8 passes `collection_path` (pathlib.Path); the old `path`
    # (py.path.local) argument was removed.
    if 'setup' in str(collection_path):
        return True

def pytest_addoption(parser):
    parser.addoption(
        "--disable-numba", action="store", default="1", help="my option: 0 or 1"
    )

def pytest_configure(config):
    os.environ["NUMBA_DISABLE_JIT"] = config.getoption("--disable-numba")
    os.environ["DISABLE_PREFERENCES"] = "1"
    os.environ["FILTER_WARNINGS"] = "1"
    os.environ["PY_IGNORE_IMPORTMISMATCH"] = "1"
    config.addinivalue_line(
        "markers",
        'slow: marks tests as slow (deselected in CI with -m "not slow")')
    for marker in (
            'isobutanol: tests of nskinetics.biorefineries.isobutanol '
            '(nskinetics/tests/isobutanol)',
            'simulation: spawns a fresh Python process that builds/simulates the '
            'isobutanol biorefinery (always also slow; strictly sequential, never '
            'under pytest-xdist)',
            'smoke: the 16 isobutanol smoke cases (scenario baselines and opt_* '
            'reproductions; python -m nskinetics.biorefineries.isobutanol.smoke N)',
            'canonical: the canonical gate -- smoke cases 1 and 2, each run twice, '
            'with a cross-run 3-significant-figure MPSP comparison',
            ):
        config.addinivalue_line("markers", marker)

def pytest_collection_modifyitems(config, items):
    """Refuse to run `simulation` tests under pytest-xdist.

    Two simulating processes at once corrupt numba's on-disk caches, so every
    `simulation`-marked test is skipped in an xdist worker (run them without
    `-n`)."""
    if not (hasattr(config, 'workerinput') or os.environ.get('PYTEST_XDIST_WORKER')):
        return
    skip = pytest.mark.skip(
        reason='simulation tests must run serially; re-run without pytest-xdist (-n)')
    for item in items:
        if item.get_closest_marker('simulation') is not None:
            item.add_marker(skip)

# def pytest_collection_modifyitems(session, config, items):
#     """Modifies test items in place to ensure test modules run in a given order."""
#     for name in ('test_sugarcane', 'test_cornstover'):
#         for i in items:
#             if i.name == name:
#                 items.remove(i)
#                 items.insert(0, i)
#                 break
    