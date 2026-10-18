# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Tests that the README's Python example runs and prints what it shows."""
import re
import textwrap
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[2]
# Only a source tree (repo checkout or unpacked sdist) has this package's
# README next to setup.py; in an installed wheel parents[2] is site-packages,
# where a README.rst would belong to some other distribution.
README = _ROOT / 'README.rst' if (_ROOT / 'setup.py').is_file() else None

# A ``code-block:: python`` directly followed by "prints:" and the
# ``code-block:: text`` that shows its output.
_EXAMPLE = re.compile(
    r'^\.\. code-block:: python\n\n((?:    .*\n|\n)+?)'
    r'\nprints:\n\n\.\. code-block:: text\n\n((?:    .*\n)+)', re.M)


def _examples():
    # the README is used only from a source tree (see README above)
    if README is None or not README.is_file():
        return []
    text = README.read_text(encoding='utf-8')
    return [(textwrap.dedent(code), textwrap.dedent(out))
            for code, out in _EXAMPLE.findall(text)]


def test_readme_has_an_example():
    if README is None or not README.is_file():
        pytest.skip('not run from a source tree with README.rst and setup.py')
    assert len(_examples()) == 1


# short ids (example1, ...): the code itself would be a ~2000-character id;
# without the README the empty list collects as one skipped test
@pytest.mark.parametrize('code, expected', [
    pytest.param(c, o, id=f'example{i}')
    for i, (c, o) in enumerate(_examples(), 1)])
def test_readme_example_prints_what_it_shows(code, expected, capsys):
    exec(compile(code, str(README), 'exec'), {'__name__': '__readme__'})
    assert capsys.readouterr().out == expected
