# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Runs the code examples of the documentation page of
``nskinetics.models.pd_c_3hp_oxidation``
(``docs/source/case_studies/pd_c_3hp_oxidation.rst``).

On that page every ``code-block:: python`` is followed by a
``code-block:: text`` with what the block prints. The blocks are run here in
order, in one namespace, and what they print is compared with the page: the
text must be the same, and every number must agree to within one unit of its
last printed digit, so that a last digit rounded the other way by another
platform or solver version does not fail the test. The numbers themselves
are pinned, with their own tolerances, in the other
``test_pd_c_3hp_oxidation*`` modules; this module only keeps the page and
the package from drifting apart. The shell commands of the page (the figure
CLI) are not run here; ``test_pd_c_3hp_oxidation_figures`` tests
``save_figures``, the function the CLI calls.

The page is not part of an installed package, so the test is skipped when
the file is not found; it is found in a source checkout. The examples take
about 1 s, so the module is not marked slow. It imports pytest at top level
and is therefore not registered in ``nskinetics/tests/__init__.py``.
"""
import contextlib
import io
import os
import re
import textwrap
from decimal import Decimal

import pytest

_PAGE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir,
    'docs', 'source', 'case_studies', 'pd_c_3hp_oxidation.rst')

# a python block (indented or empty lines) directly followed by a text block
_EXAMPLE = re.compile(
    r'\.\. code-block:: python\n\n((?:(?:   .*)?\n)+?)'
    r'\.\. code-block:: text\n\n((?:(?:   .*)?\n)+)')
_NUMBER = re.compile(r'-?\d+(?:\.\d+)?(?:e[+-]?\d+)?')


def _last_digit(number):
    # one unit in the last printed digit: 1e-4 for '0.6436', 1e-8 for
    # '1.094e-05', 1 for '15'
    return float(Decimal(1).scaleb(Decimal(number).as_tuple().exponent))


def _same_output(printed, shown):
    """True if two printed outputs have the same text between their numbers
    (white space aside) and numbers that agree to within one unit of the
    last printed digit (of the one printed with fewer digits)."""
    if _NUMBER.sub('#', printed).split() != _NUMBER.sub('#', shown).split():
        return False
    return all(
        abs(float(a) - float(b)) <= 1.01*max(_last_digit(a), _last_digit(b))
        for a, b in zip(_NUMBER.findall(printed), _NUMBER.findall(shown)))


def test_output_comparison_allows_only_the_last_digit():
    shown = '   time  X_MA\n    600  0.42\n   pooled R2 = 0.6436\n'
    assert _same_output(' time   X_MA\n  600 0.4201\npooled R2 = 0.6437\n',
                        shown)
    for printed in ('time X_MA\n600 0.42\npooled R2 = 0.6438\n',   # number
                    'time X_MA\n600 0.44\npooled R2 = 0.6436\n',   # number
                    'time X_D\n600 0.42\npooled R2 = 0.6436\n',    # text
                    'time X_MA\n600 0.42\n',                       # line lost
                    'time X_MA\n600 0.42 0.1\npooled R2 = 0.6436\n'):
        assert not _same_output(printed, shown), printed


def test_page_examples_print_what_the_page_shows():
    if not os.path.isfile(_PAGE):
        pytest.skip('the documentation is not part of an installed package')
    with open(_PAGE, encoding='utf-8') as f:
        text = f.read()
    examples = _EXAMPLE.findall(text)
    # every python block of the page has its printed output next to it
    assert len(examples) == text.count('.. code-block:: python') > 0
    namespace = {}
    for number, (code, shown) in enumerate(examples, 1):
        printed = io.StringIO()
        with contextlib.redirect_stdout(printed):
            exec(compile(textwrap.dedent(code), 'python block %d' % number,
                         'exec'), namespace)
        assert _same_output(printed.getvalue(), shown), (
            'python block %d of the page prints\n%s' % (
                number, printed.getvalue()))
