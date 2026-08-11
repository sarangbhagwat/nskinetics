# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Fast tests for the ethanol/isobutanol separation train's chemical set and
beer feed fixture.

The global thermo and main flowsheet are process-wide state, so this module
does its ``set_thermo`` inside an isolating fixture (same convention as
``test_processes.py``: pair ``set_thermo`` with a dedicated flowsheet) and
restores whatever was active before, leaving the rest of the session
untouched regardless of collection order.
"""

import pytest
import biosteam as bst

from nskinetics.processes import (
    create_ethanol_isobutanol_separation_chemicals,
    create_beer_feed,
)


@pytest.fixture(scope='module', autouse=True)
def separation_thermo():
    """Activate the separation chemical set and a dedicated flowsheet for this
    module only, restoring the previous thermo and flowsheet afterwards."""
    previous_thermo = getattr(bst.settings, '_thermo', None)
    previous_flowsheet = bst.main_flowsheet.get_flowsheet()
    bst.main_flowsheet.set_flowsheet('test_separation_train')
    bst.settings.set_thermo(create_ethanol_isobutanol_separation_chemicals())
    try:
        yield bst.settings.thermo
    finally:
        bst.main_flowsheet.set_flowsheet(previous_flowsheet)
        if previous_thermo is None:
            try:
                del bst.settings._thermo
            except AttributeError:
                pass
        else:
            bst.settings._thermo = previous_thermo


def test_chemicals_cover_the_beer():
    chems = create_ethanol_isobutanol_separation_chemicals()
    for ID in ('Water', 'Ethanol', 'Isobutanol', 'AceticAcid', 'Glucose',
               'Yeast', 'Fiber', 'SolubleProtein', 'InsolubleProtein',
               'TriOlein', 'Ash', 'CaO', 'H2SO4'):
        assert ID in chems, f'{ID} missing from the chemical set'


def test_beer_feed_matches_the_scenario_B_fixture():
    beer = create_beer_feed()
    assert beer.F_mass == pytest.approx(148769.7, rel=1e-4)
    assert beer.imass['Ethanol'] == pytest.approx(8897.78, rel=1e-4)
    assert beer.imass['Isobutanol'] == pytest.approx(7072.51, rel=1e-4)
    assert beer.imass['Water'] == pytest.approx(116980.42, rel=1e-4)
    assert beer.T == pytest.approx(305.15)
    assert beer.P == pytest.approx(607950.)
    assert beer.phase == 'l'
