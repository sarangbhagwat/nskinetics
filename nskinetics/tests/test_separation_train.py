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

SOLID_IDS = ('Fiber', 'SolubleProtein', 'InsolubleProtein', 'Ash', 'CaO',
             'Yeast', 'Glucose')


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


@pytest.fixture
def unit_flowsheet(request):
    """Give a unit-building test its own flowsheet, so unit IDs never collide
    across tests, and hand the module's flowsheet back afterwards."""
    previous_flowsheet = bst.main_flowsheet.get_flowsheet()
    bst.main_flowsheet.set_flowsheet(bst.Flowsheet(request.node.name))
    try:
        yield bst.main_flowsheet.get_flowsheet()
    finally:
        bst.main_flowsheet.set_flowsheet(previous_flowsheet)


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


def test_beer_stripper_lifts_both_alcohols_off_the_solids(unit_flowsheet):
    from nskinetics.processes import create_beer_stripper
    beer = create_beer_feed()
    # Captured before simulation so the solids assertions below also prove the
    # bypass restores the feed it borrowed from.
    fed = {ID: beer.imass[ID] for ID in beer.chemicals.IDs}
    fed_F_mass = beer.F_mass
    C1 = create_beer_stripper('C1', ins=beer, outs=('concentrate', 'stillage'))
    C1.simulate()
    conc, stillage = C1.outs
    # both alcohols must come overhead, out of the column's own VLE solution
    assert conc.imass['Ethanol'] / fed['Ethanol'] > 0.99
    assert conc.imass['Isobutanol'] / fed['Isobutanol'] > 0.95
    # every solid must stay in the stillage
    for ID in SOLID_IDS:
        assert conc.imass[ID] == pytest.approx(0.0, abs=1e-6)
        assert stillage.imass[ID] == pytest.approx(fed[ID], rel=1e-6)
    # overall mass balance
    assert conc.F_mass + stillage.F_mass == pytest.approx(fed_F_mass, rel=1e-6)
    # and the shortcut design must not run away against the Gilliland
    # singularity at minimum reflux (which it does for k near 1)
    assert C1.design_results['Actual stages'] < 100
