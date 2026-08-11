# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Fast tests for the ethanol/isobutanol separation train: the chemical set,
the scenario-B beer feed fixture, and the C1 beer stripper (built, simulated,
and checked for its split, its solids bypass, its stage count, and its
idempotency under re-simulation).

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

# Everything the beer carries that is not one of C1's ``LIGHT_IDS``, i.e. the
# full set the solids bypass is contractually responsible for routing to the
# stillage. 'H2SO4' is a real VLE-active chemical rather than a ``_solid``
# pseudochemical, but the bypass holds it out of the column all the same, so
# the same all-or-nothing assertion applies to it.
SOLID_IDS = ('Fiber', 'SolubleProtein', 'InsolubleProtein', 'Ash', 'CaO',
             'Yeast', 'Glucose', 'TriOlein', 'H2SO4')


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
    # singularity at minimum reflux (which it does if Rmin drops back to
    # biosteam's 0.01 numerical floor: 890 actual stages, $379 M)
    assert C1.design_results['Actual stages'] < 60


def test_beer_stripper_is_idempotent(unit_flowsheet):
    """The solids bypass *adds* the held solids back into the bottoms, which is
    only correct because ``ShortcutColumn._run`` empties its outlets first. If
    that biosteam internal ever changes, a second pass double-counts."""
    from nskinetics.processes import create_beer_stripper
    beer = create_beer_feed()
    fed = {ID: beer.imass[ID] for ID in beer.chemicals.IDs}
    fed_F_mass = beer.F_mass
    C1 = create_beer_stripper('C1', ins=beer, outs=('concentrate', 'stillage'))
    C1.simulate()
    first = {ID: (C1.outs[0].imass[ID], C1.outs[1].imass[ID])
             for ID in beer.chemicals.IDs}
    C1.simulate()
    conc, stillage = C1.outs
    # rel=1e-4, not machine precision: the column's Wegstein fixed point warm
    # starts from the previous solution, so the trace species shift in the
    # sixth significant figure (isobutanol to the stillage moves 0.518733 ->
    # 0.518738 kg/hr). A double-count would be a factor of two.
    for ID in beer.chemicals.IDs:
        assert conc.imass[ID] == pytest.approx(first[ID][0], rel=1e-4,
                                               abs=1e-9), f'{ID} to concentrate'
        assert stillage.imass[ID] == pytest.approx(first[ID][1], rel=1e-4,
                                                   abs=1e-9), f'{ID} to stillage'
    for ID in SOLID_IDS:
        assert conc.imass[ID] == pytest.approx(0.0, abs=1e-6)
        assert stillage.imass[ID] == pytest.approx(fed[ID], rel=1e-6)
    assert conc.F_mass + stillage.F_mass == pytest.approx(fed_F_mass, rel=1e-6)
    assert beer.F_mass == pytest.approx(fed_F_mass, rel=1e-6)
