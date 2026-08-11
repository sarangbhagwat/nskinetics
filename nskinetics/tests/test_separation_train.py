# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
import pytest
import biosteam as bst

from nskinetics.processes import (
    create_ethanol_isobutanol_separation_chemicals,
    create_beer_feed,
)


def test_chemicals_cover_the_beer():
    chems = create_ethanol_isobutanol_separation_chemicals()
    for ID in ('Water', 'Ethanol', 'Isobutanol', 'AceticAcid', 'Glucose',
               'Yeast', 'Fiber', 'SolubleProtein', 'InsolubleProtein',
               'TriOlein', 'Ash', 'CaO', 'H2SO4'):
        assert ID in chems, f'{ID} missing from the chemical set'


def test_beer_feed_matches_the_scenario_B_fixture():
    bst.settings.set_thermo(create_ethanol_isobutanol_separation_chemicals())
    beer = create_beer_feed()
    assert beer.F_mass == pytest.approx(148769.7, rel=1e-4)
    assert beer.imass['Ethanol'] == pytest.approx(8897.78, rel=1e-4)
    assert beer.imass['Isobutanol'] == pytest.approx(7072.51, rel=1e-4)
    assert beer.imass['Water'] == pytest.approx(116980.42, rel=1e-4)
    assert beer.T == pytest.approx(305.15)
    assert beer.phase == 'l'
