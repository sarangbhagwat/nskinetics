# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
System factory recovering separate high-purity ethanol and isobutanol from a
whole-beer stream (the isobutanol biorefinery's ``P301-0``).
"""

import biosteam as bst
import thermosteam as tmo
from thermosteam import functional as fn

__all__ = ('create_ethanol_isobutanol_separation_chemicals',
           'create_beer_feed',
           'P301_SCENARIO_B_KGHR',
           'P301_SCENARIO_B_T',
           'P301_SCENARIO_B_P')

# Measured on the isobutanol biorefinery, scenario B, 2026-08-10.
P301_SCENARIO_B_KGHR = {
    'Water': 116980.42, 'Ethanol': 8897.78, 'Isobutanol': 7072.51,
    'Fiber': 6060.68, 'InsolubleProtein': 2800.29, 'TriOlein': 1931.24,
    'SolubleProtein': 1931.30, 'Yeast': 1866.50, 'Ash': 795.22,
    'Glucose': 244.71, 'AceticAcid': 134.82, 'H2SO4': 48.43, 'CaO': 5.81,
}
P301_SCENARIO_B_T = 305.15
P301_SCENARIO_B_P = 607950.


def _solid(ID, MW=1., rho=1540., Cp=1.):
    """Build an inert solid pseudochemical used only as a mass carrier."""
    chem = tmo.Chemical(ID, phase='s', phase_ref='s', search_db=False,
                        MW=MW, Cp=Cp, default=True)
    chem.V.add_model(fn.rho_to_V(rho=rho, MW=chem.MW), top_priority=True)
    return chem


def create_ethanol_isobutanol_separation_chemicals():
    """
    Create the chemical set for
    :func:`create_ethanol_isobutanol_separation_system`.

    Contains the volatiles the columns resolve (``Water``, ``Ethanol``,
    ``Isobutanol``, ``AceticAcid``) plus the corn-mash solids the beer
    carries, defined as inert solid pseudochemicals that act purely as mass
    carriers to the stillage.

    Returns
    -------
    thermosteam.Chemicals
        Compiled chemical set with synonym ``H2O`` (Water).
    """
    chemicals = tmo.Chemicals([
        'Water', 'Ethanol', 'Isobutanol', 'AceticAcid',
        tmo.Chemical('Glucose', phase='l'),
        tmo.Chemical('CO2', phase='g'),
    ])
    Glucose = chemicals.Glucose
    Glucose.N_solutes = 1
    Glucose.V.add_model(fn.rho_to_V(rho=1e5, MW=Glucose.MW), top_priority=True)
    Yeast = tmo.Chemical('Yeast', phase='s', phase_ref='s', search_db=False,
                         formula='CH1.61O0.56', rho=1540,
                         Cp=Glucose.Cp(298.15), default=True)
    Yeast.Hf = Glucose.Hf / Glucose.MW * Yeast.MW
    Yeast.V.add_model(fn.rho_to_V(rho=1540, MW=Yeast.MW), top_priority=True)
    chemicals.append(Yeast)
    for ID in ('Fiber', 'SolubleProtein', 'InsolubleProtein', 'Ash', 'CaO'):
        chemicals.append(_solid(ID, MW=100., Cp=Glucose.Cp(298.15) / 100.))
    # 'TriOlein' is not in the thermo database under any of its searchable
    # names, so the corn oil it stands for is carried as an inert solid
    # pseudochemical (MW of triolein, 885.4 g/mol).
    chemicals.append(_solid('TriOlein', MW=885.4))
    chemicals.append(tmo.Chemical('H2SO4', default=True))
    for chemical in chemicals:
        chemical.default()
    chemicals.compile()
    chemicals.set_synonym('Water', 'H2O')
    return chemicals


def create_beer_feed(ID='beer'):
    """
    Create the scenario-B ``P301-0`` beer stream as a standalone feed.

    Parameters
    ----------
    ID : str
        Stream ID. Defaults to ``'beer'``.

    Returns
    -------
    thermosteam.Stream
        Liquid beer at 305.15 K and 607950 Pa with the measured
        scenario-B composition.
    """
    beer = bst.Stream(ID, units='kg/hr', T=P301_SCENARIO_B_T,
                      P=P301_SCENARIO_B_P, phase='l',
                      **P301_SCENARIO_B_KGHR)
    return beer
