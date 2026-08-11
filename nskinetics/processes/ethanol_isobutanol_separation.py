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
           'create_beer_stripper',
           'LIGHT_IDS',
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
    # Cp is specific (J/g/K) both as returned by Chemical.Cp and as taken by
    # the Chemical constructor, so glucose's value is used as-is here (~1.5
    # J/g/K, the right order for corn solids), matching the Yeast definition
    # above and the sibling sugar-prep + fermentation chemical set.
    for ID in ('Fiber', 'SolubleProtein', 'InsolubleProtein', 'Ash', 'CaO'):
        chemicals.append(_solid(ID, MW=100., Cp=Glucose.Cp(298.15)))
    # 'TriOlein' is not in the thermo database under any of its searchable
    # names, so the corn oil it stands for is carried as an inert solid
    # pseudochemical (MW of triolein, 885.4 g/mol; density and specific heat
    # of corn oil, 915 kg/m3 and 2.0 J/g/K, rather than the generic
    # mineral-solid defaults).
    chemicals.append(_solid('TriOlein', MW=885.4, rho=915., Cp=2.0))
    chemicals.append(tmo.Chemical('H2SO4', default=True))
    for chemical in chemicals:
        chemical.default()
    chemicals.compile()
    chemicals.set_synonym('Water', 'H2O')
    return chemicals


def create_beer_feed(ID='beer'):
    """
    Create the scenario-B ``P301-0`` beer stream as a standalone feed.

    Requires the chemical set from
    :func:`create_ethanol_isobutanol_separation_chemicals` (or any set
    containing all of its IDs) to already be the active thermo, e.g. via
    ``bst.settings.set_thermo(...)``; otherwise :class:`bst.Stream` raises.

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


LIGHT_IDS = ('Water', 'Ethanol', 'Isobutanol', 'AceticAcid')


def create_beer_stripper(ID='C1', ins=None, outs=None, light_IDs=LIGHT_IDS,
                         ethanol_recovery=0.9995, water_to_stillage=0.60,
                         k=1.5, P=101325.):
    """
    Create the beer stripper: lifts every volatile off the corn solids.

    Modelled as a :class:`biosteam.ShortcutColumn` on the volatile subset
    only. The solids are zeroed from the column feed, then restored to the
    bottoms after the column runs -- the same zero-out/restore pattern the
    isobutanol biorefinery uses for its extractor, and necessary because a
    stage-by-stage VLE model cannot carry inert solids.

    Ethanol and water are the light and heavy keys. Isobutanol is not
    specified: although it boils above water (108 vs 100 degC), its activity
    coefficient in dilute aqueous solution makes it far more volatile than
    that, and the column's own Hengsteback-Gaddes solution -- which ranks
    non-keys by mean relative volatility taken from real bubble- and
    dew-point calculations, not by boiling point -- carries essentially all
    of it overhead on its own.

    Parameters
    ----------
    ID : str
        Unit ID.
    ins : stream
        Whole beer.
    outs : tuple
        ``(volatile concentrate, stillage)``.
    light_IDs : tuple of str
        Chemicals the column is allowed to resolve. Everything else is
        routed unchanged to the stillage.
    ethanol_recovery : float
        Fraction of the feed ethanol -- the light key -- recovered to the
        concentrate.
    water_to_stillage : float
        Fraction of the feed water -- the heavy key -- leaving in the
        stillage. The remainder goes overhead and sets how dilute the
        concentrate is.
    k : float
        Reflux ratio over the minimum reflux ratio. This, not a stage count,
        is what sizes a shortcut column: the Gilliland correlation diverges
        as ``k`` approaches 1, and because a beer column's minimum reflux is
        essentially zero (it is a stripper), ``k`` must be kept comfortably
        above 1. At the default the column comes out near 21 theoretical /
        67 actual stages; at ``k = 1.05`` it runs away past 300.
    P : float
        Column pressure [Pa].

    Returns
    -------
    biosteam.Unit
    """
    C1 = bst.ShortcutColumn(ID, ins=ins, outs=outs,
                            LHK=('Ethanol', 'Water'),
                            Lr=ethanol_recovery, Hr=water_to_stillage,
                            k=k, P=P, partial_condenser=False)
    C1.light_IDs = tuple(light_IDs)

    @C1.add_specification(run=False)
    def C1_volatiles_only():
        feed = C1.ins[0]
        held = {}
        for chem in feed.chemicals:
            if chem.ID not in C1.light_IDs:
                held[chem.ID] = feed.imol[chem.ID]
                feed.imol[chem.ID] = 0.
        try:
            C1._run()
        finally:
            # Give the feed back whatever the column borrowed, even if it
            # failed to converge, so a retry sees an intact stream.
            for ID, mol in held.items():
                feed.imol[ID] = mol
        bottoms = C1.outs[1]
        for ID, mol in held.items():
            bottoms.imol[ID] = bottoms.imol[ID] + mol
    return C1
