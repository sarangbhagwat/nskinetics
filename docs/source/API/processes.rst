Process section factories (biosteam)
====================================

.. currentmodule:: nskinetics.processes

Sugar preparation and fermentation
----------------------------------

:obj:`create_sugar_prep_and_fermentation_system` packages the sugar-solution
preparation and fermentation section of the
:doc:`isobutanol biorefinery case study </case_studies/isobutanol/index>`
(whose flowsheet builds its fermentor ``V406`` with this factory) as a
biosteam ``SystemFactory``: saccharified sugar slurry and seed culture in;
fermentation vent, fermentation effluent, and two evaporator condensates
out. The factory builds a
:class:`~nskinetics.units.FedBatchStrategySpecification` from its own units
and attaches it to the fermentor as ``V406.fbs_spec`` — construction is
side-effect-free; the strategy is imposed only when the caller invokes its
``load_specifications()``.

.. autodata:: create_sugar_prep_and_fermentation_system
   :no-value:

:func:`create_sugar_prep_and_fermentation_chemicals` is the factory's own
thermosteam chemical set, activated only when the factory is called with
``set_thermo=True`` (by default, the caller's pre-set thermo is used
unchanged).

.. autofunction:: create_sugar_prep_and_fermentation_chemicals

Beer stripper for ethanol + isobutanol beer
-------------------------------------------

``nskinetics.processes.ethanol_isobutanol_separation`` provides the first
step of recovering ethanol and isobutanol from a whole beer: a shortcut
**beer stripper** that lifts every volatile off the corn solids, plus a
matching chemical set and a standalone beer feed with the composition of
the case study's scenario-B fermentation effluent (``P301-0``). The full
solvent-free separation train of the case study is
:obj:`nskinetics.biorefineries.isobutanol.separations.create_separation_system`.

.. autofunction:: create_ethanol_isobutanol_separation_chemicals

.. autofunction:: create_beer_feed

.. autofunction:: create_beer_stripper

Module constants: ``LIGHT_IDS`` (``('Water', 'Ethanol', 'Isobutanol',
'AceticAcid')``, the chemicals the stripper resolves by default) and the
beer composition :func:`create_beer_feed` uses, ``P301_SCENARIO_B_KGHR``
(component flows in kg/h), ``P301_SCENARIO_B_T`` (305.15 K) and
``P301_SCENARIO_B_P`` (607,950 Pa).
