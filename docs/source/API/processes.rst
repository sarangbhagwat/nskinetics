Process section factories (biosteam)
====================================

.. currentmodule:: nskinetics.processes

:obj:`create_sugar_prep_and_fermentation_system` packages the isobutanol
biorefinery's sugar-solution preparation and fermentation section as a
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
