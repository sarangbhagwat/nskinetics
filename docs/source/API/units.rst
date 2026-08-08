Process bridge (biosteam)
=========================

.. currentmodule:: nskinetics.units

:class:`NSKBatchReactor` drives any :class:`~nskinetics.KineticModel`
inside a biosteam ``BatchBioreactor``, and :class:`FermentationSaccharomycesEthanolIsobutanol` is its
fed-batch fermentation subclass.

.. autoclass:: NSKBatchReactor
   :members:
   :show-inheritance:

.. autoclass:: FermentationSaccharomycesEthanolIsobutanol
   :members:
   :show-inheritance:

.. autoclass:: AerationSpec
   :members:

.. autoclass:: SpikeReduceRetry
   :members:

:class:`FedBatchStrategySpecification` orchestrates a fed-batch feeding
strategy across the upstream flowsheet, configured by
:class:`SpikeControlVariables` (kinetic-model coupling) and
:class:`ConcentrationActuator` (upstream unit knobs).

.. autoclass:: FedBatchStrategySpecification
   :members:

.. autoclass:: SpikeControlVariables
   :members:

.. autoclass:: ConcentrationActuator
   :members:

.. autofunction:: select_tau_index
