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

.. autofunction:: select_tau_index
