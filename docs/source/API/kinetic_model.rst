KineticModel
============

.. currentmodule:: nskinetics

:class:`KineticModel` is a wrapper around a Tellurium extended
RoadRunner model. It adds unit-aware value access and a Python event API.
:meth:`~nskinetics.KineticModel.simulate` is the primary entry
point: it integrates the model and stores the trajectory on the object
(``km.results_df`` / ``.results_array`` / ``.results_dict`` /
``.results_col_names``).
:meth:`~nskinetics.KineticModel.plot_simulation_results` (aliases
:meth:`~nskinetics.KineticModel.plot_time_course` and
:meth:`~nskinetics.KineticModel.plot_trajectory`) plots the most
recent run. The underlying RoadRunner object remains directly accessible as
``km._te`` for lower-level use; the :doc:`process bridge <units>` can also
drive it.

.. autoclass:: KineticModel
   :members:
   :show-inheritance:
