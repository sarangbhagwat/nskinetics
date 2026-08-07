TelluriumReactionSystem
=======================

.. currentmodule:: nskinetics

:class:`TelluriumReactionSystem` is a wrapper around a Tellurium extended
RoadRunner model. It adds unit-aware value access and a Python event API.
:meth:`~nskinetics.TelluriumReactionSystem.simulate` is the primary entry
point: it integrates the model and stores the trajectory on the object
(``trs.results_df`` / ``.results_array`` / ``.results_dict`` /
``.results_col_names``).
:meth:`~nskinetics.TelluriumReactionSystem.plot_simulation_results` (aliases
:meth:`~nskinetics.TelluriumReactionSystem.plot_time_course` and
:meth:`~nskinetics.TelluriumReactionSystem.plot_trajectory`) plots the most
recent run. The underlying RoadRunner object remains directly accessible as
``trs._te`` for lower-level use; the :doc:`process bridge <units>` can also
drive it.

.. autoclass:: TelluriumReactionSystem
   :members:
   :show-inheritance:
