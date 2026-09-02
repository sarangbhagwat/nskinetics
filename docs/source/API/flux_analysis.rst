Flux analysis and mapping
=========================

.. currentmodule:: nskinetics

:func:`compute_flux_summary` re-evaluates a simulated model's own rate laws
along its stored trajectory to obtain per-reaction cumulative flux and
per-inhibitor fraction-lost, returned as a :class:`FluxSummary`.
:func:`nskinetics.visualization.draw_flux_map` renders one or more summaries
side by side on a :class:`nskinetics.visualization.FluxMapSpec` layout.

.. autoclass:: FluxSummary
   :members:

.. autofunction:: compute_flux_summary

.. currentmodule:: nskinetics.visualization

.. autoclass:: FluxMapSpec
   :members:

.. autofunction:: draw_flux_map
