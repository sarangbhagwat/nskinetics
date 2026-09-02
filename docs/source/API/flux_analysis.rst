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

Shipped model: scenarios and flux map
-------------------------------------

.. currentmodule:: nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo

The shipped *S. cerevisiae* ethanol/isobutanol model carries a ready-made
layout and a one-call end-to-end helper. ``FLUX_MAP_SPEC.reactions`` is the
reaction list to hand to :func:`~nskinetics.compute_flux_summary`: the drawn
edges plus mapped-but-undrawn reactions such as biomass decay (``r10``).

.. autodata:: FLUX_MAP_SPEC
   :annotation:

.. autofunction:: draw_scenario_flux_map

.. autofunction:: apply_scenario_A

.. autofunction:: apply_scenario_B
