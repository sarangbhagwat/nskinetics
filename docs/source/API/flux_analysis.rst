Flux analysis and mapping
=========================

.. currentmodule:: nskinetics

:func:`compute_flux_summary` re-evaluates a simulated model's own rate laws
along its stored trajectory to obtain per-reaction cumulative flux and
per-inhibitor fraction-lost, returned as a :class:`FluxSummary`.
:func:`nskinetics.visualization.draw_flux_map` renders one or more summaries
side by side on a :class:`nskinetics.visualization.FluxMapSpec` layout.

Cumulative flux is reported in grams of *each step's own substrate* per litre
of final broth. Two reactions are therefore directly comparable where they
compete for the same substrate -- a branch point -- but not along a chain,
because the stoichiometry, and so the mass carried per unit of flux, differs
from step to step.

.. autoclass:: FluxSummary
   :members:

.. autofunction:: compute_flux_summary

.. currentmodule:: nskinetics.visualization

.. autoclass:: FluxMapSpec
   :members:

.. autofunction:: draw_flux_map

The shipped *S. cerevisiae* ethanol/isobutanol model ships its own layout
(``FLUX_MAP_SPEC``), scenario presets and an end-to-end
``draw_scenario_flux_map`` helper; see :doc:`models`.
