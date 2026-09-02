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

Shipped model: scenarios and flux map
-------------------------------------

.. currentmodule:: nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo

The shipped *S. cerevisiae* ethanol/isobutanol model carries a ready-made
layout and a one-call end-to-end helper. ``FLUX_MAP_SPEC.reactions`` is the
reaction list to hand to :func:`~nskinetics.compute_flux_summary`: the drawn
edges plus mapped-but-undrawn reactions such as biomass decay (``r10``).

Its ``inhibition_map`` includes, alongside the product-inhibition
coefficients, the thermodynamic reverse-reaction (product) terms of the two
reversible steps (``k_6r`` on ``r6``, ``k_16r`` on ``r16``). A strip row
therefore reports everything the product does to that step: ADH's "fraction
lost to ethanol" counts both the inhibition of the forward rate and the
ethanol-driven reverse flux.

:func:`draw_scenario_flux_map` changes only the Ehrlich rate constants between
its two panels, never the fed-batch feeding strategy, so panel **b** is
scenario A's process with the Ehrlich branch switched on -- not the isobutanol
biorefinery's full scenario-B configuration.

.. autodata:: FLUX_MAP_SPEC
   :annotation:

.. autofunction:: draw_scenario_flux_map

.. autofunction:: apply_scenario_A

.. autofunction:: apply_scenario_B
