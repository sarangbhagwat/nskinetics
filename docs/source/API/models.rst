Shipped model
=============

.. currentmodule:: nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo

The shipped *S. cerevisiae* ethanol/isobutanol model (``te_r``) ships with
kinetic-only scenario presets, a flux-map layout, and a classification of
its kinetic parameters that names a parameter change in a few words.

Every name documented on this page is re-exported from the subpackage, so
``from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import
describe_parameter_change`` works without naming the defining submodule.

Scenarios and flux map
----------------------

The shipped model carries a ready-made layout and a one-call end-to-end
helper. ``FLUX_MAP_SPEC.reactions`` is the reaction list to hand to
:func:`~nskinetics.compute_flux_summary`: the drawn edges plus
mapped-but-undrawn reactions such as biomass decay (``r10``).

Its ``inhibition_map`` includes, alongside the product-inhibition
coefficients, the thermodynamic reverse-reaction (product) terms of the two
reversible steps (``k_6r`` on ``r6``, ``k_16r`` on ``r16``). A strip row
therefore reports everything the product does to that step: ADH's "fraction
lost to ethanol" counts both the inhibition of the forward rate and the
ethanol-driven reverse flux.

:func:`~nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo.flux_map_spec.draw_scenario_flux_map`
changes only the Ehrlich rate constants between its two panels, never the
fed-batch feeding strategy, so panel **b** is scenario A's process with the
Ehrlich branch switched on -- not the isobutanol biorefinery's full
scenario-B configuration.

.. currentmodule:: nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo.flux_map_spec

.. autodata:: FLUX_MAP_SPEC
   :annotation:

.. autofunction:: draw_scenario_flux_map

.. currentmodule:: nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo.scenarios

.. autofunction:: apply_scenario_A

.. autofunction:: apply_scenario_B

Parameter categories
--------------------

.. currentmodule:: nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo.parameter_categories

Every kinetic parameter of the shipped model has a *role* -- what kind of
term it is in its rate law -- and belongs to one or more pathway *modules*
through the reactions it appears in. :func:`describe_parameter_change`
turns a baseline-versus-current diff into direction + role + effector +
modules, e.g. ``stronger isobutanol inhibition of glycolysis/fermentation,
overflow/acetate and growth`` or ``Ehrlich branch on``.

Roles, in the order their clauses appear: ``capacity`` (multiplies the whole
rate term), ``affinity`` (a ``S/(S + K)`` saturation constant),
``substrate_regulation`` (glucose repression or the acetaldehyde signal),
``product_inhibition`` (an ``exp(-k*P)`` term), ``product_self_inhibition``
(the ``K_6e``/``k_6r``-style denominator and reverse terms), ``lethality``
and ``lethality_threshold`` (the steepness and onset of the product-enhanced
biomass decay on ``r10``), and ``initial_state`` (``X_a``, ``X_AcDH``).
Modules: glycolysis/fermentation (r1, r3, r6), respiration (r2, r5),
overflow/acetate (r4), growth (r7, r8), physiological state (r9--r11), and
the Ehrlich branch (r13--r16).

Operation parameters -- fed-batch feeding, aeration staging, dilution, and
the vessel's oxygen-transfer capacity -- are not kinetics and are skipped by
the helpers. The last of these is the pair ``kLa`` (volumetric O2
mass-transfer coefficient, default 200 /h) and ``C_O2_sat`` (dissolved-O2
saturation, default 0.20 mmol/L): the model caps the realized respiratory
O2 uptake of its aerobic fluxes at ``OTR_max = kLa*C_O2_sat`` (40 mmol
O2/L/h at the defaults) through the assignment-rule factor ``f_O2``, which
scales r2, r5, r8 and the aerobic share of r7 so that uptake equals
``OTR_max`` exactly while the bound binds, and reads 1 wherever it does not
(and 0 while ``is_aerobic`` is off). Only the growth-associated O2 of the
ungated ``anaerobic_growth_mult`` share of r7 (1.0 by default) is committed
and never throttled. The shipped fed-batch process (aeration off at 5 g/L
biomass, peak demand ~20 mmol/L/h) is never bound; a glucose-lean, fully
aerobic batch that grows past ~10 g/L biomass is, while glucose-rich runs
stay under the cap because r2, r5 and r8 are glucose-inhibited. Set either
parameter on the model (``te_r._te.kLa = 100``) to describe a different vessel;
``f_O2`` is a default result column of the process factory so a bound run
is visible in ``V406.nsk_results_df``.

.. autodata:: ROLES
   :annotation:

.. autodata:: MODULES
   :annotation:

.. autodata:: MODULE_LABELS
   :annotation:

.. autodata:: KINETIC_PARAMETERS
   :annotation:

.. autodata:: OPERATION_PARAMETERS
   :annotation:

.. autodata:: REACTION_MODULES
   :annotation:

.. autoclass:: ParameterInfo
   :members:

.. autofunction:: snapshot_parameters

.. autofunction:: diff_parameters

.. autofunction:: describe_parameter_change

.. autofunction:: categorize
