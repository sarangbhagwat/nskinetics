Shipped models
==============

This page documents the two model packages under ``nskinetics.models``, each
a subpackage that holds the Antimony source and the SBML of its models next
to the Python helpers documented here:

- `S. cerevisiae ethanol/isobutanol model`_
  (``nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo``): the fed-batch
  fermentation model behind :doc:`/tutorial/05_real_model` and the
  :doc:`isobutanol biorefinery demonstration </case_studies/isobutanol/index>`.
- `Pd/C oxidation of 3-HP to malonic acid`_
  (``nskinetics.models.pd_c_3hp_oxidation``): four kinetic models that
  reproduce the model curves of a published kinetic study, presented on
  :doc:`/case_studies/pd_c_3hp_oxidation`.

``import nskinetics`` imports neither package.

*S. cerevisiae* ethanol/isobutanol model
----------------------------------------

.. currentmodule:: nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo

The shipped *S. cerevisiae* ethanol/isobutanol model (``te_r``) ships with
kinetic-only scenario presets, a flux-map layout, and a classification of
its kinetic parameters that names a parameter change in a few words.
Importing the subpackage builds and configures the model.

Every name documented in this section is re-exported from the subpackage, so
``from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import
describe_parameter_change`` works without naming the defining submodule.

Scenarios and flux map
~~~~~~~~~~~~~~~~~~~~~~

The shipped model carries a ready-made layout and a one-call end-to-end
helper. ``FLUX_MAP_SPEC.reactions`` is the reaction list to hand to
:func:`~nskinetics.compute_flux_summary`: the drawn edges plus
mapped-but-undrawn reactions such as biomass decay (``r10``).

Its ``inhibition_map`` includes, alongside the product-inhibition
coefficients, the thermodynamic reverse-reaction (product) terms of the two
reversible steps (``k_6r`` on ``r6``, ``k_17r`` on ``r17``). A strip row
therefore reports everything the product does to that step: ADH's "fraction
lost to ethanol" counts both the inhibition of the forward rate and the
ethanol-driven reverse flux.

The Ehrlich branch ends in two steps. ``r16`` is the Aro10 decarboxylase --
irreversible Michaelis--Menten in KIV, releasing isobutyraldehyde and CO2,
with no product or cross-product term of any kind, so it carries no strip at
all. ``r17`` is the Adh6 reductase, a structural mirror of ``r6``, and it
carries the ethanol and acetate strips together with its own isobutanol
terms ``K_17e``/``k_17r``. The lumped step's ``K_16i``/``k_16r`` and
``k_16ia``/``k_16ie`` stay declared at 0 for the isobutanol demonstration's
parameter workbooks that set them, but do nothing; the live knobs with those meanings are ``K_17e``
and ``k_17ia``/``k_17ie``.

:func:`~nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo.flux_map_spec.draw_scenario_flux_map`
changes only the Ehrlich rate constants between its two panels, never the
fed-batch feeding strategy, so panel **b** is scenario A's process with the
Ehrlich branch switched on -- not the full scenario-B configuration of the
:doc:`isobutanol biorefinery demonstration </case_studies/isobutanol/index>`.

.. currentmodule:: nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo.flux_map_spec

.. autodata:: FLUX_MAP_SPEC
   :annotation:

.. autofunction:: draw_scenario_flux_map

.. currentmodule:: nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo.scenarios

.. autofunction:: apply_scenario_A

.. autofunction:: apply_scenario_B

Parameter categories
~~~~~~~~~~~~~~~~~~~~

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
the Ehrlich branch (r13--r17).

Operation parameters -- fed-batch feeding, aeration staging, dilution, and
the vessel's oxygen-transfer capacity -- are not kinetics and are skipped by
the helpers. The last of these is the pair ``kLa`` (volumetric O2
mass-transfer coefficient, default 200 /h) and ``C_O2_sat`` (dissolved-O2
saturation, default 0.20 mmol/L): the model caps the realized respiratory
and growth-associated O2 uptake of its aerobic fluxes -- the quantity it
reports as ``qO2_TCA_growth_only``, times biomass -- at ``OTR_max =
kLa*C_O2_sat`` (40 mmol O2/L/h at the defaults) through the
assignment-rule factor ``f_O2``, which scales r2, r5, r8 and the aerobic
share of r7 so that uptake equals ``OTR_max`` exactly while the bound
binds, and reads 1 wherever it does not (and 0 while ``is_aerobic`` is
off). Only the
growth-associated O2 of the ungated ``anaerobic_growth_mult`` share of r7
(0.75 by default) is committed and never throttled. The cap does *not*
cover the model's total ``qO2`` -- the quantity ``AerationSpec`` sizes the
air stream from -- which additionally carries the aerobic NADH
re-oxidation of r1 and r4, net of the NAD(P)H that r6, r14 and r17 consume. That
term is not throttled, so ``qO2*x`` can exceed ``OTR_max`` (by ~20 % in a
bound 100 g/L batch); extending the cap to total ``qO2`` is a follow-up.
The shipped fed-batch process (aeration off at 5 g/L biomass, total
``qO2*x`` peaking at ~20 mmol/L/h, with the ``OUR_demand`` that the cap
compares against lower still while aerated) is never bound; a
glucose-lean, fully aerobic batch that grows past ~10 g/L biomass is,
while glucose-rich runs stay under the cap because r2, r5 and r8 are
glucose-inhibited. Set either
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

Pd/C oxidation of 3-HP to malonic acid
--------------------------------------

.. currentmodule:: nskinetics.models.pd_c_3hp_oxidation

``nskinetics.models.pd_c_3hp_oxidation`` holds four kinetic models of the
selective oxidation of 3-hydroxypropionic acid (3-HP) to malonic acid (MA)
over 10 % Pd/C, after S. Y. Kim, J. E. Chavarrio, J. Install and G. W. Huber,
*Appl. Catal. B: Environ. Energy* 386 (2026) 126403, together with the
authors' deposited data, curves digitized from the published figures and the
functions that redraw the paper's Figs. 9, 10 and 11 and SI Fig. S12.
:doc:`/case_studies/pd_c_3hp_oxidation` presents the models, the reproduced
figures, how to change the operating conditions, and what does not match the
paper; ``NOTES.md`` in the package folder is the full record.

Unlike the *S. cerevisiae* package, importing this subpackage loads no model
and runs no simulation: each loader builds a new
:class:`~nskinetics.KineticModel` when called (time in min, concentrations
in M). Temperatures are in K where an argument is named ``T`` and in
degrees Celsius where its name ends in ``_C``. It needs no optional
dependencies. Every name documented in this section except the figure
functions is re-exported from the subpackage, so
``from nskinetics.models.pd_c_3hp_oxidation import load_o2_model`` works; the
figure module is imported explicitly
(``from nskinetics.models.pd_c_3hp_oxidation import figures``).

Loaders
~~~~~~~

Each of the four models is an Antimony file with a generated SBML file next
to it; the module attributes ``h2o2_antimony_filename``,
``o2_antimony_filename``, ``o2_reversible_antimony_filename`` and
``o2_simplified_reversible_antimony_filename`` hold the file names. What a
loader selects (temperature, constant set, initial condition) is part of the
model it builds and survives ``reset()`` and ``compile_events()``; each
docstring names the parameters meant to be changed with ``set_value``.

.. currentmodule:: nskinetics.models.pd_c_3hp_oxidation.pd_c_3hp_oxidation

.. autofunction:: load_h2o2_model

.. autofunction:: load_o2_model

.. autofunction:: load_o2_reversible_model

.. autofunction:: load_o2_simplified_reversible_model

.. autofunction:: simulate_at_times

Species, plotted quantities and constant sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: SPECIES

.. autodata:: NORMALIZED

.. autodata:: H2O2_NODE_PARAMETER_SETS
   :annotation:

.. autodata:: O2_NODE_PARAMETER_SETS
   :annotation:

.. autodata:: O2_REVERSIBLE_PARAMETER_SETS
   :annotation:

.. autodata:: SIMPLIFIED_REVERSIBLE_PARAMETER_SETS
   :annotation:

.. autodata:: SIMPLIFIED_REVERSIBLE_INITIAL_CONDITIONS

.. autodata:: TABLE_S4
   :annotation:

Reference data
~~~~~~~~~~~~~~

.. currentmodule:: nskinetics.models.pd_c_3hp_oxidation.reference

The readers return the reference data shipped in ``reference_data/`` as
tidy tables whose columns carry the models' selection names. The deposited
files (Dryad, https://doi.org/10.5061/dryad.jsxksn0q4, CC0) are verbatim
copies; the readers undo their known quirks. Curves marked *digitized* were
read from the published figures, and ``reference_data/README.md`` gives the
provenance of every file. ``reference_data_dir`` is the path of the folder.

.. autofunction:: read_fig9a_h2o2

.. autofunction:: read_fig9b_o2

.. autofunction:: read_figS12a_h2o2_inset_digitized

.. autofunction:: read_figS12b_opa_inset_digitized

.. autofunction:: read_fig10_o2

.. autofunction:: read_fig10d_arrhenius

.. autofunction:: read_fig11a_selectivity_digitized

.. autofunction:: read_fig11b_yield_digitized

Heat maps (Fig. 11)
~~~~~~~~~~~~~~~~~~~

.. currentmodule:: nskinetics.models.pd_c_3hp_oxidation.heat_maps

Temperature--time maps of the two quantities of Fig. 11, MA "selectivity"
(``X_MA``) and MA "yield" (``(1 - X_HP3)*X_MA``), and the measures that
compare a computed map with the contour lines digitized from the published
figure.

.. autodata:: FIG11_T_C
   :annotation:

.. autodata:: FIG11_TIMES
   :annotation:

.. autodata:: FIG11_TIME_OFFSET

.. autodata:: FIG11_LEVELS

.. autodata:: FIG11_MIN_PER_PX

.. autodata:: FIG11_C_PER_PX

.. autofunction:: simulate_heat_maps

.. autofunction:: fig11_panel_heat_map

.. autofunction:: fig11_heat_maps

.. autofunction:: table3_o2_heat_maps

.. autofunction:: heat_map_maximum

.. autofunction:: heat_map_ridge

.. autofunction:: contour_offsets

.. autofunction:: contour_shifts

Temperature series and Arrhenius plot (Fig. 10)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: nskinetics.models.pd_c_3hp_oxidation.temperature_series

The measures that compare a model with the experimental points of
Fig. 10a--c, and the points and straight lines of the Arrhenius plot
(Fig. 10d).

.. autodata:: FIG10_TEMPERATURES_C

.. autodata:: FIG10D_REGRESSION_TEMPERATURES_C

.. autodata:: FIG10D_LINE_ENDS_C

.. autodata:: TEXT_ACTIVATION_ENERGIES

.. autofunction:: pooled_r2

.. autofunction:: rms_difference

.. autofunction:: fig10d_points

.. autofunction:: fig10d_regression_lines

.. autofunction:: fig10d_text_lines

Emulation of the authors' solver run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: nskinetics.models.pd_c_3hp_oxidation.solver_emulation

The model curves the authors deposited for Fig. 9 and Fig. 10a–c are solver
steps that lie up to 2.7e-3 from the converged solutions; they are
reproduced step for step by the method of MATLAB ``ode23s`` at default
tolerances (emulated; MATLAB was not run). ``solver_emulation`` repeats
that solver run with the rates of a shipped model, so that the models can
be compared with the deposited nodes to the digits the files print. It is
used for verification only (the package tests and
:doc:`/case_studies/pd_c_3hp_oxidation`); it is not a general-purpose
solver.

.. autofunction:: ode23s

.. autofunction:: emulate_authors_run

.. autofunction:: emulate_fig9a

.. autofunction:: emulate_fig9b

.. autofunction:: emulate_fig10abc

.. autofunction:: node_deviations

Figures
~~~~~~~

.. currentmodule:: nskinetics.models.pd_c_3hp_oxidation.figures

``figures`` regenerates Fig. 9, Fig. 10, SI Fig. S12 and Fig. 11 of the
paper from the models, and draws three comparison figures and a
verification figure. Each ``plot_*`` function returns a new matplotlib
figure; each ``draw_*`` function draws one panel on an existing axes. Run
the module as a script to write the figures to a folder (``--figures``,
``--formats`` and ``--dpi`` select the figures, the file formats and the
resolution of the PNG files)::

    python -m nskinetics.models.pd_c_3hp_oxidation.figures OUTPUT_DIR

.. autofunction:: save_figures

.. autofunction:: plot_fig9

.. autofunction:: plot_figS12

.. autofunction:: plot_fig10

.. autofunction:: plot_fig10_table3_o2

.. autofunction:: plot_fig11

.. autofunction:: plot_fig11_verification

.. autofunction:: plot_fig11_table3_o2

.. autofunction:: plot_solver_emulation

.. autofunction:: draw_fig9a

.. autofunction:: draw_fig9b

.. autofunction:: draw_figS12a

.. autofunction:: draw_figS12b

.. autofunction:: draw_fig10abc

.. autofunction:: draw_fig10d

.. autofunction:: draw_fig11a

.. autofunction:: draw_fig11b
