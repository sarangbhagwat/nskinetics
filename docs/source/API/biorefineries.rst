Isobutanol biorefinery (case study)
===================================

``nskinetics.biorefineries.isobutanol`` is the corn-to-ethanol and
isobutanol biorefinery of the :doc:`/case_studies/isobutanol/index`: a
BioSTEAM flowsheet whose fermentor runs the shipped *S. cerevisiae* kinetic
model, with its scenarios, TEA and the engines behind the analyses. It
needs the optional ``isobutanol`` dependencies and an NSKinetics version
that ships ``nskinetics.biorefineries``; until such a release is on PyPI,
install it from source (see :ref:`isobutanol-installation`).

Importing the package is **build-free**: it imports the modules below but
builds and simulates nothing, and ``import nskinetics`` never imports it.
:func:`~nskinetics.biorefineries.isobutanol.load` builds the flowsheet once
per process. This page documents the curated public API. The analysis
scripts (``analyses/``) and the figure scripts (``plots/plot_*.py`` and the
schematics) are files, not modules, and are catalogued on
:doc:`/case_studies/isobutanol/running_analyses`; the two plotting libraries
in ``plots/`` are importable and documented under
`Plotting library (plots)`_.

.. contents:: On this page
   :local:
   :depth: 1
   :backlinks: none


Package
-------

.. py:module:: nskinetics.biorefineries.isobutanol

.. autofunction:: nskinetics.biorefineries.isobutanol.load

After :func:`load`, the package namespace holds every published build
object, among them ``corn_EtOH_IBO_sys`` and its TEA
``corn_EtOH_IBO_sys_tea``, the fermentor ``V406``, the fed-batch strategy
specification ``fbs_spec``, the flowsheet ``f``, the separation units by
factory ID ``sep_udct``, the cost and utility groups ``unit_groups_dict``,
the uncertainty ``model`` and its workbook namespace ``namespace_dict``.
:func:`~nskinetics.biorefineries.isobutanol.system.solve_TEA`,
:func:`~nskinetics.biorefineries.isobutanol.scenarios.load_scenario` and
:func:`~nskinetics.biorefineries.isobutanol.utils.generate_save_kinetic_parameter_distributions`
are re-exported at package level.


Scenarios
---------

.. automodule:: nskinetics.biorefineries.isobutanol.scenarios
   :no-members:

.. autodata:: nskinetics.biorefineries.isobutanol.scenarios.SCENARIOS
   :no-value:

.. autoclass:: nskinetics.biorefineries.isobutanol.scenarios.ScenarioSpec

.. autofunction:: nskinetics.biorefineries.isobutanol.scenarios.load_scenario

.. autofunction:: nskinetics.biorefineries.isobutanol.scenarios.assert_MPSPs_pinned

.. autofunction:: nskinetics.biorefineries.isobutanol.scenarios.assert_MPSPs_stable

.. autofunction:: nskinetics.biorefineries.isobutanol.scenarios.assert_spike_feed_residual

.. autofunction:: nskinetics.biorefineries.isobutanol.scenarios.assert_objective_reproduced


Smoke cases
-----------

.. automodule:: nskinetics.biorefineries.isobutanol.smoke
   :no-members:

.. autodata:: nskinetics.biorefineries.isobutanol.smoke.SMOKE_CASES
   :no-value:

.. autodata:: nskinetics.biorefineries.isobutanol.smoke.ETHANOL_PRIMARY_PINS

.. autoclass:: nskinetics.biorefineries.isobutanol.smoke.SmokeCase

.. autofunction:: nskinetics.biorefineries.isobutanol.smoke.get_case

.. autofunction:: nskinetics.biorefineries.isobutanol.smoke.load_simulate_baseline

.. autofunction:: nskinetics.biorefineries.isobutanol.smoke.main


Flowsheet, simulation and TEA (``system``)
------------------------------------------

.. automodule:: nskinetics.biorefineries.isobutanol.system
   :no-members:

Use the package-level :func:`nskinetics.biorefineries.isobutanol.load`
rather than ``system.load``: it also builds the uncertainty ``Model``.

.. autofunction:: nskinetics.biorefineries.isobutanol.system.model_specification

.. autofunction:: nskinetics.biorefineries.isobutanol.system.load_simulate

.. autofunction:: nskinetics.biorefineries.isobutanol.system.solve_TEA

.. autofunction:: nskinetics.biorefineries.isobutanol.system.solve_TEA_at_IRR

.. autofunction:: nskinetics.biorefineries.isobutanol.system.set_active_burden

.. autofunction:: nskinetics.biorefineries.isobutanol.system.get_active_burden

.. autofunction:: nskinetics.biorefineries.isobutanol.system.snapshot_flowsheet_state

.. autofunction:: nskinetics.biorefineries.isobutanol.system.restore_flowsheet_state

.. autofunction:: nskinetics.biorefineries.isobutanol.system.plot_kinetic_results

Feeding-strategy optimizers used by the analyses:

.. autofunction:: nskinetics.biorefineries.isobutanol.system.optimize_max_n_glu_spikes_for_MPSP

.. autofunction:: nskinetics.biorefineries.isobutanol.system.optimize_1D_feeding_strategy_for_MPSP

.. autofunction:: nskinetics.biorefineries.isobutanol.system.optimize_stage_1_time_and_max_n_glu_spikes_for_MPSP


Separation trains
-----------------

.. automodule:: nskinetics.biorefineries.isobutanol.separations
   :no-members:

.. autodata:: nskinetics.biorefineries.isobutanol.separations.create_separation_system
   :no-value:

.. autodata:: nskinetics.biorefineries.isobutanol.separations.create_IBO_EtOH_separation_system
   :no-value:

.. autodata:: nskinetics.biorefineries.isobutanol.separations.create_EtOH_primary_separation_system
   :no-value:

.. autofunction:: nskinetics.biorefineries.isobutanol.separations.create_separation_chemicals

.. autofunction:: nskinetics.biorefineries.isobutanol.separations.create_scenario_B_feed


Cost index and price year
-------------------------

.. automodule:: nskinetics.biorefineries.isobutanol.process_settings
   :no-members:

Capital is costed at ``CEPCI = 797.9`` (2023 annual average) and every price
is expressed in ``PRICE_YEAR = 2023`` dollars, indexed with the BLS
chemicals producer price index ``chem_index`` (series WPU06).

.. autofunction:: nskinetics.biorefineries.isobutanol.process_settings.load_process_settings

.. autofunction:: nskinetics.biorefineries.isobutanol.process_settings.price_index_ratio


Uncertainty model
-----------------

.. automodule:: nskinetics.biorefineries.isobutanol.models.models_EtOH_IBO_corn
   :no-members:

.. autofunction:: nskinetics.biorefineries.isobutanol.models.models_EtOH_IBO_corn.create_model


Workbook helpers (``utils``)
----------------------------

.. automodule:: nskinetics.biorefineries.isobutanol.utils
   :no-members:

.. autofunction:: nskinetics.biorefineries.isobutanol.utils.generate_save_kinetic_parameter_distributions


Filesystem locations (``paths``)
--------------------------------

.. automodule:: nskinetics.biorefineries.isobutanol.paths
   :no-members:

.. autodata:: nskinetics.biorefineries.isobutanol.paths.RESULTS_DIR_ENV

.. autofunction:: nskinetics.biorefineries.isobutanol.paths.results_dir

.. autofunction:: nskinetics.biorefineries.isobutanol.paths.publication_dir

.. autofunction:: nskinetics.biorefineries.isobutanol.paths.analyses_dir

.. autofunction:: nskinetics.biorefineries.isobutanol.paths.workbook_dir

.. autofunction:: nskinetics.biorefineries.isobutanol.paths.workbook_path

.. autofunction:: nskinetics.biorefineries.isobutanol.paths.require_input

.. autofunction:: nskinetics.biorefineries.isobutanol.paths.load_by_path


Enzyme burden
-------------

.. py:module:: nskinetics.biorefineries.isobutanol.enzyme_burden

The proteome-allocation constraint described on
:doc:`/case_studies/isobutanol/enzyme_burden`. The module is stdlib-only and
can also be loaded by file path. Its sector constants are
``PROTEIN_CONTENT = 0.49`` g protein/gDCW, ``HOUSEKEEPING_FRACTION = 0.50``,
``TRANSLATION_FRACTION_WT = 0.225`` and ``SIGMA_EFF = 0.50``, giving the
flexible sector ``F_FLEX = 0.245`` g/gDCW and the wild-type translation
sector ``PHI_T_WT = 0.11025`` g/gDCW.

.. autoclass:: nskinetics.biorefineries.isobutanol.enzyme_burden.BurdenModel
   :members:

.. autoclass:: nskinetics.biorefineries.isobutanol.enzyme_burden.BurdenResult
   :no-members:

.. autoexception:: nskinetics.biorefineries.isobutanol.enzyme_burden.BurdenInfeasibleError

.. autofunction:: nskinetics.biorefineries.isobutanol.enzyme_burden.ehrlich_unit_cost


Kinetic optimization
--------------------

.. py:module:: nskinetics.biorefineries.isobutanol.kinetic_optimization

The engine behind the :doc:`/case_studies/isobutanol/kinetic_optimization`
campaigns: Bayesian optimization (Optuna TPE or Gaussian process) and dual
annealing over kinetic parameters and the fed-batch feeding variables, under
the enzyme-burden and fed-batch volume constraints, with one trajectory-CSV
row per trial. The module has no top-level import of NSKinetics or BioSTEAM
(it can be loaded by file path); the engines themselves need a loaded
biorefinery. Production campaigns are launched with
``analyses/kinetic_optimization/optimize_kinetics_BO_supervised.py``. Only
the main entry points are listed here.

.. autodata:: nskinetics.biorefineries.isobutanol.kinetic_optimization.OBJECTIVE_REGISTRY
   :no-value:

Running a campaign
~~~~~~~~~~~~~~~~~~

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.run_kinetic_optimization

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.run_kinetic_dual_annealing

.. autoclass:: nskinetics.biorefineries.isobutanol.kinetic_optimization.AnnealingResult
   :no-members:

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.resolve_study_preset

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.default_study_name

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.get_handles

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.default_results_dir

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.parameter_distributions_workbook

Reproducing a recorded trial
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.reproduce_split12d_trial

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.split12d_trajectory_path

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.read_trajectory_row

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.reconstruct_trial_kinetics

.. autofunction:: nskinetics.biorefineries.isobutanol.kinetic_optimization.compare_tracked_metrics


Global sensitivity analysis
---------------------------

.. automodule:: nskinetics.biorefineries.isobutanol.sensitivity_analysis
   :no-members:

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.campaign_engine_kwargs

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.screening_search_space

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.custom_search_space

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.design_record

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.feasible_sobol_stream

.. autoclass:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.VectorizedFeasibility

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.sample_feasible

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.all_closed_indices

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.first_order

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.total_order

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.shapley_effects

.. autofunction:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.fit_surrogates

.. autoclass:: nskinetics.biorefineries.isobutanol.sensitivity_analysis.Surrogate
   :no-members:


.. _ibo-api-plots:

Plotting library (``plots``)
----------------------------

.. automodule:: nskinetics.biorefineries.isobutanol.plots
   :no-members:

Two modules of ``plots/`` are importable libraries; both need the optional
``contourplots`` package and run no simulation. Import them explicitly (or
through the lazy names ``plots.plot_spearman_matrix``,
``plots.plot_metric_with_trajectories`` and ``plots.JBEI_UCB_colormap``).
Every other file in the folder is a figure script, run as a file and listed
on :doc:`/case_studies/isobutanol/running_analyses`.

Spearman correlation matrix (``plots.spearman_matrix``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. py:module:: nskinetics.biorefineries.isobutanol.plots.spearman_matrix

Ellipse correlation matrix of the Spearman rank correlations in an
uncertainty-analysis workbook (``*_1_full_evaluation.xlsx``, sheets
``Spearman`` and ``Spearman p-values``). Run as a script
(``python plots/spearman_matrix.py [WORKBOOK] [--out-dir DIR]``), it draws
the newest scenario-A workbook in the results root.

.. autofunction:: nskinetics.biorefineries.isobutanol.plots.spearman_matrix.plot_spearman_matrix

Sweep contours with greedy trajectories (``plots.trajectory_contourplot``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: nskinetics.biorefineries.isobutanol.plots.trajectory_contourplot
   :no-members:

.. autofunction:: nskinetics.biorefineries.isobutanol.plots.trajectory_contourplot.plot_metric_with_trajectories

.. autofunction:: nskinetics.biorefineries.isobutanol.plots.trajectory_contourplot.JBEI_UCB_colormap
