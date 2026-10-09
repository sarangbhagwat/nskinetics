.. _isobutanol-running-analyses:

Running the analyses
====================

This page is the user guide to ``nskinetics.biorefineries.isobutanol``: how
to install it, build the biorefinery and load a scenario, check an
installation with the pinned smoke cases, run its test suite, and run,
supervise and re-plot every analysis behind the case-study pages. The API of the modules named
here is on :doc:`/API/biorefineries`.

.. admonition:: Three rules for every simulating run
   :class: important

   1. **One simulation process at a time.** BioSTEAM, ThermoSTEAM and
      NSKinetics keep compiled numba kernels in an on-disk cache; two
      simulating processes writing it at once can corrupt it. Run analyses
      strictly one after another, including smoke cases and test runs.
   2. **One** ``load()`` **per process.** :func:`isobutanol.load()
      <nskinetics.biorefineries.isobutanol.load>` builds the flowsheet once
      per Python process; a second call is a no-op that ignores a different
      separation configuration. Run a different configuration in a fresh
      process.
   3. **One scenario per fresh process to reproduce a pin.** Loading a
      scenario mutates global fermentation state, so a later scenario in the
      same process can land slightly differently from its pin.

   If a run crashes at import with ``ReferenceError: ... underlying object
   has vanished`` (from ``numba/core/caching.py``) or with ``OSError: access
   violation`` from ``llvmlite``, the numba cache is corrupted. Delete the
   ``*.nbi`` and ``*.nbc`` cache files (and stale ``*.pyc``) under the
   installed packages and run again; the first run then recompiles and is
   slower.


.. _isobutanol-installation:

Installation
------------

.. admonition:: Which NSKinetics version has the case study
   :class: important

   The case study needs the first NSKinetics release that ships
   ``nskinetics.biorefineries``. The releases published on PyPI so far do
   not include it: with them, ``pip install "nskinetics[isobutanol]"`` only
   warns that the ``isobutanol`` extra does not exist, and the smoke command
   below fails with ``ModuleNotFoundError``. Until that release is out,
   install from source. These pages are built from the ``dev`` branch,
   which carries the case study; the ``main`` branch, also numbered 0.6.0,
   does not.

**From source (works now).** Clone the repository, check out the ``dev``
branch, which carries the case study, and install it in editable mode with
the ``isobutanol`` extra:

.. code-block:: bash

   $ git clone https://github.com/sarangbhagwat/nskinetics
   $ cd nskinetics
   $ git checkout dev
   $ pip install -e ".[isobutanol]"

Instead of the extra you can install ``requirements_isobutanol.txt`` on top
of ``requirements_install.txt``.

**From PyPI (once that release is out).** ``pip install
"nskinetics[isobutanol]"`` will install the same package and optional
dependencies.

The ``isobutanol`` extra adds the dependencies that plain NSKinetics does not
need: the Bioindustrial-Park ``biorefineries`` package (the flowsheet reuses
its ``corn``, ``cellulosic`` and ``ethanol`` modules), ``hensmith``
(heat-exchanger network synthesis), ``contourplots``, ``optuna``, ``SALib``
and ``openpyxl``. Check the installation with smoke case 1 (about 2–4
minutes; see `Check an installation: the smoke cases`_):

.. code-block:: bash

   $ python -m nskinetics.biorefineries.isobutanol.smoke 1

.. note::

   **Validated stack.** The pinned results were computed with development
   versions installed from git rather than PyPI: BioSTEAM and ThermoSTEAM
   clones (biosteam 2.54), a Bioindustrial-Park clone (``biorefineries``
   2.35.1) and a hensmith clone (0.1.2, which requires
   ``biosteam>=2.54.0``). The PyPI releases resolve the same imports but
   can move TEA results by a few percent; heat-exchanger-network synthesis
   in particular shifts with the hensmith version. The *build tags* in
   :ref:`isobutanol-conventions` (``2b5b27d``, ``6d4776f``, ``95e9362``) are
   commits of the `hensmith repository
   <https://github.com/BioSTEAMDevelopmentGroup/hensmith>`__, and the
   validated hensmith clone is at ``95e9362``. To reproduce a number
   recorded on another build exactly, switch hensmith to that commit
   (`Re-running on a recorded build`_).

Importing the package is build-free: ``from nskinetics.biorefineries import
isobutanol`` imports its modules (and with them BioSTEAM and the optional
dependencies) but builds and simulates nothing. ``import nskinetics`` alone
never imports the case study.

.. _isobutanol-switch-build:

Re-running on a recorded build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every analysis runs on the validated stack (hensmith ``95e9362``, the build
of the ``opt_*`` pins), but a number recorded on another build then differs
by the TEA shift between builds, up to a few percent. To reproduce it
exactly, check out the build's tag in a hensmith clone that is installed in
editable mode (or placed first on ``PYTHONPATH``), run the analysis in a
fresh process, and restore the clone afterwards:

.. code-block:: bash

   $ git -C <hensmith clone> checkout 2b5b27d   # e.g. the campaigns' build
   $ python <analysis script>                   # fresh process
   $ git -C <hensmith clone> checkout -         # back to the previous commit

- **Required for the two 1-D slice scripts** of :doc:`objective_landscape`,
  ``evaluate_axis_try_informed_optimum.py`` and
  ``evaluate_k13_try_informed_optimum.py``. Each first re-simulates the
  flagship design #1912 and aborts with ``anchor PI off by … (> 0.02)``
  unless the live PI is within 2 % of the recorded 0.808. That holds only
  on 2b5b27d; on 95e9362 the same design reads PI 0.8559 (+5.9 %).
- **Optional elsewhere.** The campaigns, the TEA breakdowns, the sweeps and
  the Sobol' runs complete on any build; only their TEA numbers shift.
- **Restore before running the smoke cases.** The ``opt_*`` pins were
  computed on 95e9362, and on another build the TEA shift can exceed the
  smoke cases' 1 % tolerance.


Build, load a scenario, solve the TEA
-------------------------------------

.. code-block:: python

   from nskinetics.biorefineries import isobutanol

   isobutanol.load()                                # build + baseline simulation
   bundle = isobutanol.load_scenario('A')           # workbook kinetics + feeding
   print(bundle['solve_TEA']())                     # {'IRR': ..., 'MPSPs': {...}}

:func:`~nskinetics.biorefineries.isobutanol.load` builds the corn
dry-grind flowsheet with the NSKinetics fermentation section (``V406``),
simulates it to its baseline state and builds the BioSTEAM uncertainty
``Model``. It publishes the built objects into the package namespace:
``isobutanol.V406``, ``isobutanol.fbs_spec`` (the fed-batch strategy
specification), ``isobutanol.corn_EtOH_IBO_sys`` and its TEA,
``isobutanol.unit_groups_dict``, ``isobutanol.sep_udct`` (the separation
units by their factory IDs), ``isobutanol.model`` and more. Its
``separation_processes`` argument chooses the separation trains to build,
a non-empty subset of ``('IBO_EtOH', 'ethanol')`` (default: both, with all
broth gated to the ethanol + isobutanol train). Expect a load to take
tens of seconds.

**The scenario registry.** Every scenario is a
:class:`~nskinetics.biorefineries.isobutanol.scenarios.ScenarioSpec` in
:data:`~nskinetics.biorefineries.isobutanol.scenarios.SCENARIOS`: a
parameter-distribution workbook (baseline kinetics plus uncertainty
distributions, shipped under the package's
``data/parameter_distributions``), a baseline feeding strategy, an
enzyme-burden policy and the pinned MPSPs. The ten registered scenarios are
listed on :doc:`index`.
:func:`~nskinetics.biorefineries.isobutanol.scenarios.load_scenario` applies
one to the built model: the workbook's baseline column is applied by one
baseline ``Model`` evaluation, and one baseline ``model_specification``
follows. The first call in a process with the enzyme burden on also
evaluates scenario A once, to calibrate the burden, so a first
``load_scenario('A')`` simulates the flowsheet three times. It returns a
*bundle* dict with ``model_specification`` (the entry point to simulate at a
feeding strategy), ``solve_TEA``, ``V406``, ``fbs_spec``, ``tea``,
``model``, ``namespace_dict`` (the namespace the workbook's load
statements are executed in), ``feeding_kwargs``, ``expected`` (the pins),
``burden_model``, ``burden_on`` and ``spec``.

**Enzyme burden per scenario.** ``load_scenario`` installs the
enzyme-burden constraint (:doc:`enzyme_burden`) according to the scenario's
``burden_default``: on for ``A`` (where it is inert) and for the eight
``opt_*`` scenarios, off for ``B`` (whose pathway constants exceed the
proteome cap). Pass ``burden=True`` or ``burden=False`` to override it. A
bare :func:`~nskinetics.biorefineries.isobutanol.load` installs no burden.

**Simulate, then solve.** Every change goes through
``model_specification`` (which loads the feeding specifications and
simulates the flowsheet to convergence, with recovery from solver
failures) followed by ``solve_TEA`` (which reads the current state without
re-simulating):

.. code-block:: python

   bundle['model_specification'](threshold_conc=226.1, target_conc=230.,
                                 max_n_spikes=0)    # a plain batch at 230 g/L
   print(bundle['solve_TEA']())

``solve_TEA`` returns the purity-adjusted MPSP of each product at a 15 %
IRR (NaN for a product that is not made) and the IRR at market prices
(−inf when the plant loses money outright).


Check an installation: the smoke cases
--------------------------------------

The module :mod:`nskinetics.biorefineries.isobutanol.smoke` holds 16 pinned
cases. Each builds the biorefinery in one separation configuration, applies
a scenario and simulates it three times, checking every simulation: MPSPs
within 1 % of their pins, within 0.5 % of the first simulation, a closed
spike-feed balance in ``V406`` and, for the ``opt_*`` scenarios, the
campaign objective reproduced within 2 %. Run one case per process:

.. code-block:: bash

   $ python -m nskinetics.biorefineries.isobutanol.smoke 1
   $ python -m nskinetics.biorefineries.isobutanol.smoke 10 --json

It prints one line per simulation (``sim1: IRR=... MPSPs={...}``) and the
elapsed time, and exits non-zero on any failed check; ``--json`` adds a
machine-readable final line and ``--n-sims K`` changes the number of
simulations. A case takes about 2–4 minutes.

.. list-table:: Smoke cases
   :header-rows: 1
   :widths: 10 26 34 30

   * - Case
     - Scenario
     - Separation trains built
     - Pins
   * - 1, 2
     - A, B
     - both (all broth to the ethanol + isobutanol train)
     - registry
   * - 3, 4
     - A, B
     - ethanol + isobutanol train only (reproduces 1, 2)
     - registry
   * - 5, 6
     - A, B
     - ethanol-primary train only
     - own (``ETHANOL_PRIMARY_PINS``)
   * - 7, 8
     - A, B
     - both, all broth re-gated to the ethanol-primary train
     - own (track 5, 6)
   * - 9–16
     - the eight ``opt_*`` scenarios
     - ethanol + isobutanol train only
     - registry + campaign objective


Testing the case study
----------------------

The case study's test suite is in the repository, under
``nskinetics/tests/isobutanol``: 13 pytest modules, a ``fixtures/`` folder
with a small golden trajectory CSV and a ``scripts/`` folder of
script-style checks that the tests run in fresh processes. Run it from the
root of a source install (see `Installation`_). Its markers (``isobutanol``,
``slow``, ``simulation``, ``smoke``, ``canonical``) are registered in the
repository's ``conftest.py``.

.. code-block:: bash

   # offline tests: no biorefinery build, seconds to a minute
   $ python -m pytest nskinetics/tests/isobutanol -m "not slow" -v

   # slow offline checks: kinetic-optimization checks, surrogate fits
   $ python -m pytest nskinetics/tests/isobutanol -m "slow and not simulation" -v

   # simulation tests: flowsheet regressions and the smoke cases (an hour or more)
   $ python -m pytest nskinetics/tests/isobutanol -m simulation -v

   # the canonical gate: smoke cases 1 and 2, each run twice
   $ python -m pytest nskinetics/tests/isobutanol/test_smoke.py -m canonical -v

- **Offline tests** load the stdlib-only and sim-safe modules (``paths``,
  the enzyme burden, the Sobol' estimators, the sweep helpers, the plotting
  libraries and figure scripts) by file path and never build the
  biorefinery. A few skip when an optional module such as ``contourplots``
  or ``openpyxl`` is missing.
- **Simulation tests** (marker ``simulation``, always also ``slow``) run
  each check in a fresh Python process that builds and simulates the
  biorefinery: the flowsheet regressions of ``scripts/`` through
  ``test_simulations.py``, and the smoke cases through ``test_smoke.py``
  (marker ``smoke``). They need the ``isobutanol`` extra and are skipped
  when its simulation dependencies are missing.
- **The canonical gate** (marker ``canonical``) runs smoke cases 1 and 2
  twice each, in fresh processes, and requires the two runs' MPSPs to agree
  to three significant figures.

The one-process rule holds for tests too: run simulation tests strictly one
at a time and never under pytest-xdist (``-n``), which the suite refuses by
skipping them in an xdist worker. Each simulation subprocess also waits a
few seconds after it exits so that numba's cache writes settle. The first
run after a numba-cache clear recompiles the kernels and can hit a test
timeout; set ``NSKINETICS_TEST_TIMEOUT_SCALE=2`` to allow more time.


Where results go
----------------

Every analysis writes under one flat **results root**:
``analyses/results`` inside the installed package by default, or the folder
named by the ``NSKINETICS_IBO_RESULTS_DIR`` environment variable. Child
processes inherit the variable, so supervisors and the scripts they launch
agree on it. Publication figures go to ``<results root>/publication/<theme>``
(``Feed-strat``, ``Kinetic-sweeps``, ``Objective-landscape``,
``Optimization-figures``, ``Uncertainty``, ``Schematics``).

**No analysis results are distributed.** The PyPI package and the git
repository ship the scripts, the parameter workbooks and the Sobol' band
file, but no outputs. The figures on
these pages come from the authors' own runs of these scripts, kept locally
under ``analyses/results/publication`` and excluded from version control. A
script whose input is another analysis's output stops at once with a
``FileNotFoundError`` that names the analysis producing it; run that
analysis first, or point ``NSKINETICS_IBO_RESULTS_DIR`` at a results folder
that already holds its output.

Caches count as outputs too: for example,
``plots/plot_feeding_strategy_baseline_A.py`` keeps its simulated time course
as ``feeding_strategy_baseline_A_trajectory.csv`` and ``.json`` directly
under the results root.

.. tip::

   Keep the results root short on Windows. Kinetic-optimization file names
   are long, and Windows paths are limited to 260 characters; the default
   location in a typical checkout leaves only a few characters to spare.
   For a regular (non-editable) install, set ``NSKINETICS_IBO_RESULTS_DIR``
   to a short, writable folder in any case: the default then sits inside
   ``site-packages``, which may not be writable, is replaced on upgrade and,
   in a conda environment, makes the path longer still.


Running the scripts
-------------------

The analysis scripts are plain files, not importable modules. They live in
thematic folders under the package's ``analyses/`` directory; the figure
scripts (``plot_*.py`` and the schematics) live in ``plots/`` and are run
as files too. The two plotting libraries in that folder,
``plots.spearman_matrix`` and ``plots.trajectory_contourplot``, are the
exception: they are importable modules (see
:ref:`the API page <ibo-api-plots>`). Every path on these pages is relative
to the package folder, which you can find without importing anything:

.. code-block:: python

   import importlib.util, pathlib
   spec = importlib.util.find_spec('nskinetics')
   print(pathlib.Path(spec.origin).parent / 'biorefineries' / 'isobutanol')

Run a script from the ``analyses/`` folder (or give its full path), with
``MPLBACKEND=Agg`` for unattended runs:

.. tab-set::

   .. tab-item:: Bash

      .. code-block:: bash

         cd <package folder>/analyses
         export MPLBACKEND=Agg
         IBO_SWEEP_SCENARIO=A IBO_SWEEP_STEPS=40 \
             python feeding_strategies/evaluate_feeding_strategies.py

   .. tab-item:: PowerShell

      .. code-block:: powershell

         cd <package folder>\analyses
         $env:MPLBACKEND = "Agg"
         $env:IBO_SWEEP_SCENARIO = "A"; $env:IBO_SWEEP_STEPS = "40"
         python feeding_strategies\evaluate_feeding_strategies.py

**Plotting scripts are sim-safe.** Except for the first stage of
``plots/plot_feeding_strategy_baseline_A.py``, they read CSV and Excel
outputs only, load their helpers by file path and never import BioSTEAM or
build the biorefinery, so they can run while a simulation is in progress.


Supervising long runs
---------------------

A native integrator crash (a CVODE segfault: exit code 5, no traceback) can
kill a long run outright. Two stdlib-only supervisors relaunch it in a fresh
process, strictly one process at a time:

- ``analyses/supervise_sweep.py`` runs a **checkpointed** sweep. The sweep
  logs every evaluated point to ``<results root>/<script stem>_checkpoint.csv``
  and writes ``<script stem>_inflight.json`` just before each simulation; a
  relaunch logs the point the previous process died in as ``LOST`` and
  resumes after it. The supervisor kills a child whose in-flight sidecar is
  older than ``--stall-timeout-min`` (default 10) and gives up after
  ``--max-empty-attempts`` (default 3) attempts that add no row or after
  ``--max-attempts`` (default 25). ``--results-dir`` overrides the results
  root it polls.

  .. code-block:: bash

     python supervise_sweep.py \
         kinetic_sweeps/evaluate_EtOH_glycolysis_inhib_ethanol.py
     python supervise_sweep.py sensitivity/evaluate_sobol_split12d.py \
         --stem <STUDY_NAME> --checkpoint-suffix _trajectory.csv

  A script configured by environment variables keeps them under the
  supervisor, e.g. ``IBO_AXIS_SWEEP=screen python supervise_sweep.py
  objective_landscape/evaluate_axis_try_informed_optimum.py --stem
  evaluate_axis_try_informed_optimum_screen``.
- ``analyses/kinetic_optimization/optimize_kinetics_BO_supervised.py`` is
  the launcher for Bayesian-optimization campaigns. It relaunches after a
  crash, kills a child that logs no trial for ``--stall-timeout-min``
  minutes (default 3), logs the killed trial as a ``LOST`` row, and resumes
  the study (``--help`` lists its options).

Sweeps that are not checkpointed (``evaluate_feeding_strategies.py``,
``evaluate_EtOH_k1e_inhib_ethanol.py``, ``evaluate_EtOH_k13_k7ii.py``,
``evaluate_EtOH_k1e_k1ie.py``, ``evaluate_EtOH_k1e_k7ie.py``) cannot be
supervised: a crash restarts them from the beginning.

.. list-table:: Environment variables
   :header-rows: 1
   :widths: 34 66

   * - Variable
     - Effect
   * - ``NSKINETICS_IBO_RESULTS_DIR``
     - Results root (default ``<package>/analyses/results``).
   * - ``IBO_SWEEP_SCENARIO``, ``IBO_SWEEP_STEPS``
     - Feeding sweep: scenario ``A`` or ``B`` (default B) and grid size
       ``<n>`` or ``<nx>x<ny>`` (default 25 × 25).
   * - ``IBO_SWEEP_REPLOT_FROM_CSV=1``
     - Grid sweeps (all but ``evaluate_EtOH_k1e_k1ie.py``): skip every
       simulation and redraw the plots from the CSVs of an earlier run.
   * - ``IBO_SWEEP_FRESH=1``
     - Checkpointed sweeps: discard a leftover checkpoint and start over.
   * - ``IBO_AXIS_SWEEP``, ``IBO_AXIS_N_POINTS``, ``IBO_AXIS_LOW``,
       ``IBO_AXIS_HIGH``
     - 1-D sweeps through the TRY-informed optimum: mode (``screen``, the
       default, or a variable name), points per variable (default 30 in
       screen mode, 160 for a single variable) and the bounds of a
       single-variable sweep. Run directly, the script also takes them as
       ``--sweep``, ``--n-points``, ``--low`` and ``--high``.
   * - ``IBO_SOBOL_MEASURE``, ``IBO_SOBOL_N_POINTS``, ``IBO_SOBOL_SEED``,
       ``IBO_SOBOL_BANDS``
     - Sobol' stage 1: sampling measure (``campaign``, ``screening`` or
       ``custom``), design size, seed and the band file of the custom
       measure (a bare file name such as ``linear_0-4x_2026-09-22.json`` is
       looked up in the package's ``data/sobol_bands``).
   * - ``VORONOI_NODE``
     - ``plots/plot_proteome_voronoi.py``: path to the Node.js executable
       that renders the treemap (``--node`` overrides it).
   * - ``NSKINETICS_TEST_TIMEOUT_SCALE``
     - Case-study tests: multiplies every simulation-subprocess timeout,
       e.g. ``2`` after clearing the numba caches (see
       `Testing the case study`_).


Catalog of analyses
-------------------

Runtime classes: **sim-free** scripts only read outputs (seconds to
minutes); **minutes** means one load plus a handful of simulations; **hours**
and **a day or more** are long, opt-in runs. A single flowsheet simulation
takes roughly 15–30 s. Outputs are under the results root.

Feeding strategies
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 36 32 12 20

   * - Script
     - What it does
     - Runtime
     - Outputs · page
   * - ``analyses/feeding_strategies/evaluate_feeding_strategies.py``
     - Threshold × target grid with the spike cap optimized for MPSP at
       every valid point (19 metrics). Set ``IBO_SWEEP_SCENARIO=A`` and
       ``IBO_SWEEP_STEPS=40`` for the published grid. Not checkpointed.
     - a day or more (at least ~3,400 converged simulations, so at least
       ~15–30 h)
     - one CSV per metric · :doc:`feeding_strategies`
   * - ``plots/plot_feeding_strategy_baseline_A.py``
     - MESP map + baseline time course. Stage 1 simulates scenario A once
       when its trajectory cache is missing (or with ``--resimulate``);
       stage 2 is sim-free.
     - minutes, then seconds
     - ``publication/Feed-strat`` · :doc:`feeding_strategies`

Kinetic sweeps
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 36 32 12 20

   * - Script
     - What it does
     - Runtime
     - Outputs · page
   * - ``analyses/kinetic_sweeps/evaluate_EtOH_glycolysis_inhib_ethanol.py``
     - Glycolysis × ethanol-inhibition plane, scenario A, 40 × 40,
       checkpointed.
     - hours (~7–11 h)
     - one CSV per metric · :doc:`kinetic_sweeps`
   * - ``.../evaluate_EtOH_k1e_inhib_ethanol.py``
     - ``k_1e`` × ethanol inhibition, 40 × 40, not checkpointed.
     - hours (~7–11 h)
     - CSVs · :doc:`kinetic_sweeps`
   * - ``.../evaluate_EtOH_k3_inhib_ethanol.py``, ``.../evaluate_EtOH_k6_inhib_ethanol.py``
     - Pdc and Adh1 × ethanol inhibition, 40 × 40, checkpointed.
     - hours (~7–11 h each)
     - CSVs · :doc:`kinetic_sweeps`
   * - ``.../evaluate_EtOH_k13_inhib_isobutanol.py``
     - ALS × isobutanol inhibition around ``opt_PI_TRY_informed``,
       40 × 40, checkpointed.
     - hours (~7–11 h)
     - CSVs · :doc:`kinetic_sweeps`
   * - ``.../evaluate_EtOH_k13_k7ii.py``, ``.../evaluate_EtOH_k1e_k1ie.py``,
       ``.../evaluate_EtOH_k1e_k7ie.py``
     - Exploratory 20 × 20 or 25 × 25 planes on absolute parameter values;
       not checkpointed. The two ``k_1e`` planes run without the enzyme
       burden; the ``k_13`` × ``k_7ii`` plane runs with it (its anchor's
       default). No results are presented in these pages.
     - hours (~2–4 h)
     - CSVs
   * - ``.../evaluate_OAT_kinetics_B.py``
     - One-at-a-time sweep of every kinetic parameter under scenario B
       (0–10× baseline); resumable.
     - a day or more (> 1,000 simulations)
     - ``OAT_kinetics_B/``
   * - ``.../plot_OAT_kinetics_B.py``
     - Per-parameter and overview figures of the OAT sweep.
     - sim-free
     - ``OAT_kinetics_B/``
   * - ``.../probe_glycolysis_cliff_mechanism.py``
     - Counterfactual re-simulations across the diagonal cliff of the
       glycolysis plane, with carbon routing from cumulative fluxes.
     - minutes (~20 simulations)
     - CSV
   * - ``plots/plot_glycolysis_inhib_ethanol_A.py``, ``plots/plot_k1e_inhib_ethanol_A.py``,
       ``plots/plot_k13_inhib_isobutanol_opt_PI_TRY_informed.py``,
       ``plots/plot_greedy_trajectories_*.py``
     - Sweep figures (MESP and IRR twins, greedy hill-climb overlays).
     - sim-free
     - ``publication/Kinetic-sweeps`` · :doc:`kinetic_sweeps`

Kinetic optimization and the ``opt_*`` scenarios
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 36 32 12 20

   * - Script
     - What it does
     - Runtime
     - Outputs · page
   * - ``analyses/kinetic_optimization/optimize_kinetics_BO_supervised.py``
     - Launches and supervises one Bayesian-optimization campaign
       (``--objective``, ``--study-type``, ``--method``, ``--seed``,
       ``--relay-from``, ...). Most launcher defaults differ from the
       publication campaigns' settings; :ref:`ibo-opt-reproduce` gives the
       exact flags. The campaigns ran on hensmith 2b5b27d; on another build
       they still run, but the TEA numbers shift.
     - hours (~8–17 h per 2,000-trial campaign)
     - ``<study>_trajectory.csv``, optuna store · :doc:`kinetic_optimization`
   * - ``.../optimize_kinetics_BO.py``
     - The campaign driver the supervisor launches; rerunning a study name
       resumes it. Importing it builds the biorefinery.
     - as above
     - as above
   * - ``.../reproduce_split12d_trial.py``
     - Re-simulates one recorded trial and compares it metric by metric.
       Writes nothing.
     - minutes
     - — · :doc:`kinetic_optimization`
   * - ``analyses/opt_scenarios/build_opt_split12d_workbooks.py <scenario>``
     - Regenerates one ``opt_*`` workbook and its registry pin from a
       recorded trial (overwrites the shipped workbook).
     - minutes
     - workbook + reproduction JSON · :ref:`ibo-opt-scenarios`
   * - ``analyses/opt_scenarios/verify_opt_scenarios.py <scenario>``
     - Checks that a scenario's baseline equals the reproduced campaign
       trial (``MATCH``). Read-only.
     - minutes
     - — · :ref:`ibo-opt-scenarios`
   * - ``plots/optimization_figure_ideas/fig_main_scouts.py``,
       ``fig_s2_fingerprints.py``, ``check_facts.py``
     - The two campaign figures, and a self-test of every number in them.
     - sim-free
     - ``publication/Optimization-figures`` · :doc:`kinetic_optimization`
   * - ``plots/plot_kin_opt_parameter_sets.py``, ``plots/plot_proteome_voronoi.py``,
       ``plots/plot_relay_campaign_progress.py``
     - Parameter-set and proteome-treemap figures of selected trials, and
       the progress of a relay campaign (``--study``). The treemap's second
       stage renders with Node.js 22: run ``npm --prefix plots/voronoi
       install`` once (``node_modules`` is not shipped), and point
       ``--node`` or ``VORONOI_NODE`` at the Node executable if it is not on
       ``PATH``; ``--no-render`` writes only the allocation JSON.
     - sim-free
     - ``publication/Optimization-figures``; results root

Objective landscape
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 36 32 12 20

   * - Script
     - What it does
     - Runtime
     - Outputs · page
   * - ``analyses/objective_landscape/evaluate_axis_try_informed_optimum.py``
     - 1-D sweeps through the TRY-informed optimum: a coarse screen of 11
       variables or one variable at full resolution (``IBO_AXIS_SWEEP``).
       Checkpointed. **Requires hensmith 2b5b27d**; on any other build the
       anchor check aborts with ``anchor PI off by … (> 0.02)`` (see
       `Re-running on a recorded build`_).
     - hours
     - ``evaluate_axis_try_informed_optimum_<name>.csv`` · :doc:`objective_landscape`
   * - ``.../evaluate_k13_try_informed_optimum.py``
     - 161-point ``k_13`` sweep through the same optimum. Checkpointed.
       **Requires hensmith 2b5b27d**, with the same anchor check.
     - minutes (~15–25 min)
     - CSV · :doc:`objective_landscape`
   * - ``plots/plot_objective_landscape_{fitness_distance,needle_plateau,embedding}.py``,
       ``plots/plot_axis_try_informed_sweep.py``, ``plots/plot_k13_try_informed_sweep.py``
     - Landscape statistics of the pooled campaign trials, and the 1-D
       slices.
     - sim-free (the t-SNE fit takes minutes once, then is cached)
     - ``publication/Objective-landscape`` · :doc:`objective_landscape`

Global sensitivity
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 36 32 12 20

   * - Script
     - What it does
     - Runtime
     - Outputs · page
   * - ``analyses/sensitivity/evaluate_sobol_split12d.py``
     - Stage 1: simulates 4,096 feasible Sobol' designs of the campaign
       space (``IBO_SOBOL_*``). Resumable and nested.
     - a day or more (~20 h)
     - ``<study>_trajectory.csv``, ``_design.json`` · :doc:`sensitivity`
   * - ``analyses/sensitivity/analyze_sobol_split12d.py``
     - Stage 2: surrogates, closed and Shapley indices, tables and figures
       (``--study-name``; ``--self-test`` needs no data).
     - sim-free (minutes; hours with GP surrogates)
     - ``<study>_sobol_*`` · :doc:`sensitivity`

TEA breakdowns
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 36 32 12 20

   * - Script
     - What it does
     - Runtime
     - Outputs · page
   * - ``analyses/tea_breakdowns/collect_tea_breakdowns_split12d.py``
     - Re-simulates scenario A and the eight campaign optima and records
       their unit-group cost and utility breakdowns, checking the closures.
       The published breakdowns are on hensmith 2b5b27d; on another build
       it still runs, but the TEA numbers shift (`Re-running on a recorded
       build`_).
     - minutes (~4 min)
     - ``tea_breakdowns_split12d_<stamp>.json`` · :doc:`process_and_tea`
   * - ``plots/plot_tea_breakdowns_split12d.py``
     - 3 × 3 breakdown figure from the newest JSON.
     - sim-free
     - ``publication/Optimization-figures`` · :doc:`process_and_tea`

Uncertainty
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 36 32 12 20

   * - Script
     - What it does
     - Runtime
     - Outputs · page
   * - ``analyses/uncertainty/uncertainties_IBO_EtOH.py``
     - 6,000-sample Latin-hypercube Monte Carlo of scenario A with Spearman
       analysis. Autosaves every 10 samples and resumes.
     - a day or more (~30 h)
     - ``_IBO_<stamp>_..._1_full_evaluation.xlsx`` · :doc:`uncertainty`
   * - ``analyses/uncertainty/plot_uncertainty_scatters.py``
     - Per-metric scatter grids against every uncertain input.
     - sim-free
     - next to the workbook
   * - ``plots/plot_uncertainty_MPSP_vs_TCI.py``,
       ``plots/plot_uncertainty_spearman_tornado.py``, ``plots/spearman_matrix.py``
     - Bivariate MESP figure, Spearman tornado (+ CSV), correlation matrix.
     - sim-free
     - ``publication/Uncertainty`` · :doc:`uncertainty`

Schematics
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 36 32 12 20

   * - Script
     - What it does
     - Runtime
     - Outputs · page
   * - ``plots/plot_block_flow_diagram.py``,
       ``plots/conceptual_feeding_aeration_sugar_prep.py``
     - The block flow diagram and the sugar-preparation, feeding and
       aeration concept figure. Matplotlib only; need no results.
     - seconds
     - ``publication/Schematics`` · :doc:`process_and_tea`
