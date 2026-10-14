.. _isobutanol-running-analyses:

Running the analyses
====================

How to install the case study, check an installation and run its scripts.
The API is on :doc:`/API/biorefineries`.

.. admonition:: Rules for every simulating run
   :class: important

   - **One simulation process at a time.** Simulating processes share an
     on-disk numba cache, and two writers can corrupt it. If a run then
     fails with ``ReferenceError: ... underlying object has vanished`` or
     ``OSError: access violation``, delete the ``*.nbi``, ``*.nbc`` and
     ``*.pyc`` files under the installed packages and run again.
   - **One** ``load()`` **per process**, and one scenario per fresh process to
     reproduce a pin.

.. _isobutanol-installation:

Installation
------------

The case study needs NSKinetics 0.7.0 or later. Until 0.7.0 is on PyPI,
install from the ``dev`` branch with the ``isobutanol`` extra (or
``requirements_isobutanol.txt``), then check it with smoke case 1:

.. code-block:: bash

   git clone https://github.com/sarangbhagwat/nskinetics
   cd nskinetics
   git checkout dev
   pip install -e ".[isobutanol]"
   python -m nskinetics.biorefineries.isobutanol.smoke 1

The extra adds the Bioindustrial-Park ``biorefineries`` package, ``hensmith``,
``contourplots``, ``optuna``, ``SALib`` and ``openpyxl``. The pins were
computed with git clones of BioSTEAM/ThermoSTEAM (biosteam 2.54),
Bioindustrial-Park (``biorefineries`` 2.35.1) and hensmith (``95e9362``); PyPI
releases can move TEA results by a few percent.

.. _isobutanol-switch-build:

Re-running on a recorded build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The build tags on these pages are hensmith commits. To reproduce a number
recorded on another build, check that commit out in an editable hensmith
clone, run the analysis in a fresh process and switch back (the smoke pins
need ``95e9362``). The two 1-D slice scripts of :doc:`objective_landscape`
*require* ``2b5b27d``.

.. code-block:: bash

   git -C <hensmith clone> checkout 2b5b27d
   python <analysis script>
   git -C <hensmith clone> checkout -

Smoke cases and tests
---------------------

``python -m nskinetics.biorefineries.isobutanol.smoke N`` loads one pinned
case, simulates it three times and exits non-zero if an MPSP is more than 1 %
off its pin (about 2–4 min per case). Cases 1–2 are scenarios A and B, 3–8
the same in other separation configurations, 9–16 the ``opt_*`` scenarios.

.. code-block:: bash

   python -m pytest nskinetics/tests/isobutanol -m "not slow" -v                 # offline, ~1 min
   python -m pytest nskinetics/tests/isobutanol/test_smoke.py -m canonical -v   # cases 1 and 2, twice each

Never run the simulation tests under pytest-xdist.

Running the scripts
-------------------

The scripts are plain files under the package folder,
``nskinetics/biorefineries/isobutanol``: analyses in ``analyses/<topic>/`` and
figure scripts in ``plots/``. The commands on these pages run from that folder
with ``MPLBACKEND=Agg``; in PowerShell set variables with
``$env:NAME = "value"`` first. To find the folder:

.. code-block:: bash

   python -c "import importlib.util, pathlib; print(pathlib.Path(importlib.util.find_spec('nskinetics').origin).parent / 'biorefineries' / 'isobutanol')"

- **Results.** Every script writes under one results root,
  ``analyses/results`` by default or ``NSKINETICS_IBO_RESULTS_DIR``. No
  results ship with NSKinetics: a script whose input is missing raises a
  ``FileNotFoundError`` naming the analysis that produces it. On Windows, keep
  the root short (MAX_PATH).
- **Plotting scripts are sim-safe**: they read CSV and Excel outputs only and
  can run beside a simulation (except the first stage of
  ``plots/plot_feeding_strategy_baseline_A.py``).
- **Crash recovery.** ``analyses/supervise_sweep.py <script>`` relaunches a
  checkpointed sweep after a native integrator crash, and
  ``analyses/kinetic_optimization/optimize_kinetics_BO_supervised.py`` does
  the same for optimization campaigns (``--help`` lists the options).

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
     - Grid sweeps: redraw the plots from an earlier run's CSVs without
       simulating.
   * - ``IBO_SWEEP_FRESH=1``
     - Checkpointed sweeps: discard a leftover checkpoint.
   * - ``IBO_AXIS_SWEEP``, ``IBO_AXIS_N_POINTS``, ``IBO_AXIS_LOW``,
       ``IBO_AXIS_HIGH``
     - 1-D sweeps through #1912 (also ``--sweep``, ``--n-points``,
       ``--low``, ``--high``).
   * - ``IBO_SOBOL_MEASURE``, ``IBO_SOBOL_N_POINTS``, ``IBO_SOBOL_SEED``,
       ``IBO_SOBOL_BANDS``
     - Sobol' stage 1: sampling measure, design size, seed and band file.
   * - ``VORONOI_NODE``
     - Node.js executable for ``plots/plot_proteome_voronoi.py``.
   * - ``NSKINETICS_TEST_TIMEOUT_SCALE``
     - Multiplies the simulation-test timeouts (e.g. ``2`` after clearing the
       numba caches).
