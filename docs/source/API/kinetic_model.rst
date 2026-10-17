KineticModel
============

.. currentmodule:: nskinetics

:class:`KineticModel` is a wrapper around a Tellurium extended
RoadRunner model. It adds unit-aware value access and a Python event API.

* **Building.** :meth:`~nskinetics.KineticModel.from_antimony` and
  :meth:`~nskinetics.KineticModel.from_sbml` load a model;
  :meth:`~nskinetics.KineticModel.copy` rebuilds an independent one from the
  model's origin SBML (``RoadRunner.getSBML()``, read at the time of the
  copy).
* **Running.** :meth:`~nskinetics.KineticModel.run` simulates a
  :class:`Condition` and returns a :class:`SimulationResult`; this is the
  entry point for reproducible work.
  :meth:`~nskinetics.KineticModel.simulate` integrates the model directly
  (over ``[t0, t_end]`` or an explicit ``times`` grid). Both store the
  trajectory on the object (``km.results_df`` / ``.results_array`` /
  ``.results_dict`` / ``.results_col_names``).
* **Plotting.** :meth:`~nskinetics.KineticModel.plot_simulation_results`
  (aliases :meth:`~nskinetics.KineticModel.plot_time_course` and
  :meth:`~nskinetics.KineticModel.plot_trajectory`) plots the most recent
  run.
* **State.** :meth:`~nskinetics.KineticModel.reset` restores the model's
  baseline (or runs a custom reset function; see its ``mode`` and
  ``preserve`` arguments and :attr:`~nskinetics.KineticModel.baseline_variables`).
  :meth:`~nskinetics.KineticModel.overrides` changes values temporarily.
  :meth:`~nskinetics.KineticModel.copy` returns an independent model.
* **Introspection.** :attr:`~nskinetics.KineticModel.parameters` and
  :attr:`~nskinetics.KineticModel.species` describe the model as tables.

:doc:`../concepts` explains resets and baselines. The underlying RoadRunner
object remains directly accessible as ``km._te`` for lower-level use; the
:doc:`process bridge <units>` can also drive it.

.. autoclass:: KineticModel
   :members:
   :show-inheritance:

Conditions and results
----------------------

A :class:`Condition` is the immutable description of one run (initial values,
parameter overrides, output grid, reset behaviour);
:meth:`KineticModel.run` simulates it and returns a :class:`SimulationResult`.

.. autoclass:: Condition
   :members:

.. autoclass:: SimulationResult
   :members:

Batch runs and scans
--------------------

:func:`run_many` runs many conditions on a private copy of a model (serially,
or on a process pool with ``workers > 1``) and returns a :class:`RunBatch`;
:func:`scan` runs a condition over a grid of values and returns a
:class:`ScanResult`. Failed runs are recorded, not raised, by default.
Every condition must have ``reset=True`` (and not ``mode='roadrunner'`` in
its ``reset_kwargs``), so that results do not depend on run order or
worker; :func:`scan` puts each grid key in the conditions'
``initial`` (state) or ``parameters`` (settings), as described in
:doc:`../concepts`. :meth:`RunBatch.to_csv` saves a batch or a scan (metrics,
conditions and grid; not the trajectories) and :meth:`RunBatch.from_csv`
reads it back.

.. autofunction:: run_many

.. autofunction:: scan

.. autoclass:: RunBatch
   :members:

.. autoclass:: ScanResult
   :members:
   :show-inheritance:
