Changelog
=========

.. currentmodule:: nskinetics

Unreleased
----------

Behaviour changes
~~~~~~~~~~~~~~~~~

* **The default reset now restores the model's baseline.** This applies to
  a :class:`KineticModel` without a custom ``reset_func``, and to
  everything that resets it: :meth:`KineticModel.reset`,
  ``simulate(reset=True)``, :meth:`KineticModel.run` and
  :class:`~nskinetics.units.NSKBatchReactor`.

  - *What changed.* The reset now also restores compartments and
    event-assigned variables to their ``init(x)`` values: working volumes,
    :class:`FeedSpike` counters, stage flags, and the targets of events
    written in the SBML or Antimony source (RoadRunner's reset already
    restored species and rate-rule variables). They are set
    before RoadRunner's reset as well as after it, so initial
    concentrations are computed with the origin volumes, and a reset is
    bit-identical to loading the model afresh.
  - *Why.* Previously only RoadRunner's ``reset()`` ran. The second run of
    an event-driven model therefore started silently from the first run's
    final volume and spike counter. On a three-spike fed-batch toy model,
    the final product titer dropped from 161.9 g/L to 37.6 g/L, with no
    error.
  - *What stays the same.* Global parameter values you set yourself (for
    example ``km.set_value('k', 0.4)``) are still kept by a reset.
    For models with a ``reset_func``, :meth:`KineticModel.reset` calls it
    exactly as before (but see the shipped model's reset function below).
  - *Getting the old behaviour.* ``km.reset(mode='roadrunner')``.
  - *Migration.* A compartment or an event-assigned variable that you set
    by hand after ``compile_events()`` is now put back by the next reset.
    Pass it for one run in ``Condition.initial``, or keep it with
    ``reset(preserve=(...))``.

  :doc:`concepts` (section "Resets, baselines and reproducible runs") has
  the details.

* **Repeated runs of the shipped** *S. cerevisiae* **model no longer
  alternate.** Its reset function called ``RoadRunner.reset()`` while the
  working volume was still at the previous run's end value, so a run that
  did not set ``[s_glu]`` itself started at 15 g/L and 229.65 g/L in turn.
  It now sets the volume and its bookkeeping variables before
  ``RoadRunner.reset()`` as well as after it, and every run starts at
  15 g/L. A species whose starting value a run sets explicitly is
  unaffected; every other species now starts from the same value on every
  run, whatever the previous run's end volume.

* :class:`~nskinetics.units.NSKBatchReactor` keeps its spike cap across
  the resets of its spike-count retry and its spike-feed reconciliation
  passes (``reset(preserve=...)``). A baseline reset therefore cannot undo
  a reduced or frozen cap.

New
~~~

* :class:`Condition`, an immutable description of one run. It holds
  initial values, temporary parameter values, an output grid (``t_end`` and
  ``n_steps``, or explicit ``times``), reset options and a label. Derive
  variants with ``replace``, ``with_parameters`` and ``with_initial``.
  ``to_dict`` and ``from_dict`` give a JSON-safe form. The order of the
  ``initial`` keys does not matter: volumes are set before concentrations.
* :meth:`KineticModel.run` runs a condition: it validates every id first,
  restores the overridden parameters afterwards, and returns a
  :class:`SimulationResult`.
* :class:`SimulationResult` holds a trajectory with its units and
  condition. It offers ``time``, ``array``, ``final``, ``at`` (linear
  interpolation) and ``plot``. ``to_csv`` and ``from_csv`` round-trip
  exactly.
* :func:`run_many` returns a :class:`RunBatch`: one row of metrics per
  condition, with failed runs recorded rather than raised. It runs on a
  private copy of the model, serially, or on a process pool with
  ``workers > 1``.
* :func:`scan` returns a :class:`ScanResult`, with ``pivot`` and
  ``contourf``, for product or zipped grids of parameter and initial
  values. ``RunBatch.to_csv`` and ``RunBatch.from_csv`` save a batch or a
  scan as a long-format CSV. The metrics table, the conditions and the grid
  round-trip exactly; trajectories are not stored.
* :meth:`KineticModel.copy` returns an independent model rebuilt from the
  origin SBML. :meth:`KineticModel.overrides` is a context manager for
  temporary value changes.
* :attr:`KineticModel.parameters` and :attr:`KineticModel.species` are
  introspection tables: values, origin values, role, event targets and the
  reactions that use each id.
* ``KineticModel.reset(mode=..., preserve=...)`` and
  :attr:`KineticModel.baseline_variables`.
* ``KineticModel.simulate(times=...)`` returns rows at an explicit output
  grid.
* :meth:`KineticModel.from_antimony` loads Antimony text or a file.
  :meth:`KineticModel.from_sbml` now takes ``units`` and ``reset``.
* Top-level exports: ``nsk.Condition``, ``nsk.SimulationResult``,
  ``nsk.run_many``, ``nsk.scan``, ``nsk.RunBatch`` and ``nsk.ScanResult``.

Fixed
~~~~~

* :meth:`KineticModel.state_selections` now includes the targets of events
  written in the SBML or Antimony source, not only those of :class:`Event`
  objects. For such a model, :func:`compute_flux_summary` therefore also
  needs those targets in the trajectory: record
  ``km.state_selections()`` (an
  :class:`~nskinetics.units.NSKBatchReactor` already does). A trajectory
  without them now raises :class:`~nskinetics.exceptions.KineticSimulationError`
  instead of being replayed with stale values.
* The shipped *S. cerevisiae* model module no longer imports the unused
  ``simplesbml`` package. Without that package installed, the module failed
  to import.
* The :attr:`KineticModel.time_factor` docstring now says what the property
  returns: model time units per hour. Its behaviour is unchanged.
* The shipped model's reset function and
  :class:`~nskinetics.units.NSKBatchReactor`'s spike-cap writes use item
  access (``r['env'] = 1.``) instead of attribute access. libroadrunner
  removes a model's ``r.<id>`` properties whenever another model is loaded,
  after which an attribute write silently missed the model.
