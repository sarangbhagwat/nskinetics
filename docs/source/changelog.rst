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

* ``nskinetics.utils.hybrid_global_optimize`` now runs serially by
  default (``parallel='none'``; it was ``'threads'``). One RoadRunner
  must never be used from several threads at once, and its docstring now
  warns against threads for any objective that runs a model. Its
  ``success`` flag is no longer always ``True`` (see "Fixed").

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
* ``nskinetics.analysis``, a new subpackage that asks which parameters
  matter and which ones measurements can determine (:doc:`API/analysis`,
  :doc:`tutorial/07_sensitivity_identifiability`, and the
  :doc:`concepts` section "Sensitivities across events"). Import it as
  ``from nskinetics.analysis import ...``; it is not re-exported at the top
  level.

  - Metrics with an explicit basis: :class:`~nskinetics.analysis.Final`,
    :class:`~nskinetics.analysis.At`, :class:`~nskinetics.analysis.Max`,
    :class:`~nskinetics.analysis.Min`,
    :class:`~nskinetics.analysis.TimeAtMax`,
    :class:`~nskinetics.analysis.Integral`,
    :class:`~nskinetics.analysis.EventCount`,
    :class:`~nskinetics.analysis.Productivity`,
    :class:`~nskinetics.analysis.Yield` and
    :class:`~nskinetics.analysis.CallableMetric`. They also work as
    ``metrics=`` in :func:`run_many` and :func:`scan`, which then record the
    columns they read. A zero denominator gives ``NaN``.
  - :func:`~nskinetics.analysis.local_sensitivity`: sensitivities of output
    time courses and metrics by central finite differences in
    log-parameter space. Every perturbation is a run from the restored
    baseline on a private copy of the model, whose integrator is tightened
    to at most ``(1e-12, 1e-10)`` by default (tighter model values are
    kept), with ``atol`` scaled to the peak of outputs and metric columns
    that start at 0, so results do not depend on the concentration unit.
    When a perturbation changes how often an event fires, or moves a
    firing across an output time, the step is refined, and parameters
    whose event counters still differ are flagged
    (``event_count_changed``, ``reliable``). A parameter at 0 is
    differenced in linear space; without an ``abs_step`` its step is
    ``rel_step`` in its own units, unrelated to the model's scale, and a
    warning says so.
    :class:`~nskinetics.analysis.LocalSensitivityResult` offers ``rank``,
    ``time_course``, ``metric_sensitivities``, two plots and an exact CSV
    round trip.
  - :func:`~nskinetics.analysis.identifiability`: the Fisher information
    matrix of weighted log-sensitivities, its singular values, numerical
    rank and null directions, Cramer-Rao standard errors, a practical
    ``identifiable`` verdict per parameter, and the Brun et al. (2001)
    importance and collinearity measures with
    :meth:`~nskinetics.analysis.IdentifiabilityResult.best_subset`. The
    noise is absolute, or relative with a floor at ``relative_floor``
    (default 1e-3) of each output's peak, per output. A warning, recorded
    as ``rank_test_doubtful``, says when the estimated finite-difference
    error could hide an exact non-identifiability; the documented
    step-convergence check (recompute at ``rel_step / 3``) settles a full
    rank.
  - :func:`~nskinetics.analysis.morris` and
    :func:`~nskinetics.analysis.sobol`: global screening over parameter
    ranges (:class:`~nskinetics.analysis.ParamRange`, linear or log). Runs
    go through :func:`run_many`, failed runs and non-finite metric values
    are excluded and counted, and the same seed gives identical results.
  - Every analysis result's ``from_csv`` refuses an incomplete file (a
    missing ``meta`` key, or a missing or duplicated entry) with a
    ``ValueError`` instead of loading it.
  - RoadRunner's native forward sensitivities are deliberately not used.
    With libroadrunner 2.10.0 they drop event jumps, mislabel parameter
    subsets, and return zeros, an error or wrong values on the shipped
    model (see :doc:`concepts`).

* Parameter estimation for fed-batch data in ``nskinetics.analysis``
  (the :doc:`calibration tutorial <tutorial/08_calibration>`, and the
  :doc:`concepts` section "Calibration and identifiability").

  - :class:`~nskinetics.analysis.Dataset` holds measured time courses
    (from long or wide tables or a CSV file, with optional per-point
    standard deviations). :class:`~nskinetics.analysis.Observable` maps a
    model selection to a measured quantity, with an optional scale,
    offset and transform, and an additive, proportional or combined
    noise model. The standard deviations are computed once from the data
    and held fixed, so the fit is a weighted least-squares problem; the
    proportional model has the same floor as ``identifiability``.
  - :class:`~nskinetics.analysis.Experiment` pairs a :class:`Condition`,
    with its own feeding parameters and initial values, with its data.
    :class:`~nskinetics.analysis.FitParameter` estimates a global
    parameter, an initial value (also one per experiment,
    ``scope='local'``) or an observable's scale or offset, in log space
    by default; by default each starts from the value its run already
    starts from (the condition's or the model's).
  - :class:`~nskinetics.analysis.CalibrationProblem` validates
    everything before the first simulation and evaluates weighted
    residuals on a private copy of the model. A failed integration
    gives a penalty residual and is counted, never raised.
  - :func:`~nskinetics.analysis.fit_parameters` runs
    ``scipy.optimize.least_squares`` from several reproducible starts
    (``method='multistart'``, optionally on a process pool, and the
    first optionally from a previous fit, ``initial=``) or after a
    global stage (``method='hybrid'``), with the event-aware
    finite-difference Jacobian of ``local_sensitivity`` and optional
    robust losses (with absolute standard deviations). It raises
    :class:`FittingError` only when every start fails.
  - :class:`~nskinetics.analysis.FitResult` reports the estimates, the
    chi-square, AIC, AICc and BIC, the Wald covariance, standard errors,
    asymmetric 95% intervals (Student-t with ``absolute_sigma=False``)
    and correlations, and flags (lack of fit, rank deficiency,
    ill-conditioning, high relative standard errors or correlations,
    intervals wider than the bounds of a zero-capable linear parameter,
    estimates at a bound, and spike counts that change within a
    finite-difference step). It offers ``summary`` (warnings first),
    ``plot`` (data against the model, with its spikes marked), ``apply``
    (writes the global estimates into a model), ``sample`` and
    ``predict`` (Laplace parameter draws and prediction bands through
    :func:`run_many`), and an exact CSV round trip.
  - :func:`~nskinetics.analysis.profile_likelihood` returns a
    :class:`~nskinetics.analysis.ProfileResult` per parameter, with its
    profile-likelihood interval, ``is_identifiable``, a plot and an exact
    CSV round trip.
  - :class:`FittingError`, a new :class:`NSKError` subclass.

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
* ``nskinetics.utils.hybrid_global_optimize`` reported ``success=True``
  whatever happened. It now reports success only when differential
  evolution converged or a local polish did, and its message says which.
  Its result class, ``GlobalOptResult``, is exported from
  ``nskinetics.utils``.
* The shipped model's reset function and
  :class:`~nskinetics.units.NSKBatchReactor`'s spike-cap writes use item
  access (``r['env'] = 1.``) instead of attribute access. libroadrunner
  removes a model's ``r.<id>`` properties whenever another model is loaded,
  after which an attribute write silently missed the model.
