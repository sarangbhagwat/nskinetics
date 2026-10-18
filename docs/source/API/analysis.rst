Analysis
========

.. currentmodule:: nskinetics.analysis

``nskinetics.analysis`` holds the tools that run a model many times to
answer a question about it. Each one works on a private copy of the model
(:meth:`nskinetics.KineticModel.copy`), so the caller's model is never
modified, and runs through :meth:`nskinetics.KineticModel.run` and
:func:`nskinetics.run_many`. It is a subpackage: import from
``nskinetics.analysis``, e.g. ``from nskinetics.analysis import Final,
local_sensitivity``. :doc:`../tutorial/07_sensitivity_identifiability`
walks through the sensitivity, identifiability and screening tools on a
fed-batch model, and the :doc:`calibration tutorial
<../tutorial/08_calibration>` fits a fed-batch model to the data of two
runs. The :doc:`../concepts` sections
"Sensitivities across events" and "Calibration and identifiability"
explain how the derivatives are computed, what their flags mean, and the
choices behind the fits.

Metrics
-------

A metric turns one :class:`~nskinetics.SimulationResult` into a float and
names the columns it reads (``name`` and ``selections``), so the same object
works in :func:`nskinetics.run_many` / :func:`nskinetics.scan` (as a
``metrics=`` entry), in :func:`local_sensitivity` and in screening. Each
metric states its basis explicitly; a zero denominator gives ``NaN``.

.. autoclass:: Metric
   :members: from_function
   :special-members: __call__

.. autoclass:: Final

.. autoclass:: At

.. autoclass:: Max

.. autoclass:: Min

.. autoclass:: TimeAtMax

.. autoclass:: Integral

.. autoclass:: EventCount

.. autoclass:: Productivity

.. autoclass:: Yield

.. autoclass:: CallableMetric

Local sensitivity
-----------------

:func:`local_sensitivity` differentiates output trajectories and metrics
with respect to parameters and initial values by finite differences in
log-parameter space, re-running the condition from the restored baseline for
each perturbation, and flags the parameters whose perturbation changed how
often, or before which output times, an event fired. Its result's
:attr:`~LocalSensitivityResult.metric_sensitivities` carry those flags on
every row.

.. autofunction:: local_sensitivity

.. autoclass:: LocalSensitivityResult
   :members: metric_sensitivities, time_course, rank, plot_time_course,
             plot_heatmap, to_csv, from_csv

Identifiability
---------------

:func:`identifiability` weights the log-scaled sensitivities of a
:class:`LocalSensitivityResult` by a noise model (absolute, or relative and
bounded below at a small fraction of each output's peak, per output) and
reports which parameter combinations the measured outputs can determine:
the Fisher information matrix, its singular values, numerical rank and
null-space directions, Cramer-Rao standard errors, and the Brun et al.
(2001) importance and collinearity measures behind
:meth:`IdentifiabilityResult.best_subset`. :ref:`identifiability-noise`
explains the noise models, what limits the rank test and how to check it.

.. autofunction:: identifiability

.. autoclass:: IdentifiabilityResult
   :members: correlation, table, collinearity, gamma, best_subset, plot,
             to_csv, from_csv

Global screening
----------------

:func:`morris` (elementary effects) and :func:`sobol` (first-order and total
variance-based indices) sample parameter ranges, run every sample through
:func:`nskinetics.run_many` and evaluate metrics; failed runs, and
non-finite metric values, are excluded and counted.

.. autofunction:: morris

.. autoclass:: MorrisResult
   :no-show-inheritance:
   :inherited-members:
   :members: n_runs, n_failed, n_used, design, samples, elementary_effects,
             rank, plot, to_csv, from_csv

.. autofunction:: sobol

.. autoclass:: SobolResult
   :no-show-inheritance:
   :inherited-members:
   :members: n_runs, n_failed, n_used, design, samples, rank, plot, to_csv,
             from_csv

Calibration
-----------

:class:`CalibrationProblem` compares :class:`Observable` values of a model
with measured :class:`Dataset` time courses in one or more
:class:`Experiment` conditions (each with its own feeding parameters and
initial state), and :func:`fit_parameters` estimates the
:class:`FitParameter` values by weighted least squares from several starts,
returning a :class:`FitResult` with Wald uncertainty, correlations and
identifiability and goodness-of-fit flags. Failed integrations inside the
optimizer are penalized and counted, never raised;
:exc:`nskinetics.FittingError` is raised only when every start fails.
:func:`profile_likelihood` re-optimizes the other parameters along a grid
of each profiled one, and :meth:`FitResult.sample` /
:meth:`FitResult.predict` propagate the Gaussian (Laplace) approximation to
parameter draws and prediction bands. Standard
deviations come from each observable's noise model (absolute,
proportional to the data with the floor of :func:`identifiability`, or
both) or from the data's ``sigma`` column, and are fixed before the fit.
Data times are in the model's time unit.

.. autoclass:: Dataset
   :members: from_frame, from_long_csv, from_wide, frame, experiments,
             observables

.. autoclass:: Observable
   :members: parameters, values, noise_sigma

.. autoclass:: Experiment

.. autoclass:: FitParameter
   :members: range, initial_for

.. autoclass:: CalibrationProblem
   :members: residuals, cost, to_natural, data, bounds, x0, log,
             n_failed_evals

.. autofunction:: fit_parameters

.. autoclass:: FitResult
   :members: summary, table, warnings, apply, plot, sample, predict, to_csv,
             from_csv

.. autofunction:: profile_likelihood

.. autoclass:: ProfileResult
   :members: estimate, plot, to_csv, from_csv

The global stage of ``fit_parameters(method='hybrid')`` is
:func:`nskinetics.utils.hybrid_global_optimize`, which also works on its
own for any scalar objective:

.. autofunction:: nskinetics.utils.hybrid_global_optimize

.. autoclass:: nskinetics.utils.GlobalOptResult
   :exclude-members: x, fun, success, message, nfev, history, meta

Parameter ranges
----------------

.. autoclass:: ParamRange
   :members: bounds, transformed_bounds, from_unit, to_unit,
             to_transformed, from_transformed, contains

References
----------

The identifiability measures, the screening designs and estimators, the
profile-likelihood intervals, the model-selection criteria, the
least-squares uncertainty and the column equilibration of its rank test
follow:

* Brun, R., Reichert, P., Kuensch, H. R. (2001). Practical identifiability
  analysis of large environmental simulation models. *Water Resources
  Research* 37(4), 1015-1030.
* Morris, M. D. (1991). Factorial sampling plans for preliminary
  computational experiments. *Technometrics* 33(2), 161-174.
* Campolongo, F., Cariboni, J., Saltelli, A. (2007). An effective screening
  design for sensitivity analysis of large models. *Environmental Modelling
  & Software* 22(10), 1509-1518.
* Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M.,
  Tarantola, S. (2010). Variance based sensitivity analysis of model output.
  Design and estimator for the total sensitivity index. *Computer Physics
  Communications* 181(2), 259-270.
* Jansen, M. J. W. (1999). Analysis of variance designs for model output.
  *Computer Physics Communications* 117(1-2), 35-43.
* Sobol', I. M. (2001). Global sensitivity indices for nonlinear
  mathematical models and their Monte Carlo estimates. *Mathematics and
  Computers in Simulation* 55(1-3), 271-280.
* Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M.,
  Klingmueller, U., Timmer, J. (2009). Structural and practical
  identifiability analysis of partially observed dynamical models by
  exploiting the profile likelihood. *Bioinformatics* 25(15), 1923-1929.
* Akaike, H. (1974). A new look at the statistical model identification.
  *IEEE Transactions on Automatic Control* 19(6), 716-723.
* Hurvich, C. M., Tsai, C.-L. (1989). Regression and time series model
  selection in small samples. *Biometrika* 76(2), 297-307.
* Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of
  Statistics* 6(2), 461-464.
* Bates, D. M., Watts, D. G. (1988). *Nonlinear Regression Analysis and Its
  Applications*. Wiley, New York.
* van der Sluis, A. (1969). Condition numbers and equilibration of
  matrices. *Numerische Mathematik* 14(1), 14-23.
