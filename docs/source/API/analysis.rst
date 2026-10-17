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
walks through every tool on a fed-batch model, and the :doc:`../concepts`
section "Sensitivities across events" explains how the derivatives are
computed and what their flags mean.

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

Parameter ranges
----------------

.. autoclass:: ParamRange
   :members: bounds, transformed_bounds, from_unit, to_unit,
             to_transformed, from_transformed, contains

References
----------

The identifiability measures and the screening designs and estimators
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
