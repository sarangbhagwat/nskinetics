Calibration
===========

Fitting rate constants to measured time courses is the first thing a
model of a new strain needs, and in fed-batch operation it comes with a
twist: every run has its own feeding settings and its own, slightly
uncertain, initial charge. ``nskinetics.analysis`` fits a model to several
such runs at once, each described by its own :class:`~nskinetics.Condition`,
and reports how well the data determine each estimate. This page fits two
synthetic fed-batch runs, diagnoses a fit that went wrong, profiles the
estimates, propagates their uncertainty to predictions, and writes the
result back into the model. Every fit works on a private copy of the
model, so ``km`` changes only when the result is applied at the end. The
whole page runs in about 10 s, plus the import.

The model and the data
----------------------

The model is the fed-batch model of :doc:`../concepts`: substrate ``S`` is
converted to product ``P`` at the rate ``k*S``, and a
:class:`~nskinetics.FeedSpike` tops ``S`` back up to ``target`` (50 g/L)
whenever it falls to ``threshold``, at most ``n_max`` (3) times, with a feed
of ``feed_conc`` (500 g/L) that dilutes the broth:

.. code-block:: python

   import numpy as np
   import pandas as pd
   import nskinetics as nsk
   from nskinetics.analysis import (Dataset, Observable, Experiment,
                                    FitParameter, CalibrationProblem,
                                    fit_parameters, profile_likelihood,
                                    FitResult)

   km = nsk.KineticModel.from_antimony("""
   model fb()
     compartment env = 1; species S in env, P in env;
     S = 50; P = 0; k = 0.5; threshold = 5; target = 50; feed_conc = 500;
     n_max = 3; n_spk = 0; last_vol = 0; tot_vol = 0;
     J: S => P; k*S*env;
   end
   """, units={'time': 'h', 'conc': 'g/L'})
   spike = nsk.FeedSpike('S', 'S <= threshold', 'target', 'feed_conc',
                         volume_var='env', max_count='n_max', count_var='n_spk',
                         last_vol_var='last_vol', tot_vol_var='tot_vol')
   for e in spike.expand():
       km.add_event(e)
   km.compile_events()

Two fed-batch runs differ in their feeding: ``FB1`` starts from 50 g/L of
substrate and is fed when ``S`` falls to 5 g/L, ``FB2`` starts from 30 g/L
and is fed at 10 g/L. Each recipe is a :class:`~nskinetics.Condition`.
The synthetic "truth" behind the measurements has ``k = 0.4`` 1/h (the
model says 0.5), and the initial charges missed their recipes: 48.5 g/L
in ``FB1`` and 27 g/L in ``FB2``. ``S`` and ``P`` are measured every hour
for 20 h, ``S`` with a standard deviation of 1 g/L and ``P`` with one of
5%:

.. code-block:: python

   recipes = {
       'FB1': nsk.Condition(initial={'[S]': 50}, parameters={'threshold': 5},
                            t_end=20),
       'FB2': nsk.Condition(initial={'[S]': 30}, parameters={'threshold': 10},
                            t_end=20)}
   true_s0 = {'FB1': 48.5, 'FB2': 27.0}
   t = np.arange(1.0, 21.0)
   rng = np.random.default_rng(0)
   tables = []
   for name, recipe in recipes.items():
       truth = (recipe.with_parameters(k=0.4)
                .with_initial({'[S]': true_s0[name]})
                .replace(times=np.r_[0.0, t], selections=['[S]', '[P]']))
       res = km.run(truth)
       tables.append(pd.DataFrame({
           'run': name, 'time_h': t,
           'S': res['[S]'][1:] + rng.normal(0.0, 1.0, t.size),
           'P': res['[P]'][1:] * (1.0 + rng.normal(0.0, 0.05, t.size))}))
   wide = pd.concat(tables, ignore_index=True).round(2)
   print(wide.head(7))

prints:

.. code-block:: text

      run  time_h    S    P
   0  FB1       1 32.6 15.9
   1  FB1       2 21.7 28.5
   2  FB1       3 15.2 32.8
   3  FB1       4  9.9 39.4
   4  FB1       5 6.03 43.8
   5  FB1       6 44.4 45.8
   6  FB1       7 30.8 57.8

The first spike of ``FB1`` came between 5 and 6 h, when ``S`` was topped
back up to 50 g/L.

The data are a wide table, as a lab would keep it. A
:class:`~nskinetics.analysis.Dataset` turns it into the long format the
fit uses, one row per measurement
(:meth:`~nskinetics.analysis.Dataset.from_long_csv` and
:meth:`~nskinetics.analysis.Dataset.from_frame` read long tables
directly):

.. code-block:: python

   data = Dataset.from_wide(wide, experiment='run', time='time_h',
                            columns=['S', 'P'])
   print(data)
   print(data['FB2'].frame.head(3))

prints:

.. code-block:: text

   Dataset(80 rows, experiments=['FB1', 'FB2'], observables=['S', 'P'])
     experiment observable  time  value  sigma
   0        FB2          S     1   16.8    NaN
   1        FB2          S     2   13.7    NaN
   2        FB2          S     3     42    NaN

Times must be in the model's time unit, here hours: no unit conversion is
applied. The ``sigma`` column is empty, so each point's standard deviation
comes from its observable's noise model.

The calibration problem
-----------------------

An :class:`~nskinetics.analysis.Observable` maps a model selection to a
measured quantity and gives its noise model: an absolute standard
deviation of 1 g/L for ``S``, and 5% of the measured value for ``P``
(never less than 5% of 0.1% of the run's largest ``P``, so that values
near 0 do not get an unbounded weight). An
:class:`~nskinetics.analysis.Experiment` pairs a run's condition with its
data. A :class:`~nskinetics.analysis.FitParameter` names what to estimate:
the rate constant ``k``, shared by both runs, and the initial substrate
``[S]`` with ``scope='local'``, one free value per run. Both are fitted in
log space (the default) within their bounds. Without an ``initial`` value,
each starts from what its run already starts from: the model's ``k``
(0.5) and each recipe's ``[S]`` (50 and 30 g/L):

.. code-block:: python

   observables = [Observable('S', '[S]', sigma=1.0),
                  Observable('P', '[P]', noise='proportional', rel_sigma=0.05)]
   experiments = [Experiment(name, recipe, data[name])
                  for name, recipe in recipes.items()]
   parameters = [FitParameter('k', bounds=(0.05, 2.0)),
                 FitParameter('[S]', bounds=(20, 80), scope='local')]
   problem = CalibrationProblem(km, experiments, observables, parameters)
   print(problem.names, problem.n_data, problem.event_counters)

prints:

.. code-block:: text

   ['k', '[S]@FB1', '[S]@FB2'] 80 ('n_spk',)

The local initial values are named ``'[S]@FB1'`` and ``'[S]@FB2'``. Each
model evaluation runs every experiment's own condition, with the
candidate values merged in, on the output times ``{t0}`` plus its data
times; a value set for one run never reaches the other. Everything is
checked before the first simulation: unknown targets, data outside a run,
missing standard deviations, more free parameters than data points. The
problem found the spike counter ``n_spk``, which the finite differences
watch (see "Sensitivities across events" in :doc:`../concepts`).

A first fit, and how to tell that it failed
-------------------------------------------

:func:`~nskinetics.analysis.fit_parameters` minimizes the weighted sum of
squares from several starts (the first at the initial values, the others
spread over the bounds) and keeps the best:

.. code-block:: python

   first = fit_parameters(problem, n_starts=5, seed=0)
   print(first.summary())

prints:

.. code-block:: text

   Calibration fit (multistart, loss='linear')
   3 free parameter(s), 80 data points, 5/5 start(s) evaluated

   Warnings (3):
     - [lack_of_fit] chi2 = 9738 on 77 degrees of freedom (reduced chi2 126) exceeds the
       99.9% chi-square quantile 121.1: the model misses the data by more than the given
       standard deviations, and the Wald uncertainties assume it does not. Causes: a wrong
       local minimum (compare the costs in the starts table and look at plot(); try more
       starts or method='hybrid'), a model that cannot describe the data, or standard
       deviations smaller than the real noise or outliers (if the sigmas are only relative
       weights, use absolute_sigma=False).
     - [event_count_changed] k: a perturbation changes how many times (or between which
       output times) an event fires, so the objective is not smooth there; its derivative
       is one-sided and its Wald uncertainty approximate (check it with a profile
       likelihood).
     - [event_count_changed] [S]@FB2: a perturbation changes how many times (or between
       which output times) an event fires, so the objective is not smooth there; its
       derivative is one-sided and its Wald uncertainty approximate (check it with a
       profile likelihood).

   Parameters (natural units; SE first order; 95% Wald CI):
              estimate         SE     rel SE  CI lower  CI upper  lower bound  upper bound space
     k         0.48064 0.00131621 0.00273845  0.478067  0.483226         0.05            2   log
     [S]@FB1   45.0204   0.587364  0.0130466   43.8837   46.1864           20           80   log
     [S]@FB2   30.2491   0.380579  0.0125815   29.5123   31.0043           20           80   log

   Fit statistics:
     chi2 (weighted SSR)      9738.25
     cost (0.5 sum rho)       4869.12
     n_data                   80
     n_free                   3
     dof                      77
     reduced chi2             126.471
     AIC                      9744.25
     AICc                     9744.56
     BIC                      9751.39
     rank / cond(J D)         3 / 2.45
     sigma                    absolute
     starts at the best cost  1 of 5
     failed evaluations       0
     success                  True
     message                  `xtol` termination condition is satisfied.

   Model units: time [h], conc [g/L]; data times are in model time units.

The summary puts its warnings first. The first is ``lack_of_fit``: the
standard deviations are the true ones (``sigma absolute``), so a good fit
has a chi-square near its 77 degrees of freedom (a reduced chi-square
near 1), and 9738 is far beyond anything the noise explains. The other
two flag ``event_count_changed`` for two parameters. At this estimate the
third spike of ``FB2`` falls on the 9-h sample, within a fraction of a
second, and a change of ``k`` or of ``FB2``'s initial substrate as small
as 0.001% moves it to the other side of that sample: the objective jumps
right next to the optimum. The estimates do not contain the truth, and
the narrow Wald intervals around them mean nothing. The optimizer still
reports ``success``: it converged, to the wrong place. The statistics
give a hint: only one of the five starts reached the best cost
(``starts at the best cost  1 of 5``). The starts show why:

.. code-block:: python

   print(first.starts[['start', 'cost', 'estimate:k', 'estimate:[S]@FB1',
                       'estimate:[S]@FB2']])

prints:

.. code-block:: text

      start     cost  estimate:k  estimate:[S]@FB1  estimate:[S]@FB2
   0      0 5.32e+03       0.515              40.2              28.9
   1      1 7.68e+03       0.194              35.9              54.1
   2      2 4.87e+03       0.481                45              30.2
   3      3 1.02e+04       0.794                20                20
   4      4 7.24e+03       0.202              62.3              22.9

Every start converged, each to a different point with a different cost;
start 3 even ended on the lower bounds of both initial values. A plot of
the best one against the data shows what went wrong:

.. code-block:: python

   fig, axes = first.plot()

.. figure:: /_static/images/examples/tutorial_08_local_minimum.png
   :width: 600

   The stuck fit. Dots are the data with their standard deviations, lines
   the model at the estimates, dotted lines the model's spikes. The
   model's spikes come up to 3 h too early in ``FB1``, and its third
   spike in ``FB2`` sits on the 9-h sample, while the data spiked at
   about 10.5 h.

A spike makes ``S`` jump by 40 to 45 g/L, 40 to 45 standard deviations.
When a parameter change moves a spike across a sample time, the model
value at that sample jumps by as much, and the chi-square by 1600 or
more: the objective is smooth only between such jumps, a patchwork of
pieces, one per arrangement of spikes and samples. A local optimizer
follows the smooth piece it starts on and stops at the best point of that
piece, often at its edge, next to a jump, which is where the flags caught
this one.
:doc:`../concepts` (section "Calibration and identifiability") explains
what the flags can and cannot see.

A global stage
--------------

``method='hybrid'`` first searches the whole box of bounds with a
population-based global optimizer (differential evolution) on the same
objective, which needs no derivatives and is not held up by a jump, and
then polishes the best point it found, last with the same least-squares
solver:

.. code-block:: python

   fit = fit_parameters(problem, method='hybrid', seed=0)
   print(fit.summary())

prints:

.. code-block:: text

   Calibration fit (hybrid, loss='linear')
   3 free parameter(s), 80 data points, global stage + least-squares polish

   Warnings: none

   Parameters (natural units; SE first order; 95% Wald CI):
              estimate         SE     rel SE  CI lower  CI upper  lower bound  upper bound space
     k        0.400585 0.00103974 0.00259556  0.398553  0.402628         0.05            2   log
     [S]@FB1   48.7094   0.580066  0.0119087   47.5856   49.8596           20           80   log
     [S]@FB2   27.4932   0.325325  0.0118329    26.863   28.1383           20           80   log

   Fit statistics:
     chi2 (weighted SSR)  71.8435
     cost (0.5 sum rho)   35.9218
     n_data               80
     n_free               3
     dof                  77
     reduced chi2         0.933033
     AIC                  77.8435
     AICc                 78.1593
     BIC                  84.9896
     rank / cond(J D)     3 / 2.68
     sigma                absolute
     failed evaluations   0
     success              True
     message              Global stage (1026 objective evaluations): DE did not converge
                          (Maximum number of iterations has been exceeded.); 1/1 L-BFGS-B
                          polish(es) converged; polish: `ftol` termination condition is
                          satisfied.

   Model units: time [h], conc [g/L]; data times are in model time units.

The reduced chi-square is now 0.93, as it should be for the true standard
deviations, and every 95% interval contains the truth (0.4 1/h, 48.5 and
27 g/L). The message reports the three stages: the differential evolution
used up its budget of 30 generations without meeting its own tolerance,
which is normal, since it only has to land on the right piece; a
bounded quasi-Newton polish and the final ``least_squares`` polish then
converged. ``k`` is determined to 0.26% and each initial charge to 1.2%.
The SE column is first order in natural units; the intervals are computed
in log space, so they are slightly asymmetric.

``method='hybrid'`` is one remedy; the others are more starts, better
initial values, and a robust loss, which caps the pull of a residual that
sits across a jump (see :doc:`../concepts`). On these data, one of twenty
starts reached this optimum with the linear loss; with
``loss='soft_l1', f_scale=3``, five of ten reached the same piece.

The flags are empty, and the parameter table says why each parameter
passed:

.. code-block:: python

   print(fit.flags.empty)
   print(fit.table[['estimate', 'rel_se', 'identifiable', 'reliable',
                    'event_count_changed']])
   print(fit.correlation.round(2))

prints:

.. code-block:: text

   True
              estimate  rel_se  identifiable  reliable  event_count_changed
   parameter
   k             0.401  0.0026          True      True                False
   [S]@FB1        48.7  0.0119          True      True                False
   [S]@FB2        27.5  0.0118          True      True                False
                k  [S]@FB1  [S]@FB2
   parameter
   k            1     0.67     0.58
   [S]@FB1   0.67        1     0.39
   [S]@FB2   0.58     0.39        1

Each estimate is ``identifiable`` (finite standard error, relative
standard error at most 1; a ``log=False`` parameter whose bounds include 0
needs a 95% interval no wider than its bounds instead), its finite
differences are ``reliable``, and
no perturbation moved a spike across a sample. The correlations are
positive but moderate: a faster rate brings the first spike forward and a
larger initial charge delays it, so raising both together keeps it in
place, but the later spikes and the product data tell them apart.
Correlations above 0.95, a rank-deficient or ill-conditioned Jacobian, a
relative standard error above 1 or an estimate at a bound would each add
a flag (:class:`~nskinetics.analysis.FitResult` lists them all).

.. code-block:: python

   fig, axes = fit.plot()

.. figure:: /_static/images/examples/tutorial_08_fit.png
   :width: 600

   The fit: the model's spikes now line up with the data's.

Profile likelihoods
-------------------

The Wald intervals come from the curvature of the objective at the
optimum. :func:`~nskinetics.analysis.profile_likelihood` checks them
against the objective itself: it fixes one parameter at a grid of values
outward from the estimate (here 5 per side, over 3 standard errors),
re-optimizes the others at each, and reads the 95% interval where the
chi-square has risen by 3.84:

.. code-block:: python

   import matplotlib.pyplot as plt

   profiles = profile_likelihood(fit, fit.names, n_points=11)
   for name, prof in profiles.items():
       lo, hi = prof.confidence_interval
       wald = fit.ci.loc[name]
       print(f"{name:8s} profile ({lo:.4g}, {hi:.4g})  "
             f"Wald ({wald['lower']:.4g}, {wald['upper']:.4g})  "
             f"identifiable: {prof.is_identifiable}")
   fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
   for ax, prof in zip(axes, profiles.values()):
       prof.plot(ax=ax)

prints:

.. code-block:: text

   k        profile (0.3986, 0.4026)  Wald (0.3986, 0.4026)  identifiable: True
   [S]@FB1  profile (47.57, 49.85)  Wald (47.59, 49.86)  identifiable: True
   [S]@FB2  profile (26.86, 28.13)  Wald (26.86, 28.14)  identifiable: True

.. figure:: /_static/images/examples/tutorial_08_profiles.png
   :width: 700

   Profiles (solid) against the Wald parabolas (dotted), with the 3.84
   threshold (black) and the interval ends (red).

Profile and Wald intervals agree to within 0.05%: over these intervals the
objective is close to quadratic in the log parameters, and no spike
crosses a sample. They part ways when a parameter is poorly determined
(a skewed or flat profile, ``is_identifiable=False`` when a side never
crosses the threshold within the bounds), and when a spike sits within
the uncertainty of a sample time, where the profile shows a jump. Each
grid value is a re-optimization of the other parameters, so profiles cost
many fits: profile the parameters the flags point at, and the ones a
decision rests on.

From parameters to predictions
------------------------------

:meth:`~nskinetics.analysis.FitResult.sample` draws parameter sets from
the same Gaussian approximation (the Laplace approximation: normal in log
space around the estimates, with the fit's covariance, clipped to the
bounds), as a table for any downstream use:

.. code-block:: python

   draws = fit.sample(500, seed=1)
   print(draws.describe().loc[['mean', 'std']])

prints:

.. code-block:: text

              k  [S]@FB1  [S]@FB2
   mean     0.4     48.7     27.5
   std  0.00106     0.57     0.34

The spread of the draws matches the standard errors.
:meth:`~nskinetics.analysis.FitResult.predict` runs such draws through
:func:`~nskinetics.run_many` and returns percentile bands at every output
time. Given an experiment's name, it runs that experiment's condition with
its own local values; here, around the third spike of ``FB1``:

.. code-block:: python

   bands = fit.predict('FB1', n_samples=200, seed=1)
   s = bands[bands['output'] == '[S]']
   print(s[s['time'].between(17.1, 17.26)].to_string(
       index=False, float_format=lambda v: f'{v:.3f}'))

prints:

.. code-block:: text

     time output  nominal     q5    q50    q95
   17.117    [S]    5.125  5.035  5.134  5.248
   17.137    [S]    5.084  5.018  5.098 49.776
   17.157    [S]    5.044  5.010  5.084 49.835
   17.177    [S]    5.003  5.010  5.117 49.888
   17.197    [S]   49.634  5.008 49.128 49.907
   17.217    [S]   49.238  5.015 49.062 49.833
   17.237    [S]   48.844 47.614 48.856 49.742
   17.257    [S]   48.454 47.529 48.490 49.521

Away from the spikes the 90% band of ``[S]`` is less than 2 g/L wide.
Around a spike it covers the whole jump: at least 5% of the draws have
spiked by 17.137 h, and at least 5% have not by 17.217 h, while the
nominal run (``nominal``, at the estimates) spikes between 17.177 and
17.197 h (the condition's grid has 1000 output times over 20 h). In a
fed-batch run, the uncertainty of a rate becomes an uncertainty of *when*
the feed comes.

A plain :class:`~nskinetics.Condition` predicts a run that was not
fitted, such as a third run fed already at 20 g/L. It gets the global
estimates of each draw (``k``); the local initial charges belong to their
own runs, so the condition's own ``[S]`` is used:

.. code-block:: python

   fb3 = nsk.Condition(initial={'[S]': 50}, parameters={'threshold': 20},
                       times=[0, 4, 8, 12, 24])
   plan = fit.predict(fb3, n_samples=200, seed=1, outputs=['[P]'])
   print(plan.attrs['parameters'], plan.attrs['n_failed'])
   print(plan)

prints:

.. code-block:: text

   ['k'] 0
      time output  nominal   q5  q50  q95
   0     0    [P]        0    0    0    0
   1     4    [P]     52.9 52.8 52.9 53.1
   2     8    [P]     97.5 97.1 97.5 97.8
   3    12    [P]      123  123  123  123
   4    24    [P]      129  129  129  129

None of the 200 runs failed (``n_failed``). The early titers carry the
uncertainty of ``k``; the final one is set by the substrate fed, not by
the rate. The bands hold the parameter uncertainty only, not the
measurement noise, and the Gaussian approximation behind them is only as
good as the agreement of the profiles with the Wald intervals above.
Run-to-run variation of the initial charge, which missed its recipe by
1.5 and 3 g/L in the fitted runs, is not in them either: draw it
yourself into the conditions of :func:`~nskinetics.run_many` to include
it.

Saving and applying the result
------------------------------

A fit saves to one CSV file and loads back exactly, with its starts,
residuals, Jacobian and flags. A loaded result can still be summarized,
sampled and, given a model (``predict(..., model=km)``), used for
predictions; profiles need the problem, so they need the original fit.
:meth:`~nskinetics.analysis.FitResult.apply` writes the global estimates
into a model:

.. code-block:: python

   fit.to_csv('fed_batch_fit.csv')
   saved = FitResult.from_csv('fed_batch_fit.csv')
   print(saved.summary() == fit.summary())

   leftover = fit.apply(km)
   print({name: round(value, 2) for name, value in leftover.items()})
   print(f"k = {km.parameters.loc['k', 'value']:.4f} 1/h")

prints:

.. code-block:: text

   True
   {'[S]@FB1': 48.71, '[S]@FB2': 27.49}
   k = 0.4006 1/h

``k`` is now 0.4006 in ``km``, and stays so through resets and in
copies. The initial charges are returned instead of applied: they belong
to their runs, not to the model, so ``km`` still starts from 50 g/L.
Use them in those runs' conditions, e.g.
``recipes['FB1'].with_initial({'[S]': 48.71})``.

The :doc:`../concepts` section "Calibration and identifiability"
explains the noise models and weights, why the fit works in log space,
what the flags can and cannot see, and how the Wald, profile and Laplace
uncertainties differ. Before fitting a larger model, ask
:doc:`07_sensitivity_identifiability` which parameters the planned
measurements can determine, and fit only those. The :doc:`API reference
<../API/analysis>` lists every option.
