Sensitivity and identifiability
===============================

A model with a dozen parameters raises three questions before anyone fits
or engineers it. Which parameters control titer and yield, and are
therefore worth changing in the strain or the process? Which of them
could measurements of the broth actually determine? And does the answer
still hold over the whole range a parameter might take? This page answers
them with ``nskinetics.analysis`` on a small fed-batch model, in that
order. Local sensitivities give the first answer, the identifiability
analysis builds on them for the second, and global screening (Morris,
then Sobol') gives the third. Every analysis works on a private copy of
the model, so ``km`` itself is never changed. The analyses take a few
seconds in all.

The model and the metrics
-------------------------

The model is the fed-batch model from :doc:`../concepts`, with its rate
constant split into a catalytic constant ``kcat`` and an enzyme level
``E``, the knob a strain engineer would turn. Substrate ``S`` is converted
to product ``P``. A :class:`~nskinetics.FeedSpike` tops ``S`` back up to
``target`` (50 g/L) whenever it falls to ``threshold`` (5 g/L), at most
``n_max`` (3) times, with a feed of ``feed_conc`` (500 g/L) that dilutes
the broth:

.. code-block:: python

   import numpy as np
   import nskinetics as nsk
   from nskinetics.analysis import (Final, Yield, Productivity, EventCount,
                                    ParamRange, local_sensitivity,
                                    identifiability, morris, sobol)

   km = nsk.KineticModel.from_antimony("""
   model fb()
     compartment env = 1; species S in env, P in env;
     S = 50; P = 0; kcat = 0.5; E = 1; threshold = 5; target = 50;
     feed_conc = 500; n_max = 3; n_spk = 0; last_vol = 0; tot_vol = 0;
     J: S => P; kcat*E*S*env;
   end
   """, units={'time': 'h', 'conc': 'g/L'})
   spike = nsk.FeedSpike('S', 'S <= threshold', 'target', 'feed_conc',
                         volume_var='env', max_count='n_max', count_var='n_spk',
                         last_vol_var='last_vol', tot_vol_var='tot_vol')
   for e in spike.expand():
       km.add_event(e)
   km.compile_events()

The questions are asked of *metrics*: scalar summaries of a run, each with
an explicit basis. :class:`~nskinetics.analysis.Final` is the titer at the
end of the run. :class:`~nskinetics.analysis.Yield` with ``basis='added'``
divides the product formed by all the substrate that went in: the initial
charge plus ``tot_vol * feed_conc``, both as amounts (concentration times
the volume ``env``). :class:`~nskinetics.analysis.Productivity` divides the
titer gained by the run time plus 12 h of turnaround.
:class:`~nskinetics.analysis.EventCount` counts the spikes. The same
objects work as ``metrics=`` in :func:`~nskinetics.run_many`, which gives
their nominal values:

.. code-block:: python

   metrics = {'P_titer': Final('[P]'),
              'P_yield': Yield('[P]', '[S]', volume='env',
                               fed=('tot_vol', 'feed_conc'), basis='added'),
              'P_prod': Productivity('[P]', turnaround=12.0),
              'n_spikes': EventCount('n_spk')}
   cond = nsk.Condition(initial={'[S]': 50}, t_end=48, n_steps=481)
   print(nsk.run_many(km, [cond], metrics=metrics).metrics)

prints:

.. code-block:: text

     label status error  P_titer  P_yield  P_prod  n_spikes
   0  None     ok  None      162        1     2.7         3

All three spikes fire, and by 48 h practically all the substrate that went
in has become product (a yield of 1 on the added basis).

Ranking targets
---------------

:func:`~nskinetics.analysis.local_sensitivity` perturbs each parameter by
0.1% up and down and re-runs the condition from the model's baseline for
each perturbation. It returns the derivatives of the output time courses
and of the metrics. Initial values such as ``'[S]'`` can be perturbed like
parameters:

.. code-block:: python

   params = ['kcat', 'E', 'threshold', 'target', 'feed_conc', '[S]']
   ls = local_sensitivity(km, cond, params, outputs=['[S]', '[P]'],
                          metrics=metrics)
   print(ls.rank('P_titer'))

prints:

.. code-block:: text

      parameter  value  derivative  semi_elasticity  elasticity  reliable  event_count_changed  n_refinements
   0     target     50        2.25              113       0.696      True                False              0
   1        [S]     50       0.751             37.6       0.232      True                False              0
   2  feed_conc    500      0.0438             21.9       0.135      True                False              0
   3  threshold      5       -2.05            -10.2     -0.0633      True                False              0
   4       kcat    0.5    9.06e-05         4.53e-05     2.8e-07      True                False              0
   5          E      1    4.53e-05         4.53e-05     2.8e-07      True                False              0

``rank`` sorts by the *elasticity* ``d ln m / d ln p``: the percentage
change of the titer per percentage change of the parameter, which compares
parameters of different units and sizes. A 1% higher ``target`` raises the
48-h titer by 0.7%, because each spike then adds more substrate. A more
concentrated feed (``feed_conc``) dilutes the broth less. A higher
``threshold`` makes each spike top up from a higher level, so it adds less
substrate.
The enzyme level ``E`` and ``kcat`` barely matter: by 48 h the substrate is
used up whatever the rate, so the titer is set by how much substrate went
in.

That conclusion belongs to this condition, not to the model. Ask about the
titer at 12 h instead, and the rate decides:

.. code-block:: python

   at12 = ls.rank('[P]', at=12)
   print(at12[['parameter', 'elasticity', 'reliable', 'event_count_changed']])

prints:

.. code-block:: text

      parameter  elasticity  reliable  event_count_changed
   0       kcat       0.643      True                False
   1          E       0.643      True                False
   2     target         0.5      True                False
   3        [S]        0.25      True                False
   4  threshold       0.149      True                False
   5  feed_conc       0.101      True                False

``kcat`` and ``E`` have identical effects, since only their product enters
the rate. :meth:`~nskinetics.analysis.LocalSensitivityResult.plot_time_course`
draws the elasticities of ``[P]`` over the whole run:

.. code-block:: python

   fig, ax = ls.plot_time_course('[P]')

.. figure:: /_static/images/examples/tutorial_07_time_course.png
   :width: 500

   Elasticities of ``[P]`` (the ``kcat`` and ``E`` curves coincide).
   Before the first spike, ``[P]`` is proportional to the initial substrate
   (elasticity 1 for ``[S]``). The rate matters most around the spikes and
   fades once the substrate is used up, after which the feed settings set
   the titer. The peaks follow the three spikes: a faster rate brings each
   spike forward, so just after a spike the faster run has had longer to
   convert the fresh substrate.

``ls.metric_sensitivities`` lists the derivative, semi-elasticity and
elasticity of every metric with respect to every parameter, with the same
flags as ``rank`` (``reliable``, ``event_count_changed``,
``n_refinements``) on every row, and ``ls.time_course('[S]')`` the
elasticities of an output over time as a table. That table carries no
flags, so check ``ls.flags`` before using a column of it;
``plot_time_course`` draws unreliable parameters dashed.

Reading the flags
-----------------

Every row of the ranking carries the finite-difference diagnostics, also
in ``ls.flags``. ``ls.counters`` names the event counters that were
watched: here ``('n_spk',)``, the spike counter, found automatically. They
matter when a perturbation changes how often an event fires before the end
of the run: the output then jumps rather than changing smoothly, and the
derivative does not exist. They are compared at every output time, so a
spike moved across an output time, where the difference would measure the
jump, is caught as well. The third spike fires at
:math:`t_3 = 3 \ln(10)/(k_{cat} E) = 13.8155` h. Harvest at 13.8156 h,
just after it:

.. code-block:: python

   edge = cond.replace(t_end=13.8156, n_steps=101)
   ls_edge = local_sensitivity(km, edge, params, outputs=['[P]'],
                               metrics=metrics)
   print(ls_edge.flags[['step', 'n_refinements', 'event_count_changed',
                        'reliable', 'scheme']])

prints:

.. code-block:: text

              step  n_refinements  event_count_changed  reliable   scheme
   parameter
   kcat      0.001              2                 True     False  forward
   E         0.001              2                 True     False  forward
   threshold 1e-05              2                False      True  central
   target    1e-05              2                False      True  central
   feed_conc 0.001              0                False      True  central
   [S]       1e-05              2                False      True  central

A 0.1% change of most parameters moves the third spike past the harvest
time, so the spike count changes. NSKinetics then shrinks the step tenfold,
up to twice (``n_refinements``). For ``threshold``, ``target`` and
``[S]``, a step of 0.001% keeps the spike on the same side of the harvest,
and their derivatives are reliable. ``kcat`` and ``E`` shift the spike the
most: the count still changes at the smallest step, so they are flagged
``event_count_changed=True`` and ``reliable=False``. Their values are
one-sided differences (``scheme``) from the side on which the spike count
matches the nominal run, taken at the largest step at which it matched
(``step``), where integration error matters least. Do not read such a
value as a derivative.
Move the harvest time, or use the global screening below, which does not
need derivatives. ``feed_conc`` changes how much volume a spike adds but
not when it fires, so it needs no refinement. :doc:`../concepts` (section
"Sensitivities across events") explains the mechanism.

Which parameters can the data determine?
----------------------------------------

Suppose ``[S]`` and ``[P]`` are measured every 4 h with a standard
deviation of 0.5 g/L. :func:`~nskinetics.analysis.identifiability` weights
the log-scaled sensitivities at those times by the noise and analyzes the
resulting Fisher information matrix:

.. code-block:: python

   ident = identifiability(ls, sigma={'[S]': 0.5, '[P]': 0.5},
                           times=np.arange(0, 49, 4.0))
   print(ident.numerical_rank, 'of', len(ident.parameters))
   print(ident.table[['delta_msqr', 'se', 'identifiable']])
   print(ident.null_directions.round(3) + 0.0)   # + 0.0 shows -0.0 as 0

prints:

.. code-block:: text

   5 of 6
              delta_msqr      se  identifiable
   parameter
   kcat             88.7     inf         False
   E                88.7     inf         False
   threshold        30.2  0.0724          True
   target            131   0.023          True
   feed_conc        26.2  0.0941          True
   [S]              53.5 0.00766          True
                 v6
   parameter
   kcat       0.707
   E         -0.707
   threshold      0
   target         0
   feed_conc      0
   [S]            0

The matrix has rank 5, not 6. Its null direction (the overall sign of a
direction is arbitrary) moves ``ln kcat`` and ``ln E`` by the same amount
in opposite directions, and leaves every measurement as it was: raise
``kcat`` by some factor, lower ``E`` by the same factor, and nothing
changes. The data determine the product ``kcat * E`` but neither factor,
so their Cramer-Rao standard errors are infinite. The others are well
determined: ``se`` is the lower bound on the standard error of ``ln p``,
so ``target`` could be estimated to about 2%. The ``identifiable``
column asks for both: no loading on a null direction, and a relative
standard error (``rel_se`` in the table) of at most 1, so a parameter that
no exact compensation hides but that the data pin down only to within a
factor of several still reads ``False``. ``delta_msqr`` (Brun et al.,
2001) is each parameter's root-mean-square effect on the weighted
measurements, a measure of its importance.

The rank rests on the precision of the finite differences: a combination
the data cannot see shows up as a singular value at the level of the
integration noise and truncation error, not at 0. ``local_sensitivity``
tightens the integrator of its private copy to at most
``(1e-12, 1e-10)`` by default and uses central differences at
``rel_step=1e-3``; keep these defaults when the sensitivities feed
``identifiability``. Here ``kcat * E`` cancels exactly, but a
compensation such as ``(k1 + k2) * S`` does not: at CVODE's own
tolerances (``tolerances=(1e-12, 1e-6)``), with ``method='forward'`` or
with a large ``rel_step``, it can come out as full rank with finite
standard errors. ``identifiability`` warns when its estimate of that error
could hide such a direction, and records it as
``ident.rank_test_doubtful=True``; ``False`` is no proof, notably at a
depletion front (Monod growth, Michaelis-Menten kinetics with a small
``Km``). Before trusting a full rank, recompute the sensitivities with
``rel_step=1e-3 / 3`` and compare ``1 / ident.condition_number``: it falls
about 9 times for a null direction lifted by truncation error (3 times
with one-sided differences; integration noise makes it rise about 3
times), and barely moves for a genuine weak direction
(:ref:`identifiability-noise` in :doc:`../concepts` explains the rules and
shows the check in code).

The collinearity index :math:`\gamma` of a group of parameters measures how
well a change in one of them can be compensated by changes in the others
(1: not at all; above about 10: so well that they cannot be estimated
together). :meth:`~nskinetics.analysis.IdentifiabilityResult.best_subset`
picks, in order of importance, the parameters that can be estimated
together: with :math:`\gamma \le` ``gamma_max``, full rank, and a relative
standard error of at most 1 each (the others fixed):

.. code-block:: python

   print(ident.collinearity(subset_size=2).head(3))
   print(ident.best_subset(gamma_max=10))
   print(round(ident.gamma(['target', 'kcat', '[S]', 'threshold']), 1))

prints:

.. code-block:: text

               parameters  gamma  rho
   0            (kcat, E)    inf    0
   1  (target, feed_conc)    6.2 27.9
   2    (kcat, threshold)   3.29 33.5
   ['target', 'kcat', '[S]']
   13.1

``kcat`` and ``E`` cannot be told apart at all, so a fit should estimate
one of them with the other fixed (:doc:`08_calibration` shows how to fit
a chosen set of parameters). ``threshold`` and ``feed_conc`` are each
identifiable, but added to the three most important parameters they push
:math:`\gamma` above 10. Their estimates would then be strongly correlated
with the others'. Fix them at known values (they are process settings
here), or measure something that separates them. The analysis is local: it
holds at these parameter values, for these measurement times and this
noise level.

Measurement errors are often relative instead, a few percent of the
measured value. ``sigma=('relative', 0.05)`` gives each row a standard
deviation of 5% of the output's value there, but never less than 5% of
``relative_floor`` (default ``1e-3``) times the output's peak over the run:

.. code-block:: python

   rel = identifiability(ls, sigma=('relative', 0.05),
                         times=np.arange(0, 49, 4.0))
   print(rel.numerical_rank, 'of', len(rel.parameters))
   print(rel.table[['delta_msqr', 'se', 'identifiable']])
   print(rel.best_subset(gamma_max=10))

prints:

.. code-block:: text

   5 of 6
              delta_msqr     se  identifiable
   parameter
   kcat             88.1    inf         False
   E                88.1    inf         False
   threshold        24.3  0.167          True
   target           25.7  0.163          True
   feed_conc        1.67  0.791          True
   [S]                12 0.0271          True
   ['kcat', 'target', 'threshold', '[S]']

The floor acts like a detection limit. From 28 h on, the substrate is below
0.05 g/L, 1e-3 of its 50 g/L peak, and so is the product at 0 h; those rows
count as measurements at that limit. A pure relative error would give them
a weight that grows without bound as the value approaches 0, while the
integration error of the finite differences there does not shrink with it.
The floor is a fraction of each output's own peak, so it does not depend on
the units, and no row is left out (:ref:`identifiability-noise` in
:doc:`../concepts` explains the choice). Only the product ``kcat * E`` is
determined, as before. The subset differs from the one above because the
two noise models weigh the measurement times differently: 5% of a value is
less than 0.5 g/L only while the value is below 10 g/L. A mapping can also
give each output its own model, e.g.
``sigma={'[S]': 0.5, '[P]': ('relative', 0.05)}``.

Screening over ranges: Morris, then Sobol'
------------------------------------------

Derivatives describe small changes. Over wide ranges, and for a discrete
setting such as the spike cap ``n_max``, use global screening. At
``n_max = 3``, any increase allows a fourth spike (the titer jumps from 162
to 193 g/L) and a small decrease changes nothing, so there is no
derivative: :func:`~nskinetics.analysis.local_sensitivity` would flag
``n_max`` with ``event_count_changed=True``. The ranges are given as
``(lo, hi)`` or as :class:`~nskinetics.analysis.ParamRange`, sampled in log
space with ``log=True``. :func:`~nskinetics.analysis.morris` is the cheap
first pass: ``r`` random one-at-a-time trajectories, ``r * (d + 1)`` runs
for ``d`` parameters:

.. code-block:: python

   bounds = {'kcat': ParamRange(0.05, 5, log=True), 'threshold': (1, 20),
             'target': (30, 80), 'n_max': (0, 5)}
   g = morris(km, cond, bounds, metrics, r=20, seed=0)
   print(g.n_runs, g.n_failed)
   print(g.rank('P_titer'))

prints:

.. code-block:: text

   100 0
      parameter    mu  mu_star  sigma  mu_star_conf
   0      n_max   116      116   93.7          40.1
   1     target  95.6     95.6   60.8          25.6
   2       kcat  67.4     67.4    103          42.3
   3  threshold -4.62     37.9   55.8          17.1

``mu_star`` is the mean absolute change of the titer (g/L) per full range
of a parameter (of ``ln p`` for a log range), ``mu`` the signed mean,
``sigma`` the spread of the changes, which comes from nonlinearity and
interactions, and ``mu_star_conf`` a 95% bootstrap half-width. Over these
ranges the spike cap matters most, which the local analysis could not
show. ``kcat``, which barely mattered at the nominal point, now matters: at
its low end the substrate is no longer used up within 48 h. ``threshold``
has a small mean effect ``mu`` but a large ``mu_star``, because its effects
have both signs. A higher threshold makes each spike add less substrate,
but when the rate is low it also lets more spikes fire within 48 h (at
``kcat = 0.1``, raising it from 5 to 15 g/L allows a third spike and raises
the titer from 87 to 125 g/L).

:func:`~nskinetics.analysis.sobol` then apportions the variance of each
metric. ``S1`` is the share a parameter explains alone and ``ST`` its share
including all its interactions, from ``n * (d + 2)`` runs (``n`` a power of
2):

.. code-block:: python

   s = sobol(km, cond, bounds, metrics, n=256, seed=0)
   print(s.n_runs, s.n_failed)
   print(s.rank('P_titer'))
   fig, ax = s.plot('P_titer')

prints:

.. code-block:: text

   1536 0
      parameter    S1     ST  S1_conf  ST_conf
   0      n_max 0.386  0.494    0.127    0.112
   1       kcat  0.25  0.378    0.131    0.121
   2     target 0.245  0.247    0.134   0.0515
   3  threshold 0.042 0.0557    0.132   0.0171

.. figure:: /_static/images/examples/tutorial_07_sobol.png
   :width: 450

   First-order (``S1``) and total (``ST``) Sobol' indices of the 48-h titer,
   with 95% bootstrap intervals.

The spike cap and the rate explain most of the titer's variance, and the
gaps between ``ST`` and ``S1`` (and ``S1`` summing to about 0.92 rather
than 1) suggest that they interact, although at ``n=256`` the gaps are
within the intervals: a higher cap only adds titer if the rate is high
enough for the substrate to fall to the threshold again within 48 h.
``target`` acts almost alone. The intervals are wide at ``n=256``. Each
doubling of ``n`` narrows them by about :math:`\sqrt{2}` and doubles the
runs, so increase ``n`` before trusting the order of close parameters.
``P_prod`` has the same indices as ``P_titer`` here, because with a fixed
run time it is the titer divided by a constant.

Each sample is the base condition with the sampled values merged in, so
every value not screened keeps its value from ``cond``. Runs that fail are
recorded, excluded from the indices and counted in ``n_failed``, with a
warning when more than 5% fail. The same ``seed`` gives identical
results. ``workers=4`` spreads the runs over four processes (see
:func:`~nskinetics.run_many`).

Saving results
--------------

Every result saves to one CSV file and loads back exactly:

.. code-block:: python

   from nskinetics.analysis import LocalSensitivityResult, SobolResult

   ls.to_csv('sensitivity.csv')
   s.to_csv('sobol.csv')
   back = SobolResult.from_csv('sobol.csv')
   print(back.ST.equals(s.ST),
         LocalSensitivityResult.from_csv('sensitivity.csv').flags.equals(ls.flags))

prints:

.. code-block:: text

   True True

:doc:`../concepts` (section "Sensitivities across events") explains why the
derivatives are computed by finite differences, what the flags guard
against, and the basis of each metric. The :doc:`API reference
<../API/analysis>` lists every option.

Next: :doc:`08_calibration` fits the identifiable parameters to fed-batch
data.
