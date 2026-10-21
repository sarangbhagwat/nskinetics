A real fed-batch model
=======================

:doc:`03_events` and :doc:`04_fed_batch_feeding` built ``Event`` and
``FeedSpike`` from small, hand-written Antimony models. This page runs the
same machinery on a **shipped, production-scale** model: the *S. cerevisiae*
fed-batch fermentation SBML that ships inside the NSKinetics package and
drives the fermentor of the
:doc:`isobutanol biorefinery case study <../case_studies/isobutanol/index>`. It has
ten floating species (glucose, ethanol, isobutanol, biomass, and six
pathway intermediates) and dozens of kinetic parameters — the same shape of
model as :doc:`04_fed_batch_feeding`'s toy example, just larger.

Loading the shipped SBML
-------------------------

The file lives under ``nskinetics/models/s_cerevisiae_ferm_fb_inhib_mod_ibo/``,
next to the model's Antimony source and builder module:

.. code-block:: python

   import os, tellurium as te, nskinetics as nsk

   p = os.path.join(os.path.dirname(nsk.__file__),
                     'models', 's_cerevisiae_ferm_fb_inhib_mod_ibo',
                     's_cerevisiae_ferm_fb_inhib_mod_ibo_sbml.xml')
   r = te.loadSBMLModel(p)
   km = nsk.KineticModel(r, units={'time': 'h', 'conc': 'g/L'})

Inspecting the loaded model before doing anything else confirms what it
exposes:

.. code-block:: text

   species: ['s_glu', 's_pyr', 's_acetate', 's_acetald', 's_EtOH', 'x',
             's_AL', 's_DHI', 's_KIV', 's_IBO']
   n events: 0

``s_glu`` (glucose), ``s_EtOH`` (ethanol), ``s_IBO`` (isobutanol), and ``x``
(biomass) are the four species this page tracks. ``n events: 0`` matters: the
shipped SBML ships **without** any baked-in events — the fed-batch glucose
feed and the aerobic/anaerobic stage switch are not part of the file itself,
and must be added through the same :class:`~nskinetics.FeedSpike` /
:class:`~nskinetics.Event` API from the last two pages, exactly as the
isobutanol biorefinery's own fermentor does.

Adding the events conditionally
---------------------------------

Because a *different* copy of this same model might ship with events already
compiled in (for example, if it were re-exported after a prior
``compile_events()`` call), the robust pattern checks
``r.getNumEvents()`` first and only adds events if there are none. This
makes the snippet correct either way, and is exactly what runs here since
the shipped file has zero events:

.. code-block:: python

   # --- add events via the API ONLY if the SBML has none ---
   if r.getNumEvents() == 0:
       spike = nsk.FeedSpike(
           species='s_glu', when='s_glu <= threshold_conc_glu_spike',
           target='target_conc_glu_spike', feed_conc='conc_glu_feed_spike',
           volume_var='env', max_count='max_n_glu_spikes', count_var='n_glu_spikes',
           last_vol_var='last_vol_glu_feed_added', tot_vol_var='tot_vol_glu_feed_added',
           delay='glucose_feed_spikeDelay', priority=5, name='glucose_feed_spike')
       for e in spike.expand():
           km.add_event(e)
       km.add_event(nsk.Event(when='time >= stage_1_max_time',
                               do={'is_aerobic': '0'}, name='stage_1_complete_max_time'))
       km.add_event(nsk.Event(when='x >= stage_1_max_x',
                               do={'is_aerobic': '0'}, name='stage_1_complete_x_target'))
       km.compile_events()

The ``FeedSpike`` watches ``s_glu`` and tops it back up to
``target_conc_glu_spike`` (in the growing ``env`` compartment) whenever it
drops to ``threshold_conc_glu_spike`` — the same pattern from
:doc:`04_fed_batch_feeding`, just wired to this model's own parameter names.
The two plain :class:`~nskinetics.Event`\ s implement a **stage switch**:
whichever happens first — elapsed time reaching ``stage_1_max_time``, or
biomass ``x`` reaching ``stage_1_max_x`` — flips ``is_aerobic`` to ``0``,
switching the model from its aerobic growth phase into its anaerobic
production phase. Both ``stage_1_max_time`` and ``stage_1_max_x`` default to
infinity in the shipped SBML, so unless a caller sets them to finite values
the stage switch simply never fires and the run stays in the aerobic phase
throughout — that is the case exercised below.

Tolerances and initial conditions — strictly after compiling
----------------------------------------------------------------

As :doc:`03_events` warns, ``compile_events()`` regenerates the underlying
RoadRunner model and resets it to its SBML-defined origin values. Integrator
tolerances and initial conditions must therefore be set **after** it, never
before. Rather than poking the raw RoadRunner object directly, this page
attaches a custom ``reset_func`` to ``km`` — the same mechanism
``KineticModel(..., reset=...)`` accepts at construction time — so
that every ``km.reset()`` call, including ones a process unit triggers
internally, re-seeds the fed-batch bookkeeping consistently.

Which variables make up that bookkeeping? There is no need to guess from
their names. :attr:`~nskinetics.KineticModel.parameters` lists every global
parameter with its value, origin value and role, read from the model's own
structure: its rules, events and rate laws. The parameters the events
assign are the run-to-run state that a reset has to re-seed:

.. code-block:: python

   print(km.parameters.query('assigned_by_events').index.tolist())

.. code-block:: text

   ['is_aerobic', 'tot_vol_glu_feed_added', 'n_glu_spikes', 'last_vol_glu_feed_added']

These event-assigned parameters, together with the working volume ``env``
(a compartment) and the model's two rate-rule variables, are exactly what a
baseline reset restores to their origin values;
:attr:`~nskinetics.KineticModel.baseline_variables` lists them:

.. code-block:: python

   print(km.baseline_variables)

.. code-block:: text

   ('env', 'X_a', 'X_AcDH', 'last_vol_glu_feed_added', 'tot_vol_glu_feed_added', 'n_glu_spikes', 'is_aerobic')

So the reset below starts from the baseline reset instead of listing those
ids by hand, and adds only what the baseline does not cover:

.. code-block:: python

   integ = r.getIntegrator()
   integ.absolute_tolerance = 1e-10
   integ.relative_tolerance = 1e-9

   # Give the kinetic model its own reset: the baseline reset plus the spike
   # cap. Attaching it to `km` means every `km.reset()` (including the ones a
   # process unit calls internally) runs it.
   def reset_model(km, reset_spike_cap=True):
       km.reset(mode='baseline')   # volume, counters, stage flag, species
       if reset_spike_cap:
           km._te['max_n_glu_spikes'] = 20

   km.reset_func = reset_model
   km.reset()
   r['[s_glu]'] = 140; r['[x]'] = 1
   r['threshold_conc_glu_spike'] = 100   # refill trigger (SBML default: 10 g/L)
   r['target_conc_glu_spike'] = 140      # refill target  (SBML default: 100 g/L)

A custom ``reset_func`` **replaces** the default reset rather than layering
on top of it, so it must restore the model's state itself — otherwise the
model would start wherever the previous run left it. There are two ways to
do that. Start with ``km.reset(mode='baseline')``, as ``reset_model`` does.
Or call ``km._te.reset()`` with the volumes and event-assigned variables set
*before* that call as well as after it: RoadRunner's reset converts initial
concentrations to amounts at the *current* ``env`` volume, so resetting at
the previous run's end volume would start the species from values that
depend on that run. The shipped model's builder module
(``reset_nsk_kinetic_model``) uses that hand-written form; the two are
equivalent for this model (see :doc:`../concepts`, "Writing a reset
function").

Every write here uses item access, ``km._te['max_n_glu_spikes'] = 20`` (or
``km.set_value('max_n_glu_spikes', 20)``), never attribute access such as
``km._te.max_n_glu_spikes = 20``. libroadrunner removes the ``r.<id>``
attribute properties whenever another model is loaded, and from then on an
attribute write silently creates a Python attribute the model never sees
(see :doc:`../concepts`, "Writing a reset function").

The default reset of a model *without* a ``reset_func`` is that same
baseline reset: it already restores the working volume and every
event-assigned variable to its origin value (see :doc:`../concepts`).
``reset_model`` exists for what it adds on top: the spike cap
``max_n_glu_spikes``, which is not event-assigned, set to 20 (the SBML says
5). Keyword arguments passed to ``km.reset(...)`` reach the ``reset_func``.
``reset_spike_cap=False`` keeps whatever cap the model holds; a process
unit's spike-count retry relies on this to keep a reduced cap while it
re-simulates.

The initial concentrations ``s_glu`` and ``x`` are set *after*
``km.reset()`` returns, not inside ``reset_model`` and not before it:
``reset()`` restores every floating species to its SBML-defined origin value,
so setting them any earlier would simply be overwritten.

The two feed parameters are overridden the same way, immediately after
``km.reset()``: this run widens the glucose band from the shipped SBML's
defaults (refill at ``threshold_conc_glu_spike`` = 10 g/L, top up to
``target_conc_glu_spike`` = 100 g/L) to **refill at 100 g/L, top up to
140 g/L**. Note that ``s_glu`` also *starts* at the 140 g/L target rather
than at the 100 g/L threshold — this matters because the ``FeedSpike``
trigger fires only on a *downward* crossing of the threshold. Starting glucose
exactly at (or below) the threshold means the trigger is already satisfied at
``t=0`` and never transitions false-to-true, so no spike would ever fire;
starting above it lets glucose fall through 100 g/L and engage the feed.

The tightened tolerances (``1e-10``/``1e-9`` versus RoadRunner's looser
defaults) matter here because the model's rate laws are numerically stiff
around the repeated glucose-spike discontinuities; the isobutanol
biorefinery's own setup uses the same values.

Simulating 200 hours
----------------------

Requesting the bracketed (concentration) selectors — ``'[s_glu]'`` rather
than plain ``'s_glu'`` — matters here just as it did in
:doc:`04_fed_batch_feeding`: a plain species name returns the raw *amount*
(mass in the current, growing ``env`` compartment), not the g/L
concentration this page has been narrating throughout. Since this run's
``env`` grows from ``1`` to well over ``1`` as fed-batch feed is added,
plain ``s_glu``/``s_EtOH``/``s_IBO``/``x`` would silently mix amount and
concentration — bracket every species selector to stay in g/L:

.. code-block:: python

   cols = ['time', '[s_glu]', '[s_EtOH]', '[s_IBO]', '[x]', 'n_glu_spikes', 'f_O2']
   res = km.simulate(0, 200, 2001, cols)
   print('final:', dict(zip(cols, res[-1])))

(``n_glu_spikes`` stays unbracketed — it is a bookkeeping counter, not a
species, so it has no amount/concentration distinction.)

(``f_O2`` is the model's oxygen-transfer factor — 1 while the vessel's
``kLa``-bounded O2 transfer keeps up with the aerobic fluxes' demand, below 1
while it throttles them; see
:doc:`the shipped-model API page </API/models>`.)

Running the full snippet above end to end prints:

.. code-block:: text

   final: {'time': 200.0, '[s_glu]': 124.84688365540983, '[s_EtOH]': 106.52919734291154, '[s_IBO]': 0.0, '[x]': 12.899820003343597, 'n_glu_spikes': 8.0, 'f_O2': 1.0}

.. code-block:: python

   km.plot_simulation_results(
       variables=['[s_glu]', '[s_EtOH]', '[s_IBO]', '[x]'],
       labels=['glucose', 'ethanol', 'isobutanol', 'biomass'])

.. figure:: /_static/images/examples/tutorial_05_real_model.png
   :width: 450

   Eight glucose spikes hold glucose in the 100–140 g/L band until feeding
   tapers off past ~47 h; ethanol and biomass climb while the isobutanol
   branch stays at 0 (its rate constants are 0 in the shipped baseline).

Reading the final state
-------------------------

``n_glu_spikes: 8.0`` confirms the ``FeedSpike`` fired eight times over the
200 h run. Each spike is a discrete, exact jump: it snaps ``[s_glu]`` from
the ``threshold_conc_glu_spike`` trigger (100 g/L) back up toward
``target_conc_glu_spike`` (140 g/L) the instant it fires, and between spikes
— with this model's compartment volume not diluted by any outflow —
``[s_glu]`` can only decrease as it is consumed, never exceeding the 140 g/L
target it was just set to. Tracing the eight spikes confirms this pattern:
they land at roughly ``t=5.5``, ``7.2``, ``9.1``, ``11.4``, ``14.5``,
``19.2``, ``27.1``, and ``46.9`` h, each jumping ``[s_glu]`` from the 100 g/L
threshold back up to essentially the 140 g/L target. The spikes cluster early
and spread further apart later, as the volumetric glucose draw-down between
spikes slows; after the eighth spike, ``[s_glu]`` decays to a ~125 g/L
plateau within roughly 50 h of that spike and holds there for the rest of
the run, ending at ``124.8`` g/L — still well clear of the 100 g/L
threshold, so no ninth spike fires. ``[s_EtOH]`` climbs to ``106.5`` g/L
and biomass ``[x]`` grows almost thirteenfold (``1 → 12.9``), both driven
by the aerobic growth/fermentation rate laws active throughout this run
(recall from the previous section that ``stage_1_max_time``/``stage_1_max_x``
default to infinity, so the run never switches to the anaerobic stage).

This run also shows how the shipped model's oxygen-transfer bound behaves.
Aeration stays on throughout and biomass climbs past 10 g/L, yet ``f_O2``
reads 1 for all 200 h: the respiratory O2 demand of the aerobic fluxes peaks
at only ``17.5`` mmol O2/L/h (at ``t ≈ 9.0 h``, biomass ``11.6`` g/L), well
under the ``OTR_max`` = 40 mmol O2/L/h that the default ``kLa`` = 200 /h and
``C_O2_sat`` = 0.20 mmol/L allow, because the respiration and acetate-growth
rate laws are glucose-inhibited and this feeding strategy holds glucose at
100–140 g/L. The bound bites in glucose-lean, fully aerobic runs — a
100 g/L batch under the model's default 10 → 100 g/L refills, glucose-lean
while it binds, exceeds the cap near peak biomass (unbounded demand
~69 mmol O2/L/h; realized respiratory uptake held at the 40 mmol O2/L/h
cap, ``f_O2`` down to ~0.58), which is what the package's
``test_o2_transfer_bound`` tests pin — while the shipped fed-batch process
in :doc:`06_process_tea_bridge` is never bound (its aeration stops at 5 g/L
biomass). Lower ``r.kLa`` to describe a poorly aerated vessel and watch
``f_O2`` dip below 1 — in this run nothing happens until ``kLa`` falls
below about 88 /h (demand 17.5 mmol/L/h versus a cap of 0.2 × ``kLa``);
:doc:`the shipped-model API page </API/models>` documents the parameters.

``[s_IBO]`` stays at exactly ``0.0``. This is not a bug — it is the shipped
model's own default: the first three rate constants of the
``s_pyr -> s_AL -> s_DHI -> s_KIV -> s_isobutald -> s_IBO`` branch (``k_13``, ``k_14``,
``k_15``) are ``0.0`` in the SBML as shipped; only the constants of the last
two steps (``k_16``, ``k_17``) are nonzero. With no ALS (``k_13``) flux into
the branch, it carries zero flux regardless of how long the model runs. This
mirrors the :doc:`isobutanol biorefinery's <../case_studies/isobutanol/index>`
*baseline* scenario (scenario A), which runs this exact model with the
isobutanol branch switched off and activates it only in co-production
scenarios: scenario B, with literature-derived constants, and the ``opt_*``
strain-design scenarios, with optimized ones (see
:doc:`../case_studies/isobutanol/kinetic_optimization`). Nothing about the
fed-batch feeding or stage-switch machinery above changes between the two
cases — only the kinetic parameters do.

The same run as a ``Condition``
--------------------------------

The run above depends on four values set by hand after ``km.reset()``: the
two initial concentrations and the two feed parameters. The two parameters
also stay changed on ``km`` afterwards. A :class:`~nskinetics.Condition`
(see :doc:`02_build_and_simulate`) records all four with the run. Its
``reset_kwargs`` are passed to ``reset_model``. Its ``parameters`` are put
back when the run ends:

.. code-block:: python

   import numpy as np

   band = nsk.Condition(
       initial={'[s_glu]': 140, '[x]': 1},
       parameters={'threshold_conc_glu_spike': 100, 'target_conc_glu_spike': 140},
       t_end=200, n_steps=2001, selections=cols,
       reset_kwargs={'reset_spike_cap': True}, label='refill 100 -> 140 g/L')
   run1 = km.run(band)
   run2 = km.run(band)
   print(np.array_equal(run1.array, res), np.array_equal(run1.array, run2.array))

.. code-block:: text

   True True

Running the condition reproduces the hand-driven run bit for bit, and
running it twice gives the same result twice. To see the ``reset_kwargs``
at work, lower the cap to three spikes, as a spike-count retry would, and
run once with ``reset_spike_cap=False`` and once with the condition as
written:

.. code-block:: python

   km.set_value('max_n_glu_spikes', 3)
   kept = km.run(band.replace(reset_kwargs={'reset_spike_cap': False}))
   reseeded = km.run(band)
   print(kept.final['n_glu_spikes'], reseeded.final['n_glu_spikes'],
         km.get_value('max_n_glu_spikes'))

.. code-block:: text

   3.0 8.0 20.0

With ``reset_spike_cap=False``, ``reset_model`` leaves the lowered cap in
place, and the feed stops after three spikes. The condition as written
re-seeds the cap to 20 and gets all eight spikes back. A cap meant for
one run only belongs in the condition instead, for example
``band.with_parameters(max_n_glu_spikes=3)``, which puts the cap back when
the run ends.

Describing a parameter change
-----------------------------

When a scenario or an uncertainty run changes a set of kinetic parameters,
the shipped model can name the change in a few words. Every kinetic
parameter is classified by *role* (what kind of term it is in its rate law:
a capacity, a saturation constant, a product-inhibition coefficient, ...)
and by pathway *module* (glycolysis/fermentation, growth, the Ehrlich
branch, ...); :func:`~nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo.describe_parameter_change`
diffs a baseline snapshot against the current values and phrases each group
of changes as direction + role + effector + modules. Snapshot the baseline
first, then apply the change:

.. code-block:: python

   from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import (
       snapshot_parameters, describe_parameter_change, SCENARIO_B_EHRLICH)

   baseline = snapshot_parameters(km)          # km from the top of this page
   for name, value in SCENARIO_B_EHRLICH.items():
       km.set_value(name, value)
   print(describe_parameter_change(baseline, km))

.. code-block:: text

   Ehrlich branch on; faster Ehrlich branch

The preset turns on the engineered upstream capacities ``k_13`` .. ``k_15``
(0 -> nonzero -- a single ``on`` clause, since they share one role, effector
and module) and further up-regulates the native Aro10 capacity ``k_16`` from
its constitutive default to the overexpressed scenario-B level, a separate
``faster`` clause. (Aro10, like Adh6/``k_17``, is a native enzyme with a
nonzero default, so it reads as sped up rather than switched on; applying
scenario A first, which re-zeros ``k_16``, would instead collapse the whole
change to a single ``Ehrlich branch on``.) Note that the snippet mutates the
model in place: ``km`` keeps scenario B applied -- the Ehrlich branch on --
for the rest of the session. A partial override dict can also be passed
directly as the
"current" side, and ``verbose=True`` appends the parameter names and values
to each clause:

.. code-block:: python

   print(describe_parameter_change(
       baseline, {'k_7ie': 0.08, 'k_1ie': 0.04}, verbose=True))

.. code-block:: text

   stronger ethanol inhibition of glycolysis/fermentation and growth (k_1ie 0.02 -> 0.04, k_7ie 0.04 -> 0.08)

Operation parameters -- the fed-batch feeding thresholds, the stage-switch
cutoffs, the dilution rate -- are not kinetics and are ignored by the
helper; describe those changes in words yourself.

Next: :doc:`06_process_tea_bridge` takes a kinetic model like this one and
drives a biosteam process unit with it, connecting kinetics to process
design and TEA.
