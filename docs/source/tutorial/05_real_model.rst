A real fed-batch model
=======================

:doc:`03_events` and :doc:`04_fed_batch_feeding` built ``Event`` and
``FeedSpike`` from small, hand-written Antimony models. This page runs the
same machinery on a **shipped, production-scale** model: the *S. cerevisiae*
fed-batch fermentation SBML that ships inside the NSKinetics package and
backs the downstream ``isobutanol`` biorefinery's fermentation unit. It has
ten floating species (glucose, ethanol, isobutanol, biomass, and six
pathway intermediates) and dozens of kinetic parameters — the same shape of
model as :doc:`04_fed_batch_feeding`'s toy example, just larger.

Loading the shipped SBML
-------------------------

The file lives under ``nskinetics/examples/``, next to the package's other
example models:

.. code-block:: python

   import os, numpy as np, tellurium as te, nskinetics as nsk

   p = os.path.join(os.path.dirname(nsk.__file__),
                     'examples', 's_cerevisiae_ferm_fb_inhib_mod_ibo_sbml.xml')
   r = te.loadSBMLModel(p)
   trs = nsk.TelluriumReactionSystem(r, units={'time': 'h', 'conc': 'g/L'})

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
downstream biorefinery's own setup script does.

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
           trs.add_event(e)
       trs.add_event(nsk.Event(when='time >= stage_1_max_time',
                               do={'is_aerobic': '0'}, name='stage_1_complete_max_time'))
       trs.add_event(nsk.Event(when='x >= stage_1_max_x',
                               do={'is_aerobic': '0'}, name='stage_1_complete_x_target'))
       trs.compile_events()

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
before:

.. code-block:: python

   integ = r.getIntegrator()
   integ.absolute_tolerance = 1e-10
   integ.relative_tolerance = 1e-9
   trs.reset()
   r.n_glu_spikes = 0; r.last_vol_glu_feed_added = 0.; r.tot_vol_glu_feed_added = 0.
   r.env = 1.; r.is_aerobic = 1; r.max_n_glu_spikes = 20; r.s_glu = 100; r.x = 1

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

   cols = ['time', '[s_glu]', '[s_EtOH]', '[s_IBO]', '[x]', 'n_glu_spikes']
   res = trs.simulate(0, 200, 2001, cols)
   print('final:', dict(zip(cols, res[-1])))

(``n_glu_spikes`` stays unbracketed — it is a bookkeeping counter, not a
species, so it has no amount/concentration distinction.)

Running the full snippet above end-to-end in the ``HP_2024`` environment
prints:

.. code-block:: text

   final: {'time': 200.0, '[s_glu]': 98.63213916292953, '[s_EtOH]': 105.91411061058339,
           '[s_IBO]': 0.0, '[x]': 10.553862807704611, 'n_glu_spikes': 4.0}

.. code-block:: python

   trs.plot_simulation_results(
       variables=['[s_glu]', '[s_EtOH]', '[s_IBO]', '[x]'],
       labels=['glucose', 'ethanol', 'isobutanol', 'biomass'],
       markers=[7.5, 11.9, 21.0, 94.0])

.. figure:: /_static/images/examples/tutorial_05_real_model.png
   :width: 450

   Four glucose spikes (dashed lines) over 200 h; ethanol and biomass climb
   while the isobutanol branch stays at 0 (its rate constants are 0 in the
   shipped baseline).

Reading the final state
-------------------------

``n_glu_spikes: 4.0`` confirms the ``FeedSpike`` fired four times over the
200 h run. Each spike is a discrete, exact jump: it snaps ``[s_glu]`` from
the ``threshold_conc_glu_spike`` trigger (10 g/L) back up toward
``target_conc_glu_spike`` (100 g/L) the instant it fires, and between spikes
— with this model's compartment volume not diluted by any outflow —
``[s_glu]`` can only decrease as it is consumed, never exceeding the 100 g/L
target it was just set to. Tracing the four spikes confirms this pattern:
they land at roughly ``t=7.4``, ``11.9``, ``21.0``, and ``94.0`` h, each
jumping ``[s_glu]`` from just above the 10 g/L threshold up to
essentially the 100 g/L target; after the fourth spike, ``[s_glu]`` simply
decays for the remaining ~106 h, reaching the final ``98.63`` g/L — just
below the 100 g/L target, consistent with slow consumption rather than any
further feeding. ``[s_EtOH]`` climbs to ``105.9`` g/L and biomass ``[x]``
grows tenfold (``1 → 10.6``), both driven by the aerobic growth/fermentation
rate laws active throughout this run (recall from the previous section that
``stage_1_max_time``/``stage_1_max_x`` default to infinity, so the run never
switches to the anaerobic stage).

``[s_IBO]`` stays at exactly ``0.0``. This is not a bug — it is the shipped
model's own default: the four isobutanol-pathway rate constants
(``k_13``, ``k_14``, ``k_15``, ``k_16``, governing the
``s_pyr -> s_AL -> s_DHI -> s_KIV -> s_IBO`` branch) are all ``0.0`` in the
SBML as shipped, so that branch carries zero flux regardless of how long the
model runs. This mirrors the downstream ``isobutanol`` biorefinery's own
*baseline* scenario, which runs this exact model with the isobutanol branch
switched off and only activates it (by setting ``k_13`` and friends to
literature-derived values) in scenarios explicitly exploring isobutanol
co-production. Nothing about the fed-batch feeding or stage-switch machinery
above changes between the two cases — only the kinetic parameters do.

Next: :doc:`06_process_tea_bridge` takes a kinetic model like this one and
drives a biosteam process unit with it, connecting kinetics to process
design and TEA.
