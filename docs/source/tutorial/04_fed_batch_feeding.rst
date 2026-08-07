Fed-batch feeding
=================

:doc:`03_events` covered :class:`~nskinetics.Event`, a single trigger +
assignment pair. Fed-batch feeding — topping a species back up whenever it
runs low, while tracking the volume added — needs three events kept in
sync, not one. :class:`~nskinetics.FeedSpike` is that pattern, pre-built.

What ``FeedSpike`` represents
--------------------------------

A ``FeedSpike`` watches one species and, whenever it drops to a trigger
condition, adds enough feed (at a known feed concentration) to bring it
back up to a target concentration — while growing the working volume by
the amount of feed added. The constructor parameters:

- ``species`` — the model variable being spiked (e.g. ``'s_glu'``).
- ``when`` — the base trigger expression (e.g. ``'s_glu <= threshold'``).
- ``target`` — target concentration after the spike; a number or
  model-parameter name.
- ``feed_conc`` — concentration of the feed being added; a number or
  parameter name.
- ``volume_var`` — the working-volume variable the reactor reads (default
  ``'env'``).
- ``max_count`` — optional cap on the number of spikes; a number or
  parameter name. When given, ``&& (count_var < max_count)`` is folded into
  every trigger, so once the count is reached the spike simply stops firing
  and the species is left to evolve under whatever rate law governs it.
- ``count_var``, ``last_vol_var``, ``tot_vol_var`` — bookkeeping variables
  for the spike count, the volume added by the last spike, and the running
  total volume added. Each defaults to a name derived from ``species``
  (e.g. ``f'{species}_spike_count'``) if not given explicitly.
- ``delay``, ``priority`` — as for :class:`~nskinetics.Event`; ``priority``
  is the priority of the *first* of the three expanded events (default
  ``5``) — the other two are offset below it to enforce firing order (see
  next section).
- ``name`` — base event name; the three expanded events are named
  ``{name}_a``, ``{name}_b``, ``{name}_c``.

The three-event expansion
----------------------------

``FeedSpike`` itself is not an SBML event — :meth:`~nskinetics.FeedSpike.expand`
turns it into a list of three coordinated :class:`~nskinetics.Event`
objects, all sharing the same trigger, that must fire in a fixed order.
Reading directly off ``FeedSpike.expand()``:

- **a** (highest priority, fires first) computes the volume of feed needed
  to hit the target and stores it in ``last_vol_var``, using the
  species/volume values *at trigger time*
  (``use_trigger_time_values=True``):
  ``last_vol := volume_var*(target - species)/(feed_conc - target)``.
- **b** (priority − 1) applies that volume: it adds ``last_vol`` to
  ``volume_var`` (the working volume grows) and to ``tot_vol_var`` (the
  running total).
- **c** (priority − 2, fires last) sets ``species`` to ``target`` and
  increments ``count_var`` by 1.

The ordering matters: **b** must update the volume *before* **c** sets the
species value, since the species is a concentration inside the
now-larger volume. SBML fires simultaneously-triggered events in priority
order, which is exactly why ``a``/``b``/``c`` are given descending
priorities (``priority``, ``priority-1``, ``priority-2``) rather than being
left to fire in an arbitrary order.

``max_count`` works by folding a count comparison into the *shared*
trigger passed to all three events — once ``count_var`` reaches
``max_count``, none of ``a``/``b``/``c`` fire again, so the species is left
to whatever the model's own rate law does next.

Compiling a ``FeedSpike``
----------------------------

Because ``expand()`` returns plain ``Event`` objects, a ``FeedSpike`` plugs
directly into the same ``add_event``/``compile_events`` workflow from
:doc:`03_events` — add all three expanded events, then compile once:

.. code-block:: python

   import tellurium as te, nskinetics as nsk

   model = """
   model spiker()
     compartment env; species s_glu in env;
     s_glu = 100; env = 1; k = 1;
     threshold = 10; target = 100; feed_conc = 600;
     n_spk = 0; n_max = 2; last_vol = 0; tot_vol = 0; dly = 0;
     n_spk has dimensionless; n_max has dimensionless;
     s_glu' = -k*s_glu;
   end
   """
   r = te.loadAntimonyModel(model)
   trs = nsk.TelluriumReactionSystem(r, units={'time': 'h', 'conc': 'g/L'})
   fs = nsk.FeedSpike(species='s_glu', when='s_glu <= threshold',
                      target='target', feed_conc='feed_conc', volume_var='env',
                      max_count='n_max', count_var='n_spk',
                      last_vol_var='last_vol', tot_vol_var='tot_vol',
                      delay='dly', priority=5, name='spk')
   for e in fs.expand():
       trs.add_event(e)
   trs.compile_events()
   trs.reset()
   res = trs.simulate(0, 40, 401, ['time', '[s_glu]', 'n_spk'])
   print('max spikes:', res[:, 2].max(), 'final s_glu:', res[-1, 1])
   trs.plot_simulation_results(labels=['s_glu'],
                               markers=[(2.3, 'spike'), (4.6, 'spike')])

The model starts ``s_glu`` at 100 g/L, decaying at ``k=1``; it will cross
the ``threshold`` of 10 g/L repeatedly over 40 h, but ``max_count=2``
(``n_max``) caps it at two spikes. Note that ``max_count``/``count_var``
are passed as the *model-parameter names* ``'n_max'``/``'n_spk'`` here, not
Python numbers — ``FeedSpike`` folds them into the trigger expression as
text, so a plain float works equally well.

Running this in the ``HP_2024`` environment prints:

.. code-block:: text

   max spikes: 2.0 final s_glu: -6.19487843009339e-12

.. figure:: /_static/images/examples/tutorial_04_fed_batch_i.png
   :width: 400

   Two feed spikes (dashed lines) refill ``s_glu`` to ~100 g/L; after the
   ``n_max=2`` cap it decays to ~0.

``max spikes: 2.0`` confirms the cap was reached — exactly two spikes
fired, matching ``n_max=2``. After the second spike, no further trigger
fires (the ``count_var < max_count`` guard is false), so ``s_glu`` decays
freely under ``s_glu' = -k*s_glu`` for the rest of the 40 h window; by
``t=40`` it has decayed to ``-6.19e-12`` g/L — zero within floating-point
noise, and well below the 10 g/L ``threshold`` it would otherwise have
re-triggered a spike at, had the cap not been reached.

Seeing ``a``/``b``/``c`` fire
--------------------------------

Requesting the bracketed (concentration) selector ``'[s_glu]'`` alongside
the working volume ``'env'`` at a finer time resolution shows the
mechanism directly — each spike snaps the concentration back up near
``target`` while ``env`` grows by the volume ``a`` computed. This rebuilds
the same setup from scratch and simulates it **once**, rather than
re-simulating the ``r``/``trs`` from the block above (see the note below
for why):

.. code-block:: python

   import numpy as np, tellurium as te, nskinetics as nsk

   model = """
   model spiker()
     compartment env; species s_glu in env;
     s_glu = 100; env = 1; k = 1;
     threshold = 10; target = 100; feed_conc = 600;
     n_spk = 0; n_max = 2; last_vol = 0; tot_vol = 0; dly = 0;
     n_spk has dimensionless; n_max has dimensionless;
     s_glu' = -k*s_glu;
   end
   """
   r = te.loadAntimonyModel(model)
   trs = nsk.TelluriumReactionSystem(r, units={'time': 'h', 'conc': 'g/L'})
   fs = nsk.FeedSpike(species='s_glu', when='s_glu <= threshold',
                      target='target', feed_conc='feed_conc', volume_var='env',
                      max_count='n_max', count_var='n_spk',
                      last_vol_var='last_vol', tot_vol_var='tot_vol',
                      delay='dly', priority=5, name='spk')
   for e in fs.expand():
       trs.add_event(e)
   trs.compile_events()
   trs.reset()
   res = trs.simulate(0, 40, 4001, ['time', '[s_glu]', 'env', 'n_spk'])
   changes = np.where(np.diff(res[:, 3]) > 0)[0]
   for i in changes:
       print(res[i], '->', res[i + 1])

prints:

.. code-block:: text

   [2.3   10.026  1.     0.   ] -> [2.31 99.26  1.18  1.  ]
   [4.6   10.052  1.18   1.   ] -> [4.61  99.515  1.392  2.   ]

.. code-block:: python

   trs.plot_simulation_results(variables=['[s_glu]', 'env'],
                               labels=['[s_glu]', 'env'],
                               markers=[(2.3, 'spike'), (4.6, 'spike')])

.. figure:: /_static/images/examples/tutorial_04_fed_batch_ii.png
   :width: 400

   Each spike snaps ``[s_glu]`` back up while the working volume ``env`` steps
   up by the feed added (1.0 → 1.18 → 1.39).

.. warning::

   Reusing the *same* compiled ``r``/``trs`` for a second ``simulate()``
   call — even after calling ``reset()`` in between — is not reliable, and
   this is independent of point count.
   :meth:`~nskinetics.TelluriumReactionSystem.reset` calls a bare
   ``self._te.reset()``, which only resets RoadRunner's default selection
   (time, floating species, and rate-rule-governed variables) to their
   origin values. It does **not** restore compartments (like ``env`` here)
   or plain, event-only bookkeeping parameters that carry no rate rule
   (``n_spk``, ``last_vol``, ``tot_vol`` here) — once an event has changed
   them, ``reset()`` leaves them exactly as the event left them. Re-running
   :meth:`~nskinetics.TelluriumReactionSystem.compile_events` does not help
   either: it is guarded (``if self._events_compiled: return``, see
   :doc:`03_events`), so a second call after the first succeeded is a
   no-op. The reliable way to run an independent simulation is to build a
   fresh model and a fresh ``TelluriumReactionSystem``, as the standalone
   script above does, rather than re-simulating an already-fired
   event-bearing model in place.

At the first spike (between ``t=2.3`` and ``t=2.31`` h), ``[s_glu]`` jumps
from just under the 10 g/L threshold back up to about ``99.26`` g/L —
close to ``target=100`` (the small gap is because the trigger actually
fires a moment before ``t=2.31``, and ``s_glu`` had already decayed a
little further by the time this row was recorded) — while ``env`` grows
from ``1.0`` to ``1.18``, matching ``a``'s formula:
``1*(100 - 10.026)/(600 - 100) ≈ 0.18``. The second spike repeats the same
pattern from a slightly different starting concentration, growing ``env``
from ``1.18`` to about ``1.392``. No third spike occurs, consistent with
``max_count=2``.

Next: :doc:`05_real_model` puts this together in a real, shipped fed-batch
*S. cerevisiae* fermentation model.
