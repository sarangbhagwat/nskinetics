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
     compartment env; species s_glu in env, p_eth in env;
     s_glu = 100; p_eth = 0; env = 1; k = 1;
     threshold = 10; target = 100; feed_conc = 600;
     n_spk = 0; n_max = 2; last_vol = 0; tot_vol = 0; dly = 0;
     n_spk has dimensionless; n_max has dimensionless;
     J: s_glu => p_eth; k*s_glu*env;
   end
   """
   r = te.loadAntimonyModel(model)
   km = nsk.KineticModel(r, units={'time': 'h', 'conc': 'g/L'})
   fs = nsk.FeedSpike(species='s_glu', when='s_glu <= threshold',
                      target='target', feed_conc='feed_conc', volume_var='env',
                      max_count='n_max', count_var='n_spk',
                      last_vol_var='last_vol', tot_vol_var='tot_vol',
                      delay='dly', priority=5, name='spk')
   for e in fs.expand():
       km.add_event(e)
   km.compile_events()
   km.reset()
   res = km.simulate(0, 40, 401, ['time', '[s_glu]', '[p_eth]', 'n_spk'])
   print('max spikes:', res[:, 3].max(), 'final s_glu:', res[-1, 1],
         'final p_eth:', res[-1, 2])
   km.plot_simulation_results(variables=['[s_glu]', '[p_eth]'],
                               labels=['s_glu', 'p_eth'],
                               markers=[(2.3, 'spike'), (4.6, 'spike')])

The model starts ``s_glu`` at 100 g/L and converts it to product ``p_eth``
through the mass-action reaction ``J: s_glu => p_eth`` at ``k=1``. Writing
it as a *reaction* rather than a rate rule is what makes the 1:1
stoichiometry explicit: every gram of substrate consumed appears as a gram
of product, so ``p_eth`` is a running record of everything the spikes have
fed in. ``s_glu`` will cross the ``threshold`` of 10 g/L repeatedly over
40 h, but ``max_count=2`` (``n_max``) caps it at two spikes. Note that
``max_count``/``count_var`` are passed as the *model-parameter names*
``'n_max'``/``'n_spk'`` here, not Python numbers — ``FeedSpike`` folds them
into the trigger expression as text, so a plain float works equally well.

Running this in the ``HP_2024`` environment prints:

.. code-block:: text

   max spikes: 2.0 final s_glu: -5.060339805135531e-11 final p_eth: 240.90778511926948

.. figure:: /_static/images/examples/tutorial_04_fed_batch_i.png
   :width: 400

   Two feed spikes (dashed lines) refill ``s_glu`` to ~100 g/L, each one
   feeding another batch of substrate into ``p_eth``; after the ``n_max=2``
   cap the remaining ``s_glu`` is converted and ``p_eth`` levels off at
   ~241 g/L. The small drop in ``p_eth`` at each spike is dilution by the
   feed volume, not product disappearing.

``max spikes: 2.0`` confirms the cap was reached — exactly two spikes
fired, matching ``n_max=2``. After the second spike, no further trigger
fires (the ``count_var < max_count`` guard is false), so the remaining
``s_glu`` is simply converted to ``p_eth`` over the rest of the 40 h
window; by ``t=40`` the substrate is down to ``-5.06e-11`` g/L — zero
within floating-point noise, and well below the 10 g/L ``threshold`` it
would otherwise have re-triggered a spike at, had the cap not been reached.

The final ``p_eth`` of ``~240.9`` g/L is the whole fed-batch mass balance
in one number. Three charges of substrate enter the reactor: the initial
100 g in 1 L, then ``0.18 L × 600 g/L = 108 g`` at the first spike and
``0.212 L × 600 g/L = 127.4 g`` at the second — 335.4 g of substrate
total, ending in a working volume of 1.392 L. With 1:1 stoichiometry and
essentially complete conversion, that is
``335.4 / 1.392 ≈ 240.9`` g/L of product, exactly what the simulation
reports. Note that this is *higher* than the 100 g/L ``target``
concentration ever reached: product accumulates across spikes while
substrate is held near a fixed setpoint, which is the entire reason
fed-batch operation is used.

Seeing ``a``/``b``/``c`` fire
--------------------------------

Requesting the bracketed (concentration) selectors ``'[s_glu]'`` and
``'[p_eth]'`` alongside the working volume ``'env'`` at a finer time
resolution shows the mechanism directly — each spike snaps the substrate
concentration back up near ``target`` while ``env`` grows by the volume
``a`` computed, and ``[p_eth]`` dips by exactly the corresponding dilution
factor. This rebuilds the same setup from scratch and simulates it
**once**, rather than re-simulating the ``r``/``km`` from the block above
(see the note below for why):

.. code-block:: python

   import numpy as np, tellurium as te, nskinetics as nsk

   model = """
   model spiker()
     compartment env; species s_glu in env, p_eth in env;
     s_glu = 100; p_eth = 0; env = 1; k = 1;
     threshold = 10; target = 100; feed_conc = 600;
     n_spk = 0; n_max = 2; last_vol = 0; tot_vol = 0; dly = 0;
     n_spk has dimensionless; n_max has dimensionless;
     J: s_glu => p_eth; k*s_glu*env;
   end
   """
   r = te.loadAntimonyModel(model)
   km = nsk.KineticModel(r, units={'time': 'h', 'conc': 'g/L'})
   fs = nsk.FeedSpike(species='s_glu', when='s_glu <= threshold',
                      target='target', feed_conc='feed_conc', volume_var='env',
                      max_count='n_max', count_var='n_spk',
                      last_vol_var='last_vol', tot_vol_var='tot_vol',
                      delay='dly', priority=5, name='spk')
   for e in fs.expand():
       km.add_event(e)
   km.compile_events()
   km.reset()
   res = km.simulate(0, 40, 4001, ['time', '[s_glu]', '[p_eth]', 'env', 'n_spk'])
   changes = np.where(np.diff(res[:, 4]) > 0)[0]
   for i in changes:
       print(res[i], '->', res[i + 1])

prints:

.. code-block:: text

   [ 2.3   10.026 89.974  1.     0.   ] -> [ 2.31  99.262 77.01   1.18   1.   ]
   [  4.6    10.052 166.219   1.18    1.   ] -> [  4.61   99.519 141.389   1.392   2.   ]

The two concentrations go on one plot, since they share a y-axis:

.. code-block:: python

   km.plot_simulation_results(variables=['[s_glu]', '[p_eth]'],
                               labels=['[s_glu]', '[p_eth]'],
                               markers=[(2.3, 'spike'), (4.6, 'spike')])

.. figure:: /_static/images/examples/tutorial_04_fed_batch_ii.png
   :width: 400

   Each spike snaps ``[s_glu]`` back up to ~100 g/L and knocks ``[p_eth]``
   down by the dilution factor; ``[p_eth]`` then climbs to ~241 g/L as the
   last charge of substrate is converted.

``env`` needs its own plot. It is a *volume*, not a concentration, so it
does not belong on the axis above — and on a 0–250 g/L scale its 1.0 → 1.39
range would be flattened into an unreadable line at the axis floor.
``ylabel``/``ylabel_unit`` relabel the y-axis for exactly this case:

.. code-block:: python

   km.plot_simulation_results(variables=['env'], labels=['env'],
                               ylabel='Working volume', ylabel_unit='L',
                               markers=[(2.3, 'spike'), (4.6, 'spike')])

.. figure:: /_static/images/examples/tutorial_04_fed_batch_iii.png
   :width: 400

   On its own axis the working volume is a clean staircase: ``env`` holds at
   1.0 L, then steps to 1.18 L and 1.392 L as each spike's event **b** adds
   the feed volume that event **a** computed. It is flat everywhere else —
   volume changes only at a spike.

.. warning::

   Reusing the *same* compiled ``r``/``km`` for a second ``simulate()``
   call — even after calling ``reset()`` in between — is not reliable, and
   this is independent of point count.
   :meth:`~nskinetics.KineticModel.reset` calls a bare
   ``self._te.reset()``, which only resets RoadRunner's default selection
   (time, floating species, and rate-rule-governed variables) to their
   origin values. It does **not** restore compartments (like ``env`` here)
   or plain, event-only bookkeeping parameters that carry no rate rule
   (``n_spk``, ``last_vol``, ``tot_vol`` here) — once an event has changed
   them, ``reset()`` leaves them exactly as the event left them. Re-running
   :meth:`~nskinetics.KineticModel.compile_events` does not help
   either: it is guarded (``if self._events_compiled: return``, see
   :doc:`03_events`), so a second call after the first succeeded is a
   no-op. The reliable way to run an independent simulation is to build a
   fresh model and a fresh ``KineticModel``, as the standalone
   script above does, rather than re-simulating an already-fired
   event-bearing model in place.

At the first spike (between ``t=2.3`` and ``t=2.31`` h), ``[s_glu]`` jumps
from just under the 10 g/L threshold back up to about ``99.26`` g/L —
close to ``target=100`` (the small gap is because the trigger actually
fires a moment before ``t=2.31``, and ``s_glu`` had already been consumed a
little further by the time this row was recorded) — while ``env`` grows
from ``1.0`` to ``1.18``, matching ``a``'s formula:
``1*(100 - 10.026)/(600 - 100) ≈ 0.18``. The second spike repeats the same
pattern from a slightly different starting concentration, growing ``env``
from ``1.18`` to about ``1.392``. No third spike occurs, consistent with
``max_count=2``.

``[p_eth]`` makes the volume update visible from the other side. Event
**c** assigns only ``s_glu``; ``p_eth`` is never touched by any of the
three events, so its *amount* carries straight through the spike and only
its *concentration* changes — ``89.974`` g/L drops to ``77.01`` g/L across
the first spike as the same product mass is spread over the now-larger
working volume. This is the clearest confirmation that **b** really did
grow ``volume_var`` before **c** ran: had the events fired in the other
order, ``s_glu`` would have been set to ``target`` in the *old* volume and
then diluted along with ``p_eth``, landing well short of 100 g/L. Checking
total mass at ``t=2.31`` confirms the bookkeeping —
``(99.262 + 77.01) g/L × 1.18 L ≈ 208 g``, exactly the 100 g charged
initially plus the ``0.18 L × 600 g/L = 108 g`` the spike added.

Next: :doc:`05_real_model` puts this together in a real, shipped fed-batch
*S. cerevisiae* fermentation model.
