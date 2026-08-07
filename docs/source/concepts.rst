Key Concepts
============

.. currentmodule:: nskinetics

This page explains the ideas behind NSKinetics' current, recommended
workflow: build a kinetic model in SBML/Antimony, drive it with Tellurium's
RoadRunner engine, and (optionally) plug it into a biosteam process model.
It assumes you have already worked through the :doc:`tutorial <tutorial/index>`.

SBML and Antimony
------------------

Kinetic models are declared as SBML, either loaded from an existing file or
authored in the more readable Antimony syntax and compiled to SBML under the
hood. To load a model from disk, call
:meth:`TelluriumReactionSystem.from_sbml(path) <TelluriumReactionSystem.from_sbml>`.
To build one inline, write Antimony text, compile it with
``tellurium.loadAntimonyModel(text)``, and wrap the resulting RoadRunner
object with :class:`nsk.TelluriumReactionSystem(r, units=...) <TelluriumReactionSystem>`.
Either route produces the same wrapper, so the rest of the workflow (units,
events, simulation) is identical regardless of how the model was authored.

Simulating and the RoadRunner engine
-------------------------------------

:class:`TelluriumReactionSystem` provides :meth:`~TelluriumReactionSystem.simulate`,
which integrates the model over ``[t0, t_end]`` and stores the trajectory on
the reaction system (``trs.results`` / ``.results_array`` / ``.results_dict``) —
the object is the source of truth for full-trajectory results. The underlying
RoadRunner engine remains directly accessible as ``trs._te`` for lower-level use.
:meth:`~TelluriumReactionSystem.plot_simulation_results` (aliases
``plot_time_course`` / ``plot_trajectory``) plots the most recent run. Use
bracketed ``'[S]'`` selections to record concentrations, and plain variable
names to record amounts.

Units
-----

Every :class:`TelluriumReactionSystem` carries a ``units`` dict with
``'time'`` and ``'conc'`` keys (default ``{'time': 'min', 'conc': 'M'}``).
Recognized time tokens are ``h``, ``hr``, ``min``, ``m``, ``s``, and ``sec``
(case-insensitive); recognized concentration tokens are ``M`` and ``mol/L``
for molar concentrations, and ``g/L``, ``kg/m3``, and ``kg/m^3`` for mass
concentrations — note that ``kg/m3`` is a *mass* unit, not molar, since it is
numerically equivalent to g/L. The :attr:`~TelluriumReactionSystem.time_factor`
property converts model time to hours, and
:attr:`~TelluriumReactionSystem.material_indexer` resolves the configured
concentration unit to the matching biosteam stream indexer (``'imol'`` for
molar, ``'imass'`` for mass) for use by the process bridge.
:meth:`~TelluriumReactionSystem.validate_units` raises a
``KineticSimulationError`` if either token is unrecognized.

The event lifecycle (important)
--------------------------------

The typical event workflow is: ``add_event(...)`` (repeatable) →
``compile_events()`` → set any custom initial conditions or parameters →
``reset()`` → simulate. **The ordering matters.** Compiling events
regenerates the underlying RoadRunner model and, as a side effect, resets all
model state to its SBML-defined origin values. Any custom initial conditions
or parameter overrides must therefore be applied *after* ``compile_events()``,
never before — values set beforehand are silently discarded when the model is
regenerated. This is documented directly on
:meth:`TelluriumReactionSystem.compile_events`, :meth:`Event.compile`, and
:meth:`FeedSpike.compile`; when in doubt, compile events first and set state
last.

``Event`` vs ``FeedSpike``
----------------------------

:class:`Event` represents a single SBML event: a trigger expression paired
with one or more assignments. :class:`FeedSpike` is a higher-level
convenience that expands into three coordinated :class:`Event` objects fired
in sequence: one computes the volume of feed needed to hit a target
concentration, one updates the working volume (and a running total), and one
sets the species to its target concentration and increments a spike counter.
An optional ``max_count`` caps how many times a given ``FeedSpike`` may fire,
by folding a count comparison into every trigger it expands to.

The process bridge
-------------------

:class:`~nskinetics.units.NSKBatchReactor` drives any
:class:`TelluriumReactionSystem` inside a biosteam ``BatchBioreactor``,
handling the kinetic simulate loop, reaction-time selection, species-to-
chemical mapping, and effluent construction so the model can participate in
design, costing, and TEA. :class:`~nskinetics.units.NSKFermentation` is a
fed-batch fermentation subclass of it that adds aeration, sucrose hydrolysis,
fed-batch spike-count retry, and yield/mass-balance validators on top of the
shared machinery.

.. note::
   The string-based ``ReactionSystem`` / ``SpeciesSystem`` / ``Reaction`` API
   is legacy and intentionally undocumented; new work should use the
   Tellurium/SBML workflow described here.
