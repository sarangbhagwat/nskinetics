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
:meth:`KineticModel.from_sbml(path) <KineticModel.from_sbml>`.
To build one inline, write Antimony text, compile it with
``tellurium.loadAntimonyModel(text)``, and wrap the resulting RoadRunner
object with :class:`nsk.KineticModel(r, units=...) <KineticModel>`.
Either route produces the same wrapper, so the rest of the workflow (units,
events, simulation) is identical regardless of how the model was authored.

Simulating and the RoadRunner engine
-------------------------------------

:class:`KineticModel` provides :meth:`~KineticModel.simulate`,
which integrates the model over ``[t0, t_end]`` and stores the trajectory on
the kinetic model (``km.results_df`` / ``.results_array`` / ``.results_dict``) —
the object is the source of truth for full-trajectory results. The underlying
RoadRunner engine remains directly accessible as ``km._te`` for lower-level use.
:meth:`~KineticModel.plot_simulation_results` (aliases
``plot_time_course`` / ``plot_trajectory``) plots the most recent run. Use
bracketed ``'[S]'`` selections to record concentrations, and plain variable
names to record amounts.

Units
-----

Every :class:`KineticModel` carries a ``units`` dict with
``'time'`` and ``'conc'`` keys (default ``{'time': 'min', 'conc': 'M'}``).
Recognized time tokens are ``h``, ``hr``, ``min``, ``m``, ``s``, and ``sec``
(case-insensitive); recognized concentration tokens are ``M`` and ``mol/L``
for molar concentrations, and ``g/L``, ``kg/m3``, and ``kg/m^3`` for mass
concentrations — note that ``kg/m3`` is a *mass* unit, not molar, since it is
numerically equivalent to g/L. The :attr:`~KineticModel.time_factor`
property converts model time to hours, and
:attr:`~KineticModel.material_indexer` resolves the configured
concentration unit to the matching biosteam stream indexer (``'imol'`` for
molar, ``'imass'`` for mass) for use by the process bridge.
:meth:`~KineticModel.validate_units` raises a
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
:meth:`KineticModel.compile_events`, :meth:`Event.compile`, and
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
:class:`KineticModel` inside a biosteam ``BatchBioreactor``,
handling the kinetic simulate loop, reaction-time selection, species-to-
chemical mapping, and effluent construction so the model can participate in
design, costing, and TEA. :class:`~nskinetics.units.FermentationSaccharomycesEthanolIsobutanol` is a
fed-batch fermentation subclass of it that adds aeration, sucrose hydrolysis,
fed-batch spike-count retry, and yield/mass-balance validators on top of the
shared machinery.

Process section factories
--------------------------

Beyond single units, ``nskinetics.processes`` ships whole flowsheet sections
as biosteam system factories.
:obj:`~nskinetics.processes.create_sugar_prep_and_fermentation_system`
builds the isobutanol biorefinery's sugar-solution preparation and
fermentation section — a splitter, two evaporator/dilution/heat-exchange
conditioning trains (initial feed and spike feed), the
:class:`~nskinetics.units.FermentationSaccharomycesEthanolIsobutanol`
fermentor, and its aeration loop — with defaults that reproduce the inline
configuration of that biorefinery's ``system.py`` exactly. The factory
constructs a :class:`~nskinetics.units.FedBatchStrategySpecification` from
its own units and attaches it to the fermentor
(``V406.fed_batch_strategy_specification``, short alias ``V406.fbs_spec``);
construction is side-effect-free — the strategy is imposed only when the
caller invokes ``load_specifications()``. The factory also ships its own
chemical set,
:func:`~nskinetics.processes.create_sugar_prep_and_fermentation_chemicals`,
activated only by calling the factory with ``set_thermo=True``; by default
the caller's pre-set thermo is used unchanged. A few couplings are
deliberately left caller-side — fermentor specifications that reach into
upstream flowsheet streams, and the docking of the vent and effluent outlets
downstream; the factory's docstring lists them. See the
:doc:`quickstart tutorial <tutorial/01_quickstart>` for the factory in
action and :doc:`the API page <API/processes>` for every knob.
