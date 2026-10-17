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
:meth:`KineticModel.from_antimony(text_or_path, units=...)
<KineticModel.from_antimony>` builds one from Antimony text (or an Antimony
file) in one step. To wrap a RoadRunner object you already have, for
example one compiled with ``tellurium.loadAntimonyModel(text)``, use
:class:`nsk.KineticModel(r, units=...) <KineticModel>`.
Every route produces the same wrapper, so the rest of the workflow (units,
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
names to record amounts. For reproducible work, describe each run as a
:class:`Condition` and call :meth:`~KineticModel.run`, which returns the
trajectory as a :class:`SimulationResult` (see
`Resets, baselines and reproducible runs`_).

Units
-----

Every :class:`KineticModel` carries a ``units`` dict with
``'time'`` and ``'conc'`` keys (default ``{'time': 'min', 'conc': 'M'}``).
Recognized time tokens are ``h``, ``hr``, ``min``, ``m``, ``s``, and ``sec``
(case-insensitive); recognized concentration tokens are ``M`` and ``mol/L``
for molar concentrations, and ``g/L``, ``kg/m3``, and ``kg/m^3`` for mass
concentrations — note that ``kg/m3`` is a *mass* unit, not molar, since it is
numerically equivalent to g/L. The :attr:`~KineticModel.time_factor`
property gives model time units per hour (``60.0`` for a model in minutes),
and
:attr:`~KineticModel.material_indexer` resolves the configured
concentration unit to the matching biosteam stream indexer (``'imol'`` for
molar, ``'imass'`` for mass) for use by the process bridge.
:meth:`~KineticModel.validate_units` raises a
``KineticSimulationError`` if either token is unrecognized.

The event lifecycle (important)
--------------------------------

The typical event workflow is: ``add_event(...)`` (repeatable) →
``compile_events()`` → run. **The ordering matters.** Compiling events
regenerates the underlying RoadRunner model and, as a side effect, resets all
model state to its SBML-defined origin values. Any custom initial conditions
or parameter overrides must therefore be applied *after* ``compile_events()``,
never before — values set beforehand are silently discarded when the model is
regenerated. This is documented directly on
:meth:`KineticModel.compile_events`, :meth:`Event.compile`, and
:meth:`FeedSpike.compile`; when in doubt, compile events first and set state
last.

The best place for those overrides is a :class:`Condition`:
:meth:`km.run(cond) <KineticModel.run>` resets the model, sets
``cond.parameters``, then ``cond.initial``, simulates, and finally restores
the overridden parameters, so the order is right by construction and nothing
leaks into the next run. If you drive the model by hand instead, the order
is ``compile_events()`` → ``reset()`` → set values → ``simulate()``. The
default reset restores floating species, compartments and event-assigned
variables to their origin values, so state you set *before* it is
overwritten (see the next section).

Resets, baselines and reproducible runs
----------------------------------------

A study simulates the same model many times: once per condition, scan
point, perturbation or fit iteration. Every run must start from the same,
well-defined state, or results quietly depend on what ran before.
NSKinetics calls that state the model's *baseline*. Here it is on a small
fed-batch model whose :class:`FeedSpike` tops ``S`` back up to 50 g/L up to
three times:

.. code-block:: python

   import nskinetics as nsk

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
   print(km.baseline_variables)

   cond = nsk.Condition(t_end=48, n_steps=481, selections=['[P]', 'env', 'n_spk'])
   for mode in ('baseline', 'baseline', 'roadrunner'):
       res = km.run(cond.replace(reset_kwargs={'mode': mode}))
       print(f"{mode:<10} start: env={res['env'][0]:g}, n_spk={res['n_spk'][0]:g}"
             f"  end: [P]={res.final['[P]']:.1f} g/L")

prints:

.. code-block:: text

   ('env', 'last_vol', 'tot_vol', 'n_spk')
   baseline   start: env=1, n_spk=0  end: [P]=161.9 g/L
   baseline   start: env=1, n_spk=0  end: [P]=161.9 g/L
   roadrunner start: env=1.331, n_spk=3  end: [P]=37.6 g/L

The two baseline runs are identical, down to the last bit. The plain
RoadRunner reset in the third run restores the species but not the working
volume or the spike counter. That run therefore starts in the grown vessel
with its three spikes already used up, and ends at 37.6 g/L instead of
161.9 g/L. It raises no error.

Reset modes
~~~~~~~~~~~

:meth:`KineticModel.reset` takes a ``mode``. The default, ``mode=None``,
depends on whether the model has a custom :attr:`~KineticModel.reset_func`:

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - ``mode``
     - ``reset_func``
     - What happens
   * - ``None`` (default)
     - not set
     - ``'baseline'``
   * - ``None`` (default)
     - set
     - ``'custom'``, as before this change
   * - ``'baseline'``
     - either
     - Baseline restore; any ``reset_func`` is bypassed
   * - ``'roadrunner'``
     - either
     - Plain ``RoadRunner.reset()``, the default before this change
   * - ``'custom'``
     - set
     - ``reset_func(km, **kwargs)``, which is fully responsible for
       restoring state (see `Writing a reset function`_)
   * - ``'custom'``
     - not set
     - ``ValueError``, raised before anything changes

``simulate(reset=True)``, :meth:`KineticModel.run` and
:class:`~nskinetics.units.NSKBatchReactor` all reset with the default mode.
Other keyword arguments of ``reset()``, and a condition's ``reset_kwargs``,
are forwarded to a custom ``reset_func``; the other modes ignore them:

.. code-block:: python

   km.reset()                     # baseline, or the reset_func if one is set
   km.reset(mode='roadrunner')    # the plain RoadRunner reset (the old default)
   km.reset(preserve=('n_max',))  # keep a modified spike cap across the reset

What the baseline restores, and what it keeps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RoadRunner's own ``reset()`` restores time, floating species and rate-rule
variables. It leaves compartments and parameters wherever the last run left
them, which breaks as soon as an event changes them, as the third run above
shows. A baseline reset also sets every entry of
:attr:`KineticModel.baseline_variables` back to its ``init(x)`` value. These
are the compartments, rate-rule variables and event-assigned non-species
variables, such as working volumes, spike counters and stage flags. They
include the targets of :class:`Event` / :class:`FeedSpike` objects and of
events written in the SBML or Antimony source.

The baseline variables are set before RoadRunner's reset as well as after
it. RoadRunner converts initial concentrations to amounts with the
*current* compartment volumes. Under a grown volume, an initial-assigned
species would start from the wrong value, and the others would be off by
rounding errors that the integrator amplifies across events. Restoring the
volumes first makes a baseline reset bit-identical to loading the model
afresh.

Global parameters that no rule or event changes are *not* restored, and
this is deliberate. ``km.set_value('k', 0.4)`` changes the current value of
``k`` but leaves ``init(k)`` at the model's ``0.5``. A baseline reset
therefore keeps your edit, and the common workflow "set ``k``, then
``simulate(reset=True)``" works as it always has. Parameters are settings;
the baseline is about state. To change a parameter for one run only, put it
in a condition's ``parameters`` or use :meth:`KineticModel.overrides`.

The flip side is that a compartment or an event-assigned variable you set by
hand *is* state, so the next baseline reset puts it back. Pass such values
for a run in ``Condition.initial``, or keep them across resets with
``preserve``.

``preserve``
~~~~~~~~~~~~

``reset(preserve=(...))`` reads the named selections just before the reset
and writes them back right after it, in every mode. It is never forwarded
to a ``reset_func``. It is for values a caller changed on purpose and must
keep across resets. For example,
:class:`~nskinetics.units.NSKBatchReactor` keeps a frozen or reduced spike
cap with it while it re-simulates. A ``preserve`` id that is not a model
selection raises ``ValueError`` before anything is reset, and so does an
``init(...)`` selection: a reset never changes ``init(...)`` values, and
writing one back would re-initialize the model and revert your parameter
edits. In the conditions
of :func:`run_many` and :func:`scan`, ``preserve`` may name settings only:
preserving a state variable carries the previous run's end state into the
next run, so it is rejected there.

Writing a reset function
~~~~~~~~~~~~~~~~~~~~~~~~

A custom ``reset_func`` replaces the baseline reset, so it must restore
every state variable itself. Two rules keep its runs reproducible:

* Set compartment volumes and event-assigned variables *before* calling
  ``km._te.reset()`` as well as after it, or simply start with
  ``km.reset(mode='baseline')`` and then do the extra work. RoadRunner's
  reset converts initial concentrations to amounts at the *current*
  volumes, so a reset at the previous run's end volume starts the species
  from values that depend on that run.
* Write with item access (``r['env'] = 1.``) or
  :meth:`~KineticModel.set_value`, not attribute access (``r.env = 1.``).
  libroadrunner removes the attribute properties whenever another model is
  loaded, and from then on an attribute write silently creates a Python
  attribute that the model never sees.

For the model above, a reset function that also restores a spike cap reads:

.. code-block:: python

   def reset_model(km, reset_cap=True, **kwargs):
       km.reset(mode='baseline')    # volumes, counters and species
       if reset_cap:
           km._te['n_max'] = 3      # the extra work: the cap is a setting

   km.reset_func = reset_model

The shipped yeast model's reset function follows the first rule by setting
its bookkeeping variables before and after ``km._te.reset()``.

``Condition.parameters`` vs ``Condition.initial``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :class:`Condition` is the immutable description of one run.
:meth:`KineticModel.run` executes it in a fixed order:

1. Validate every id. Unknown, rule-defined or read-only ids, and the
   ``init(x)`` forms that would rewrite the model (see below), raise
   ``ValueError`` before anything changes.
2. Set any ``init(x)`` selections (see below).
3. Reset the model, if ``reset`` is true, with ``reset_kwargs``.
4. Set ``parameters``.
5. Set ``initial``.
6. Simulate, then restore.

The two mappings mean different things:

* ``parameters`` are settings for this run: rate constants, thresholds,
  caps. They are restored to their pre-run values when the run ends, even if
  it fails.
* ``initial`` is the starting state of the run: ``'[S]'`` (a
  concentration), ``'S'`` (an amount), a compartment volume, a counter. Like
  after :meth:`~KineticModel.simulate`, state is left where the run ends,
  and the next reset restores it. An ``initial`` entry that a baseline
  reset would not restore by itself (a global parameter, a boundary species
  or an ``init(x)`` selection) is restored after the run, like a parameter;
  a boundary species gets its pre-run amount back.
  The order of the keys does not matter: ``init(x)`` selections are set
  first, then compartment volumes, then everything else, so every value you
  give holds at ``t0``. Setting a volume keeps the species *amounts*, so
  ``{'env': 2}`` alone halves every concentration of a 1 L start, while
  ``{'env': 2, '[S]': 50}`` (in either order) starts 2 L at 50. A species
  may be given as an amount (``'S'``) or a concentration (``'[S]'``), not
  both. ``parameters`` follow the same fixed order, so their key order does
  not matter either.

``parameters`` are applied after the reset, so they do not re-evaluate SBML
initial assignments. In a model with ``S = 20*k``,
``parameters={'k': 1.0}`` still starts ``S`` at the value computed from the
model's own ``k``. Set such starting values through ``initial``. An id may
appear in only one of the two mappings.

To run the model as if it were *defined* with another parameter value, set
``initial={'init(k)': v}``. RoadRunner re-initializes the whole model
whenever an ``init(...)`` value is written, re-evaluating every initial
assignment, so :meth:`~KineticModel.run` writes these entries first, before
the reset. Your other settings (global parameters and boundary species not
defined by an initial assignment) are written back right after, and every
setting is restored when the run ends. Because that restore re-initializes
the model too, the model is left at its initial state rather than its
end-of-run state; read the final state from the result. The three ways of
changing ``k`` differ when initial assignments read it. For a model with
``k = 0.5; k2 = 2*k; S = 20*k`` and a new value ``v``:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * - Route
     - ``k`` in the run
     - ``S`` starts at
     - ``k2``
   * - ``Condition(parameters={'k': v})``
     - ``v``
     - ``20*0.5``
     - ``2*0.5``
   * - ``km.set_value('k', v)``, then a reset
     - ``v``
     - ``20*v``
     - ``2*0.5``
   * - ``Condition(initial={'init(k)': v})``
     - ``v``
     - ``20*v``
     - ``2*v``

Only the last route starts exactly like a freshly loaded model with
``k = v``: a reset re-evaluates the initial values of species and
compartments, but not the initial assignments of other parameters. It is
also the slowest on a large model: every ``init(...)`` write makes
RoadRunner re-initialize the whole model, which takes about 0.8 s on the
shipped yeast model, and a run with such an entry writes it twice (once to apply it, once
to restore it). Prefer ``parameters`` when no initial assignment reads the
parameter. The re-initialization also resets some integrators' settings
(e.g. rk45's ``epsilon``, gillespie's ``seed``), so NSKinetics re-applies
the integrator's settings after each such write.

Two ``init(...)`` forms are rejected with ``ValueError`` before anything
changes, because RoadRunner would rewrite the model's SBML for good and
restoring the old value does not undo it (later runs, :meth:`~KineticModel.copy`
and every batch row would see the change):

* ``init()`` of an id that is itself defined by an initial assignment, such
  as ``init(S)``, ``init([S])`` or ``init(k2)`` above: the write deletes the
  assignment. Set ``init()`` of its inputs instead (``init(k)``).
* ``init()`` of a species in the other form than its SBML declares. An
  Antimony ``S = 50`` declares an initial *concentration*, so use
  ``init([S])``; ``init(S)`` would redeclare it as an initial amount, and a
  later run that changes the volume would start ``S`` elsewhere. A species
  declared with an initial amount (e.g. ``substanceOnly``) takes ``init(S)``
  and rejects ``init([S])``.

The same checks apply to :meth:`~KineticModel.overrides`.

:func:`scan` makes the choice for you. A grid key the base condition
already sets stays in that mapping, and its value is replaced. Any other
key goes into ``initial`` if it is state: a bracketed selection, a species
or compartment id, or a rate-rule or event-assigned variable. The rest go
into ``parameters``. Override this with ``placement=`` for keys the base
condition does not set. :func:`run_many` and :func:`scan` require
``reset=True`` in every condition, and reject ``reset_kwargs`` with
``mode='roadrunner'``, which resumes from the previous run's events and
volumes, or with a ``preserve`` of state variables. Every row also starts
from your model's settings (global parameters and boundary species), which
are written back before each condition, so a ``reset_func`` that changes a
setting, such as the shipped model's ``reset_spike_cap``, cannot leak into
later rows. A result therefore never depends on which run came before it or
which worker ran it, as long as a custom ``reset_func`` restores every state
variable, setting the volumes before RoadRunner's reset (see `Writing a
reset function`_).

Why ``copy()`` rebuilds from the origin SBML
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Batch runs and analyses work on a private :meth:`KineticModel.copy`, so they
never change your model. A copy is the model *as defined*, not as last run.
It reloads the model's origin SBML, which is ``RoadRunner.getSBML()`` read
when the copy is made: the document as loaded plus every structural edit
made on the RoadRunner -- events compiled by
:meth:`~KineticModel.compile_events` or directly by ``Event.compile`` /
``FeedSpike.compile``, ``r.addEvent(...)`` -- and values written to
``init(x)`` selections, with origin values and initial assignments intact.
Run state and parameter edits such as ``r.k = 0.4`` never enter this
document. The copy then re-applies the settings that make up the model:

* user-edited values of global parameters and boundary species (a reset
  keeps these, so the copy must too);
* the integrator and the settings that differ from its defaults;
* units, events and the ``reset_func``;
* extra public attributes.

The copy finally applies a baseline reset (never the ``reset_func``), so it
starts at its baseline, whatever the original ran before.

Exporting with ``RoadRunner.getCurrentSBML()`` instead would copy the wrong
thing. It bakes the post-run state into the initial values: ``S`` would
start at the 1.9e-6 g/L it ended at. It also drops initial assignments, so
a later parameter change would no longer flow into the initial values that
depend on it.

Pickling or deep-copying the model reproduces the RoadRunner's *current*
state: an unpickled copy of the model above, taken after a run, sits in the
grown vessel (``env`` 1.331) with three spikes already counted.
:meth:`~KineticModel.run` and ``simulate(reset=True)`` reset first, so they
still give the right trajectory on such a copy, but ``simulate()`` without a
reset would continue from that state. A :meth:`~KineticModel.copy` is
already at its baseline. Its recipe is plain data: SBML text, numbers, the
event objects and the reset function, with integrator settings and user
edits carried explicitly. It is what :func:`run_many` sends to each worker
process when ``workers > 1``, and each worker builds its model from it
once. With ``workers > 1``, an unpicklable ``reset_func`` (a lambda or a
closure) raises a ``ValueError`` that names it before any process starts.
On macOS and Windows, workers start with ``spawn`` and re-import the main
module: call :func:`run_many` and :func:`scan` under
``if __name__ == '__main__':`` in scripts, and define metrics and the
``reset_func`` at module level.

libroadrunner keeps attribute access (``r.k = 0.4``) as properties on the
RoadRunner *class* and replaces them with the newest model's on every load.
The loads that :meth:`~KineticModel.copy`, :func:`run_many` and
:func:`scan` make on your behalf keep the other live models' properties in
place; explicit Tellurium loads still follow libroadrunner's behaviour.
Item access, ``r['k'] = 0.4``, is always safe.

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
builds the sugar-solution preparation and fermentation section of the
:doc:`isobutanol biorefinery case study <case_studies/isobutanol/index>` — a
splitter, two evaporator/dilution/heat-exchange conditioning trains (initial
feed and spike feed), the
:class:`~nskinetics.units.FermentationSaccharomycesEthanolIsobutanol`
fermentor, and its aeration loop — with defaults that match the
configuration the biorefinery itself uses (its flowsheet,
``nskinetics.biorefineries.isobutanol.system``, builds this section by
calling this factory). The factory
constructs a :class:`~nskinetics.units.FedBatchStrategySpecification` from
its own units and attaches it to the fermentor
(``V406.fed_batch_strategy_specification``, short alias ``V406.fbs_spec``);
construction is side-effect-free — the strategy is imposed only when the
caller invokes ``load_specifications()``. The same specification is also
attached as the fermentor's ``spike_feed_reconciler``: every fermentor run
ends by comparing the spike-feed sugar its kinetics added with what the
spike inlet delivers, and, if they disagree, re-derives the splitter split
from that run and re-simulates the feed trains (with the spike count
frozen) until they agree — so the effluent is conserved at the reactor
boundary on every run, not only when the last ``load_specifications()``
happened to match it. The factory also ships its own
chemical set,
:func:`~nskinetics.processes.create_sugar_prep_and_fermentation_chemicals`,
activated only by calling the factory with ``set_thermo=True``; by default
the caller's pre-set thermo is used unchanged. A few couplings are
deliberately left caller-side — fermentor specifications that reach into
upstream flowsheet streams, and the docking of the vent and effluent outlets
downstream; the factory's docstring lists them. See the
:doc:`quickstart tutorial <tutorial/01_quickstart>` for the factory in
action and :doc:`the API page <API/processes>` for every knob.
