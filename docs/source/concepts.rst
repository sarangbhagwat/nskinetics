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

Sensitivities across events
---------------------------

A local sensitivity is the derivative of an output with respect to a
parameter or an initial value, for one condition: ``dy/dp``, or the
elasticity ``d ln y / d ln p``, the percentage change of ``y`` per
percentage change of ``p``. Events make it delicate. An event makes a state
jump, and a parameter change also moves *when* the event fires; a correct
derivative must carry both effects. In the smallest example, ``S`` decays
at rate ``k*S`` and is topped up by 5 whenever it falls to 5, and ``n``
counts the top-ups:

.. code-block:: python

   import numpy as np
   import nskinetics as nsk
   from nskinetics.analysis import local_sensitivity

   toy = nsk.KineticModel.from_antimony("""
   model ev()
     S = 10; n = 0; k = 0.5;
     J: S => ; k*S;
     E1: at S < 5: S = S + 5, n = n + 1;
   end
   """)
   cond = nsk.Condition(t_end=6, n_steps=61, selections=['S', 'n'])
   ls = local_sensitivity(toy, cond, ['k'], outputs=['S'])
   print(ls.counters)
   print(f'dS/dk = {ls.raw[-1, 0, 0]:.4f}, '
         f'd ln S/d ln k = {ls.scaled[-1, 0, 0]:.5f}')

prints:

.. code-block:: text

   ('n',)
   dS/dk = -47.7956, d ln S/d ln k = -3.00000

Each top-up restores ``S`` to 10, at :math:`t_n = n \ln 2 / k`, so
:math:`S(t) = 10 \cdot 2^{n} e^{-kt}` after ``n`` top-ups. Four of them
happen before :math:`t = 6` at :math:`k = 0.5`, so :math:`S(6) = 160
e^{-3}`, :math:`dS/dk = -6\,S(6) = -960 e^{-3} = -47.7956` and
:math:`d \ln S / d \ln k = -kt = -3`, as computed.

Central differences from a restored baseline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`~nskinetics.analysis.local_sensitivity` sets each parameter to
:math:`p\,e^{\pm h}` (``rel_step`` :math:`h`, default ``1e-3``) and divides
the change of the output by :math:`2h`. That is a central difference in
:math:`\ln p`: accurate to :math:`O(h^2)`, with the same relative step for
parameters of any size. Each perturbation is a complete run of the
condition:

* Every run starts from the restored baseline (see "Resets, baselines and
  reproducible runs" above), so a perturbed run differs from the nominal
  one only in the perturbed value. The difference therefore contains
  everything the parameter does, jumps and shifted firing times included,
  without any event-specific bookkeeping. Without the baseline restore, a
  perturbed run would start from the previous run's spike counter and
  volume, and the difference would measure that instead.
* Every run, the nominal one included, sets each perturbed id explicitly.
  The values are merged into the condition with ``with_parameters`` and
  ``with_initial``, so the condition's own values stay in force. Each
  column is therefore a partial derivative with the *other listed ids*
  held fixed, and it can change when another id is added. Changing a
  compartment keeps the species amounts unless a concentration ``'[S]'``
  is also listed, in which case the concentration is held: with a rate
  ``k*S*env``, ``d[P]/d env`` is -9.50 with ``params=['env']`` and 0 with
  ``['env', '[S]']``. Morris and Sobol' screening set their sampled ids the
  same way.
* The ids are set after the reset, so SBML initial assignments are not
  re-evaluated. If ``S = S_init`` is an initial assignment, perturbing
  ``S_init`` leaves ``S`` at its nominal start and its effect through
  ``S`` is missing (a sensitivity of 0, an "inert" parameter in
  screening); NSKinetics warns about such ids. Perturb ``'init(S_init)'``
  instead, which re-evaluates the assignment.
* The runs use a private :meth:`~KineticModel.copy`. Your model is never
  changed, and only the copy's integrator is tightened: to at most
  ``(atol, rtol) = (1e-12, 1e-10)`` by default, or to
  ``tolerances=(atol, rtol)``. Integration error enters the difference
  divided by the step. At CVODE's own ``rtol = 1e-6`` and
  :math:`h = 10^{-3}` the derivative of even a plain decay
  ``S => P; (k1 + k2)*S`` is off by about 1%, and refined steps (below)
  amplify the error further; at the default ``rtol = 1e-10`` it stays
  below 1e-6. A ``UserWarning`` is raised when ``rtol`` exceeds ``rel_step /
  100``, or 1/100 of the relative size of a smaller step that a reported
  value used (a refined step, or one set by ``abs_step``). It flags only
  gross mismatches: at CVODE's ``rtol = 1e-6``, the 1% error above
  passes without a warning. Tighter tolerances need more integrator steps,
  so the copy's step budget (``maximum_num_steps``) is raised tenfold with
  them; if a run still fails where your model's own run succeeds, the error
  says so, and ``tolerances=`` loosens them.
* A parameter at 0 (or below) cannot be perturbed in :math:`\ln p`, so it is
  differenced in linear space, with ``abs_step``. Without one, a parameter
  at 0 gets ``rel_step`` in its own units, which has nothing to do with the
  model's scale: 1e-3 mol/L is ten times a 100 µM substrate, and on a
  product-inhibition toy at that scale it gave a derivative of the wrong
  sign, still flagged reliable. A warning names such parameters; give each
  an ``abs_step`` on the scale of the quantity, e.g. 1e-3 of the
  concentrations it competes with. Central differences then
  also run a parameter at 0 one step below zero. If that run fails (a
  negative value the model cannot integrate) while the other side runs,
  the other side's one-sided difference is used: ``scheme`` says
  ``'forward'`` or ``'backward'``, the value stays reliable, and a warning
  names the parameter. If the running side counted other firings than the
  nominal run, the step is refined first, as for any count change.
* Log scaling divides by the output, so it is only as good as the
  output's difference. Where an output sits at the integrator's noise
  level -- a depleted substrate, say -- its difference is integration
  error of about :math:`3\,\mathrm{atol}_y/h` and
  :math:`\mathrm{d}\ln y/\mathrm{d}\ln p` is noise divided by a tiny
  :math:`y`. Here :math:`\mathrm{atol}_y` is the absolute tolerance in the
  output's units: RoadRunner multiplies a scalar ``absolute_tolerance`` by
  each state's start value, so :math:`\mathrm{atol}_y =
  \mathrm{atol}\,|y_0|` for an output that starts at :math:`y_0 \ne 0`,
  and :math:`\mathrm{atol}` for one that starts at 0. The noise level
  is therefore relative to the start value: a 1 µM substrate in mol/L is
  as accurate as a 10 g/L one. A state that starts at 0 keeps the fixed
  tolerance :math:`\mathrm{atol}` for the whole run, so by default the
  copy's ``atol`` is also scaled down to the smallest peak below 1 of the
  outputs, and of the columns the metrics read, that start at 0 (one
  extra run; explicit ``tolerances=`` are used as given): a product at
  the nanomolar scale is then as accurate as one at 10 g/L, and the
  tolerances used are in ``ls.tolerances``.
  Log-scaled output entries below :math:`|y| = 30\,\mathrm{atol}_y/h`
  (:math:`|y|/|y_0|` below 3e-8 at the defaults, or :math:`|y|` below
  :math:`3 \cdot 10^4\,\mathrm{atol}` for an output starting at 0) are
  NaN; up to about :math:`1000\,\mathrm{atol}_y/h` (1e-6) they are kept
  but noisy, and the elasticity of a near-zero metric such as
  ``Final('[S]')`` of a residual substrate is not masked at all.
  ``reliable`` speaks about event counts and failed runs, not about such
  entries: use ``scaled='semi'`` or ``'none'`` there, or leave those
  times out.

With central differences the cost is two runs per parameter plus the
nominal one, one more run to read start values the condition does not set,
one to scale ``atol`` (default tolerances only), and two more per step
refinement (forward differences: one per parameter and per refinement), at
a few milliseconds per run for a small model.

RoadRunner's native sensitivities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

libroadrunner can integrate forward sensitivities alongside the model
(CVODES, ``RoadRunner.timeSeriesSensitivities``). NSKinetics does not use
them. With libroadrunner 2.10.0, the version this page was checked against,
they fail in ways that return plausible numbers and no error:

* **Event jumps are dropped.** Add a product to the toy above
  (``J: S => P``), so that the model has two state variables. The native
  ``dS/dk`` at :math:`t = 6` is then -13.91 instead of -47.80. It agrees
  with the finite difference up to the first top-up, at :math:`t = 1.39`,
  and departs from it from the next output time on.
* **Single-state models give zeros.** On the toy exactly as written above,
  and on a plain decay ``S => ; k*S`` without any event, every native
  sensitivity is 0.
* **Parameter subsets are mislabelled.** With ``params=[...]``, the values
  come back in the model's parameter order, under the names in the order
  asked for. Asking for ``['k', 'a']`` in a model that declares ``a``
  first returns the sensitivity to ``a`` labelled ``k``. A single parameter
  that is not the model's first gets 0. On the shipped *S. cerevisiae*
  model, the sensitivities of ``s_glu`` and ``s_IBO`` to ``k_13`` and
  ``k_16`` are all 0, although ``s_glu`` at 100 h does depend on ``k_13``.
  Ethanol's sensitivity to ``k_7ie``, one of the rate constants it depends
  on most, is also 0. Asking for all parameters (the default) raises a
  ``RuntimeError``, because some of the model's global parameters are
  defined by assignment rules. Listing the 87 other parameters in model
  order gives values that are off by factors of 2 to 13 for the three
  checked against finite differences, and NaN for two parameters.

A wrapper around them would silently report "insensitive", or a wrong
number. Finite differences cost more runs, but they are right across events
by construction, and where no derivative exists they say so.

What ``event_count_changed`` and ``reliable`` mean
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A derivative exists only where a small change of the parameter changes the
output smoothly. A change that adds or removes an event firing before the
end of the run makes the output jump. So does, at one output time, a change
that moves a firing from just after that time to just before it. The
finite difference then measures the jump divided by the step, which grows
without bound as the step shrinks. NSKinetics watches *event counters* to
catch this:

* Counters are found automatically: every variable that only events
  change and that an event sets to ``v + 1``, such as a :class:`FeedSpike`
  ``count_var`` or ``n`` above. A species that also reacts, or a rule
  variable, is not a counter even if an event adds 1 to it (a pulse feed
  ``S = S + 1``): its final value differs between runs anyway. Nor is a
  species in a compartment whose size changes: ``n = n + 1`` adds 1 to its
  concentration, so its amount is the count times the volume at each
  firing (NSKinetics warns and leaves it out; count with a parameter, as a
  :class:`FeedSpike` does). A species counter in a constant compartment is
  watched through its concentration ``'[n]'``, which rises by exactly 1
  per firing even when the compartment's size is perturbed.
  ``ls.counters`` lists them. ``event_counters=[...]`` replaces the list,
  and ``event_counters=()`` turns the check off. An event that switches a
  stage without counting (``is_aerobic = 0``) is not detected; list the
  flag it sets in ``event_counters``, since every watched variable is
  compared at every output time.
* After the perturbed runs, the firings each counter has counted by every
  output time (its value minus its start value, so perturbing a counter's
  own start value changes nothing) are compared with the nominal run's.
  If they differ anywhere -- a firing added or removed, or moved across an
  output time -- the step is divided by 10 and the parameter re-run, up to
  ``max_refinements`` times (default 2, so down to a step of ``1e-5``). If
  the counters then agree, the derivative is reliable. ``n_refinements``
  says how many refinements it took, and ``step`` gives the step of the
  reported difference.
* If they still differ, the parameter is flagged
  ``event_count_changed=True`` and ``reliable=False``. Its reported value
  is the difference on the side whose counters agree with the nominal run
  (``scheme`` is ``'forward'`` or ``'backward'``), when exactly one side
  does, taken at the largest step at which that side agreed: the smaller
  the step, the more the integration error weighs. Across the end of the
  run that is a one-sided derivative, not a derivative; across an output
  time, it is the derivative from the side that does not cross.
* ``reliable=False`` also marks a perturbed run that failed, when no
  one-sided difference could stand in for it (both sides of a central
  difference failed, or the forward run). Its column is NaN, and
  ``error`` says why.

At :math:`k = 4 \ln 2 / 6`, the fourth top-up of the toy lands exactly at
:math:`t = 6`:

.. code-block:: python

   edge = local_sensitivity(toy, cond.with_parameters(k=4 * np.log(2) / 6),
                            ['k'], outputs=['S'])
   print(edge.flags[['step', 'n_refinements', 'event_count_changed',
                     'reliable', 'scheme']])
   print(f'dS/dk = {edge.raw[-1, 0, 0]:.4f}')

prints:

.. code-block:: text

              step  n_refinements  event_count_changed  reliable    scheme
   parameter
   k         0.001              2                 True     False  backward
   dS/dk = -30.0266

A slightly larger ``k`` fires the fourth top-up before :math:`t = 6`, and a
slightly smaller one does not, so ``S(6)`` jumps from 5 to 10 at this
``k``. The reported value is the derivative from below,
:math:`-6 \cdot 5 = -30`, up to the truncation error of a one-sided
difference at the step 0.001 (0.1%). From above it is
:math:`-6 \cdot 10 = -60`, and no two-sided derivative exists.

The same check covers the rows of a time course. A perturbation that moves
a firing across an output time leaves the final count alone, but the
counter differs at that output time, so the step is refined: each
refinement narrows the window around the firing time tenfold, from 0.1% of
it for the toy (a rate constant at the default step) to 0.001%. A crossing
is flagged only when the output time lies inside that last window, and
then the reported value comes from the side that does not cross. Events
that switch a stage without counting still need ``event_counters=`` for
this.
:meth:`~nskinetics.analysis.LocalSensitivityResult.rank` shows the flags
next to every value,
:meth:`~nskinetics.analysis.LocalSensitivityResult.plot_time_course` draws
unreliable parameters dashed, and
:func:`~nskinetics.analysis.identifiability` warns when it uses one.

The basis of each metric
~~~~~~~~~~~~~~~~~~~~~~~~

A metric reduces a run to one number, for :func:`run_many`,
:func:`~nskinetics.analysis.local_sensitivity` and the screening
functions. Each one in ``nskinetics.analysis`` states exactly what it
computes, so two analyses mean the same thing by "yield":

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Metric
     - Value
   * - ``Final(sel)``
     - The value at the last output time.
   * - ``At(sel, t)``
     - The value at ``t``: linear interpolation between output times,
       exact at an output time. ``ValueError`` outside the run.
   * - ``Max(sel)``, ``Min(sel)``
     - The largest and smallest value over the output times.
   * - ``TimeAtMax(sel)``
     - The time of the maximum: the vertex of the parabola through the
       first grid maximum and its two neighbours (the grid time at the
       first or last row, or on a plateau).
   * - ``Integral(sel)``
     - The trapezoid rule over the output times, in the selection's units
       times model time units.
   * - ``EventCount(counter)``
     - The final minus the initial value of the counter.
   * - ``Productivity(P, turnaround=0)``
     - ``(P_end - P_0) / (t_end - t_0 + turnaround)``, in concentration
       per model time unit. It uses concentrations, so dilution by fed
       volume lowers it.
   * - ``Yield(P, S, volume=None, fed=None, basis='consumed')``
     - Amounts are concentrations times the ``volume`` selection (1
       without one). The fed amount is ``F = (V_fed,end - V_fed,0) *
       feed_conc``, from ``fed=(V_fed, feed_conc)``. ``'added'`` divides
       the product formed by ``A_S,0 + F``; ``'consumed'`` by
       ``A_S,0 + F - A_S,end``. No stoichiometric factor is applied.
   * - ``CallableMetric(fn, name, selections)``
     - Anything else: ``fn(result)`` must return a number, and
       ``selections`` names the columns it reads.

A zero denominator gives NaN, never an error or an infinity, although
importing NSKinetics makes numpy raise on division by zero. The metrics
over the output times (``Max``, ``Min``, ``TimeAtMax``, ``Integral``) are
only as fine as the output grid. Local sensitivities difference a metric as
computed, so one that is piecewise constant in the parameters has no
usable derivative: a plain grid time would give 0, or a jump of one grid
step over the step size, which is why ``TimeAtMax`` refines the grid
maximum by a parabola: the vertex moves continuously with the parameters,
so its derivative converges as the grid is refined, though only to first
order in the output spacing (use a dense grid around the peak).
``EventCount`` stays piecewise constant: its derivative is 0 unless a
count changes, which the flags report. Screen such metrics with Morris or
Sobol' instead.

.. _identifiability-noise:

Noise models and the rank test of identifiability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`~nskinetics.analysis.identifiability` divides each measured row by
its standard deviation: an absolute one in the output's units, or a
relative one, ``sigma=('relative', rel)``; a mapping can give each output
either, as in ``{'[S]': 0.5, '[P]': ('relative', 0.05)}``. Five rules keep
its verdicts sound.

**Keep the floor of the relative model at its default.** A pure relative
error, :math:`\sigma_i = \mathrm{rel}\,|y_i|`, gives unbounded weight to
the rows where an output is near 0, such as a depleted substrate or a
product at :math:`t_0`, while the integration and finite-difference errors
of :math:`y` there do not shrink with :math:`y`. The relative model is
therefore bounded below:

.. math::

   \sigma_i = \mathrm{rel} \cdot \max\left(|y_i|,\;
   f \max_t |y(t)|\right),

with :math:`f` = ``relative_floor`` (default ``1e-3``) and the maximum over
the whole output grid of the sensitivity run, so ``times=`` does not
change it. The floor acts like a detection limit at a small fraction of
the output's own peak, and does not depend on the output's unit. On
``S => P; (k1 + k2)*S`` with the substrate depleted (61 output times to
:math:`t = 30`, ``sigma=('relative', 0.1)``), the default floor leaves
the null singular value of ``k1`` and ``k2`` at 4e-7 of the largest,
below the default ``rank_rtol`` of 1e-6; a floor of 1e-4 already lifts it
to full rank. Smaller floors let the noise back in, so for them the rank
test scales its error estimate (below) by the extra weight they give the
rows, at most :math:`10^{-3}/f`. An output that is 0 throughout is left
out with a warning; give it an absolute sigma in a mapping to count it.

**Leave parameters at 0 out.** A parameter at 0 (or below) has a column
per model unit of :math:`p` instead of per :math:`\ln p`, so the singular
values, the numerical rank and with them the other parameters' standard
errors and ``identifiable`` verdicts depend on that unit. A warning says
so; judge the other parameters with ``params=`` leaving such parameters
out.

**Keep the differences precise.** The null singular value of an exact
non-identifiability sits at the finite-difference error level, not at 0:
integration noise of about :math:`3\,\mathrm{rtol}/h` plus the truncation
error that a sum such as ``(k1 + k2) * S`` does not cancel (a product
``k1 * k2`` does), about :math:`0.3\,h^2` for central differences and
:math:`h/2` for one-sided ones (``method='forward'``, or a column that
fell back to one side). At the defaults of
:func:`~nskinetics.analysis.local_sensitivity` (central differences,
``rel_step=1e-3``, tolerances of at most ``(1e-12, 1e-10)``) that level
stayed below the default ``rank_rtol`` on smooth toys. One-sided
differences or CVODE's own ``rtol=1e-6`` lift it to 1e-6 to 1e-3 and
beyond, where no ``rank_rtol`` separates an exact non-identifiability from
a genuine weak direction. The estimate counts ``rtol`` only: an ``atol``
above 1e-8 is warned about, and the level is then recorded as unknown
(``fd_error_level`` and ``rank_test_doubtful`` are ``None``).

**Read the warning as a screen.**
:func:`~nskinetics.analysis.identifiability` warns, for two or more
parameters, when its estimate of this error level exceeds ``rank_rtol``
and the smallest singular value the rank test retained lies within 100
times it, relative to the largest, and records the doubt as
``rank_test_doubtful=True`` (kept by ``to_csv``). A value already counted
as zero is the conservative verdict and raises no doubt. ``False`` is no
proof: the estimate assumes that the higher derivatives in :math:`\ln p`
are of the order of the first, while at a depletion front or a switch a
step in :math:`\ln p` shifts the front in time. Monod growth,
``S => X; (m1 + m2)*S/(Ks + S)*X`` with 20 g/L of substrate, 0.01 g/L of
biomass and ``Ks = 0.5``, with ``[X]`` and ``[S]`` measured on 301 output
times to :math:`t = 30` h and ``sigma=('relative', 0.1)``, depends on
``m1 + m2`` only. Yet at the defaults its null singular value sits at
about 630 times the estimate: the result reads full rank, with all three
parameters identifiable and no warning. Michaelis-Menten kinetics with
:math:`K_m` far below the substrate behaves alike. The relative standard
error check of ``identifiable`` is no substitute: it weakens as the data
get more precise or denser.

**Check step convergence.** The decisive test is to recompute
:func:`~nskinetics.analysis.local_sensitivity` at ``rel_step / 3`` with the
same tolerances and compare :math:`s_{min}/s_{max}`, which is
``1 / condition_number`` of a full-rank result (in general, the smallest
retained singular value over the largest,
``singular_values[numerical_rank - 1] / singular_values[0]``). On the
Monod example:

.. code-block:: python

   import nskinetics as nsk
   from nskinetics.analysis import identifiability, local_sensitivity

   monod = nsk.KineticModel.from_antimony("""
   model monod()
     X = 0.01; S = 20; m1 = 0.1; m2 = 0.3; Ks = 0.5;
     J1: S => X; (m1 + m2)*S/(Ks + S)*X;
   end
   """, units={'time': 'h', 'conc': 'g/L'})
   cond = nsk.Condition(t_end=30, n_steps=301)
   params, outputs = ['m1', 'm2', 'Ks'], ['[X]', '[S]']
   ls = local_sensitivity(monod, cond, params, outputs=outputs)
   fine = local_sensitivity(monod, cond, params, outputs=outputs,
                            rel_step=1e-3 / 3)
   for subset in (params, ['m2', 'Ks']):
       ident = identifiability(ls, ('relative', 0.1), params=subset)
       ident_fine = identifiability(fine, ('relative', 0.1), params=subset)
       falls = ident_fine.condition_number / ident.condition_number
       print(f'{subset}: rank {ident.numerical_rank}, falls {falls:.3g}x')

prints:

.. code-block:: text

   ['m1', 'm2', 'Ks']: rank 3, falls 9.06x
   ['m2', 'Ks']: rank 2, falls 0.975x

The ratio of the condition numbers is the factor by which
:math:`s_{min}/s_{max}` falls. About 9 means a null direction lifted by
the truncation error of central differences, which falls as :math:`h^2`
(``m1`` and ``m2``, despite the full rank); about 3 means the same for
one-sided differences, whose error falls as :math:`h`; about 1/3 means one
lifted by integration noise, which grows as :math:`1/h`; about 1 means a
genuine weak direction (``m2`` and ``Ks``). Run it before trusting a full
rank that matters.

Results depend on the condition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A sensitivity holds at the nominal parameter values, and for one
condition: its initial values, its feed settings, its run length and its
output grid. In :doc:`tutorial/07_sensitivity_identifiability`, the
elasticity of the titer with respect to the enzyme level is 3e-7 at 48 h,
when the substrate is used up, and 0.64 at 12 h. Change the condition and
the ranking can change with it. The same holds for identifiability, which
also depends on what is measured, when, and how precisely. Morris and Sobol'
screening vary only the sampled ranges; every other value comes from the
base condition. Sensitivity and screening results keep their condition
(``.condition``, also saved by ``to_csv``), so a number's context travels
with it.

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
