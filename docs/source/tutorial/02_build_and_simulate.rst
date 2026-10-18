Build and simulate a model
===========================

This page walks through the full minimal workflow: write a tiny kinetic
model in Antimony, wrap it in a :class:`~nskinetics.KineticModel`,
set and validate units, read values, and simulate. See
:doc:`../concepts` for the ideas behind each step.

The model
---------

NSKinetics models are declared as SBML. The easiest way to author one by
hand is Antimony, a compact text syntax that Tellurium compiles to SBML.
Here is a one-reaction model: species ``S`` decays into ``P`` inside a
compartment ``env``, with mass-action rate constant ``k``.

.. code-block:: python

   import tellurium as te
   import nskinetics as nsk

   model = """
   model demo()
     compartment env; species S in env, P in env;
     S = 10; P = 0; env = 1; k = 0.3;
     J: S => P; k*S*env;
   end
   """
   r = te.loadAntimonyModel(model)

``te.loadAntimonyModel`` compiles the Antimony text to SBML and returns a
Tellurium extended RoadRunner instance — the engine that actually
integrates the ODEs.

Wrap it and set units
----------------------

:class:`~nskinetics.KineticModel` wraps the RoadRunner object,
adding unit-aware value access and a Python event API on top. It carries a
``units`` dict with ``'time'`` and ``'conc'`` keys (default
``{'time': 'min', 'conc': 'M'}``); pass your own to match the model:

.. code-block:: python

   km = nsk.KineticModel(r, units={'time': 'h', 'conc': 'g/L'})
   km.validate_units()

The two steps can also be done in one:
``nsk.KineticModel.from_antimony(model, units={'time': 'h', 'conc': 'g/L'})``
(see :meth:`~nskinetics.KineticModel.from_antimony`) compiles the Antimony
text, or an Antimony file, and wraps the result. This page keeps ``r`` in
hand to show what the wrapper adds.

:meth:`~nskinetics.KineticModel.validate_units` raises a
``KineticSimulationError`` if either unit token is unrecognized — call it
right after construction to catch typos early. Recognized time tokens are
``h``, ``hr``, ``min``, ``m``, ``s``, and ``sec``; recognized concentration
tokens are ``M``/``mol/L`` (molar) and ``g/L``/``kg/m3``/``kg/m^3`` (mass).

Once units are set, two properties become meaningful:

.. code-block:: python

   print('time_factor:', km.time_factor, 'indexer:', km.material_indexer)

:attr:`~nskinetics.KineticModel.time_factor` gives model time
units per hour (here ``1.0``, since you declared the model's time unit as
hours via ``units={'time': 'h'}``; a model in minutes gives ``60.0``); :attr:`~nskinetics.KineticModel.material_indexer`
resolves the concentration unit to the matching biosteam stream indexer
(``'imass'`` for a mass unit like ``g/L``, ``'imol'`` for a molar one) —
used by the process bridge described in later tutorial pages.

Reading values: brackets vs. plain names
------------------------------------------

:meth:`~nskinetics.KineticModel.get_value` reads a RoadRunner
selection. A bracketed name (``'[S]'``) returns a *concentration*; a plain
name (``'S'``) returns an *amount*:

.. code-block:: python

   print('[S]0 =', km.get_value('[S]'))

Since the compartment volume ``env`` is ``1``, concentration and amount
coincide here — but they will not in general, so use whichever the
question calls for.

Reset, then simulate
----------------------

:meth:`~nskinetics.KineticModel.reset` restores the model's
*baseline* state without touching its *structure*. Time, floating species,
compartment volumes, rate-rule variables and any variable an event changes
go back to their SBML-defined origin values. Parameter values you set yourself (for example
``km.set_value('k', 0.4)``) are kept. :doc:`../concepts` explains exactly
what is restored and why. Call it before simulating to guarantee a clean
start, even on a freshly-constructed model:

.. code-block:: python

   km.reset()

:meth:`~nskinetics.KineticModel.simulate` integrates the model over
``[t0, t_end]`` and stores the trajectory on the kinetic model itself —
``km.results_df`` (DataFrame), ``km.results_array`` (2-D array),
``km.results_dict``, and ``km.results_col_names``. Request ``'time'`` plus
whichever bracketed concentration selections you want recorded (bracketed so the
values are concentrations, and so plotting works with no extra arguments); the
returned array is the same one stored on ``km``:

.. code-block:: python

   result = km.simulate(0, 10, 101, ['time', '[S]', '[P]'])
   print('t=10:', result[-1])

The underlying RoadRunner object is still reachable as ``km._te`` for
lower-level calls, but ``km.simulate`` is the recommended entry point.

Plotting the run
-----------------

:meth:`~nskinetics.KineticModel.plot_simulation_results` (aliases
:meth:`~nskinetics.KineticModel.plot_time_course`,
:meth:`~nskinetics.KineticModel.plot_trajectory`) plots the most
recent simulation. Pass ``labels`` for a readable legend:

.. code-block:: python

   fig, ax = km.plot_simulation_results(labels=['S', 'P'])

.. figure:: /_static/images/examples/tutorial_02_build_and_simulate.png
   :width: 400

   ``S`` → ``P`` under first-order decay (``k=0.3``).

Running it
----------

Putting the whole snippet together and running it prints:

.. code-block:: text

   time_factor: 1.0 indexer: imass
   [S]0 = 10.0
   t=10: [10.     0.498  9.502]

``time_factor`` and ``indexer`` confirm the units were recognized;
``[S]0 = 10.0`` matches the Antimony model's initial condition
(``S = 10``, ``env = 1``, so concentration equals amount); and the final
row of the simulated array shows ``S`` decayed from ``10`` to about
``0.498`` g/L by ``t=10`` h while ``P`` rose to about ``9.502`` g/L — the
stoichiometric (1:1) conversion of ``S`` into ``P`` under first-order decay,
as expected from ``k=0.3``.

Reproducible runs: ``Condition`` and ``run``
---------------------------------------------

``reset()``, then a few ``set_value`` calls, then ``simulate()`` works, but
it leaves you to remember what you changed and to change it back. A
:class:`~nskinetics.Condition` describes one run completely, and
:meth:`~nskinetics.KineticModel.run` executes it. A condition holds:

* the initial values;
* temporary parameter values;
* the output grid;
* an optional label.

Here the run starts ``S`` at 20 g/L and uses a faster ``k``:

.. code-block:: python

   cond = nsk.Condition(initial={'[S]': 20}, parameters={'k': 0.5},
                        t_end=10, n_steps=101, label='fast')
   res = km.run(cond)
   print(res)
   print(res.final)
   print('k after the run:', km.get_value('k'))

prints:

.. code-block:: text

   SimulationResult(label='fast', rows=101, t=[0.0, 10.0], columns=['time', '[S]', '[P]'], units={'time': 'h', 'conc': 'g/L'})
   time      10
   [S]    0.135
   [P]     19.9
   Name: 100, dtype: float64
   k after the run: 0.3

:meth:`~nskinetics.KineticModel.run` follows a fixed order:

1. Reset the model to its baseline.
2. Set the condition's ``parameters`` (``k = 0.5``).
3. Set its ``initial`` values (``[S] = 20``).
4. Simulate.
5. Put ``k`` back to the value it had before the run.

The final ``[S]`` matches the closed form
:math:`20\,e^{-0.5 \cdot 10} \approx 0.135` g/L. Running the same condition
again gives an identical result, whatever ran in between. An id the model
does not have raises a ``ValueError`` listing it, before anything is
changed. Conditions are immutable. Derive variants with
:meth:`~nskinetics.Condition.with_parameters` and
:meth:`~nskinetics.Condition.with_initial`, which merge new values into the
existing ones, or with :meth:`~nskinetics.Condition.replace`, which swaps
whole fields.

The returned :class:`~nskinetics.SimulationResult` holds the trajectory
(``res.df``), the model's units, the condition and its label.

* ``res['[P]']`` is a column as a numpy array, and ``res.time`` is the time
  column.
* ``res.final`` is the last row.
* :meth:`~nskinetics.SimulationResult.at` interpolates linearly between
  output points.
* :meth:`~nskinetics.SimulationResult.plot` draws the trajectory with the
  same styling as ``km.plot_simulation_results``.

.. code-block:: python

   print(res.at(2.5))
   fig, ax = res.plot(labels=['S', 'P'])

prints:

.. code-block:: text

   time    2.5
   [S]    5.73
   [P]    14.3
   Name: 2.5, dtype: float64

Output times
~~~~~~~~~~~~

``t_end`` and ``n_steps`` give an evenly spaced grid. To get rows at
specific times instead, for example the sampling times of an experiment,
pass ``times``. Here :meth:`~nskinetics.Condition.replace` swaps the grid,
and ``t_end`` is dropped:

.. code-block:: python

   coarse = km.run(cond.replace(times=[0, 1, 2, 5, 10]))
   print(coarse.df)

prints:

.. code-block:: text

      time   [S]  [P]
   0     0    20    0
   1     1  12.1 7.87
   2     2  7.36 12.6
   3     5  1.64 18.4
   4    10 0.135 19.9

The rows fall exactly on the requested times; the integrator still takes
whatever internal steps it needs, and events still fire between output
points. ``km.simulate`` accepts the same ``times=`` argument.

Saving a run
~~~~~~~~~~~~

:meth:`~nskinetics.SimulationResult.to_csv` writes the trajectory as a
plain CSV file. The file starts with a few ``#`` header lines that hold the
units, the label and the full condition.
:meth:`~nskinetics.SimulationResult.from_csv` reads it back exactly:

.. code-block:: python

   res.to_csv('fast.csv')
   back = nsk.SimulationResult.from_csv('fast.csv')
   print(back.df.equals(res.df), back.condition == cond, back.units)

prints:

.. code-block:: text

   True True {'time': 'h', 'conc': 'g/L'}

The first lines of ``fast.csv`` are:

.. code-block:: text

   # nskinetics: "SimulationResult"
   # units: {"time": "h", "conc": "g/L"}
   # label: "fast"
   # condition: {"initial": {"[S]": 20.0}, "parameters": {"k": 0.5}, "t0": 0.0, "t_end": 10.0, "n_steps": 101, "times": null, "selections": null, "reset": true, "reset_kwargs": {}, "label": "fast"}
   time,[S],[P]
   0.0,20.0,0.0

Many runs: ``run_many`` and ``scan``
-------------------------------------

:func:`~nskinetics.run_many` runs a list of conditions and evaluates
*metrics* on each result. A metric is any function that takes a
:class:`~nskinetics.SimulationResult` and returns a number. It can read only
the columns the run recorded (by default ``time`` and the species
concentrations), so to measure a counter or a volume such as ``'n_spk'`` or
``'env'``, list it in the condition's ``selections`` or give the metric
function a ``selections`` attribute (e.g. ``fn.selections = ('n_spk',)``);
otherwise the row fails with a message saying so. The batch runs
on a private copy of the model (:meth:`~nskinetics.KineticModel.copy`), so
``km`` itself is never changed:

.. code-block:: python

   P_end = {'P_end': lambda res: res.final['[P]']}
   conds = [cond.with_parameters(k=k).replace(label=f'k={k}')
            for k in (0.05, 0.1, 0.2)]
   batch = nsk.run_many(km, conds, metrics=P_end)
   print(batch.metrics)

prints:

.. code-block:: text

       label status error  P_end
   0  k=0.05     ok  None   7.87
   1   k=0.1     ok  None   12.6
   2   k=0.2     ok  None   17.3

Each value is :math:`20\,(1 - e^{-10k})`, the closed form for first-order
conversion. A run that fails, for example because the integrator gives up,
does not stop the batch. Its row gets ``status='failed'``, the error message
and ``NaN`` metrics, and ``on_error='raise'`` stops at the first failure
instead. ``batch.results`` keeps every
:class:`~nskinetics.SimulationResult`. ``workers=4`` spreads the runs over
up to four processes, each with its own copy of the model; the default
``workers=1`` runs them one after another. On macOS and Windows the
processes start with ``spawn`` and re-import your script, so call
``run_many`` under ``if __name__ == '__main__':`` and define metrics as
module-level functions: the ``lambda`` above works with ``workers > 1`` only
on Linux. The ready-made metrics of :mod:`nskinetics.analysis` (e.g.
``Final('[P]')``, ``Yield`` and ``Productivity``; see
:doc:`07_sensitivity_identifiability`) state their basis explicitly and work
with any number of workers.

:func:`~nskinetics.scan` builds the conditions for you from a grid of
values. It returns a :class:`~nskinetics.ScanResult`, which can arrange a
metric over the grid. Each grid key is set the way it should be: ``k`` is
a parameter, so it goes into the condition's ``parameters``; ``[S]`` is
state, so it goes into ``initial``. A key the base condition already sets
stays where it is, and the grid value replaces the base value; here
``cond`` sets both.

.. code-block:: python

   grid = nsk.scan(km, cond, {'k': [0.05, 0.1, 0.2], '[S]': [10, 20]},
                   metrics=P_end)
   print(grid.placement)
   print(grid.pivot('P_end'))

prints:

.. code-block:: text

   {'k': 'parameters', '[S]': 'initial'}
   [S]    10   20
   k
   0.05 3.93 7.87
   0.1  6.32 12.6
   0.2  8.65 17.3

``grid.metrics`` holds the full table, with one column per grid variable.
For a two-variable grid, :meth:`~nskinetics.ScanResult.contourf` draws the
metric as a filled contour plot.

Next: :doc:`03_events` introduces the event API used to change model state
mid-simulation (feed spikes, parameter switches, and more).
