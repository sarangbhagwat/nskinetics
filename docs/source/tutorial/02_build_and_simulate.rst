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

:meth:`~nskinetics.KineticModel.validate_units` raises a
``KineticSimulationError`` if either unit token is unrecognized — call it
right after construction to catch typos early. Recognized time tokens are
``h``, ``hr``, ``min``, ``m``, ``s``, and ``sec``; recognized concentration
tokens are ``M``/``mol/L`` (molar) and ``g/L``/``kg/m3``/``kg/m^3`` (mass).

Once units are set, two properties become meaningful:

.. code-block:: python

   print('time_factor:', km.time_factor, 'indexer:', km.material_indexer)

:attr:`~nskinetics.KineticModel.time_factor` gives hours per
model time unit (here ``1.0``, since you declared the model's time unit as
hours via ``units={'time': 'h'}``); :attr:`~nskinetics.KineticModel.material_indexer`
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

:meth:`~nskinetics.KineticModel.reset` restores model *state*
(species amounts, time) to its SBML-defined origin values, without
touching model *structure*. Call it before simulating to guarantee a clean
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

Putting the whole snippet together and running it in the ``HP_2024``
environment prints:

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

Next: :doc:`03_events` introduces the event API used to change model state
mid-simulation (feed spikes, parameter switches, and more).
