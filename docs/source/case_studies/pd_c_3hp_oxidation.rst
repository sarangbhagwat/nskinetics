.. _pd-c-3hp-oxidation:

Pd/C oxidation of 3-HP to malonic acid
======================================

.. The images on this page are our own reproductions; no image of the paper
   is in the repository. They are written, from the repository root, by
       python -m nskinetics.models.pd_c_3hp_oxidation.figures \
           docs/source/_static/images/pd_c_3hp_oxidation --formats png --dpi 100
   *.png is gitignored; .gitignore makes an exception for
   docs/source/_static/images/pd_c_3hp_oxidation/*.png.

.. Each python code block below is followed by a text block with what it
   prints. nskinetics/tests/test_pd_c_3hp_oxidation_docs.py runs the blocks
   in order and compares the two, so keep that pairing when editing.

``nskinetics.models.pd_c_3hp_oxidation`` holds kinetic models of the
selective oxidation of 3-hydroxypropionic acid (3-HP) to malonic acid (MA)
over 10 % Pd/C, with O2 or H2O2 as the oxidant, after

    S. Y. Kim, J. E. Chavarrio, J. Install, G. W. Huber, "Selective
    oxidation of 3-hydroxypropionic acid to malonic acid over Pd/C:
    Mechanistic and kinetic study", *Appl. Catal. B: Environ. Energy* 386
    (2026) 126403, https://doi.org/10.1016/j.apcatb.2026.126403

and the authors' data deposit (Dryad,
https://doi.org/10.5061/dryad.jsxksn0q4, CC0), called *the deposit* below.
The package provides:

- three models, each an Antimony file with its generated SBML file and a
  loader that returns a :class:`~nskinetics.KineticModel` (time in min,
  concentrations in M, temperatures in K);
- the fit that defines one constant of the third model, as a function;
- the deposited tables of Figs. 9 and 10 (verbatim copies), the 50 °C data
  table of the authors' heat-map script, curves digitized from SI Fig. S12
  and Fig. 11, and readers for them;
- functions that draw Fig. 9, SI Fig. S12 and Fig. 11 from the models, one
  figure that compares the computed Fig. 11 with the published one, and one
  that verifies the Table 3 models against the deposited model nodes;
- ``NOTES.md``, the full record of what is printed in the paper, what was
  inferred from the deposit and where the two disagree. This page
  summarizes it.

The package needs no optional dependencies, and importing it loads no
model. Its API is documented on :doc:`the shipped-models page </API/models>`.

The paper describes one model per oxidant, fitted at 50 °C (its Table 3).
Its deposit holds a third model, in the script with which the
temperature–time maps of its Fig. 11 were computed. The package ships these
three. Three points apply, and `What does not match the paper, and what is
inferred`_ gives the details; read it before using a model for anything
other than redrawing the figures:

- The two Table 3 models reproduce Fig. 9 and SI Fig. S12b with the
  constants printed in the paper (for O2, together with one factor,
  ``x_O2_ref``, that is stated in neither the paper nor the deposit). They
  are models of the 50 °C experiments and have no temperature parameter.
- The third model, the model of the authors' deposited heat-map script
  (called here the *simplified reversible model*; the name is this
  package's), is the only one with a temperature. The constants its
  script's run fitted are not in the deposit. Here every constant the
  deposit gives is used as given (nine of ten; where the deposit gives one
  twice, the full-precision value: ``k_1``, ``k_2`` and ``k_5`` from its
  Arrhenius table, not the rounded values of the script's ``k0``), and the
  one it does not give, ``k_3``, is fitted to the 50 °C data with the
  authors' fitting set-up (`The constants of the simplified reversible
  model`_). No constant is fitted to a published figure. The model comes close to the
  published Fig. 11 without being fitted to it, and it does **not** describe
  the paper's experiments at 25, 75 and 100 °C (`Limits of the temperature
  dependence`_). The model curves of Fig. 10a–c are not reproduced by this
  package.
- Several statements on this page rest on analyses that are not in the
  repository. They are marked "(not shipped)"; `What is not in the
  repository`_ lists them. The repeat of the authors' solver run and the
  fit of ``k_3`` are in the repository.

The three models
----------------

Each model is an Antimony file with the SBML file generated from it next to
it (same name, ending in ``_sbml.xml``).

.. list-table::
   :header-rows: 1
   :widths: 15 37 18 30

   * - Model
     - Loader and Antimony file
     - Compared with
     - Where it comes from
   * - H2O2 model
     - | ``load_h2o2_model()``
       | ``pd_c_3hp_oxidation_h2o2_antimony.txt``
     - Fig. 9a (50 °C), which it reproduces. Its predicted H2O2 is compared
       with SI Fig. S12a, whose H2O2 curve lies up to 0.38 equivalents lower
     - **Paper**: Table 3, H2O2 column, as printed, in the rate form of SI
       Fig. S2. How H2O2 enters the rates was **inferred** from the
       deposited model curve (the authors posted only the O2 version of
       their script).
   * - O2 model
     - | ``load_o2_model()``
       | ``pd_c_3hp_oxidation_o2_antimony.txt``
     - Fig. 9b and SI Fig. S12b (50 °C), which it reproduces
     - **Paper**: Table 3, O2 column, as printed, in the rate form of SI
       Fig. S2 with one **inferred** factor (``x_O2_ref``) in place of the
       O2 activity. It is a 50 °C model without a temperature parameter.
   * - Simplified reversible O2 model
     - | ``load_o2_simplified_reversible_model(T)``
       | ``pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt``
     - Fig. 11 (contour lines 1.05 px rms from the published ones in panel a
       and 0.30 px in panel b); the 50 °C experiment; the experiments at
       25, 75 and 100 °C, which it does not describe
     - **Deposit**: the authors' heat-map script. It is not in the paper.
       Nine of its ten 50 °C constants are values the deposit gives; the
       tenth, ``k_3``, is **fitted** to the 50 °C data of the script.

The two Table 3 models describe the batch of the Fig. 9 caption: 60 mM 3-HP,
NaOH/3-HP = 1:1, 0.2 g Pd/C, 50 °C, and 3 bar O2 or 15 equivalents of H2O2
under 3 bar Ar. Their species are ``HP3``, ``OPA`` (3-oxopropanoic acid, the
aldehyde intermediate, which the paper reports as not detected), ``MA``,
``AAc`` (acetic acid), ``OXA`` (oxalic acid) and ``D`` (lumped other
by-products). The quantities plotted in the paper are assignment rules of
the models and are selected by name: ``X_HP3``, ``X_MA``, ``X_AAc``,
``X_OXA`` and ``X_D`` are each species over the sum of the five measured
species (OPA is left out of the sum, as in the authors' post-processing).

**Units.** Time is in min and concentrations are in M. Temperatures are in
K wherever an argument or a model parameter is named ``T``, and in °C
wherever a name ends in ``_C`` (``T_C``, ``temperature_C``). The rate
constants ``k_1`` to ``k_6`` of the two Table 3 models are in
1/(g\ :sub:`cat` min), as printed in Table 3, and ``k_7`` of the H2O2 model
is in 1/min. The rate constants of the simplified model are lumped for the
paper's batch and are in 1/min.

Using the models
----------------

Load a model and simulate it. Every loader call builds a new, independent
model:

.. code-block:: python

   from nskinetics.models.pd_c_3hp_oxidation import load_o2_model

   km = load_o2_model()                 # Table 3 O2 model at 50 °C
   km.simulate(0, 1440, 1441)           # 0 to 1440 min, 1441 output points
   print(km.results_df.iloc[[0, 60, 240, 600]].to_string(
       index=False, float_format='%.4g'))

.. code-block:: text

    time     [HP3]     [OPA]    [MA]     [AAc]     [OXA]      [D]
       0      0.06         0       0         0         0        0
      60   0.03389  0.001866 0.01742  3.25e-05 3.064e-05 0.006766
     240  0.006109 0.0003365 0.03374 0.0003618 0.0003411  0.01911
     600 0.0001986 1.094e-05 0.02519  0.001051 0.0009906  0.03256

``km.simulate`` continues from the current state of the model: a second
call starts where the first one ended, although its time column starts
again at the start time given. To start a new run, call ``km.reset()``
first, pass ``reset=True`` to ``simulate``, or use ``simulate_at_times``.

``simulate_at_times`` resets a model and returns the plotted quantities (or
any other selection) at exactly the times asked for, for example the
sampling times of an experiment:

.. code-block:: python

   from nskinetics.models.pd_c_3hp_oxidation import (
       NORMALIZED, pooled_r2, read_fig9b_o2, simulate_at_times)

   sim = simulate_at_times(km, [0, 60, 240, 600],
                           list(NORMALIZED) + ['[OPA]'])
   print(sim.to_string(index=False, float_format='%.4g'))

   experiment = read_fig9b_o2()['experiment']     # deposited points, Fig. 9b
   model = simulate_at_times(km, experiment['time'])
   print(f'pooled R2 = {pooled_r2(experiment, model):.4f}')

.. code-block:: text

    time   X_HP3   X_MA     X_AAc     X_OXA    X_D     [OPA]
       0       1      0         0         0      0         0
      60   0.583 0.2996 0.0005591 0.0005271 0.1164  0.001866
     240  0.1024 0.5655  0.006065  0.005717 0.3204 0.0003365
     600 0.00331   0.42   0.01752   0.01651 0.5427 1.094e-05
   pooled R2 = 0.9945

The H2O2 model adds the state ``H2O2`` and the plotted quantity ``eq_H2O2``,
the H2O2 left in equivalents of the initial 3-HP (SI Fig. S12a):

.. code-block:: python

   from nskinetics.models.pd_c_3hp_oxidation import load_h2o2_model

   sim = simulate_at_times(load_h2o2_model(), [0, 60, 600, 1440],
                           ['X_HP3', 'X_MA', 'eq_H2O2'])
   print(sim.to_string(index=False, float_format='%.4g'))

.. code-block:: text

    time  X_HP3   X_MA   eq_H2O2
       0      1      0        15
      60 0.8331 0.1286      9.88
     600 0.3511  0.397    0.2433
    1440 0.2078 0.3551 0.0008299

The two Table 3 models are at 50 °C and take no temperature. The simplified
reversible model takes one **in K**
(``load_o2_simplified_reversible_model(T=...)``, or ``km.set_value('T',
...)`` on a loaded model), which scales its ten constants by the Arrhenius
factors of the authors' script. It also takes one of the two initial
conditions of that script: ``'heatmap'`` (the default, as for Fig. 11) or
``'fitting'`` (as when the script fits the 50 °C data). The example
compares it with the 50 °C experiment it was fitted to and with the 75 °C
experiment of Fig. 10b. Read `Limits of the temperature dependence`_ before
using it away from 50 °C:

.. code-block:: python

   from nskinetics.models.pd_c_3hp_oxidation import (
       load_o2_simplified_reversible_model, read_fig10_o2)

   km_75 = load_o2_simplified_reversible_model(
       T=75 + 273.15, initial_condition='fitting')         # kelvin
   factor = km_75.get_value('k_1')/km_75.get_value('k_1_ref')
   print(f'k_1(75 C)/k_1(50 C) = {factor:.2f}')

   for label, model, experiment in (
           ('50 C', load_o2_simplified_reversible_model(
               initial_condition='fitting'), read_fig9b_o2()['experiment']),
           ('75 C', km_75, read_fig10_o2(75)['experiment'])):
       sim = simulate_at_times(model, experiment['time'])
       print(f'{label}: pooled R2 = {pooled_r2(experiment, sim):.3f}')

.. code-block:: text

   k_1(75 C)/k_1(50 C) = 4.94
   50 C: pooled R2 = 0.994
   75 C: pooled R2 = 0.641

A temperature given to the loader in the wrong unit is refused:
``load_o2_simplified_reversible_model(T=75)`` raises a ``ValueError`` (``T``
must lie within 250 to 450 K), and ``simulate_heat_maps`` refuses
temperatures outside −50 to 200 °C. These ranges only catch that mistake;
they are not ranges of validity. ``set_value('T', ...)`` is not checked.

Set and read model values with ``set_value`` and ``get_value``, not through
attributes of the underlying RoadRunner object: attribute access is
unreliable once more than one model has been loaded in a process.

What a loader selects (a temperature, a constant set, an initial condition)
is part of the model it builds, so it survives both ``reset()`` and
``compile_events()``. Values that *you* set afterwards with ``set_value``
behave as for any :class:`~nskinetics.KineticModel`: ``reset()`` restores
the species and keeps a changed parameter such as ``T`` or ``m_cat``, and
``compile_events()`` regenerates the model and discards them, so set values
after compiling events. A parameter changed by an event during a simulation
also stays changed after ``reset()``: set it back with ``set_value`` before
the next run.

Changing the operating conditions
---------------------------------

The parameters below are meant to be changed with ``set_value``. Everything
else in a model (the constants, ``x_per_M``, ``x_O2_ref``) is what the
paper or the deposit gives, or what was fitted to the 50 °C data, and
should be left alone unless you know why you are changing it.

.. list-table::
   :header-rows: 1
   :widths: 20 12 10 10 48

   * - Model
     - Parameter
     - Unit
     - Default
     - Meaning
   * - H2O2 model, O2 model
     - ``m_cat``
     - g
     - 0.2
     - Catalyst mass (10 % Pd/C) **of the paper's batch**, not a mass per
       litre. The rates of reactions 1 to 6 are proportional to it.
   * - all three
     - ``C0``
     - M
     - 0.06
     - Initial 3-HP concentration. Takes effect at the next ``reset()``
       (``simulate_at_times`` resets).
   * - H2O2 model
     - ``eq_H2O2_0``
     - –
     - 15
     - Initial H2O2 in equivalents of ``C0`` (0.9 M). The initial H2O2
       concentration is ``eq_H2O2_0*C0``. Takes effect at the next
       ``reset()``, like ``C0``.
   * - simplified reversible model
     - ``T``
     - K
     - 323.15
     - Temperature (50 °C by default). The model was fitted at 50 °C only
       and does not describe the experiments at 25, 75 and 100 °C (`Limits
       of the temperature dependence`_).
   * - simplified reversible model
     - ``P``
     - bar
     - 3.5
     - Pressure of the authors' heat-map script (the captions of the paper
       give 3 bar of O2). Do not use it to explore pressure: the model was
       fitted at no other pressure, and the paper's Fig. 3a contradicts its
       trend above 3 bar (see below).

.. code-block:: python

   km = load_o2_model()
   for m_cat in (0.1, 0.2, 0.4):               # g in the paper's batch
       km.set_value('m_cat', m_cat)
       sim = simulate_at_times(km, [60, 240], ['X_HP3', 'X_MA'])
       print(f'm_cat = {m_cat} g: X_HP3 = {sim["X_HP3"][0]:.3f}, '
             f'{sim["X_HP3"][1]:.3f}; X_MA = {sim["X_MA"][0]:.3f}, '
             f'{sim["X_MA"][1]:.3f}')

   km = load_h2o2_model()
   for equivalents in (5, 15, 30):
       km.set_value('eq_H2O2_0', equivalents)
       sim = simulate_at_times(km, [600], ['X_HP3', 'X_MA'])
       print(f'{equivalents} equivalents of H2O2: X_HP3 = '
             f'{sim["X_HP3"][0]:.3f}, X_MA = {sim["X_MA"][0]:.3f} at 600 min')

.. code-block:: text

   m_cat = 0.1 g: X_HP3 = 0.783, 0.325; X_MA = 0.153, 0.471
   m_cat = 0.2 g: X_HP3 = 0.583, 0.102; X_MA = 0.300, 0.565
   m_cat = 0.4 g: X_HP3 = 0.325, 0.010; X_MA = 0.471, 0.481
   5 equivalents of H2O2: X_HP3 = 0.600, X_MA = 0.110 at 600 min
   15 equivalents of H2O2: X_HP3 = 0.351, X_MA = 0.397 at 600 min
   30 equivalents of H2O2: X_HP3 = 0.153, X_MA = 0.571 at 600 min

Read this before changing them:

- **Only one condition per model was ever fitted.** The Table 3 constants
  were fitted by the authors to one experiment each: 60 mM 3-HP, 0.2 g of
  catalyst, 50 °C, with 15 equivalents of H2O2 (samples up to 1440 min) or
  3 bar O2 (samples up to 1680 min). Any other value of ``m_cat``, ``C0``,
  ``eq_H2O2_0``, ``T`` or ``P``, and any time beyond the last sample, is an
  extrapolation of the fitted rate laws. The example above shows what the
  models compute, not what was measured. For ``T`` the comparison with the
  time courses of Fig. 10 is poor (`Limits of the temperature
  dependence`_). For ``m_cat`` and ``C0`` no experiment was found in the
  paper or the deposit to compare with.
- **The paper's screening experiments contradict some of these
  extrapolations.** The time courses behind the models are those of Fig. 9;
  the O2 one is at 3 bar of O2. The paper also has experiments
  analysed once, after 8 h, in which the H2O2 charge, the temperature or
  the O2 pressure was varied (its Figs. 2 and 3; SI Table S2 has the
  temperature series with O2, and SI Table S1 other catalysts), and
  experiments with malonic acid alone (its Fig. 4). The deposit holds the
  values of Figs. 2 to 4. The models were not built on these experiments
  and do not describe them (``NOTES.md`` has the table; the model values
  are pinned by the package tests). The 3-HP conversions and MA yields of
  the models below are ``1 - [HP3]/C0`` and ``[MA]/C0``, the quantities of
  the paper's Eqs. 1 and 3:

  - *O2 pressure* (Fig. 3a: 1, 3, 5, 10 and 20 bar at 50 °C). The measured
    conversion is 21.7, 72.8, 63.2, 25.0 and 25.9 %, and the MA yield 2.7,
    41.1, 20.5, 1.3 and 0.6 %: both are highest at 3 bar, and the paper
    calls the decline "likely due to catalyst poisoning by excess O2". (The
    deposited 5 bar values equal an entry of SI Table S1 for another
    catalyst under 3 bar, see ``NOTES.md``; the 10 and 20 bar values show
    the decline as well.) The simplified reversible model was fitted at no
    other pressure, and its rate form only speeds up with ``P``: after
    480 min it gives 87.7, 98.2, 99.7, 99.9 and 100.0 % conversion with an
    MA fraction of 0.30, 0.48, 0.51, 0.53 and 0.51 (and at 240 min ``X_MA``
    is 0.27, 0.57 and 0.70 for 1, 3.5 and 10 bar). Its pressure dependence
    is therefore contradicted by the paper, not merely untested: **do not
    use** ``P`` **to explore pressure.** The O2 model has no pressure
    parameter, and ``x_O2_ref`` is not one.
  - *H2O2 charge* (Fig. 2a). The measured conversion and MA yield are
    highest at 15 equivalents (78.4 and 58.5 %) and fall beyond (56.4 and
    28.4 % in the next experiment, at 20 equivalents by the axis of the
    figure and 30 by the deposit). In the H2O2 model the conversion rises
    with every equivalent added, and the MA yield after 8 h rises up to
    42 equivalents (59.0 %) and falls beyond (58.2 % at 50): the model's
    optimum is far beyond the measured one.
  - *Temperature with O2* (Fig. 3c). After 8 h at 25 °C the simplified
    reversible model has converted 53 % where 8.5 % was measured; at 75 °C
    it has an MA yield of 6.0 % where 35.7 % was measured.
  - *Malonic acid alone* (Fig. 4). Without 3-HP, 3.8 % (Fig. 4a) or 1.2 %
    (the time course of Fig. 4b) of the MA was converted after 1300 min at
    50 °C; the O2 model started from 60 mM MA converts 80 %.
  - *Catalyst reuse* (Fig. 6b: three cycles at 50 °C and 3 bar O2). This
    is not a screening experiment: up to 120 min the deposited 3-HP and MA
    values of its first cycle are those of Fig. 9b. At 120 min the
    deposited 3-HP and MA values (normalized concentrations in the figure)
    are 0.255 and 0.476 in the first cycle, 0.652 and 0.092 in the second
    and 0.769 and 0.002 in the third, where the O2 model has ``X_HP3`` =
    0.325 and ``X_MA`` = 0.471: the models have no term for this loss of
    activity and do not describe a reused catalyst.

  The screening experiments are not on the same footing as the time
  courses either. At the same nominal conditions they give 72.8 % (O2) and
  78.4 % (H2O2) conversion after 480 min, where the time courses of Fig. 9
  are at 96.1 % after 450 min and 60.6 % after 420 min. No explanation was
  found in the paper; it states that for time-resolved runs substrate,
  NaOH and water were tripled and that at every sample 1 mL was discarded
  and 0.5 mL sampled (1.5 mL in all).
- **The catalyst mass is per batch of the paper, not per litre.** The
  compartment of the model files is a nominal 1 L that cancels from every
  concentration, and no rate has a volume term. The rates scale with
  ``m_cat/0.2``, the catalyst loading relative to the paper's:
  ``m_cat = 0.4`` is twice the paper's loading whatever the liquid volume
  (in the example, 0.4 g reaches at 60 min what 0.1 g reaches at 240 min),
  and ``m_cat = 20`` is not "20 g/L" but a hundred times the paper's
  loading. The paper gives its batch as 0.2 g in 10 mL of liquid (20 g/L)
  and adds that for time-resolved runs substrate, NaOH and water were
  tripled, without saying whether the catalyst was, so the loading of the
  kinetic runs is either 20 or 6.7 g/L.
- **What depends on the initial concentration.** In the O2 model the
  plotted fractions hardly depend on ``C0`` (halving or doubling it changes
  them by less than 1e-4), because every rate is first order in its liquid
  species. In the H2O2 model, changing ``C0`` alone also changes the H2O2
  concentration and with it the rates. In the simplified reversible model
  ``x_per_M`` has to be changed together with ``C0`` (see the loader's
  docstring).
- **The simplified reversible model has lumped constants.** It has no
  catalyst-mass parameter; its rate constants are effective ones for the
  paper's batch.
- **Inputs are not validated.** ``set_value`` is the generic NSKinetics
  call and checks nothing: zero or negative ``m_cat``, ``C0`` or
  ``eq_H2O2_0`` run without an error and give meaningless results.

Heat maps and figures
---------------------

The heat maps of Fig. 11 are computed by ``fig11_heat_maps``, which runs the
simplified reversible model as the authors' script does, once per
temperature (71 integrations), with its default constants. Its time axis is
the *plotted* time of the published figure, which is the reaction time plus
10 min (``FIG11_TIME_OFFSET``).
``simulate_heat_maps`` does the same for any model with a temperature
parameter on any grid; its time axis is the reaction time unless a
``time_offset`` is given, and its temperatures are in °C:

.. code-block:: python

   from nskinetics.models.pd_c_3hp_oxidation import (
       FIG11_TIME_OFFSET, fig11_heat_maps, heat_map_maximum,
       load_o2_simplified_reversible_model, simulate_heat_maps)

   T_C, times, selectivity, MA_yield = fig11_heat_maps()   # 71 x 300 grids
   for matrix in (selectivity, MA_yield):
       value, T_max, t_max = heat_map_maximum(T_C, times, matrix)
       print(f'{value:.4f} at {T_max:g} C, plotted time {t_max:.1f} min '
             f'(reaction time {t_max - FIG11_TIME_OFFSET:.1f} min)')

   # rows: 50, 75 and 100 C; columns: reaction times of 30, 60, 120, 240 min
   T_C, times, selectivity, MA_yield = simulate_heat_maps(
       load_o2_simplified_reversible_model(), T_C=[50, 75, 100],
       times=[30, 60, 120, 240])
   print(selectivity.round(3))

.. code-block:: text

   0.5931 at 95 C, plotted time 26.8 min (reaction time 16.8 min)
   0.5484 at 84 C, plotted time 49.1 min (reaction time 39.1 min)
   [[0.189 0.324 0.483 0.571]
    [0.535 0.584 0.449 0.226]
    [0.441 0.207 0.045 0.004]]

The figure module writes every figure of this page, including the
verification figure ``solver_emulation`` (see `Repeating the authors' solver
run`_), as PNG (600 dpi) and PDF, to a folder; ``--figures``, ``--formats``
and ``--dpi`` select a subset, the formats and the resolution of the PNG
files:

.. code-block:: bash

   python -m nskinetics.models.pd_c_3hp_oxidation.figures OUTPUT_DIR
   python -m nskinetics.models.pd_c_3hp_oxidation.figures OUTPUT_DIR \
       --figures fig9 fig11 --formats png --dpi 150

In Python, ``figures.plot_fig9()`` and its siblings return the matplotlib
figure, and ``figures.draw_fig9a(ax)`` and its siblings draw one panel on an
axes of your own.

Figures drawn from the models
-----------------------------

All figures below are drawn by the package from the models. Deposited and
digitized data are overlaid where the package does that; no image of the
paper is reproduced here. "Normalized concentration" is the plotted fraction
``X_i``, and all differences are in that unit.

Fig. 9 and SI Fig. S12: H2O2 and O2 at 50 °C
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig9.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig9.png
   :alt: Two panels of normalized concentration against time from 0 to 1440
         min, each next to its reaction list. Panel a, H2O2: 3-HP falls to
         about 0.2, malonic acid rises to a plateau near 0.4 and the other
         by-products rise to about 0.33. Panel b, O2: 3-HP is consumed within
         about 500 min, malonic acid peaks near 0.57 at about 250 min and
         falls to 0.15, and the other by-products rise to about 0.79. In
         both panels the solid model lines lie on the black dotted deposited
         curves; circles are the experimental points.
   :class: pub-figure
   :width: 100%

   Fig. 9, drawn by ``plot_fig9``: experiment (circles), the NSKinetics
   models (solid lines) and the model curves the authors deposited (black
   dotted lines), with (a) H2O2 and (b) O2.

At the nodes of the deposited curves after t = 0 the models differ from
them by at most 9.6e-4 (H2O2, 23 nodes) and 2.0e-3 (O2, 19 nodes). These
residuals are the solver error of the deposited curves (see `Published lines
carry a solver error`_) and, for O2, the rounding of Table 3; for H2O2 the
constants behind the deposited curve differ from Table 3 by up to 0.6 %,
more than rounding allows for five of them. With the authors' solver run
repeated, the models come within 4.4e-4 (H2O2) and 4.0e-4 (O2) of the nodes
with Table 3 as printed, and within 6.8e-6 and 5.0e-6, the print rounding
of the deposited files, with constants 0.01 to 0.6 % off Table 3 (see
`Repeating the authors' solver run`_); the rate laws are the same. Against
the experimental points the pooled R2 is 0.9945 for O2 (16 rows) and, for
H2O2, 0.945 on the 11 rows plotted in SI Fig. S11b and 0.941 on all 13 (see
`Parity R2 values (SI Fig. S11)`_).

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_figS12.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_figS12.png
   :alt: The two profile panels of Fig. 9 with insets. Panel a inset: the
         predicted H2O2 falls from 15 equivalents to near zero within about
         700 min; the dashed model line lies slightly above the digitized
         published curve and the deposited trace. Panel b inset: the
         predicted OPA mole fraction peaks at 1.4e-5 near 20 min and decays
         to zero by 600 min; the dashed model line passes through the
         digitized vertices of the published curve.
   :class: pub-figure
   :width: 100%

   SI Fig. S12, drawn by ``plot_figS12``: the predicted (a) H2O2 and
   (b) OPA as dashed lines, in the main axes and in the insets. The insets
   add the published inset curves, digitized from the SI figure, and in (a)
   the H2O2 trace the authors deposited.

The predicted OPA of the O2 model agrees with the 13 digitized vertices of
the published inset to 3.0e-7 on a peak of 1.4e-5. The predicted H2O2 of the
Table 3 model does not reproduce the published trace: it lies up to 0.38
equivalents (rms 0.22) above it (see `The H2O2 curve of SI Fig. S12a is
matched only with other constants`_).

Fig. 11: heat maps over temperature and time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11.png
   :alt: Two heat maps over time (10 to 1680 min) and temperature (30 to 100
         degrees C), panel a for malonic acid selectivity and panel b for
         yield, on a colour scale from 0 to 1. In both, a ridge of values
         above 0.5 runs from short times at high temperature to long times
         at low temperature, marked by a red dashed line, with contour lines
         at 0.1, 0.3 and 0.5 on both sides. A red star marks the largest
         cell, 0.593 at 95 degrees and 27 min in panel a and 0.548 at 84
         degrees and 49 min in panel b; a purple star marks the maximum of
         the 50 degree row. An inset in each panel zooms in on 47 to 57
         degrees.
   :class: pub-figure
   :width: 100%

   Fig. 11, drawn by ``plot_fig11``: (a) MA "selectivity" (``X_MA``) and
   (b) MA "yield" (``(1 - X_HP3)*X_MA``), computed with the simplified
   reversible model as the authors' heat-map script runs it, both panels
   with the model's default constants. These are not the published maps:
   the constants that the script's run fitted are not in the deposit, and
   the contour lines lie 1.05 px (a) and 0.30 px (b) rms from the published
   ones (next figure). Two things differ from the published figure
   on purpose: the largest cell of each map is annotated with its value,
   and the colour-bar label of panel a is spelled "MA Selectivity"
   (published: "MA Selectvity"). The stars sit on the computed maxima, and
   the insets are zooms of the same matrices (``NOTES.md`` gives the
   differences from the published stars and insets).

The deposit holds no table of the heat maps, so the reference is the
published figure: its contour lines at 0.1, 0.3 and 0.5, digitized. The
model was not fitted to them. Its contours lie within an rms of 0.30 px of
those of panel b and 1.05 px of those of panel a, in pixels of the published
raster (1 px = 3.07 min or 0.148 °C); the digitization itself is good to
about 0.3 px. So panel b is matched to about the reading uncertainty and
panel a is not: its late 0.1 and 0.5 contours are 1.26 and 1.28 px rms off,
single samples up to 3.0 px. The red dashed line, which joins the row
maxima and is therefore defined to one column of the time grid (1.82 px),
lies 0.25 px (panel a) and 0.29 px (panel b) from the digitized one in the
mean along the time axis. The largest cell of the selectivity map is the
one the text names (95 °C, 27 min), with 0.593 where the text has 59.1 %.
The largest cell of the yield map is also the one the text names (84 °C,
49 min), with 0.548 where the text has 55.2 %.

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11_verification.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11_verification.png
   :alt: Top row: the two computed heat maps with the digitized published
         contour lines, ridge line and stars drawn on top as fine black
         dots. In panel b the dots lie on the wide computed lines; in panel a
         they leave them slightly on the late side of the ridge. Bottom row:
         the offset of every digitized contour sample from the computed
         contour against the temperature of the sample. In panel b all
         offsets lie within 1.7 px, with an rms of 0.30 px. In panel a the
         late contours run from about minus 3.0 px to plus 2.2 px, with an
         rms of 1.05 px.
   :class: pub-figure
   :width: 100%

   Drawn by ``plot_fig11_verification``. Top: the computed maps with the
   lines and stars digitized from the published panels in black; the title
   of each panel says whether the published contours are matched. Bottom:
   the offset of each digitized contour sample from the computed contour,
   positive where the computed map is above the level of the contour at
   the sample: on the late side of the red dashed line the published line
   then lies earlier than the computed one, on the early side later. The
   misfit of panel a (1.05 px rms) cannot be removed for both panels at
   once by any single set of constants: the two published panels do not
   come from one run (see `The constants of the simplified reversible
   model`_).

Limits of the temperature dependence
------------------------------------

.. warning::

   The simplified reversible model is the only model of the package with a
   temperature, and its temperature dependence is an extrapolation from
   50 °C that the paper's own experiments at 25, 75 and 100 °C do not
   support. Use ``load_o2_simplified_reversible_model`` away from 50 °C
   only as what it is: the model of the authors' heat-map script, which is
   what Fig. 11 shows.

The model was fitted to the 50 °C experiment only. Its ten constants at
another temperature are the 50 °C constants times the Arrhenius factors of
the script. All numbers below are for the model as shipped and are pinned by
the package tests.

At 50 °C the model describes the experiment it was fitted to: pooled R2
0.9939 against the 16 deposited rows of Fig. 9b in the script's fitting
convention (the Table 3 O2 model has 0.9945). Against the deposited
experimental points of Fig. 10a–c (all rows of each file) it does not:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * -
     - 25 °C
     - 75 °C
     - 100 °C
   * - Pooled R2 / rms difference, fitting convention
     - 0.653 / 0.220
     - 0.641 / 0.175
     - 0.757 / 0.147
   * - Pooled R2 / rms difference, heat-map convention
     - 0.561 / 0.248
     - 0.637 / 0.176
     - 0.758 / 0.146
   * - Time at which ``X_HP3`` falls to 0.5: model / experiment
     - 475 min / not reached (0.876 left at 1890 min)
     - 15 min / between 60 and 90 min
     - 3.9 min / between 20 and 40 min
   * - Maximum MA fraction: model / experiment
     - 0.427 at 1203 min / 0.023 at 1890 min
     - 0.586 at 51 min / 0.436 at 180 min
     - 0.594 at 13 min / 0.568 at 40 min

The model is too fast at all three temperatures. At 25 °C it leaves 0.16 of
the 3-HP after 1890 min where 0.876 was measured. At 75 °C the measured
3-HP conversion is no faster than at 50 °C (``X_HP3`` falls to 0.5 between
60 and 90 min at both temperatures), whereas the Arrhenius factor of the
first step from 50 to 75 °C is 4.9. The paper's own Arrhenius plot
(Fig. 10d) shows the same: its 50 °C point of 3-HP → MA lies a factor of 4.3
above the line through the 25, 75 and 100 °C points, so a model anchored at
50 °C with the activation energy of the text overshoots at the other three
temperatures (`Fig. 10d and the activation energies of the text`_).

The 8 h screening experiments at other temperatures and pressures are not
described either (`Changing the operating conditions`_).

**No model curve of Fig. 10a–c is reproduced by this package.** The paper
does not say which model those curves are and prints no constants for
them. Neither script of the deposit produces them as posted: the
heat-map script fits the 50 °C data, and the deposited 75 and 100 °C
curves hold 3-HP constant at a value above zero for over 500 min, which
the irreversible 3-HP steps of the posted O2 script cannot do (``NOTES.md``
has the numbers). The deposited tables of Fig. 10 are shipped as reference
data (``read_fig10_o2``, ``read_fig10d_arrhenius``); the package uses their
experimental points for the comparison above and the 50 °C row of the
Arrhenius table for three constants of the simplified reversible model.

What this means for Fig. 11: its 50 °C row rests on data; its other rows
are the extrapolation described here. That holds for the published maps as
for the maps of this package.

What does not match the paper, and what is inferred
---------------------------------------------------

``NOTES.md`` in the package folder gives the evidence for each point below
and lists further, smaller ones. How the authors produced a figure is known
only where the paper or the deposit says so; elsewhere this page says what
the numbers show.

SI Fig. S2 and Eq. S1 (rate expressions)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The paper refers to SI Fig. S2 for its rate expressions. Fig. S2 is a flow
chart of the fitting algorithm that states: the rate form
``r_j = k_jf*prod(a_i^n_ij)``; the activities ``a_i = gamma_i*x_i`` for a
liquid species and ``a_i = phi_i*y_i*P/f0`` for a gas species; the bounds
``k >= 0`` and ``0 <= n <= 3``; and that the reactions are irreversible
because ``Keq >> 1``. Its caption gives the mole balance, Eq. S1:
``dN_i/dt = N dx_i/dt = sum_j v_ij r_j``, in mole fractions. Section 2.6
adds that gas and liquid are taken as ideal (unit activity coefficients).
Fig. S2 also defines the R2 that ``pooled_r2`` computes.

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * - Model
     - Where it is the Fig. S2 form
     - Where it departs
   * - H2O2 model
     - Irreversible ``k*prod(a^n)`` with ``n = 1`` for every species (the
       fitted orders of Table S4 are not used, see below) and Table 3
       constants (times ``m_cat``, except ``k_7``); liquid activities are
       mole fractions, H2O2 included.
     - The mole fractions are on the authors' scale, 0.000339 for the 60 mM
       feed (in water it would be 1.08e-3).
   * - O2 model
     - The same for the liquid species.
     - The same, and **the O2 activity**: Fig. S2 gives ``y*P/f0`` = 2.96
       for 3 bar of O2; the model has the inferred constant ``x_O2_ref`` =
       0.005085 (times ``yO2``, which stays above 0.999).
   * - Simplified reversible model
     - Liquid activities are the script's mole fractions and the gas
       activities are ``y*P/f0``, with the script's ``P`` = 3.5 bar.
     - Every step has a reverse term with an "equilibrium constant" of 0.6
       to 82 (not the ``Keq`` of Table S4); no OPA; lumped rate constants.

With the literal O2 activity of Fig. S2 the Table 3 constants would make
reactions 1, 2 and 4 of the O2 model 582 times faster: 3-HP would have a
half-life of 0.16 min instead of 73 min, and the model would miss the
deposited Fig. 9b curve by up to 0.97. The authors' posted O2 script does
not apply ``P/f0`` either (it defines ``P`` and ``f0`` and uses neither).

``x_O2_ref`` is inferred
~~~~~~~~~~~~~~~~~~~~~~~~

The Table 3 O2 constants of R1, R2 and R4 reproduce the deposited Fig. 9b
curve only when multiplied by the catalyst mass (0.2 g) and by a further
factor of about 0.0051, taken as 0.005085. That value equals the H2O2 mole
fraction of the 15-equivalent H2O2 run on the authors' scale
(15 × 0.000339). It is not in the paper, the SI, the README of the deposit
or the posted script, and how it entered the authors' code is unknown. Only
the products ``k_i*m_cat*x_O2_ref`` are identifiable. The posted fitting
script, run as posted, misses the deposited curve by an rms of 0.42 to 0.76
for 3-HP, MA, OXA and D. Solved accurately, the model cannot tell 0.005085
from a rounded 0.0051 (both are within the solver error of the deposited
curve). The choice rests on the repeated solver run (`Repeating the authors'
solver run`_): the constants with which it returns the deposited nodes,
divided by Table 3 as printed and by 0.2 g, give 0.0050854, 0.0050733 and
0.0050857 for R1, R2 and R4, and since a printed constant stands for half a
unit of its third digit either way, the three allow 0.0050818 to 0.0050889
and not 0.0051. With 0.0051 the repeated run misses the nodes by 9.4e-4,
and still by 7.0e-4 when the three constants are set to the lowest values
Table 3 can stand for and the other three keep their node-set values
(5.8e-4 with those as printed in Table 3; limit 1e-5). The package tests
compute all three. This holds if Table 3 is the rounded set of the
constants behind the deposited curve; the fit that gave those constants is
not shipped.

With this value the oxidant factor of the O2 model equals that of the H2O2
model at t = 0. If the factor is read as an O2 mole fraction, the O2 model
is the hypothetical condition of equal oxidant mole fractions that the
paper discusses, and for it the paper states a rate ratio that the models
do not give (see `The initial-rate ratio of 0.204 is not reproduced`_).

Table S4 is not used
~~~~~~~~~~~~~~~~~~~~

SI Table S4 lists fitted reaction orders and DFT equilibrium constants.
Neither enters the rate laws. That the steps are irreversible is printed in
the paper (Section 3.4: the equilibrium constants are "high enough to
consider all reactions as irreversible"). That the orders are exactly 1 is
inferred: the posted script has no orders, and the deposited curves are
reproduced by first-order steps and not with the Table S4 orders. Applied to
the organic reactant of each step in an independent integration (computed
outside the package, not shipped), those orders move the curves off the
deposited ones by 0.010 to 0.17 (O2) and 0.037 to 0.16 (H2O2), depending on
the scale on which ``x^n`` is taken: the smallest values are for fractions
of the feed, and 0.17 and 0.14 are for the authors' mole fractions, the
scale of Fig. S2. On the feed-normalized scale the Table S4 orders fit the
experiments about as well as first order does (R2 0.994 against 0.9945 for
O2, 0.941 for both with H2O2), so it is the deposited curves, not the
experiments, that exclude them: 0.010 is five times the solver error of the
deposited O2 curve. The values are kept as ``TABLE_S4`` for documentation.

Published lines carry a solver error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The deposited model curves of Fig. 9 have 19 (O2) and 23 (H2O2) unevenly
spaced solver steps (20 and 24 nodes with the one at t = 0), and the
published lines are straight segments between the nodes. A Python port of
MATLAB ``ode23s`` at its default tolerances returns those nodes to their
printed digits (next section), so they are consistent with being that
solver's raw steps; MATLAB itself was not run. They are not converged
solutions: the nodes are up to 2.0e-3 (Fig. 9b) away from the accurately
solved model, and the segments between them up to 0.010. The models here
are solved to a relative tolerance of 1e-9, so they do not lie exactly on
the published polylines, and the package tests that solve the models
accurately allow for that solver error.

Repeating the authors' solver run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``solver_emulation`` compares the models with the deposited nodes more
closely than an accurate solution can. Its ``ode23s`` is a port of the
MATLAB solver (the modified Rosenbrock (2,3) pair of Shampine and Reichelt,
*SIAM J. Sci. Comput.* 18 (1997) 1–22) with the default options of the
authors' posted script, written from the paper and the documented behaviour
of the MATLAB function; MATLAB itself was not run. ``emulate_authors_run``
integrates the state vector of the posted script with it and takes every
rate of change from a loaded model, so that the rate laws, the constants
and the stoichiometry are those of the shipped Antimony file. It is for
verification only: to simulate a model, use the model.

.. code-block:: python

   from nskinetics.models.pd_c_3hp_oxidation import (
       emulate_fig9a, emulate_fig9b, node_deviations, read_fig9a_h2o2)

   def report(label, emulated, nodes):
       d = node_deviations(emulated, nodes)
       print(f"{label}: {d['steps']} steps (deposited {d['deposited_steps']}),"
             f" node times within {d['time']:.0e} min,"
             f" values within {d['max']:.1e}")

   nodes = read_fig9b_o2()['model']
   report('Fig. 9b, constants of the nodes', emulate_fig9b(), nodes)
   report('Fig. 9b, Table 3 as printed', emulate_fig9b('table3'), nodes)
   nodes = read_fig9a_h2o2()['model']
   report('Fig. 9a, constants of the nodes', emulate_fig9a(), nodes)
   report('Fig. 9a, Table 3 as printed', emulate_fig9a('table3'), nodes)

.. code-block:: text

   Fig. 9b, constants of the nodes: 19 steps (deposited 19), node times within 6e-06 min, values within 5.0e-06
   Fig. 9b, Table 3 as printed: 19 steps (deposited 19), node times within 4e-02 min, values within 4.0e-04
   Fig. 9a, constants of the nodes: 23 steps (deposited 23), node times within 8e-06 min, values within 6.8e-06
   Fig. 9a, Table 3 as printed: 23 steps (deposited 23), node times within 5e-01 min, values within 4.4e-04

The deposited files print values and node times to five decimals, so the
differences of the first and third run are their print rounding. The same
holds for the H2O2 trace of SI Fig. S12a (24 steps, of which the file holds
all but the last; within 4.4e-6 equivalents). The package tests require the
deposited number of steps, node times within 2e-5 min and values within
1e-5; the solver error that the accurate-solver tests have to allow for is
85 to 200 times larger for the plotted fractions (8.5e-4 for Fig. 9a and
2.0e-3 for Fig. 9b; for the H2O2 trace it is 0.012 equivalents).

The constants matter at this level. They are sets inferred from the nodes
(``O2_NODE_PARAMETER_SETS``, ``H2O2_NODE_PARAMETER_SETS``; selected with
``parameter_set`` in the loaders), given to eight digits. They are not the
published constants and are for this comparison only: the O2 set differs
from Table 3 by 0.01 to 0.26 % and prints as Table 3 prints it; the H2O2
set differs by 0.1 to 0.6 %, more than rounding allows for five of its
constants. The models' defaults are Table 3 as printed, with which the
repeated run comes within 4.0e-4 (O2) and 4.4e-4 (H2O2) of the nodes.

The tests also change one inferred choice at a time, with the constants
kept, and check that the repeated run then misses the nodes:

.. list-table::
   :header-rows: 1
   :widths: 46 54

   * - Change
     - Largest difference from the nodes (limit 1e-5)
   * - ``x_O2_ref`` = 0.0051 instead of 0.005085 (Fig. 9b)
     - 9.4e-4; 7.0e-4 with the most favourable constants Table 3 can stand
       for
   * - O2 fraction held at 1 in the rates, as in the posted script (Fig. 9b)
     - 7.7e-5
   * - ``k_7`` multiplied by the catalyst mass (Fig. 9a)
     - 0.39 (16 steps instead of 23)
   * - OPA included in the sum of the plotted fractions
     - 0.038 (Fig. 9b), 0.094 (Fig. 9a)

What this shows is limited. The shipped constants return the deposited
nodes, and most of them are pinned on their own (to between 1e-6 and 1e-3
of their value; ``NOTES.md`` has the table). The exception is ``k_6`` of
the H2O2 model, which is not identifiable: a change of 10 % leaves the
nodes as they are. A control shows that the shipped constants need the
shipped choice. It does not show that another choice could not do as well
with other constants: that takes a refit of the alternative, and the refits
that were made are not shipped (`What is not in the repository`_). For the
O2 fraction, for example, a refit with the fraction held at 1 still left a
difference of 1.5e-5, three times that of the shipped model.

``figures.plot_solver_emulation`` (``--figures solver_emulation`` on the
command line) draws the three deposited curves with the emulated steps
marked on their nodes and, below each, the differences of the emulated
steps and of the accurate solution from the nodes on a logarithmic axis.

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_solver_emulation.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_solver_emulation.png
   :alt: Three panels, one per deposited curve (Fig. 9a, Fig. 9b and the
         H2O2 trace of SI Fig. S12a). In each, the emulated solver steps
         fall on the deposited nodes. The lower axes show the emulated
         steps within about 7e-6 of the nodes, and the accurate solution of
         the same model up to about 2e-3 from them in the two panels of
         plotted fractions and up to 1.2e-2 equivalents for the H2O2 trace.
   :class: pub-figure
   :width: 100%

   Drawn by ``plot_solver_emulation``: the deposited model nodes (circles),
   the emulated solver steps of the shipped models (plus signs) and the
   accurate solution (lines). Lower axes: largest difference from each node
   of the emulated step and of the accurate solution; the dotted line is
   the last printed decimal of the deposited files. The models carry the
   constant sets inferred from these same nodes
   (``H2O2_NODE_PARAMETER_SETS``, ``O2_NODE_PARAMETER_SETS``), not Table 3
   as printed: with Table 3 the emulated steps are up to 4.4e-4 (Fig. 9a)
   and 4.0e-4 (Fig. 9b) from the nodes, and the run of the H2O2 trace takes
   23 steps instead of 24.

Parity R2 values (SI Fig. S11)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The text states parity R2 values above 0.95 for both oxidants. SI Fig. S11b
prints 0.945 for H2O2, which disagrees with the text. It plots 11 of the 13
sampling times (all but t = 20 and 1205 min), and on those rows the H2O2
model reproduces the printed value (0.9452); on all 13 rows it gives 0.941.
SI Fig. S11a prints 0.997 for O2. The O2 model gives 0.9945 against the 16
deposited experimental rows and 0.9957 without the 120 min row. That row is
singled out because its 3-HP point does not appear in S11a, as far as the
image can be read, and because a refit of the six constants to the
experiment returns the Table 3 values within 0.05 % only when it is left
out (not shipped). The origin of 0.997 is not identified.

The undetected OPA pool
~~~~~~~~~~~~~~~~~~~~~~~

OPA was not detected in the experiments ("not detected with our analytical
methods", Section 3.3 of the paper) and is left out of the sum that
normalizes the plotted fractions, so those add up to 1 whatever the model
holds as OPA. The text gives the predicted OPA peak with O2 as 8.4e-4 mM
near 25 min and calls it negligible. The quantity plotted in SI Fig. S12b is
a mole fraction that peaks at 1.4e-5, and 1.4e-5 times 60 mM is 8.4e-4 mM,
which suggests that the printed value is that mole fraction multiplied by
the feed concentration. As a concentration, the OPA peak of the model is 2.5
mM, 4.2 % of the initial 3-HP, at 20.5 min. The H2O2 model holds more:
11.5 % of the initial 3-HP at 84 min (7.8 % is still there at 1440 min).
The simplified reversible model has no OPA.

The H2O2 curve of SI Fig. S12a is matched only with other constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The predicted H2O2 of SI Fig. S12a (the deposited trace and the published
inset, which agree with each other) is not the H2O2 of the Table 3 model,
which lies up to 0.38 equivalents above it. The trace is reproduced as the
H2O2 of another run of the same model: with the constant set
``'figS12a_h2o2_nodes'`` (from a fit to the trace that is not shipped) the
repeated solver run returns its 24 nodes within 4.4e-6 equivalents. That
set has ``k_7`` = 6.94e-3 1/min instead of 6.48e-3 and a ``k_5`` 22 % lower,
so the trace is consistent with a parameter set other than that of the
Fig. 9a curve and Table 3; how the authors computed it is not stated. The
species curves of the same figure are the Fig. 9a curve, which is
reproduced with constants within 0.6 % of Table 3 (the ``'fig9a_nodes'``
set; see `Repeating the authors' solver run`_). In the deposited file of
this figure the blocks labelled ``exp`` and ``model`` are also swapped; the
reader ``read_fig9a_h2o2`` corrects that.

Fig. 10d and the activation energies of the text
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The text of the paper says that its activation energies "were obtained from
Arrhenius plot in Fig. 10d" and gives 59.8, 91.6 and 59.5 kJ/mol; the
simplified reversible model has these, with the digits of the heat-map
script, for its steps 1, 2 and 5. The deposit holds the points and the
lines of that plot (``read_fig10d_arrhenius``, ``fig10d_points``): the
logarithms of three rate constants at 25, 50, 75 and 100 °C, which the
figure labels 3-HP → MA, MA → AAc and MA → Others. The paper does not say
of which model they are constants, and neither the paper nor the deposit
gives their unit; the package takes it to be 1/min, the time unit of the
authors' scripts. The slopes of the lines do not depend on it. The package
names the constants ``k_1``, ``k_2`` and ``k_5``, as the authors' heat-map
script numbers these reactions; these are not the numbers of the Table 3 O2
model.

The points are taken here as three rate constants of the simplified
reversible model (3-HP → MA, MA → AAc and MA → D), fitted to the
experiment of each temperature. Three things support this:

- the figure names a one-step 3-HP → MA reaction, which only that model has;
- fitted to the 50 °C data of the script with the equilibrium constants
  written in the script, that model returns the 50 °C row: the three
  constants come out 0.09, 1.5 and 0.37 % above it, inside their standard
  errors (see `The constants of the simplified reversible model`_);
- converged fits of that model to the experiments at 50, 75 and 100 °C come
  within 9 % of the points (the fit at 25 °C is ill-determined).

Taking the points as constants of that model means that they follow the
rate law of the heat-map script. The three reasons support that and do not
establish it: neither the paper nor the deposit says it. This reading is
what makes the 50 °C row the source of three constants of the shipped
model. In 1/min, with the O2 activity of the heat-map script
(3.5/1.01325 = 3.45) for 3-HP → MA, the three 50 °C values also lie within
the 95 % confidence intervals of the corresponding Table 3 constants
(``NOTES.md`` has the comparison; it is pinned by a test).

Each dashed line of the published figure is a regression through three of
the four points: it leaves out the 50 °C point for 3-HP → MA and the 25 °C
point for the other two. The lines give activation energies of 59.93, 68.96
and 52.09 kJ/mol (``fig10d_regression_lines``); the text gives 59.8, 91.6
and 59.5 kJ/mol. No subset of the plotted points of MA → AAc gives
91.6 kJ/mol, and none of MA → Others gives 59.5 kJ/mol. The text values are
those of the heat-map script; where they come from is not known.

The lines that the text values give through the 50 °C points
(``fig10d_text_lines``) are the constants ``k_1(T)``, ``k_2(T)`` and
``k_5(T)`` of the simplified reversible model, i.e. how the script
extrapolates from 50 °C. They miss the other points: the model's ``k_1`` is
4.3, 4.3 and 4.2 times the plotted value at 25, 75 and 100 °C (see `Limits
of the temperature dependence`_). The package draws no figure of the
Arrhenius plot.

Fig. 11 is computed with the simplified model of the deposited script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The text presents the heat maps as the O2 model extended with the Arrhenius
relation. The caption says that "Fig. 9b data was used as a reference kref
of the Arrhenius equation", which matches the deposited heat-map script:
the script fits a model of its own to the Fig. 9b data and extrapolates
that model. It is a different model from that of Table 3 (five reactions,
no OPA, a reverse term in every step) and is not written out in the paper.
(The paper numbers two equations (6); its "Arrhenius relationship (Eq. 6)"
is the second one, ``k(T) = k_ref*exp[-(Ea/R)*(1/T - 1/T_ref)]``, printed at
the end of its Section 3.4.) The published figure shows three properties of
that script, which the model keeps:

- **Arrhenius factors on the equilibrium constants.** The script scales each
  equilibrium constant with the Arrhenius factor of its rate constant, which
  the paper does not mention.
- **A start from unit mole fractions.** The heat maps start from a 3-HP mole
  fraction of 1, not from the dilute solution the same script fits at 50 °C.
  The 50 °C row of Fig. 11 is therefore not the 50 °C curve that the script
  fits: the MA fraction peaks at 0.571 in the first and at 0.549 in the
  second.
- **A 10 min time offset.** The time axis of Fig. 11 is the reaction time
  plus 10 min: the "27 min" and "49 min" of the text are 16.8 and 39.1 min of
  reaction.

Without the time offset, or with the dilute start in place of the unit
start, the same constants miss the published contours of panel b by far
more than the 0.30 px of the model as the script runs it: by 3.8 and 14 px
on the linearized measure of ``contour_offsets``, clipped at 30 px per
sample (6.6 and 91 px without the clip); offsets of this size are
indicative only (tested).

The script also uses a pressure of 3.5 bar (the captions give 3 bar of O2).
What Fig. 11 calls "MA selectivity" and "MA yield" are ``X_MA`` and
``(1 - X_HP3)*X_MA``.

The constants of the simplified reversible model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The model has five reactions, numbered and written as in the comments of
the authors' script: R1 3-HP + O2 → MA + H2O (rate constant ``k_1``,
equilibrium constant ``K_6``), R2 MA → AAc + CO2 (``k_2``, ``K_7``), R3
MA + O2 → OXA + H2 + CO2 (``k_3``, ``K_8``), R4 3-HP → D (``k_4``,
``K_9``) and R5 MA → D (``k_5``, ``K_10``); every step has a reverse
term. The constants that the script's run fitted are not in the
deposit. The script writes, to four decimals, the vector ``k0`` =
[0.0021 0.0001 0.0000 0.0023 0.0010 4.7840 0.5918 3.5116 81.6870 52.0472],
from which its search takes its bounds. The model's ten constants at
50 °C follow one rule: every constant the deposit gives is used as given,
the full-precision value where it gives one twice, and the one constant it
does not give is fitted to the 50 °C data with the authors' fitting set-up.

.. list-table::
   :header-rows: 1
   :widths: 22 24 54

   * - Constant
     - Value
     - Source
   * - ``k_1``, ``k_2``, ``k_5``
     - 2.139999e-3, 5.199984e-5, 1.062501e-3 1/min
     - **Deposit**: ``exp()`` of the 50 °C row of the deposited Arrhenius
       table of Fig. 10d (``11-4_Arrhenius.csv``). ``k0`` has them to four
       decimals as 0.0021, 0.0001 and 0.0010; the last disagrees with the
       table's 1.0625e-3, which would round to 0.0011.
   * - ``k_4``
     - 0.0023 1/min
     - **Deposit**: ``k0``, its only source (two significant digits)
   * - ``K_6`` to ``K_10``
     - 4.7840, 0.5918, 3.5116, 81.6870, 52.0472
     - **Deposit**: ``k0``
   * - ``k_3``
     - 1.4215e-5 1/min
     - **Fitted** to the 50 °C data table of the script, the other nine
       held; standard error 31 %. ``k0`` writes 0.0000 (below 5e-5).
   * - ``Ea_1`` to ``Ea_5``
     - 59781.2, 91594.6, 10000, 51788.73, 59549.8 J/mol
     - **Deposit**: as written in the heat-map script. ``Ea_i`` scales
       ``k_i`` and ``K_(i+5)`` with temperature (reference 50 °C).

In the model these 50 °C constants are the parameters ``k_1_ref`` to
``K_10_ref``; ``k_1`` to ``K_10`` are their values at ``T`` (assignment
rules, which ``set_value`` cannot change). The fit below accepts either
name and returns the ``_ref`` names.

``fit_simplified_reversible_constants`` makes the fit. It runs the model as
the script runs it when it fits (an initial 3-HP "mole fraction" of
0.000339) and minimizes the squared differences between the five liquid
fractions of the model and of the script's table at its 16 sample times, 80
values in all, unweighted. The table is shipped
(``read_heatmap_script_data``); it is the experimental block of the
deposited Fig. 9b file with its last two rows rounded to four decimals, as
the script has them. The model is always run at 50 °C, so the fitted values
are the constants at the temperature of the data.

.. code-block:: python

   from nskinetics.models.pd_c_3hp_oxidation import (
       fit_simplified_reversible_constants)

   fit = fit_simplified_reversible_constants()        # k_3
   for name, value in fit['fitted'].items():
       error = fit['standard_errors'][name]
       print(f'{name} = {value:.4e} 1/min (standard error {error:.1e})')
   print(f"sum of squares = {fit['sum_of_squares']:.5f} "
         f"over {fit['n_residuals']} residuals")

   # k_4 free as well: it comes out as k0 has it, 0.0023
   fit = fit_simplified_reversible_constants(free=('k_3', 'k_4'))
   print(f"k_4 = {fit['fitted']['k_4_ref']:.4e} 1/min")

   # all five rate constants free: k_1, k_2, k_5 against the Arrhenius table
   fit = fit_simplified_reversible_constants(
       free=('k_1', 'k_2', 'k_3', 'k_4', 'k_5'))
   print(', '.join(f'{name} = {fit["fitted"][name + "_ref"]:.3e}'
                   for name in ('k_1', 'k_2', 'k_4', 'k_5')))

.. code-block:: text

   k_3_ref = 1.4215e-05 1/min (standard error 4.4e-06)
   sum of squares = 0.03380 over 80 residuals
   k_4 = 2.3204e-03 1/min
   k_1 = 2.142e-03, k_2 = 5.277e-05, k_4 = 2.312e-03, k_5 = 1.066e-03

The shipped ``k_3`` is the first result, to five significant digits; a test
checks that, from four starting points. The data determine ``k_3`` poorly
(its product, oxalic acid, is at most 3.2 % of the liquid in the data).
Fitting ``k_4`` as well gives 2.3204e-3, which rounds to the 0.0023 of
``k0``, and lowers the sum of squares by only 0.07 % (``NOTES.md`` gives
this alternative).

This is not the authors' procedure, which cannot be repeated. Their script
fits all ten constants with a particle swarm of three iterations followed by
a least-squares refinement limited to five iterations, integrating with
``ode23s`` at default tolerances. The swarm searches 0 to 1.5 times ``k0``,
which for ``k_3`` is 0 to 0. The swarm is random and neither stage is run
to convergence, so another run gives other constants. The fit here is
carried to convergence with an accurate solver, and only the constant
without a deposited value is free. What supports holding the other nine
(``NOTES.md`` has the numbers; "not shipped" marks analyses that are not in
the repository):

- **The data alone do not give the equilibrium constants.** A converged fit
  of all ten constants determines the five rate constants and ``K_6``, but
  ``K_7`` to ``K_10`` end anywhere or run off to very large values. Its
  heat maps lie 2.7 to 2.8 px (panel a) and 1.8 to 1.9 px (panel b) from the
  published contours, further than those of the shipped constants.
- **With the equilibrium constants of k0, the data return the Arrhenius row
  and the k_4 of k0.** The third fit of the example frees all five rate
  constants: ``k_1``, ``k_2`` and ``k_5`` come out 0.09, 1.5 and 0.37 %
  above the 50 °C row of the deposited table (2.1400e-3, 5.200e-5,
  1.0625e-3), inside their standard errors (2.5, 31 and 3.9 %), and
  ``k_4`` 0.53 % above 0.0023 (standard error 5.4 %). This supports reading
  that row as constants of this model; it is an inference, stated in
  neither the paper nor the deposit.
- **Fits to the published contour lines agree with the data fit** (not
  shipped): ``k_3`` fitted to the contour lines of Fig. 11b instead of the
  data comes out 1.9 % above the shipped value, and 0.3 % below it when
  the lines of both panels are fitted together.
- **The authors' procedure scatters by more than these differences** (not
  shipped; an approximate emulation, since MATLAB's optimizers were not
  run): over 20 random seeds the rms offsets from the published contours
  ended between 1.7 and 3.9 px (panel a) and between 0.6 and 4.8 px
  (panel b), and the largest "yield" between 0.548 and 0.556.
- **The solver of the script does not matter for the fit** (not shipped):
  integrating as the script does changes the fitted ``k_3`` by 0.012 %.

**The two published panels cannot come from one set of constants.** Since
the "yield" is never larger than the "selectivity", the late side of a
yield contour must lie earlier than that of the selectivity contour of the
same level; in the published figure the late 0.1 and 0.3 contours of
panel b lie *later* than those of panel a above 62 °C, by 2.2 % in the
median and up to 4.4 % (tested on the digitized lines). No single run of
the script produces both panels; why they differ is not known. The shipped
set is 0.30 px from panel b and 1.05 px from panel a. Fits to the contour
lines (not shipped) show that no set does much better on both: with all ten
constants free and both panels fitted together, 0.60 and 0.48 px remain.
The package ships one set, defined without the figure, and shows the misfit
of panel a.

The two quoted maxima, 59.1 % and 55.2 %
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cell the text names for the maximum selectivity (95 °C, 27 min) is the
largest of the model's map, with 0.593 where the text has 59.1 %. The cell
it names for the maximum yield (84 °C, 49 min) is the largest of the
model's map as well, with 0.548 where the text has 55.2 %. The crest of the
yield map is flat (13 cells between 76 and 90 °C lie within 0.001 of the
maximum), so not much should be read into which cell is the largest. The
colour of the published yield map along its crest reads 0.549 ± 0.006 (not
shipped). The converged fit of all ten constants has 0.554 in that cell,
with contour lines 1.9 px from the published ones.

The text sets two "observed" maxima against these: 56.9 % for the
selectivity ("the observed maximum of 56.9 % at 50 °C (purple star on
plot)") and 50.5 % for the yield. The abstract places both at 3 bar of O2
and 50 °C, the condition of the Fig. 9b experiment, and the paper does not
say how either was computed. The 50.5 % equals the quantity of panel b,
``(1 - X_HP3)*X_MA``, of that experiment at 360 min (0.505), which is an
inference from the match. The 56.9 % is not identified: the quantity of
panel a, ``X_MA``, is at most 0.548 in that experiment (240 min), and MA
per 3-HP converted, ``X_MA/(1 - X_HP3)``, is 0.621, 0.577 and 0.513 at 240,
360 and 450 min and 0.569 at no sampling time. The 50 °C row of the model's
selectivity map peaks at 0.571 at 261 min, near the purple star, but that
is a model value. The match of the 50.5 % therefore stands alone and does
not settle how the paper defines the two quantities.

The initial-rate ratio of 0.204 is not reproduced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first of the paper's two equations numbered 6 is the ratio of the
initial rates of reaction 1, ``r0(H2O2)/r0(O2) = k(H2O2)*x_H2O2/(k(O2)*P_O2)``,
for which the text gives 0.204 under 3 bar O2 and 15 equivalents of H2O2 (O2
about five times faster). In the shipped models the two oxidant factors are
equal at t = 0 (0.005085), so the ratio of the initial rates is that of the
Table 3 constants, 5.31/7.14 = 0.744. With those constants the paper's value
corresponds to ``x_H2O2/P_O2`` = 0.274; the quantities stated in the paper
give other rate ratios (the H2O2 concentration of 0.9 M over 3 bar, 0.3 in
place of ``x_H2O2/P_O2``, gives 0.223; the mole fraction of the models,
0.005085, over 3 bar gives 0.0013). How 0.204 was obtained is not known.

The paper goes on to the hypothetical condition of equal oxidant mole
fractions, ``x_H2O2 = x_O2``, which it puts at about 1000 bar of O2 by
Henry's law, and states that O2 would then be about 1300 times faster than
H2O2. This is not reproduced either. At t = 0 the two shipped models have
the same oxidant factor (``x_O2_ref`` equals the initial H2O2 mole fraction,
0.005085), which is that condition if ``x_O2_ref`` is read as a mole
fraction, and their ratio is 7.14/5.31 = 1.34 in favour of O2. In the
paper's equation the O2 constant multiplies a pressure; in the O2 model the
Table 3 constant multiplies ``x_O2_ref``, which has the size of a mole
fraction, because that is what reproduces the deposited curve (see
`x_O2_ref is inferred`_).

What is not in the repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The package and its tests establish that the Table 3 models, with the
constants as documented, agree with the deposited curves and the
experiments to the numbers given on this page, that with the authors'
solver run repeated they return the deposited model nodes to their printed
digits, that the constants of the simplified reversible model are the
values of the deposit and the result of the shipped fit, how far its heat
maps are from the published contour lines, and that each model implements
its documented rate laws. The following were obtained with scripts that are
not shipped, so they cannot be recomputed from the repository:

- the least-squares fits made through the repeated solver run: the constant
  sets of the Fig. 9 nodes, and the refit with the O2 fraction held
  constant. The tests show that the shipped constants return the nodes, not
  that no other constants or rate laws would. MATLAB itself was never run;
- the fits to the digitized contour lines of Fig. 11 quoted above and in
  ``NOTES.md``. No shipped constant comes from them;
- the emulation of the authors' particle swarm and truncated least squares,
  and the fit made with the ``ode23s`` port in place of the accurate
  solver;
- the refit of the O2 constants to the experiment without the 120 min row;
- the independent integrations behind the comparison with the Table S4
  orders, the repeat of the heat-map script's solver runs for the rows of
  Fig. 11 other than 100 °C, and the reading of the colours of the
  published maps (``NOTES.md``);
- the digitization of the published figures (the digitized curves
  themselves are shipped).

The fits of the simplified reversible model to the experiments (all ten
constants free, at 50 °C and at the other temperatures) were run with the
shipped ``fit_simplified_reversible_constants`` and are not pinned by a
test; the fits of ``k_3``, of ``k_3`` and ``k_4``, and of the five rate
constants are.
