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

- four models, each an Antimony file with its generated SBML file and a
  loader that returns a :class:`~nskinetics.KineticModel` (time in min,
  concentrations in M, temperatures in K);
- the deposited tables of Figs. 9 and 10 (verbatim copies), curves digitized
  from SI Fig. S12 and Fig. 11, and readers for them;
- functions that redraw Fig. 9, SI Fig. S12, Fig. 10 and Fig. 11 from the
  models, three figures that compare the models with the published ones,
  and one that verifies them against the deposited model nodes;
- ``NOTES.md``, the full record of what is printed in the paper, what was
  inferred from the deposit and where the two disagree. This page
  summarizes it.

The package needs no optional dependencies, and importing it loads no
model. Its API is documented on :doc:`the shipped-models page </API/models>`.

The paper describes one model per oxidant. Redrawing the model curves of its
Figs. 9, 10 and 11 took four models, two of which are not in the paper, and
one factor that is stated in neither the paper nor the deposit. Two things
follow, and `What does not match the paper, and what is inferred`_ gives the
details; read it before using a model for anything other than redrawing the
figures:

- Only Fig. 9 and SI Fig. S12b are reproduced with constants printed in the
  paper (for O2, together with the inferred factor ``x_O2_ref``). Fig. 10a–c
  are reproduced with constants that were all fitted to the very curves
  being reproduced. Fig. 11 is reproduced with ten constants at 50 °C, of
  which two (panel b) or four (panel a) were fitted to the contour lines
  being reproduced; the others, and the five activation energies, are taken
  from the deposited Arrhenius table and the heat-map script. That agreement
  shows consistency, not that those models are right. The points of Fig. 10d
  are deposited data, not model output.
- Several statements on this page rest on analyses that are not in the
  repository (fits to the deposited curves and to digitized lines). They
  are marked "(not shipped)"; `What is not in the repository`_ lists them.
  The repeat of the authors' solver run is in the repository: see
  `Repeating the authors' solver run`_.

.. warning::

   ``load_o2_model(T)`` extends the published 50 °C O2 model to other
   temperatures with the Arrhenius relation of the paper. Only the 50 °C
   model is a published one. Away from 50 °C it reproduces neither the
   published heat maps (Fig. 11) nor the experiments at 25, 75 and 100 °C
   (Fig. 10): see `Limits of the temperature-extended O2 model`_.

   The paper numbers two equations (6). "The Arrhenius relation of the
   paper" means the second one, ``k(T) = k_ref*exp[-(Ea/R)*(1/T - 1/T_ref)]``,
   printed at the end of its Section 3.4. The first Eq. 6 is a ratio of
   initial rates (see `The initial-rate ratio of 0.204 is not reproduced`_),
   and Fig. 6 of the paper is unrelated (XPS and reusability data).

The four models
---------------

The paper presents one kinetic model per oxidant (Table 3). Each of the four
models here is shipped because it is needed to reproduce a published item.
Each is an Antimony file with the SBML file generated from it next to it
(same name, ending in ``_sbml.xml``).

.. list-table::
   :header-rows: 1
   :widths: 15 37 18 30

   * - Model
     - Loader and Antimony file
     - Reproduces
     - Where it comes from
   * - H2O2 model
     - | ``load_h2o2_model()``
       | ``pd_c_3hp_oxidation_h2o2_antimony.txt``
     - Fig. 9a (50 °C). Its predicted H2O2 is compared with SI Fig. S12a,
       whose H2O2 curve lies up to 0.38 equivalents lower
     - **Paper**: Table 3, H2O2 column, as printed, in the rate form of SI
       Fig. S2. How H2O2 enters the rates was **inferred** from the
       deposited model curve (the authors posted only the O2 version of
       their script).
   * - O2 model
     - | ``load_o2_model(T)``
       | ``pd_c_3hp_oxidation_o2_antimony.txt``
     - Fig. 9b and SI Fig. S12b (50 °C)
     - **Paper**: Table 3, O2 column, as printed, in the rate form of SI
       Fig. S2 with one **inferred** factor (``x_O2_ref``) in place of the
       O2 activity. Its temperature dependence is an NSKinetics addition
       that reproduces no published item.
   * - Reversible O2 variant
     - | ``load_o2_reversible_model(parameter_set)``
       | ``pd_c_3hp_oxidation_o2_reversible_antimony.txt``
     - The model curves of Fig. 10a–c (25, 75 and 100 °C), one constant set
       per curve: ``'25C'``, ``'75C'``, ``'100C'``
     - **Inferred** from the deposited curves: structure and all constants.
       It is in neither the paper nor the authors' scripts.
   * - Simplified reversible O2 model
     - | ``load_o2_simplified_reversible_model()``
       | ``pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt``
     - Fig. 11; the 50 °C point of Fig. 10d supplies three of its rate
       constants
     - **Deposit**: the authors' heat-map script. It is not in the paper.
       The constants that are not in the deposit were **inferred** from the
       published contour lines.

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
is in 1/min. The constants of the reversible variant and the rate constants
of the simplified model are lumped for the paper's batch and are in 1/min.

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

The O2 model takes a temperature **in K** (``load_o2_model(T=...)``, or
``km.set_value('T', ...)`` on a loaded model), which scales its six rate
constants by Arrhenius factors. The reversible O2 variant has no
temperature: ``load_o2_reversible_model`` takes the *name* of one of its
three constant sets, ``'25C'``, ``'75C'`` or ``'100C'``, each of which
belongs to one experiment. ``read_fig10_o2`` and the other functions that
take the temperature of one of these experiments (``emulate_fig10abc``,
``figures.draw_fig10abc``) accept the number in °C (25, 75 or 100) or the
same label. The example compares the two models with the 75 °C experiment.
Read `Limits of the temperature-extended O2 model`_ before relying on the
first:

.. code-block:: python

   from nskinetics.models.pd_c_3hp_oxidation import (
       load_o2_reversible_model, read_fig10_o2)

   km_75 = load_o2_model(T=75 + 273.15)           # kelvin
   factor = km_75.get_value('k_1')/km_75.get_value('k_1_ref')
   print(f'k_1(75 C)/k_1(50 C) = {factor:.2f}')

   experiment = read_fig10_o2(75)['experiment']   # deposited points, Fig. 10b
   for name, model in (('Table 3 O2 model at 75 C', km_75),
                       ("constant set '75C'", load_o2_reversible_model('75C'))):
       sim = simulate_at_times(model, experiment['time'])
       print(f'{name}: pooled R2 = {pooled_r2(experiment, sim):.3f}')

.. code-block:: text

   k_1(75 C)/k_1(50 C) = 4.94
   Table 3 O2 model at 75 C: pooled R2 = 0.644
   constant set '75C': pooled R2 = 0.990

A temperature given to a loader in the wrong unit is refused:
``load_o2_model(T=75)`` raises a ``ValueError`` (``T`` must lie within 250
to 450 K), and ``simulate_heat_maps`` refuses temperatures outside −50 to
200 °C. These ranges only catch that mistake; they are not ranges of
validity. ``set_value('T', ...)`` is not checked.

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
else in a model (the constants, ``x_per_M``, ``x_O2_ref``) reproduces the
paper and should be left alone unless you know why you are changing it.

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
   * - all four
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
   * - O2 model, simplified reversible model
     - ``T``
     - K
     - 323.15
     - Temperature (50 °C by default).
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
  models compute, not what was measured. For ``T`` with O2 the comparison
  with the time courses of Fig. 10 is poor (`Limits of the
  temperature-extended O2 model`_). For ``m_cat`` and ``C0`` no experiment
  was found in the paper or the deposit to compare with.
- **The paper's screening experiments contradict some of these
  extrapolations.** The time courses behind the models are those of Figs. 9
  and 10; all O2 ones are at 3 bar of O2. The paper also has experiments
  analysed once, after 8 h, in which the H2O2 charge, the temperature or
  the O2 pressure was varied (its Figs. 2 and 3; SI Table S2 has the
  temperature series with O2, and SI Table S1 other catalysts), and
  experiments with malonic acid alone (its Fig. 4). The deposit holds the
  values of Figs. 2 to 4. The models were not built on these experiments
  and do not describe them (``NOTES.md`` has the table; the model values
  are pinned by the package tests):

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
    figure and 30 by the deposit). In the H2O2 model both rise with every
    equivalent added, as with 30 equivalents in the example above.
  - *Temperature with O2* (Fig. 3c). After 8 h at 25 °C the O2 model has
    converted 53 % where 8.5 % was measured; at 75 °C it has an MA yield of
    4.7 % where 35.7 % was measured (the simplified reversible model: 53.5
    and 6.0 %).
  - *Malonic acid alone* (Fig. 4). Without 3-HP, 3.8 % (Fig. 4a) or 1.2 %
    (the time course of Fig. 4b) of the MA was converted after 1300 min at
    50 °C; the O2 model started from 60 mM MA converts 80 %. At 100 °C,
    30.1 % was converted after 120 min, where the model converts 96 %; after
    1300 min the measured 89.8 % (Fig. 4a) or 99.0 % (Fig. 4b) is close to
    the model's 100 %.

  The screening experiments are not on the same footing as the time
  courses either. At the same nominal conditions they give 72.8 % (O2) and
  78.4 % (H2O2) conversion after 480 min, where the time courses of Fig. 9
  are at 96.1 % after 450 min and 60.6 % after 420 min. No explanation was
  found in the paper; it states that for time-resolved runs substrate,
  NaOH and water were tripled and that 1.5 mL was withdrawn at every
  sample.
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
- **The two reproduction models have lumped constants.** The reversible O2
  variant and the simplified reversible model have no catalyst-mass
  parameter; their constants are effective ones for the paper's batch. The
  reversible variant has no temperature or pressure either, and each of its
  constant sets belongs to one experiment (samples up to 1890, 1440 and
  1560 min at 25, 75 and 100 °C).
- **Inputs are not validated.** ``set_value`` is the generic NSKinetics
  call and checks nothing: zero or negative ``m_cat``, ``C0`` or
  ``eq_H2O2_0`` run without an error and give meaningless results.

Heat maps and figures
---------------------

The heat maps of Fig. 11 are computed by ``fig11_heat_maps``, which runs the
simplified reversible model as the authors' script does, once per
temperature, for each of the two panels with that panel's constant set
(2 × 71 integrations). Its time axis is the *plotted* time of the published
figure, which is the reaction time plus 10 min (``FIG11_TIME_OFFSET``).
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

   0.5906 at 95 C, plotted time 26.8 min (reaction time 16.8 min)
   0.5483 at 84 C, plotted time 49.1 min (reaction time 39.1 min)
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
       --figures fig10 fig10_table3_o2 --formats png --dpi 150

In Python, ``figures.plot_fig9()`` and its siblings return the matplotlib
figure, and ``figures.draw_fig9a(ax)`` and its siblings draw one panel on an
axes of your own.

Reproduced figures
------------------

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
of the deposited file, with constants 0.01 to 0.6 % off Table 3 (see
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

Fig. 10: O2 at 25, 75 and 100 °C, and the Arrhenius plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig10.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig10.png
   :alt: Four panels. Panels a to c: normalized concentration against time
         with O2 at 25, 75 and 100 degrees C. At 25 degrees 3-HP falls only
         to about 0.9; at 75 degrees malonic acid peaks near 0.4 at about 130
         min; at 100 degrees it peaks near 0.55 at about 40 min, and the
         other by-products end near 0.9. The solid model lines lie on the
         black dotted deposited curves. Panel d: ln k against 1000/T for
         three rate constants at four temperatures, with three dashed
         regression lines; each line passes close to three of its four
         points, and the fourth lies 1.4 to 4.6 above it.
   :class: pub-figure
   :width: 100%

   Fig. 10, drawn by ``plot_fig10``. (a–c) Experiment (circles), the
   reversible O2 variant with the constant set of each curve (solid lines)
   and the deposited model curves (black dotted lines); the box gives the
   pooled R2 against the experimental points to three decimals. (d) The
   deposited points of the Arrhenius plot and the three dashed lines,
   recomputed by least squares.

The reversible O2 variant differs from the deposited curves by at most
1.8e-4 (25 °C), 2.5e-3 (75 °C) and 2.7e-3 (100 °C) at their nodes, again the
solver error of those curves: with the authors' solver run repeated it
returns the nodes within 7.2e-6 (see `Repeating the authors' solver run`_).
Its pooled R2 against the experimental points is 0.9996, 0.990 and 0.993;
the published panels print 0.99. This variant was inferred from the very
curves it is compared with, so the agreement shows that the inference is
consistent, not that the model is right (see `Fig. 10a–c are reproduced by
an inferred reversible variant`_).

The points of panel d cannot be recomputed and are read from the deposit.
The lines are recomputed and agree with the deposited ones to 7e-6 in
``ln k`` at their 100 °C ends.

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
         cell, 0.591 at 95 degrees and 27 min in panel a and 0.548 at 84
         degrees and 49 min in panel b; a purple star marks the maximum of
         the 50 degree row. An inset in each panel zooms in on 47 to 57
         degrees.
   :class: pub-figure
   :width: 100%

   Fig. 11, drawn by ``plot_fig11``: (a) MA "selectivity" (``X_MA``) and
   (b) MA "yield" (``(1 - X_HP3)*X_MA``), computed with the simplified
   reversible model as the authors' heat-map script runs it, each panel with
   its own constant set. Two things differ from the published figure on
   purpose: the largest cell of each map is annotated with its value, and
   the colour-bar label of panel a is spelled "MA Selectivity" (published:
   "MA Selectvity"). The stars sit on the computed maxima, and the insets
   are zooms of the same matrices (``NOTES.md`` gives the differences from
   the published stars and insets).

The deposit holds no table of the heat maps, so the reference is the
published figure: its contour lines at 0.1, 0.3 and 0.5, digitized. The
reproduced contours lie within an rms of 0.17 px (panel a) and 0.27 px
(panel b) of them, in pixels of the published raster (1 px = 3.07 min or
0.148 °C); the digitization itself is good to about 0.3 px. The red dashed
line, which joins the row maxima and is therefore defined to one column of
the time grid (1.82 px), lies 0.6 px (panel a) and 0.15 px (panel b) from
the digitized one in the mean along the time axis. The largest cell of the
selectivity map, 0.5906 at 95 °C and 26.8 min, is the 59.1 % (95 °C, 27 min)
of the text. The largest cell of the yield map is the one the text names
(84 °C, 49 min), but its value is 0.548, not the 55.2 % of the text.

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11_verification.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11_verification.png
   :alt: Top row: the two reproduced heat maps with the digitized published
         contour lines, ridge line and stars drawn on top as fine black
         dots; the dots lie on the wide reproduced lines everywhere.
         Bottom row: the offset of every digitized contour sample from the
         reproduced contour against the temperature of the sample. All
         offsets lie within about 1 px, with an rms of 0.17 px in panel a
         and 0.27 px in panel b.
   :class: pub-figure
   :width: 100%

   Drawn by ``plot_fig11_verification``. Top: the reproduced maps with the
   lines and stars digitized from the published panels in black. Bottom: the
   offset of each digitized contour sample from the reproduced contour.

Limits of the temperature-extended O2 model
-------------------------------------------

.. warning::

   Use ``load_o2_model`` at temperatures other than 50 °C only as what it
   is: an extrapolation that reproduces no published result and misses the
   experiments at 25, 75 and 100 °C.

**What it is.** In ``load_o2_model(T)`` each of the six rate constants of
the Table 3 O2 model follows the Arrhenius relation of the paper,
``k(T) = k_ref*exp(-(Ea/R)*(1/T - 1/T_ref))`` (the second of its two
equations numbered 6), with the reference at 50 °C. At 50 °C it is the
published model by construction. Three activation energies are those quoted
in the text (59.8, 91.6 and 59.5 kJ/mol; the first is applied to both R1
and R2); the other two (51.79 kJ/mol for 3-HP → D and 10 kJ/mol for
MA → OXA, a round value) appear only in the authors' heat-map script.
Applying them to this six-reaction model is an NSKinetics addition: the
text describes Fig. 11 this way, but none of the paper's temperature
results is reproduced with it.

**It does not reproduce the published heat maps.** Computed with this model
on the grid of Fig. 11 (``table3_o2_heat_maps``), the maps have the shape of
the published ones and maxima 0.003 to 0.006 lower, but their contour lines
are elsewhere:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * -
     - Table 3 O2 model at T
     - Simplified reversible model
   * - rms offset from the published contours, panel a / b
     - 7.6 / 7.6 px
     - 0.17 / 0.27 px
   * - Largest cell, "selectivity"
     - 0.587 at 91 °C, 21 min
     - 0.591 at 95 °C, 27 min
   * - Largest cell, "yield"
     - 0.542 at 84 °C, 38 min
     - 0.548 at 84 °C, 49 min

Offsets above about 1 px are not distances. They are the linearized
estimate of ``contour_offsets`` (value difference over gradient), clipped at
30 px per sample (``max_offset=30``), and are indicative only; without the
clip, panel b of the Table 3 model reads 8.8 px. Measured on the lines
themselves (``contour_shifts``), the late side of every contour is earlier
than published by 8 to 11 % of the time in the median per line (28 to
74 min), which is the same as 1.3 to 1.8 °C lower at a given time. Shifting
its time axis by the 10 min of the published figure leaves 6.8 and 6.1 px.
The difference is 20 times the reading uncertainty of the figure or more.

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11_table3_o2.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11_table3_o2.png
   :alt: The two heat maps computed with the Table 3 O2 model, with the
         digitized published contour lines as fine black dots. On the late
         side of the ridge the wide contour lines of the Table 3 model lie
         visibly to the left of and below the published ones; the largest
         cells are 0.587 at 91 degrees C and 21 min and 0.542 at 84 degrees C
         and 38 min.
   :class: pub-figure
   :width: 100%

   Drawn by ``plot_fig11_table3_o2``: the Fig. 11 quantities computed with
   the temperature-extended Table 3 O2 model (colours and wide lines), with
   the published contour lines, ridge line and stars (digitized, black).

**It fits the experiments at 25, 75 and 100 °C poorly.** Against the
deposited experimental points of Fig. 10a–c, next to the per-temperature
constant sets that reproduce the published curves
(``load_o2_reversible_model``):

.. list-table::
   :header-rows: 1
   :widths: 46 18 18 18

   * -
     - 25 °C
     - 75 °C
     - 100 °C
   * - Table 3 O2 model at T: pooled R2 / rms
     - 0.579 / 0.243
     - 0.644 / 0.174
     - 0.746 / 0.150
   * - Per-temperature set: pooled R2 / rms
     - 0.9996 / 0.008
     - 0.990 / 0.029
     - 0.993 / 0.026
   * - Time to half conversion of 3-HP, Table 3 at T / per-temperature set
     - 458 min / not reached
     - 16.0 / 75.0 min
     - 4.2 / 22.0 min
   * - Maximum MA fraction, Table 3 at T
     - 0.522 at 1479 min
     - 0.584 at 51.7 min
     - 0.587 at 12.8 min
   * - Maximum MA fraction, experiment
     - 0.023 at 1890 min
     - 0.436 at 180 min
     - 0.568 at 40 min

The model is too fast at all three temperatures. At 25 °C it predicts 0.051
of the 3-HP left at 1890 min, where 0.876 was measured. At 100 °C the height
of the MA maximum is about right, but it comes at 13 min instead of 40 min.

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig10_table3_o2.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig10_table3_o2.png
   :alt: Four panels. Panels a to c: the experimental points at 25, 75 and
         100 degrees C with the Table 3 O2 model as solid lines and the
         per-temperature constant sets as dashed lines. The solid lines run
         far from the points: at 25 degrees the model consumes most of the
         3-HP while the measured 3-HP stays near 0.9; at 75 and 100 degrees
         the model's malonic acid peak comes much earlier than the measured
         one. The boxes give pooled R2 values of 0.579, 0.644 and 0.746 for
         the Table 3 model and 1.000, 0.990 and 0.993 for the reversible
         variant. Panel d: the Arrhenius points with the dashed regression
         lines and solid lines with the activation energies of the text
         through the 50 degree points; the solid lines miss most of the
         other points.
   :class: pub-figure
   :width: 100%

   Drawn by ``plot_fig10_table3_o2``. (a–c) The experiments of Fig. 10 with
   the temperature-extended Table 3 O2 model (solid lines) and the
   reversible variant with its per-temperature constant sets (dashed
   lines); the boxes give the pooled R2 of both, to three decimals. (d) The
   points of Fig. 10d with the published regression lines (dashed) and the
   lines that the activation energies of the text give through the 50 °C
   points (solid).

**Why.** The 50 °C constant of 3-HP conversion lies well off the Arrhenius
regression through the other temperatures. The measured 3-HP conversion at
75 °C is no faster than at 50 °C (half converted between 60 and 90 min at
both), whereas the Arrhenius factor of the model from 50 to 75 °C is 4.9 for
R1 and R2. The same shows in Fig. 10d, where the 50 °C point of the 3-HP →
MA constant lies a factor of 4.3 above the regression line through the other
three points. An extrapolation anchored at 50 °C with the activation
energies of the text therefore overshoots at 25, 75 and 100 °C. This
accounts for part of the misfit only: with the constants of R1 and R2
lowered by that factor, the pooled R2 rises from 0.58 and 0.64 to 0.82 and
0.85 at 25 and 75 °C and stays at 0.75 at 100 °C. Most of the rest is the
3-HP → D step, which is extrapolated from 50 °C with an activation energy
that appears only in the script; with its constant lowered tenfold as well,
the pooled R2 is 0.96, 0.98 and 0.93 (both are illustrations, not fits, and
are pinned by the package tests). Part of the difference from the published
heat maps is also structural: the Table 3 model has the OPA intermediate and
no reverse terms, and its O2 factor stays above 0.999 where that of the
simplified model falls to 0.83.

**What to use instead.** For the curves of Fig. 10a–c use
``load_o2_reversible_model``; its three constant sets are separate fits
with no temperature law, so it exists only for the experiments at 25, 75
and 100 °C. For Fig. 11 use ``load_o2_simplified_reversible_model``. The
authors' heat-map script extrapolates that model from 50 °C in the same way
(its 3-HP → MA constant at 75 and 100 °C is 4.3 and 4.2 times the value of
the Fig. 10d points at those temperatures); how far the maps of Fig. 11 are
from the time courses at 75 and 100 °C has not been measured here (for the
8 h screening experiments see `Changing the operating conditions`_).

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
   * - Reversible O2 variant
     - First-order terms.
     - Reactions 1, 2, 5 and 6 are reversible; the constants are lumped, with
       no activity of O2 and no catalyst mass.
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

The deposited model curves of Figs. 9 and 10 have 18 to 24 unevenly spaced
solver steps each (19 to 25 nodes with the one at t = 0), and the published
lines are straight segments between the nodes. A Python port of MATLAB
``ode23s`` at its default tolerances returns those nodes to their printed
digits (next section), so they are consistent with being that solver's raw
steps; MATLAB itself was not run. They are not converged solutions: the
nodes are up to 2.0e-3 (Fig. 9b) and 2.7e-3 (Fig. 10b and c) away from the
accurately solved model, and the segments between them up to 0.010. The
models here are solved to a relative tolerance of 1e-9, so they do not lie
exactly on the published polylines, and the package tests that solve the
models accurately allow for that solver error.

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
       emulate_fig9b, emulate_fig10abc, node_deviations, read_fig9b_o2,
       read_fig10_o2)

   def report(label, emulated, nodes):
       d = node_deviations(emulated, nodes)
       print(f"{label}: {d['steps']} steps (deposited {d['deposited_steps']}),"
             f" node times within {d['time']:.0e} min,"
             f" values within {d['max']:.1e}")

   nodes = read_fig9b_o2()['model']
   report('Fig. 9b, constants of the nodes', emulate_fig9b(), nodes)
   report('Fig. 9b, Table 3 as printed', emulate_fig9b('table3'), nodes)
   for T_C in (25, 75, 100):
       report(f'Fig. 10, {T_C} C', emulate_fig10abc(T_C),
              read_fig10_o2(T_C)['model'])

.. code-block:: text

   Fig. 9b, constants of the nodes: 19 steps (deposited 19), node times within 6e-06 min, values within 5.0e-06
   Fig. 9b, Table 3 as printed: 19 steps (deposited 19), node times within 4e-02 min, values within 4.0e-04
   Fig. 10, 25 C: 18 steps (deposited 18), node times within 5e-03 min, values within 7.2e-06
   Fig. 10, 75 C: 19 steps (deposited 19), node times within 5e-03 min, values within 6.3e-06
   Fig. 10, 100 C: 24 steps (deposited 24), node times within 5e-03 min, values within 5.7e-06

The deposited files print values to five decimals, and node times to five
decimals (Fig. 9) or two (Fig. 10), so these differences are the print
rounding of the files. The same holds for Fig. 9a (23 steps, values within
6.8e-6) and for the H2O2 trace of SI Fig. S12a (24 steps, of which the file
holds all but the last; within 4.4e-6 equivalents). The package tests
require the deposited number of steps, node times within 2e-5 min (Fig. 9)
or 0.01 min (Fig. 10) and values within 1e-5; the solver error that the
accurate-solver tests have to allow for is 18 to 270 times larger.

The constants matter at this level. For Fig. 10 they are the six-digit sets
of the reversible O2 variant. For Fig. 9 they are sets inferred from the
nodes (``O2_NODE_PARAMETER_SETS``, ``H2O2_NODE_PARAMETER_SETS``; selected
with ``parameter_set`` in the loaders), given to eight digits. They are not
the published constants and are for this comparison only: the O2 set
differs from Table 3 by 0.01 to 0.26 % and prints as Table 3 prints it; the
H2O2 set differs by 0.1 to 0.6 %, more than rounding allows for five of its
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
     - 0.038 (Fig. 9b), 0.094 (Fig. 9a), 0.029 to 0.157 (Fig. 10a–c)
   * - Reversible variant without the reverse of reaction 1
     - 0.94, 7.6e-5 and 0.011 at 25, 75 and 100 °C
   * - Reversible variant with a linear reverse of reaction 1
       (``k_1r*[OPA]``)
     - 0.10, 8.3e-5 and 0.0095
   * - Reversible variant with water formed by reaction 1 only
     - 9.1e-4 at 25 °C; not detectable at 75 °C; at 100 °C only through
       the node times, which move by 0.012 min (limit 0.01 min)

What this shows is limited. The shipped constants return the deposited
nodes, and most of them are pinned on their own (to between 1e-6 and 1e-3
of their value for the Fig. 9 sets, and to between 1e-4 and 0.1 for most
constants of the reversible variant; ``NOTES.md`` has the tables). The
exceptions are ``k_6`` of the H2O2 model, which is not identifiable, and
``k_5r`` at 25 °C and ``k_5`` at 100 °C of the reversible variant: a change
of 10 % in any of these leaves the nodes as they are. A control shows that
the shipped constants need the shipped choice. It does not show that
another choice could not do as well with other constants: that takes a
refit of the alternative, and the refits that were made are not shipped
(`What is not in the repository`_). For the O2 fraction, for example, a
refit with the fraction held at 1 still left a difference of 1.5e-5, three
times that of the shipped model.

``figures.plot_solver_emulation`` (``--figures solver_emulation`` on the
command line) draws the six deposited curves with the emulated steps marked
on their nodes and, below each, the differences of the emulated steps and
of the accurate solution from the nodes on a logarithmic axis.

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_solver_emulation.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_solver_emulation.png
   :alt: Six panels, one per deposited curve (Fig. 9a, Fig. 9b, the H2O2
         trace of SI Fig. S12a, and Fig. 10a to c). In each, the emulated
         solver steps fall on the deposited nodes. The lower axes show the
         emulated steps within about 7e-6 of the nodes, and the accurate
         solution of the same model up to about 3e-3 from them in the five
         panels of plotted fractions and up to 1.2e-2 equivalents for the
         H2O2 trace.
   :class: pub-figure
   :width: 100%

   Drawn by ``plot_solver_emulation``: the deposited model nodes (circles),
   the emulated solver steps of the shipped models (plus signs) and the
   accurate solution (lines). Lower axes: largest difference from each node
   of the emulated step and of the accurate solution; the dotted line is
   the last printed decimal of the deposited files. In the Fig. 9a, Fig. 9b
   and S12a panels the models carry the constant sets inferred from these
   same nodes (``H2O2_NODE_PARAMETER_SETS``, ``O2_NODE_PARAMETER_SETS``),
   not Table 3 as printed: with Table 3 the emulated steps are up to 4.4e-4
   (Fig. 9a) and 4.0e-4 (Fig. 9b) from the nodes, and the run of the H2O2
   trace takes 23 steps instead of 24. In the Fig. 10 panels they carry the
   constant sets of the reversible variant, which were inferred from the
   nodes too.

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
mM, 4.2 % of the initial 3-HP, at 20.5 min. The other models hold more: 11.5
% of the initial 3-HP at 84 min in the H2O2 model (7.8 % is still there at
1440 min), and peaks of 2.9, 12.6 and 20.1 % in the reversible O2 variant
with the sets of 25, 75 and 100 °C.

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

Fig. 10a–c are reproduced by an inferred reversible variant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The paper prints no constants for 25, 75 or 100 °C, and its text treats all
reactions as irreversible. Neither posted script produces the model curves
deposited for these temperatures; the README of the deposit calls the
posted codes "an example" whose details "were adjusted for each time". The
curves are reproduced by the reaction network of the Table 3 O2 model with
reactions 1, 2, 5 and 6 made reversible and a separate constant set for
each curve. Its structure and all of its constants were inferred from the
deposited curves, with fits made through a repeat of the authors' solver run
(the fits are not shipped; the repeated run is, and with these constants it
returns the deposited nodes to their printed digits, see `Repeating the
authors' solver run`_). Nothing of it is printed in the paper, the SI or the
deposit, and it is the only structure found that reproduces the curves,
which does not prove that the authors used it. The three constant sets
follow no temperature law (the forward constant of 3-HP → OPA is 0.047,
0.0092 and 0.041 1/min in the sets of 25, 75 and 100 °C), so they are fit
parameters of each curve, not kinetic quantities, and the loader takes the
name of a set, not a temperature. The model is a reproduction aid.

Fig. 10d and the activation energies of the text
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The paper does not say of which model the points of the Arrhenius plot are
constants. They are taken here as three rate constants of the simplified
reversible model (3-HP → MA, MA → AAc and MA → D), fitted to the experiment
of each temperature. Three things support this: the figure names a one-step
3-HP → MA reaction, which only that model has; the 50 °C row supplies the
three rate constants with which that model reproduces Fig. 11b (fitted
freely to the contour lines of that panel, one of them came out within 0.1 %
of the tabulated value; not shipped); and refits of that model to the
experiments at 50, 75 and 100 °C come within 9 % of the points (not shipped;
the refit at 25 °C is ill-determined). The points agree with neither the
Table 3 model nor the model that reproduces the curves of panels a to c.
Each dashed line is a regression through three of the four points: it leaves
out the 50 °C point for 3-HP → MA and the 25 °C point for the other two. The
lines give activation energies of 59.93, 68.96 and 52.09 kJ/mol; the text
gives 59.8, 91.6 and 59.5 kJ/mol. No subset of the plotted points gives 91.6
or 59.5 kJ/mol. The text values are those of the heat-map script; where they
come from is not known.

Fig. 11 is reproduced with the simplified model of the deposited script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The text presents the heat maps as the O2 model extended with the Arrhenius
relation. The caption says that "Fig. 9b data was used as a reference kref
of the Arrhenius equation", which matches the deposited heat-map script:
the script fits a model of its own to the Fig. 9b data and extrapolates
that model. It is a different model from that of Table 3 (five reactions,
no OPA, a reverse term in every step) and is not written out in the paper.
The published figure shows three properties of that script, which are
reproduced here:

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

The two panels of the published figure also need different constant vectors,
for reasons that are not known (separate fits, or other settings). The
fitted constants are not in the deposit, and no single set reproduces both
panels: with the set of panel b panel a is off by 1.07 px, with the set of
panel a panel b by 1.45 px. Two rate constants (and, for panel a, two more
constants) were inferred by fitting the digitized contour lines (not
shipped). The other constants are in the deposit: three (panel b) or two
(panel a) rate constants are the 50 °C row of the deposited Arrhenius
table, and the other equilibrium constants (all five for panel b, four for
panel a, whose ``K_10`` is fitted) and the five activation energies are
printed in the script.

The script also uses a pressure of 3.5 bar (the captions give 3 bar of O2).
What Fig. 11 calls "MA selectivity" and "MA yield" are ``X_MA`` and
``(1 - X_HP3)*X_MA``.

The maximum yield of 55.2 % is not reproduced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cell the text names for the maximum yield (84 °C, 49 min) is the largest
of the reproduced map, but its value is 0.548. The crest is flat (13 cells
between 76 and 90 °C lie within 0.001 of the maximum), so not much should
be read into which cell is the largest. None of the constant sets found
that match the contour lines gives 0.552, and the colour of the published
map along its crest reads 0.549 ± 0.006 (both not shipped). The maximum
selectivity of the text, 59.1 % (95 °C, 27 min), is reproduced (0.5906).

The text sets two "observed" maxima against these: 56.9 % for the
selectivity ("the observed maximum of 56.9 % at 50 °C (purple star on
plot)") and 50.5 % for the yield. The abstract places both at 3 bar of O2
and 50 °C, the condition of the Fig. 9b experiment, and the paper does not
say how either was computed. The 50.5 % equals the quantity of panel b,
``(1 - X_HP3)*X_MA``, of that experiment at 360 min (0.505), which is an
inference from the match. The 56.9 % is not identified: the quantity of
panel a, ``X_MA``, is at most 0.548 in that experiment (240 min), and MA
per 3-HP converted, ``X_MA/(1 - X_HP3)``, is 0.621, 0.577 and 0.513 at 240,
360 and 450 min and 0.569 at no sampling time. The 50 °C row of the
reproduced selectivity map peaks at 0.5705 at 255.8 min, next to the purple
star, but that is a model value. The match of the 50.5 % therefore stands
alone and does not settle how the paper defines the two quantities.

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

The package and its tests establish that the shipped models, with the
constants as documented, agree with the deposited curves, the digitized
lines and the experiments to the numbers given on this page, that with the
authors' solver run repeated they return the deposited model nodes to
their printed digits, and that each model implements its documented rate
laws. They do not re-derive the inferred constants. The following were
obtained with scripts that are not shipped, so they cannot be recomputed
from the repository:

- the least-squares fits made through the repeated solver run: the constant
  sets of the Fig. 9 nodes, the thirty constants of the reversible O2
  variant and their uncertainties, and the refits of the alternatives
  excluded for it (other water stoichiometries, a refitted linear reverse
  of reaction 1, the O2 fraction held constant). The tests show that the
  shipped constants return the nodes, not that no other constants or rate
  laws would. MATLAB itself was never run;
- the fits to the digitized contour lines of Fig. 11 that gave the inferred
  constants of the simplified reversible model;
- the refits to the experiments: of the O2 constants without the 120 min
  row, and of the simplified reversible model behind the reading of the
  Fig. 10d points;
- the independent integrations behind the comparison with the Table S4
  orders, and the repeat of the heat-map script's solver runs for the rows
  of Fig. 11 other than 100 °C (``NOTES.md``);
- the digitization of the published figures (the digitized curves
  themselves are shipped).
