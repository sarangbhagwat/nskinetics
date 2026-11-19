.. _pd-c-3hp-oxidation:

Pd/C oxidation of 3-HP to malonic acid
======================================

.. The images on this page are our own reproductions; no image of the paper
   is in the repository. Each is drawn by the plot_* function of
   nskinetics/models/pd_c_3hp_oxidation/figures.py named in its caption,
   exactly as save_figures() draws it (inside
   plt.rc_context(figures._RC_PARAMS)), saved with
   fig.savefig(path, dpi=100, bbox_inches='tight', facecolor='white') and
   converted to an opaque RGB PNG. *.png is gitignored; .gitignore makes an
   exception for docs/source/_static/images/pd_c_3hp_oxidation/*.png.

.. Each python code block below is followed by a text block with what it
   prints. nskinetics/tests/test_pd_c_3hp_oxidation_docs.py runs the blocks
   in order and compares the two, so keep that pairing when editing.

``nskinetics.models.pd_c_3hp_oxidation`` holds kinetic models of the
selective oxidation of 3-hydroxypropionic acid (3-HP) to malonic acid (MA)
over 10 % Pd/C, with O2 or H2O2 as the oxidant, after

    S. Y. Kim, J. Chavarrio, J. Install, G. W. Huber, "Selective oxidation
    of 3-hydroxypropionic acid to malonic acid over Pd/C: Mechanistic and
    kinetic study", *Appl. Catal. B: Environ. Energy* 386 (2026) 126403,
    https://doi.org/10.1016/j.apcatb.2026.126403

and the authors' data deposit (Dryad,
https://doi.org/10.5061/dryad.jsxksn0q4, CC0), called *the deposit* below.
The package provides:

- four models, each an Antimony file with its generated SBML file and a
  loader that returns a :class:`~nskinetics.KineticModel` (time in min,
  concentrations in M);
- the deposited tables of Figs. 9 and 10 (verbatim copies), curves digitized
  from SI Fig. S12 and Fig. 11, and readers for them;
- functions that redraw Fig. 9, SI Fig. S12, Fig. 10 and Fig. 11 from the
  models, and three figures that compare the models with the published
  ones;
- ``NOTES.md``, the full record of what is printed in the paper, what was
  inferred from the deposit and where the two disagree. This page
  summarizes it.

The package needs no optional dependencies, and importing it loads no
model. Its API is documented on :doc:`the shipped-models page </API/models>`.

The model curves of the paper's Figs. 9, 10 and 11 were not all made with
the model the paper describes. Redrawing them took four models, two of
which are not in the paper, and one factor that is stated in neither the
paper nor the deposit. `What does not match the paper, and what is inferred`_
lists these points; read it before using a model for anything other than
redrawing the figures.

.. warning::

   ``load_o2_model(T)`` extends the published 50 °C O2 model to other
   temperatures with the paper's Arrhenius equation. Only the 50 °C model is
   a published one. Away from 50 °C it reproduces neither the published heat
   maps (Fig. 11) nor the experiments at 25, 75 and 100 °C (Fig. 10): see
   `Limits of the temperature-extended O2 model`_.

The four models
---------------

The paper presents one kinetic model per oxidant (Table 3). Its figures were
made with four models, and each is shipped because each is needed to
reproduce a published item. Each model is an Antimony file with the SBML
file generated from it next to it (same name, ending in ``_sbml.xml``).

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
     - Fig. 9a and SI Fig. S12a (50 °C)
     - **Paper**: Table 3, H2O2 column, as printed. The rate laws were
       reconstructed from the deposited model curve (the authors posted only
       the O2 version of their script).
   * - O2 model
     - | ``load_o2_model(T)``
       | ``pd_c_3hp_oxidation_o2_antimony.txt``
     - Fig. 9b and SI Fig. S12b (50 °C)
     - **Paper**: Table 3, O2 column, as printed, with one **inferred**
       factor (``x_O2_ref``). Its temperature dependence is an NSKinetics
       addition that reproduces no published item.
   * - Reversible O2 variant
     - | ``load_o2_reversible_model(temperature_C)``
       | ``pd_c_3hp_oxidation_o2_reversible_antimony.txt``
     - The model curves of Fig. 10a–c (25, 75 and 100 °C), one constant set
       per temperature
     - **Inferred** from the deposited curves: structure and all constants.
       It is in neither the paper nor the authors' scripts.
   * - Simplified reversible O2 model
     - | ``load_o2_simplified_reversible_model()``
       | ``pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt``
     - Fig. 11; three of its rate constants are the points of Fig. 10d
     - **Deposit**: the authors' heat-map script. It is not in the paper.
       The constants that are not in the deposit were **inferred** from the
       published contour lines.

The two Table 3 models describe the batch of the Fig. 9 caption: 60 mM 3-HP,
NaOH/3-HP = 1:1, 0.2 g Pd/C, 50 °C, and 3 bar O2 or 15 equivalents of H2O2
under 3 bar Ar. Their species are ``HP3``, ``OPA`` (3-oxopropanoic acid, the
aldehyde intermediate, never measured), ``MA``, ``AAc`` (acetic acid),
``OXA`` (oxalic acid) and ``D`` (lumped other by-products). The quantities
plotted in the paper are assignment rules of the models and are selected by
name: ``X_HP3``, ``X_MA``, ``X_AAc``, ``X_OXA`` and ``X_D`` are each species
over the sum of the five measured species (OPA is left out of the sum, as in
the authors' post-processing).

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

The O2 model takes a temperature in K (``load_o2_model(T=...)``, or
``km.set_value('T', ...)`` on a loaded model), which scales its six rate
constants by Arrhenius factors. The example compares it with the 75 °C
experiment, next to the fit the authors plotted for that temperature. Read
`Limits of the temperature-extended O2 model`_ before relying on it:

.. code-block:: python

   from nskinetics.models.pd_c_3hp_oxidation import (
       load_o2_reversible_model, read_fig10_o2)

   km_75 = load_o2_model(T=75 + 273.15)
   factor = km_75.get_value('k_1')/km_75.get_value('k_1_ref')
   print(f'k_1(75 C)/k_1(50 C) = {factor:.2f}')

   experiment = read_fig10_o2(75)['experiment']   # deposited points, Fig. 10b
   for name, model in (('Table 3 O2 model at 75 C', km_75),
                       ('per-temperature fit', load_o2_reversible_model(75))):
       sim = simulate_at_times(model, experiment['time'])
       print(f'{name}: pooled R2 = {pooled_r2(experiment, sim):.3f}')

.. code-block:: text

   k_1(75 C)/k_1(50 C) = 4.94
   Table 3 O2 model at 75 C: pooled R2 = 0.644
   per-temperature fit: pooled R2 = 0.990

Set and read model values with ``set_value`` and ``get_value``, not through
attributes of the underlying RoadRunner object: attribute access is
unreliable once more than one model has been loaded in a process. ``reset()``
restores the species but keeps a changed parameter such as ``T``.

The heat maps of Fig. 11 are computed by ``fig11_heat_maps``, which runs the
simplified reversible model as the authors' script does (142 integrations).
``simulate_heat_maps`` does the same for any model with a temperature
parameter on any grid; its time axis is the reaction time unless a
``time_offset`` is given:

.. code-block:: python

   from nskinetics.models.pd_c_3hp_oxidation import (
       fig11_heat_maps, heat_map_maximum,
       load_o2_simplified_reversible_model, simulate_heat_maps)

   T_C, times, selectivity, MA_yield = fig11_heat_maps()   # 71 x 300 grids
   for matrix in (selectivity, MA_yield):
       value, T_max, t_max = heat_map_maximum(T_C, times, matrix)
       print(f'{value:.4f} at {T_max:g} C, {t_max:.1f} min')

   T_C, times, selectivity, MA_yield = simulate_heat_maps(
       load_o2_simplified_reversible_model(), T_C=[50, 75, 100],
       times=[30, 60, 120, 240])
   print(selectivity.round(3))

.. code-block:: text

   0.5906 at 95 C, 26.8 min
   0.5483 at 84 C, 49.1 min
   [[0.189 0.324 0.483 0.571]
    [0.535 0.584 0.449 0.226]
    [0.441 0.207 0.045 0.004]]

The figure module writes every figure of this page, as PNG (600 dpi) and
PDF, to a folder; ``--figures`` and ``--formats`` select a subset:

.. code-block:: bash

   python -m nskinetics.models.pd_c_3hp_oxidation.figures OUTPUT_DIR
   python -m nskinetics.models.pd_c_3hp_oxidation.figures OUTPUT_DIR \
       --figures fig10 fig10_table3_o2 --formats png

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

At the nodes of the deposited curves the models differ from them by at most
9.6e-4 (H2O2, 24 nodes) and 2.0e-3 (O2, 19 nodes after t = 0). These
residuals are the solver error of the deposited curves and the rounding of
Table 3, not a difference in the model (see `Published lines are
loose-tolerance solver output`_). Against the experimental points the pooled
R2 is 0.9945 for O2 (16 rows) and 0.945 for H2O2 (11 rows; 0.941 on all 13).

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
Table 3 model lies up to 0.38 equivalents (rms 0.22) above the published
trace, because that trace was generated with a different parameter set (see
`The H2O2 trace comes from a different fit`_).

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
   reversible O2 variant with the constant set of each temperature (solid
   lines) and the deposited model curves (black dotted lines); the box gives
   the pooled R2 against the experimental points to three decimals. (d) The
   deposited points of the Arrhenius plot and the three dashed lines,
   recomputed by least squares.

The reversible O2 variant differs from the deposited curves by at most
1.8e-4 (25 °C), 2.5e-3 (75 °C) and 2.7e-3 (100 °C) at their nodes, again the
solver error of those curves. Its pooled R2 against the experimental points
is 0.9996, 0.990 and 0.993; the published panels print 0.99. This variant
was inferred from the very curves it is compared with, so the agreement
shows that the reconstruction is consistent, not that the model is right
(see `Fig. 10a–c come from an unposted reversible variant`_).

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
   its own constant set.

The deposit holds no table of the heat maps, so the reference is the
published figure: its contour lines at 0.1, 0.3 and 0.5, digitized. The
reproduced contours lie within an rms of 0.17 px (panel a) and 0.27 px
(panel b) of them, in pixels of the published raster (1 px = 3.07 min or
0.148 °C); the digitization itself is good to about 0.3 px. The largest cell
of the selectivity map, 0.5906 at 95 °C and 26.8 min, is the 59.1 % (95 °C,
27 min) of the text. The largest cell of the yield map is the one the text
names (84 °C, 49 min), but its value is 0.548, not the 55.2 % of the text.

.. figure:: /_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11_verification.png
   :target: ../_static/images/pd_c_3hp_oxidation/pd_c_3hp_oxidation_fig11_verification.png
   :alt: Top row: the two reproduced heat maps with the digitized published
         contour lines, ridge line and stars drawn on top as fine black
         dots; the dots lie inside the wide reproduced lines everywhere.
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
the Table 3 O2 model follows the Arrhenius form of the paper's Eq. 6,
``k(T) = k_ref*exp(-(Ea/R)*(1/T - 1/T_ref))``, with the reference at 50 °C.
At 50 °C it is the published model by construction. Three activation
energies are those quoted in the text (59.8, 91.6 and 59.5 kJ/mol; the
first is applied to both R1 and R2); the other two (51.79 kJ/mol for
3-HP → D and 10 kJ/mol for MA → OXA, a placeholder) appear only in the
authors' heat-map script. Applying them to this six-reaction model is an
NSKinetics addition: the text describes Fig. 11 this way, but none of the
paper's temperature results was computed with it.

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

The offsets of the Table 3 model are linearized estimates. Measured on the
lines themselves, the late side of every contour is earlier than published
by 8 to 11 % of the time in the median per line (28 to 74 min), which is the
same as 1.3 to 1.8 °C lower at a given time. Shifting its time axis by the
10 min of the published figure leaves 6.8 and 6.1 px. The difference is 20
times the reading uncertainty of the figure or more.

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
fits the authors plotted (``load_o2_reversible_model``):

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
   * - Per-temperature fit: pooled R2 / rms
     - 0.9996 / 0.008
     - 0.990 / 0.029
     - 0.993 / 0.026
   * - Time to half conversion of 3-HP, Table 3 at T / per-temperature fit
     - 458 min / not reached
     - 16.5 / 75.5 min
     - 4.5 / 22 min
   * - Maximum MA fraction, Table 3 at T
     - 0.522 at 1479 min
     - 0.584 at 51.5 min
     - 0.587 at 13 min
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
         per-temperature fits as dashed lines. The solid lines run far from
         the points: at 25 degrees the model consumes most of the 3-HP while
         the measured 3-HP stays near 0.9; at 75 and 100 degrees the model's
         malonic acid peak comes much earlier than the measured one. The
         boxes give pooled R2 values of 0.579, 0.644 and 0.746 for the Table
         3 model and 1.000, 0.990 and 0.993 for the fits. Panel d: the
         Arrhenius points with the dashed regression lines and solid lines
         with the activation energies of the text through the 50 degree
         points; the solid lines miss most of the other points.
   :class: pub-figure
   :width: 100%

   Drawn by ``plot_fig10_table3_o2``. (a–c) The experiments of Fig. 10 with
   the temperature-extended Table 3 O2 model (solid lines) and the
   per-temperature fits (dashed lines); the boxes give the pooled R2 of
   both, to three decimals. (d) The points of Fig. 10d with the published
   regression lines (dashed) and the lines that the activation energies of
   the text give through the 50 °C points (solid).

**Why.** The 50 °C constant of 3-HP conversion lies well off the Arrhenius
regression through the other temperatures. The measured 3-HP conversion at
75 °C is no faster than at 50 °C (half converted between 60 and 90 min at
both), whereas the Arrhenius factor of the model from 50 to 75 °C is 4.9
for R1 and R2. The same shows in Fig. 10d, where the 50 °C point of the
3-HP → MA constant lies a factor of 4.3 above the regression line through
the other three points. An extrapolation anchored at 50 °C with the
activation energies of the text therefore overshoots at 25, 75 and 100 °C.
Part of the difference from the published heat maps is also structural: the
Table 3 model has the OPA intermediate and no reverse terms, and its O2
factor stays above 0.999 where that of the simplified model falls to 0.83.

**What to use instead.** For the curves of Fig. 10a–c use
``load_o2_reversible_model``; its three constant sets are separate fits
with no temperature law, so it exists only at 25, 75 and 100 °C. For Fig. 11
use ``load_o2_simplified_reversible_model``. The authors' heat-map script
extrapolates that model from 50 °C in the same way (its 3-HP → MA constant
at 75 and 100 °C is 4.3 and 4.2 times the value fitted to the experiments at
those temperatures); how far the maps of Fig. 11 are from the experiments
at 75 and 100 °C has not been measured here.

What does not match the paper, and what is inferred
---------------------------------------------------

``NOTES.md`` in the package folder gives the evidence for each point below
and lists further, smaller ones.

``x_O2_ref`` is inferred
~~~~~~~~~~~~~~~~~~~~~~~~

The Table 3 O2 constants of R1, R2 and R4 reproduce the deposited Fig. 9b
curve only when multiplied by the catalyst mass (0.2 g) and by a further
factor of 0.005085. That factor equals the H2O2 mole fraction of the
15-equivalent H2O2 run on the authors' scale (15 × 0.000339). It is not in
the paper, the SI, the README of the deposit or the posted script, and how
it entered the authors' code is unknown. Only the products
``k_i*m_cat*x_O2_ref`` are identifiable. The posted fitting script, run as
posted, misses the deposited curve by an rms of 0.42 to 0.76 for 3-HP, MA,
OXA and D.

Table S4 is not used
~~~~~~~~~~~~~~~~~~~~

SI Table S4 lists fitted reaction orders and DFT equilibrium constants.
Neither enters the rate laws: the deposited curves are reproduced only by
first-order, irreversible steps, and the posted script uses neither.
Applying the Table S4 orders moves the curves by up to 0.17 (O2) and 0.04
(H2O2). The values are kept as ``TABLE_S4`` for documentation.

Published lines are loose-tolerance solver output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The deposited model curves of Figs. 9 and 10 are the raw steps of MATLAB
``ode23s`` at its default tolerances, 18 to 24 per curve, and the published
lines are straight segments between them. They are not converged solutions:
the nodes are up to 2.0e-3 (Fig. 9b) and 2.7e-3 (Fig. 10b and c) away from
the accurately solved model, and the segments between them up to 0.010. The
models here are solved to a relative tolerance of 1e-9, so they do not lie
exactly on the published polylines, and the package tests allow for that
solver error.

The R2 values of the paper are not reproduced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SI Fig. S11a prints R2 = 0.997 for O2. The O2 model gives 0.9945 against
the 16 deposited experimental rows (0.9957 without the 120 min row); the
source of 0.997 could not be identified. The text states R2 > 0.95 for both
oxidants. The H2O2 model gives 0.945, and only on 11 rows (without t = 20
and 1205 min); on all 13 rows it gives 0.941.

A hidden OPA pool and a unit slip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OPA was never measured and is left out of the sum that normalizes the
plotted fractions, so those add up to 1 whatever the model holds as OPA. The
text gives the predicted OPA peak with O2 as 8.4e-4 mM near 25 min and calls
it negligible. That number is the mole fraction on the authors' scale
(1.4e-5) times 60 mM. In the model the OPA peak is 4.2 % of the initial 3-HP
(2.5 mM), at 20.5 min. The other models hold more: 11.5 % of the initial
3-HP at 84 min in the H2O2 model (7.8 % is still there at 1440 min), and
peaks of 2.9, 12.6 and 20.1 % in the reversible O2 variant at 25, 75 and
100 °C.

The H2O2 trace comes from a different fit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The predicted H2O2 of SI Fig. S12a (the deposited trace and the published
inset) was generated with a different parameter set than the Fig. 9a curve
and Table 3: recovered from the trace, its ``k_7`` is 6.94e-3 1/min instead
of 6.48e-3, and its ``k_5`` is 22 % lower. The species curves of the same
figure are the Table 3 ones. The Table 3 model therefore predicts up to 0.38
equivalents more H2O2 than the published trace. In the deposited file of
this figure the blocks labelled ``exp`` and ``model`` are also swapped; the
reader ``read_fig9a_h2o2`` corrects that.

Fig. 10a–c come from an unposted reversible variant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The paper prints no constants for 25, 75 or 100 °C, and its text treats all
reactions as irreversible. The model curves deposited for these temperatures
are produced by neither posted script. They are reproduced by the reaction
network of the Table 3 O2 model with reactions 1, 2, 5 and 6 made
reversible, fitted separately at each temperature. Its structure and all of
its constants were inferred from the deposited curves; nothing of it is
printed in the paper, the SI or the deposit. The three constant sets follow
no temperature law (the forward constant of 3-HP → OPA is 0.047, 0.0092 and
0.041 1/min at 25, 75 and 100 °C), so they are fit parameters of each
curve, not kinetic quantities. The model is a reproduction aid.

Fig. 10d and the activation energies of the text
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The points of the Arrhenius plot are three rate constants of the simplified
reversible model (3-HP → MA, MA → AAc and MA → D) fitted to the experiment
of each temperature. They are constants of neither the Table 3 model nor
the model that drew the curves of panels a to c. Each dashed line is a
regression through three of the four points: it leaves out the 50 °C point
for 3-HP → MA and the 25 °C point for the other two. The lines give
activation energies of 59.93, 68.96 and 52.09 kJ/mol; the text gives 59.8,
91.6 and 59.5 kJ/mol. No subset of the plotted points gives 91.6 or
59.5 kJ/mol. The text values are those of the heat-map script; where they
come from is not known.

Fig. 11 was computed with the simplified model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The text and the caption present the heat maps as the O2 model of Fig. 9b
extended with the Arrhenius equation. The deposited script computes them
with a different model (five reactions, no OPA, a reverse term in every
step) that is not written out in the paper. The published figure shows three
properties of that script, which are reproduced here:

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
constants) were inferred from the digitized contour lines.

The script also uses a pressure of 3.5 bar (the captions give 3 bar of O2).
What Fig. 11 calls "MA selectivity" and "MA yield" are ``X_MA`` and
``(1 - X_HP3)*X_MA``.

The maximum yield of 55.2 % is not reproduced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cell the text names for the maximum yield (84 °C, 49 min) is the largest
of the reproduced map, but its value is 0.548. No constant set that matches
the contour lines gives 0.552, and the colour of the published map along its
crest reads 0.549 ± 0.006. The maximum selectivity of the text, 59.1 %
(95 °C, 27 min), is reproduced (0.5906).
