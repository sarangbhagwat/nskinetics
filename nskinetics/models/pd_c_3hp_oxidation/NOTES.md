# Pd/C 3-HP oxidation models: provenance notes

Models of the selective oxidation of 3-hydroxypropionic acid (3-HP) to malonic
acid (MA) over 10 % Pd/C, after

> S. Y. Kim, J. E. Chavarrio, J. Install, G. W. Huber, "Selective oxidation of
> 3-hydroxypropionic acid to malonic acid over Pd/C: Mechanistic and kinetic
> study", Appl. Catal. B: Environ. Energy 386 (2026) 126403,
> https://doi.org/10.1016/j.apcatb.2026.126403

and the authors' data deposit (Dryad, https://doi.org/10.5061/dryad.jsxksn0q4,
CC0), referred to below as "the deposit".

These notes separate what is **printed** in the paper or its SI, what was
**inferred** from the deposit because the paper does not state it, and what
is an **NSKinetics addition**. They also list the places where the paper and
the deposit disagree and that matter for these models.

Two conventions are used throughout:

- How the authors produced a figure is known only where the paper or the
  deposit says so. Elsewhere these notes say what the numbers show ("is
  reproduced by", "agrees with") and not what the authors did.
- Several results quoted here come from analyses that are **not in the
  repository**: fits of constants to the deposited curves and to digitized
  contour lines, refits of alternative models, and readings of the colour
  field of a published figure. Each such result is marked "(not shipped)".
  The last section lists them. Everything else can be recomputed with the
  package; the main results are pinned by its tests. In particular the
  repeat of the authors' solver run is shipped (`solver_emulation.py`): that
  the shipped models return the deposited model nodes to their printed
  digits is tested, and so is what each inferred choice contributes to that.

## The four models

The paper presents one kinetic model per oxidant. Redrawing the model curves
of its Figs. 9, 10 and 11 takes four models, and each is shipped because each
is needed to reproduce a published item:

| Model (loader) | In the paper? | Reproduces |
|---|---|---|
| H2O2 model (`load_h2o2_model`) | Table 3, H2O2 column | Fig. 9a (50 C). Its predicted H2O2 is compared with SI Fig. S12a, whose H2O2 curve lies up to 0.38 equivalents lower and is matched only with other constants |
| O2 model (`load_o2_model`) | Table 3, O2 column, with one inferred factor | Fig. 9b, SI Fig. S12b (50 C). Its Arrhenius extension is an NSKinetics addition and reproduces no published item |
| Reversible O2 variant (`load_o2_reversible_model`) | no; inferred from the deposited curves | the model curves of Fig. 10a-c (25, 75, 100 C), one constant set per curve |
| Simplified reversible O2 model (`load_o2_simplified_reversible_model`) | no; in the authors' deposited heat-map script | Fig. 11; the 50 C point of Fig. 10d supplies three of its rate constants |

The last two models reproduce Fig. 10a-c and Fig. 11 with constants that were
fitted to the very curves and contour lines they are compared with: all
thirty constants of the reversible variant, and two (panel b) or four (panel
a) of the ten 50 C constants of the simplified model. The other constants of
the simplified model come from the deposit: three (panel b) or two (panel a)
rate constants from the 50 C row of the deposited Arrhenius table, five or
four equilibrium constants as printed in the heat-map script, and the five
activation energies of that script ("Constants at 50 C" below). That
agreement shows that the models are consistent with the published figures,
not that they are right.

## Files

| File | Content |
|---|---|
| `pd_c_3hp_oxidation_h2o2_antimony.txt` (+ `_sbml.xml`) | H2O2 model: Table 3 H2O2 column, Fig. 9a |
| `pd_c_3hp_oxidation_o2_antimony.txt` (+ `_sbml.xml`) | O2 model: Table 3 O2 column, Fig. 9b, SI Fig. S12b, plus Arrhenius temperature dependence |
| `pd_c_3hp_oxidation_o2_reversible_antimony.txt` (+ `_sbml.xml`) | reversible six-reaction O2 variant inferred from the deposited model curves of Fig. 10a-c; not in the paper or the authors' scripts |
| `pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt` (+ `_sbml.xml`) | the authors' simplified reversible O2 model (no OPA) of their heat-map script, which reproduces Fig. 11; not in the paper |
| `pd_c_3hp_oxidation.py` | loaders (`load_h2o2_model`, `load_o2_model`, `load_o2_reversible_model`, `load_o2_simplified_reversible_model`), `simulate_at_times`, `TABLE_S4`, the constant sets of the reversible variant and of the simplified model, and the constant sets inferred from the deposited Fig. 9 nodes (`H2O2_NODE_PARAMETER_SETS`, `O2_NODE_PARAMETER_SETS`) |
| `solver_emulation.py` | a port of MATLAB `ode23s` and the repeat of the authors' solver run with the rates of a shipped model, for the comparison with the deposited model nodes (verification only) |
| `heat_maps.py` | the temperature-time maps of Fig. 11 and their comparison with digitized contour lines |
| `temperature_series.py` | the lines of the Arrhenius plot (Fig. 10d) and the measures used to compare a model with the experiments of Fig. 10a-c |
| `reference.py`, `reference_data/` | deposited Fig. 9 and Fig. 10 tables, digitized Fig. S12 inset curves and Fig. 11 lines, and their readers |
| `figures.py` | regenerates Fig. 9, Fig. 10, SI Fig. S12 and Fig. 11, and draws the comparison and verification figures |

The SBML files are generated from the Antimony files (the Antimony is the
source; `//` comments exist only there).

## Documentation

The Sphinx page `docs/source/case_studies/pd_c_3hp_oxidation.rst` ("Pd/C
oxidation of 3-HP to malonic acid") summarizes these notes for users, and
`docs/source/API/models.rst` documents the API. These notes are the source:
every number on those pages is either taken from here or printed by a code
example on the page. `nskinetics/tests/test_pd_c_3hp_oxidation_docs.py` runs
the Python examples of the page and compares what they print with the output
shown there, to within one unit of the last printed digit. The shell
commands of the page (the figure CLI) were run by hand; no test runs them.

The images of that page (`docs/source/_static/images/pd_c_3hp_oxidation/`)
are the eight figures of `figures.py`, written from the repository root by

    python -m nskinetics.models.pd_c_3hp_oxidation.figures \
        docs/source/_static/images/pd_c_3hp_oxidation --formats png --dpi 100

(the default resolution of the figure CLI is 600 dpi; the comment at the top
of the page gives the same command). They show the models with the deposited
and digitized data overlaid; no image of the paper is in the repository.

## Model basis

This section and the next seven are mainly about the two Table 3 models.
The simplified reversible model and the reversible O2 variant have their own
sections further down.

Conditions of both models: 60 mM 3-HP, NaOH/3-HP = 1:1, 0.2 g Pd/C, 50 C;
3 bar O2, or 15 equivalents of H2O2 under 3 bar Ar (Fig. 9 caption).

The authors' scripts work in liquid mole fractions and plot fractions of the
initial 3-HP. Here the liquid species are concentrations in M:
`c_i = C0 * n_i`, with `n_i` the moles of species i per initial mole of 3-HP
and `C0 = 0.06` M. The mapping is exact because every rate law is first order
in the liquid species. Time is in min and temperatures are in K.

States: `HP3`, `OPA` (3-oxopropanoic acid, the aldehyde intermediate; "not
detected with our analytical methods", Section 3.3 of the paper), `MA`,
`AAc` (acetic acid), `OXA` (oxalic acid), `D` (lumped other by-products),
plus `H2O2` (H2O2 model) or the dimensionless O2 activity `yO2` (O2 model).
Co-products (H2O, H2, O2, CO2) are not tracked.

Rates, with `k_i` the Table 3 constants as printed, in 1/(g_cat min):

| | O2 model | H2O2 model |
|---|---|---|
| R1: 3-HP -> OPA | `k_1*m_cat*x_O2_ref*yO2*[HP3]` | `k_1*m_cat*x_per_M*[H2O2]*[HP3]` |
| R2: OPA -> MA | `k_2*m_cat*x_O2_ref*yO2*[OPA]` | `k_2*m_cat*x_per_M*[H2O2]*[OPA]` |
| R3: MA -> AAc + CO2 | `k_3*m_cat*[MA]` | same |
| R4: MA -> OXA | `k_4*m_cat*x_O2_ref*yO2*[MA]` | `k_4*m_cat*x_per_M*[H2O2]*[MA]` |
| R5: 3-HP -> D | `k_5*m_cat*[HP3]` | same |
| R6: MA -> D | `k_6*m_cat*[MA]` | same |
| R7: H2O2 -> H2O + 0.5 O2 | - | `k_7*[H2O2]`, `k_7` in 1/min, no `m_cat` |

O2 model: `d yO2/dt = -x_per_M*(R1 + R2 + R4)`. H2O2 model: H2O2 is consumed
1:1 by R1, R2, R4 and by R7.

Plotted quantities are assignment rules and can be selected by name:
`X_HP3, X_MA, X_AAc, X_OXA, X_D` = `c_i/([HP3]+[MA]+[AAc]+[OXA]+[D])`
(Fig. 9; OPA is left out of the sum, as in the authors' post-processing),
`x_OPA = x_per_M*[OPA]` (SI Fig. S12b inset) and `eq_H2O2 = [H2O2]/C0`
(SI Fig. S12a inset).

## Operating conditions, units and the range of the data

Parameters a user may change with `set_value`, per model (the documentation
page has the same table with a worked example):

| Model | Parameter | Unit | Default | Meaning |
|---|---|---|---|---|
| H2O2, O2 | `m_cat` | g | 0.2 | catalyst mass of the paper's batch; the rates of R1 to R6 are proportional to it |
| all four | `C0` | M | 0.06 | initial 3-HP concentration (takes effect at the next `reset()`) |
| H2O2 | `eq_H2O2_0` | - | 15 | initial H2O2 in equivalents of `C0` (0.9 M; takes effect at the next `reset()`) |
| O2, simplified | `T` | K | 323.15 | temperature |
| simplified | `P` | bar | 3.5 | pressure of the heat-map script. Not for exploring pressure: the paper's Fig. 3a contradicts the model above 3 bar (see below) |

Rate constants: `k_1` to `k_6` of the Table 3 models are in 1/(g_cat min)
and `k_7` (H2O2) in 1/min. The constants of the reversible variant and the
rate constants of the simplified model are lumped, in 1/min for the paper's
batch; those two models have no catalyst-mass parameter.

**`m_cat` is a mass per batch of the paper, not per litre.** The compartment
of the Antimony files is a nominal 1 L that cancels from every concentration;
no rate has a volume term. The rates scale with `m_cat/0.2`, the catalyst
loading relative to the paper's, whatever the liquid volume: `m_cat = 0.4`
is twice the paper's loading, and `m_cat = 20` is not "20 g/L" but a hundred
times the paper's loading.

**The catalyst loading of the kinetic runs is not stated unambiguously.**
Section 2.3 gives the batch as 60 mM 3-HP and equimolar NaOH in 10 mL of
water with 0.2 g of catalyst (20 g/L). It adds that for time-resolved
product distributions "the amounts of substrate, NaOH, and water were each
increased threefold" and does not say whether the catalyst was. The captions
of the time-resolved Figs. 1, 9 and 10 give 0.2 g. The loading behind
Table 3 is therefore 20 g/L (if the catalyst was tripled too and the caption
quotes the base recipe) or 6.7 g/L (0.2 g in 30 mL). The models do not
depend on which: they take the batch as the unit.

**Only one condition per model was fitted.** Section 3.4 of the paper:
"Kinetic models for H2O2 and O2 were built with concentration profiles
collected at the optimized temperature of 50 C and equimolar NaOH (Fig. 9)".
The Table 3 constants therefore belong to one experiment each: 60 mM 3-HP,
0.2 g of catalyst, 50 C, with 15 equivalents of H2O2 (samples up to
1440 min) or 3 bar O2 (samples up to 1680 min). Any other value of `m_cat`,
`C0`, `eq_H2O2_0` or `T` is an extrapolation of the fitted rate laws. The
paper has experiments at other H2O2 charges, temperatures and O2 pressures;
the models were not built on them and do not describe them (next
subsection). For `T` with O2 the comparison with the time courses of
Fig. 10 is poor as well (section "Temperature-extended Table 3 O2 model
versus the 25, 75 and 100 C experiments"). For `m_cat` and `C0` there is
nothing to compare with: no experiment found in the paper or the deposit
varies the 3-HP concentration or the mass of the 10 % Pd/C (where captions
state them they are 60 mM and 0.2 g; SI Table S1 adjusts the amount of the
3 % and 5 % catalysts to the same metal-to-feed ratio). The constant sets of
the reversible variant each belong to one experiment (samples up to 1890,
1440 and 1560 min at 25, 75 and 100 C). The simplified model is fitted by
its script to the 50 C experiment only (samples up to 1680 min).

### Experiments of the paper that the models were not built on

The time courses behind the four models are those of Fig. 9 (50 C) and
Fig. 10a-c (25, 75 and 100 C). All of the O2 ones are at 3 bar of O2
(captions of Figs. 9 and 10), and the heat-map script uses 3.5 bar. The
paper also reports experiments in which one condition is varied. Most are
screening experiments analysed once, after 8 h (Section 2.3: the products
were "collected after the entire reaction"):

- Fig. 2, H2O2 as the oxidant: (a) H2O2/3-HP ratio, (b) NaOH/3-HP ratio,
  (c) temperature;
- Fig. 3, O2 as the oxidant: (a) O2 pressure, 1, 3, 5, 10 and 20 bar,
  (b) NaOH/3-HP ratio, (c) temperature;
- SI Table S1 (seven combinations of catalyst and oxidant) and SI Table S2
  (the temperature series of Fig. 3c with all products);
- Fig. 4, malonic acid alone (60 mM, no 3-HP) under 3 bar of O2: (a) after
  1300 min at 50, 60, 75 and 100 C, (b) time courses at 50 and 100 C.

The values behind Figs. 2 to 4 are in the deposit (`2-1` to
`2-3_H2O2_experiments.csv`, `3-1` to `3-3_O2_experiments.csv`, `4-1` and
`4-2_MA_degradation.csv`, and the same tables in the workbook). The measured
numbers below are quoted from those files, which are not shipped with the
package; where the text of the paper gives a number it agrees (72.8 %
conversion and 56.4 % MA selectivity at 3 bar of O2, 78.4 % conversion with
15 equivalents of H2O2). None of these experiments was used to build the
models. The model numbers were measured with the package and are pinned by
its tests: conversion `1 - [HP3]/C0` and MA yield `[MA]/C0`, the quantities
of the paper's Eqs. 1 and 3, after 480 min at 50 C unless varied.

| Series (8 h) | Condition | Measured: conversion / MA yield, % | Model: conversion / MA yield, % |
|---|---|---|---|
| O2 pressure, Fig. 3a | 1 bar | 21.7 / 2.7 | simplified model, `P` set: 87.7 / 30.0 |
| | 3 bar | 72.8 / 41.1 | 98.2 / 47.6 |
| | 5 bar | 63.2 / 20.5 | 99.7 / 51.4 |
| | 10 bar | 25.0 / 1.3 | 99.9 / 52.7 |
| | 20 bar | 25.9 / 0.6 | 100.0 / 51.2 |
| Temperature with O2, Fig. 3c | 25 C | 8.5 / 1.4 | O2 model at T: 53.0 / 33.4; simplified model at T: 53.5 / 36.3 |
| | 50 C | 72.8 / 41.1 | 99.0 / 48.0; 98.9 / 49.1 |
| | 75 C | 95.5 / 35.7 | 100.0 / 4.7; 99.9 / 6.0 |
| | 100 C | 84.6 / 2.1 | 100.0 / 0.0; 100.0 / 0.2 |
| H2O2 charge, Fig. 2a | 5 equivalents | 14.0 / 4.7 | H2O2 model, `eq_H2O2_0` set: 41.8 / 9.7 |
| | 10 | 38.1 / 20.3 | 54.6 / 24.3 |
| | 15 | 78.4 / 58.5 | 64.7 / 36.5 |
| | fourth bar (see below) | 56.4 / 28.4 | 72.6 / 45.5 (20 equivalents); 83.6 / 55.7 (30) |
| | fifth bar | 20.7 / 6.6 | 83.6 / 55.7 (30 equivalents); 94.2 / 58.2 (50) |

The simplified model is as its loader builds it (`'fig11b'` set, heat-map
start); with the script's fitting start the pressure column reads 85.5,
95.9, 97.9, 98.9 and 99.5 % conversion and 27.3, 44.5, 48.9, 51.1 and 49.7 %
MA. The last two bars of Fig. 2a are labelled 20 and 30 on the axis of the
published figure and `1: 30` and `1: 50` in the deposit; the text speaks of
ratios "from 1:5 to 1:30". Both readings are given.

Two rows of the deposited Fig. 3a table differ from, or repeat, the SI:

- The 5 bar row (conversion 63.16, C balance 58.92, selectivities 32.51,
  2.25, 0.89 and 64.35 %) equals, to the printed digits, the row "3% Pd/C
  (Evonik), O2 3 bar" of SI Table S1 (63.2, 58.9, 32.5, 2.3, 0.9, 64.3),
  the entry of another catalyst under 3 bar. The published Fig. 3a draws
  these values at 5 bar (conversion near 63 %, C balance near 59 %, as read
  from the image). One of the two may be a copied row; which one is not
  known. The decline with pressure does not rest on this row: the 10 and
  20 bar rows show it as well.
- The 3 bar row agrees with SI Tables S1 and S2 in conversion and MA
  selectivity (72.8 and 56.4 %). Its C balance and "Others" are 77.43 and
  35.1 %, where SI Table S1 ("Others"), SI Table S2 ("Unknown") and the
  50 C row of the deposited Fig. 3c table ("Others") give 75.5 and 38.0 %,
  and its "Oxalic acid" is 3.93 %, which is the GLA value of SI Table S2
  (OXA is 1.0 % there). None of these enters the numbers quoted here.

What the comparison shows:

- **O2 pressure: the model is contradicted by the paper, not merely
  untested.** Measured conversion and MA selectivity (12.3, 56.4, 32.5, 5.0
  and 2.3 %) are highest at 3 bar and fall above it. The paper: "both
  declining at high pressures, likely due to catalyst poisoning by excess
  O2"; SI Fig. S6 compares the spent catalysts of the 3 and 10 bar
  experiments. The simplified model was not fitted to these experiments,
  and its rate form cannot give that decline: `P` multiplies the forward
  rates of steps 1 and 3, so the conversion only rises with `P` (tested),
  and the MA fraction after 8 h stays near 0.5 at 10 and 20 bar, where about
  1 % was measured. The response is strong at short times too: at 50 C, with
  `P` = 1, 3.5 and 10 bar, `X_HP3` is 0.77, 0.56 and 0.25 at 60 min and
  `X_MA` 0.27, 0.57 and 0.70 at 240 min. **Do not use `P` to explore
  pressure**; leave it at the script's 3.5 bar. The O2 model and the
  reversible variant have no pressure parameter. `x_O2_ref` of the O2 model
  is a factor inferred for the 3 bar experiment and is not a pressure:
  scaling it would make R1, R2 and R4 proportionally faster, which is the
  same monotonic response.
- **H2O2 charge: contradicted above 15 equivalents.** Measured conversion
  and MA yield are highest at 15 equivalents and fall beyond. The paper
  explains the lower yield by a decreased conversion, "likely due to active
  site occupation of surface peroxide", and by overoxidation to AAc and
  OXA. In the H2O2 model the conversion rises with every equivalent added,
  and the MA yield after 8 h rises up to 42 equivalents (59.0 %) and falls
  beyond (58.2 % at 50), so the model's optimum is far beyond the measured
  one at 15 (tested on the whole numbers from 1 to 100 equivalents). Below
  15 equivalents the model is too fast (41.8 against 14.0 % at
  5 equivalents), and at 15 it gives 64.7 against 78.4 %.
- **Temperature.** At 25 C the models convert 53 % in 8 h where 8.5 % was
  measured, as in the comparison with the 25 C time course. At 75 C they
  have 5 to 6 % MA left after 8 h where the screening experiment has a
  yield of 35.7 %; this holds for the simplified model, i.e. the model of
  Fig. 11, as well.
- **The screening experiments and the time courses are not on the same
  footing.** At the same nominal conditions they differ from each other,
  in both directions:

  | Condition | Screening, 480 min | Time course |
  |---|---|---|
  | O2, 50 C | 72.8 % conversion | Fig. 9b: `1 - X_HP3` = 0.961 at 450 min (0.972 at 540 min); Fig. 1a: 1.45 of 59.8 mM 3-HP left at 450 min, 97.6 % |
  | H2O2, 50 C | 78.4 % conversion | Fig. 9a: `1 - X_HP3` = 0.606 at 420 min; Fig. 1b: 22.2 of 56.4 mM left at 420 min, 60.6 % |
  | O2, 75 C | MA yield 35.7 % | Fig. 10b: `X_MA` = 0.098 at 420 min and 0.011 at 540 min |
  | O2, 100 C | 84.6 % conversion | Fig. 10c: `1 - X_HP3` = 0.967 at 360 min |
  | O2, 25 C | 8.5 % conversion | Fig. 10a: `1 - X_HP3` = 0.067 at 480 min |

  (Time-course values from the deposited files of Figs. 1, 9 and 10.) No
  explanation of the differences was found in the paper. What Section 2.3
  states about the two procedures: the screening batch is 10 mL of liquid
  with the products collected at the end; for time-resolved runs substrate,
  NaOH and water were tripled (whether the catalyst was is not stated, see
  above), samples were taken through a dip tube, and at every sample 1 mL
  was discarded and 0.5 mL sampled (1.5 mL in all). A model fitted to a
  time course therefore does not give the screening value even at its own
  condition (O2 model 99.0 against 72.8 %; H2O2 model 64.7 against 78.4 %).
- **Malonic acid alone (Fig. 4).** The two deposited files do not give the
  same MA conversion at the end:

  | | 50 C | 100 C |
  |---|---|---|
  | Fig. 4a (`4-1`, one value per temperature, "for 1300 min") | 3.8 % | 89.8 % |
  | Fig. 4b (`4-2`, time courses), 120 min | 0.8 % | 30.1 % |
  | Fig. 4b, 970 min | - | 81.0 % |
  | Fig. 4b, 1300 min | 1.2 % | 99.0 % |
  | O2 model started from 60 mM MA, 120 min | 13.7 % | 96.4 % |
  | O2 model, 1300 min | 79.7 % | 100 % |

  (Fig. 4b: 100 minus the deposited MA percentage; the 50 C series reads
  101.4 and 100.5 % MA at 1140 and 1230 min.) The published panels show the
  values of their own files, as far as the image can be read. The caption
  and both files give 1300 min; the text quotes the Fig. 4a values, 3.8 and
  89.8 %, as "after 24 h", and the 30.1 % of Fig. 4b. The model values are
  tested. The model is contradicted by the 50 C series (79.7 % against 1.2
  or 3.8 % at 1300 min; 13.7 against 0.8 % at 120 min) and by the 100 C
  series at 120 min (96.4 against 30.1 %). At 100 C after 1300 min it is
  close to the measurement (100 % against 89.8 or 99.0 %). Its steps R3, R4
  and R6 are fitted to the 3-HP run of Fig. 9b, in which the MA fraction
  falls from 0.548 at 240 min to 0.095 at 1680 min; the model has MA react
  as fast alone. It should not be used for the stability of MA on its own.
- **Catalyst reuse (Fig. 6b).** Three cycles of a "catalyst reusability
  test" at NaOH 1 eq., 3 bar of O2, 50 C and 8 h (caption); the panel shows
  the normalized 3-HP and MA concentrations of each cycle up to 120 min. In
  the deposit (`7_reusability.csv`, not shipped; columns `3HP`, `MA`, `OXA`
  and `AAc` for each cycle, starting at 1, 0, 0 and 0, with no statement of
  the normalization) the row of 120 min reads, as `3HP` / `MA`, 0.25526 /
  0.47639 (`1st`), 0.65162 / 0.09222 (`2nd`) and 0.76891 / 0.00162 (`3rd`).
  The O2 model has `X_HP3` / `X_MA` = 0.325 / 0.471 at 120 min (tested). The
  `3HP` and `MA` columns of the first cycle are the experimental points of
  Fig. 9b, the experiment the O2 model was fitted to (same values at 20 to
  120 min); the first-cycle `OXA` and `AAc` columns are not those of
  Fig. 9b (0.00254 and 0.0117 against 0.00229 and 0.00566 at 20 min). In
  the second and third cycle far more 3-HP is left than in the first cycle
  or in the model, and little MA is formed; the paper speaks of "the
  activity loss observed upon catalyst reuse" and relates it to surface
  oxidation of Pd and slight Pd leaching. The models have no term for this
  loss and do not describe a reused catalyst. The last labelled row of the
  file, `480`, is not in the published panel; in the first-cycle `3HP` and
  `MA` columns it holds the Fig. 9b values of 180 min, followed by those of
  240 to 1470 min in rows without a time, and it is not used here.

In short: the models describe the time courses they were fitted to. `P`,
`eq_H2O2_0` and `T` are in the models because the rate laws contain them,
not because the response to them was established.

What does and does not depend on `C0`, measured with the package:

- O2 model: the plotted fractions `X_i` hardly depend on `C0`, because
  every rate is first order in its liquid species. Halving or doubling `C0`
  changes them by less than 1e-4 (through the O2 depletion, which is
  proportional to `C0`).
- H2O2 model: the initial H2O2 concentration is `eq_H2O2_0*C0`, so changing
  `C0` alone changes the H2O2 concentration and with it the rates of R1, R2
  and R4 (`X_MA` at 600 min is 0.201, 0.397 and 0.564 for `C0` = 0.03, 0.06
  and 0.12 M at 15 equivalents).
- Reversible variant: the plotted fractions do not depend on `C0`
  (`k_1r` is defined per unit of `[W]/C0`).
- Simplified model: `x_per_M` holds the script's initial 3-HP mole fraction
  divided by `C0`; if `C0` is changed, `x_per_M` has to be changed with it.

`set_value` does not validate its input. Zero or negative `m_cat`, `C0` or
`eq_H2O2_0`, or a `T` set in degrees Celsius, give meaningless results
without an error. Only the loader arguments are checked (`T` must lie
within 250 to 450 K; unknown constant-set names are refused).

## SI Fig. S2 and Eq. S1 (rate expressions)

Section 3.4 says that "the rate expressions and optimization algorithm are
detailed in Figure S2". Fig. S2 is a flow chart ("Optimization algorithm for
kinetic modeling of 3-HP oxidation"); its caption introduces Eq. S1. What
they state (all printed):

- **Eq. S1 (mole balance):** `dN_i/dt = N dx_i/dt = sum_j v_ij r_j`, with
  `N` the total moles in the reactor, `x_i` (liquid) and `y_i` (gas) mole
  fractions and `v_ij` stoichiometric coefficients. "Mole fractions were
  used to simplify the data handling."
- **Rate form (box "Model"):** `r_j = k_jf * prod_i a_i^(n_ij)`.
- **Activities:** `a_i = gamma_i*x_i` for a liquid species and
  `a_i = phi_i*y_i*P/f0` for a gas species. Section 2.6: "Both the gas and
  condensed phases were assumed to behave ideally", i.e. `gamma = phi = 1`.
- **Bounds (box "Boundary conditions"):** `k >= 0` and `0 <= n <= 3`.
- **Irreversibility (box "Irreversible"):** `Keq >> 1`. The `Keq` are from
  DFT (Eq. 5; SI Table S4).
- **Parameters:** `k = {k_1, ..., k_6}` and `n = {n_1, ..., n_6}`, fitted
  by a global search (`particleswarm`, `ode23s`) followed by a local
  refinement (`lsqcurvefit`, `ode15s`), both minimizing
  `||x_exp(t) - x_mod(t)||^2`.
- **Goodness of fit (box "Residuals"):** `R2 = 1 - SSE/SST` with
  `SSE = sum_t ||x_exp - x_mod||^2` and `SST = sum_t (x_exp - mean(x_exp))^2`.
  `pooled_r2` in `temperature_series.py` is this definition, taken over all
  species and rows together (as the authors' posted script computes it).

How each shipped model relates to it, term by term:

| Term of Fig. S2 | H2O2 model | O2 model | Reversible O2 variant | Simplified reversible model |
|---|---|---|---|---|
| `r_j = k_jf*prod(a_i^n)`, irreversible | as Fig. S2 | as Fig. S2 | forward minus reverse terms in R1, R2, R5, R6 | forward minus reverse term in every step |
| orders `n` (fitted, 0 to 3; Table S4) | all exactly 1 (inferred) | all exactly 1 (inferred) | all 1 | all 1 |
| liquid activity `gamma*x` | `x_per_M*c`, `gamma = 1`: the authors' scale, 0.000339 for the 60 mM feed | the same | not separable: constants are lumped per unit of `c/C0` | the script's mole fractions `x_per_M*c` (0.000339 scale in its fits, unit scale in its heat-map run) |
| oxidant activity | H2O2 is a liquid species: `x_per_M*[H2O2]`, 0.005085 at t = 0. As Fig. S2 | **departs**: Fig. S2 gives `y*P/f0` = 3/1.01325 = 2.96 for 3 bar of O2; the model has `x_O2_ref*yO2` = 0.005085*`yO2` | none (lumped into the constants) | `yO2*P/f0` with `P` = 3.5 bar, `f0` = 1.01325 bar. As Fig. S2, with the script's pressure |
| `k_jf` | Table 3, times `m_cat` (not `k_7`) | Table 3, times `m_cat` | effective constants, 1/min | effective constants, 1/min |
| `Keq >> 1` | irreversible, as Fig. S2 | irreversible, as Fig. S2 | does not hold: four steps reversible | does not hold: "equilibrium constants" `K_6`..`K_10` of 0.6 to 82 |
| Eq. S1 | liquid species as concentrations, `c = C0*n` | the same; `yO2` loses one unit per unit of liquid mole fraction reacted in R1, R2, R4 | liquid species as concentrations | the script's states; `yH2' = r_3`, `yCO2' = r_2 + r_3` |

So the two Table 3 models **are** the Fig. S2 rate form with `n = 1`, unit
activity coefficients and the liquid mole fractions on the authors' scale.
They depart from it in two places:

1. **The O2 activity of the O2 model.** With the literal gas activity of
   Fig. S2, `P/f0` = 2.96, the Table 3 constants (times 0.2 g) make R1, R2
   and R4 582 times faster than with `x_O2_ref`: 3-HP would have a half-life
   of 0.16 min instead of 73 min, and the model would miss the deposited
   Fig. 9b curve by up to 0.97 (rms 0.42, 0.41, 0.02, 0.87 and 0.50 for
   3-HP, MA, AAc, OXA and D; computed with an independent scipy integration
   of the Table 3 rate laws). The deposited curve is reproduced with the
   constant 0.005085 in place of `P/f0` (next section). The authors' posted
   O2 script does not apply `P/f0` either: it defines `P = 3.5` and
   `f0 = 1.01325` and uses neither in its rates, which are `k*x*yO2` with
   `yO2 = 1.0`.
2. **The mole-fraction scale.** The authors' scripts map the 60 mM feed to a
   mole fraction of 0.000339. The mole fraction of 60 mM in water is about
   1.08e-3. For the first-order organic terms the scale cancels; it sets the
   size of the H2O2 factor in R1, R2 and R4 of the H2O2 model. It is taken
   as used.

The orders `n` are discussed under "Table S4 is not used" below.

## Where each number comes from

| Quantity | Value | Status |
|---|---|---|
| `k_1`..`k_7` (H2O2), `k_1_ref`..`k_6_ref` (O2) | Table 3 as printed | printed |
| `m_cat` | 0.2 g | printed (Fig. 9 caption) |
| `C0`, `eq_H2O2_0` | 0.06 M, 15 | printed (Fig. 9 caption) |
| `k_7` not multiplied by `m_cat` | - | printed (Table 3 footnote b) |
| rate form `r = k*prod(a^n)`, liquid activity `x` | - | printed (SI Fig. S2) |
| irreversible steps | - | printed (Section 3.4: "All DFT-determined Keq were found to be high enough to consider all reactions as irreversible"; SI Fig. S2) |
| reaction orders exactly 1 | - | inferred from the deposited curves and the posted script; the paper says orders were fitted and "found to be first order for most steps" (Table S4) |
| `k_1`..`k_6` multiplied by `m_cat` | - | inferred: the unit is 1/(g_cat min), the starting vector of the posted O2 script (1.43, 27.8, 6.32e-5, 1.17e-2, 2.25e-3, 1.1e-3) is the Table 3 O2 column times 0.2 to within 0.5 %, and the deposited curves are reproduced that way |
| `x_per_M` | 0.000339/0.06 per M | inferred: the authors' scripts map the 60 mM feed to a mole fraction of 0.000339 (see the previous section) |
| `x_O2_ref` | 15*0.000339 = 0.005085 | **inferred**; not in the paper, the SI, the README or the posted script. See below |
| `yO2` as a depleting state, coefficient 1 for R1, R2, R4 | - | inferred from the deposited curve (the posted script integrates the state but sets yO2 = 1.0 in its rates) |
| how H2O2 enters R1, R2, R4 (H2O2 model) | - | inferred from the deposited curve: the authors posted only the O2 version of their script |
| `Ea_1` = `Ea_2`, `Ea_3`, `Ea_6` | 59781.2, 91594.6, 59549.8 J/mol | printed to three digits in Section 3.4 (59.8, 91.6, 59.5 kJ/mol); full values from the authors' heat-map script. The text says they were obtained from the Arrhenius plot of Fig. 10d, whose lines give 59.93, 68.96 and 52.09 kJ/mol (section "Arrhenius plot (Fig. 10d)") |
| `Ea_4` | 10000 J/mol | heat-map script only (a round value there); not in the paper |
| `Ea_5` | 51788.73 J/mol | heat-map script only; not in the paper |
| `T_ref`, `R_gas` | 323.15 K, 8.314 J/(mol K) | heat-map script; 50 C reference is printed |
| Arrhenius form on all six steps of the O2 model | - | **NSKinetics addition** (the Arrhenius relation of the paper). See "Temperature dependence of the O2 model" |

### `x_O2_ref`

The Table 3 O2 entries for R1, R2 and R4 reproduce the deposited Fig. 9b model
curve only when multiplied by `m_cat` and by a further factor of about 0.0051,
taken here as 0.005085 because that is the H2O2 mole fraction of the
15-equivalent H2O2 run on the authors' scale (15*0.000339). How the factor
entered the authors' code is unknown. Only the products `k_i*m_cat*x_O2_ref`
are identifiable.

With this value the oxidant factor of the O2 model equals that of the H2O2
model at t = 0. If the factor is read as an O2 mole fraction, the O2 model
is the hypothetical condition `x_H2O2 = x_O2` that the paper discusses in
Section 3.4, and for that condition the paper states a rate ratio that the
models do not give (item 8 of "Paper versus deposit" below). This is one
more reason why the factor is recorded as inferred from the deposit only.

What fixes the value, and how firmly:

- Measured with the package (accurate solver): with 0.005085 the model is
  within 2.0e-3 of the 19 deposited nodes after t = 0 (rms 6e-4). That is the
  solver error of the deposited curve, so this comparison cannot distinguish
  0.005085 from nearby values: 0.0051 gives 1.4e-3 (rms 5e-4), 0.00508 gives
  2.3e-3. It does exclude the alternatives of the previous section: no factor
  (the posted script as posted: yO2 = 1.0 and its printed starting constants)
  misses the curve by an rms of 0.42, 0.44, 0.76 and 0.46 for 3-HP, MA, OXA
  and D, and `P/f0` by as much.
- With the authors' solver run repeated (`solver_emulation.py`; section
  "Emulation of the authors' solver run"): the constants with which the
  repeated run returns the deposited nodes to their printed digits
  (`O2_NODE_PARAMETER_SETS`; a fit that is not shipped), divided by Table 3
  as printed, give a factor of 0.0050854 (R1), 0.0050857 (R4) and 0.0050733
  (R2). A printed three-digit constant stands for half a unit of its last
  digit either way, so the three together allow 0.0050818 to 0.0050889: that
  contains 0.005085 and not 0.0051. With 0.0051 the repeated run misses the
  nodes by up to 9.4e-4 (node times by 1.8 min) with the constants of the
  node set, and still by 7.0e-4 (1.4 min) when the constants of R1, R2 and R4
  are set to the lowest values Table 3 can stand for (7.135, 138.5 and
  0.05855) and the other three keep their node-set values (5.8e-4 and 1.5 min
  with those three as printed in Table 3); the limit is 1e-5. The package
  tests compute the window and the three misses. The choice of 0.005085 over
  a rounded 0.0051 rests on this. It shows that 0.0051 is excluded if Table 3
  is the rounded set of the constants behind the deposited curve; any value
  in the window would do, and 0.005085 is taken because it is 15*0.000339.

### Table S4 is not used

SI Table S4 lists reaction orders (H2O2: 0.99, 0.99, 1.49, 1.14, 1.02, 1.01
and 1 for R7; O2: 0.99, 0.95, 0.78, 0.87, 0.98, 1.01) and DFT equilibrium
constants. Neither enters the rate laws.

**Equilibrium constants.** That the steps are irreversible is printed
(Section 3.4, SI Fig. S2); the `Keq` of Table S4 (1e4 to 1e92) are the
paper's reason for it and have no other role.

**Orders.** That the orders are exactly 1 is inferred, from two things: the
posted script has first-order rates and no orders, and the deposited curves
are reproduced with first-order rates (to 2.0e-3 for O2 and 9.6e-4 for
H2O2 when solved accurately, which is their solver error, and to their
printed digits when the solver run is repeated) and not with the Table S4
orders. How far the Table S4 orders move a curve depends on the scale of
the activity, because `x^n` with `n != 1` is not scale-free. The table
gives, for the Table S4 orders applied to the organic reactant of each step
(oxidant factor unchanged), the largest difference from the deposited nodes
and the pooled R2 against the experimental points. Its first row is
measured with the package. The other rows were computed outside the package
(not shipped) with an independent scipy integration (LSODA, relative
tolerance 1e-10) in which a species with `y` mol per mol of initial 3-HP
has the activity `x = s*y`, so that `k*x^n = k*s^(n-1)*y^n` in the units of
the first-order model:

| Basis (`s`, the activity of the feed) | O2: max. difference from the deposited curve | O2: R2, 16 rows | H2O2: max. difference from the deposited curve | H2O2: R2, 13 rows / 11 rows |
|---|---|---|---|---|
| first order (the shipped models) | 0.0020 | 0.9945 | 0.0010 | 0.9412 / 0.9452 |
| fraction of the feed (`s` = 1) | 0.010 | 0.9943 | 0.037 | 0.9413 / 0.9442 |
| authors' mole fraction (`s` = 0.000339), the scale of Fig. S2 | 0.17 | 0.943 | 0.14 | 0.925 / 0.922 |
| mole fraction in water (`s` = 1.08e-3) | 0.13 | 0.963 | 0.13 | 0.927 / 0.926 |
| M (`s` = 0.06) | 0.047 | 0.991 | 0.096 | 0.936 / 0.937 |
| mM (`s` = 60) | 0.023 | 0.993 | 0.16 | 0.906 / 0.918 |

So the Table S4 orders move the curves off the deposited ones by 0.010 to
0.17 (O2) and 0.037 to 0.16 (H2O2), depending on the basis. On the
feed-normalized basis they fit the experiments about as well as first order
does (R2 0.994 against 0.9945 for O2, 0.941 for both with H2O2): the
experiments do not discriminate. It is the deposited curves that exclude
the Table S4 orders: even the smallest difference (0.010) is five times the
solver error of the deposited O2 curve.

Table S4 also writes 0.5, 0.5 and 1.5 O2 for R1, R2 and R4, whereas the O2
depletion in the deposited run is 1 for each; with `yO2` falling only to
0.9995 this is immaterial (dropping the depletion changes the plotted values
by at most 7e-5). The values of the table are kept as `TABLE_S4` in
`pd_c_3hp_oxidation.py` for documentation, and a test pins them against the
SI.

## Agreement with the deposit

Measured with these models (CVODE, relative tolerance 1e-9), evaluated at the
deposited node or sampling times; values are normalized concentrations.

| Comparison | Result |
|---|---|
| O2 model vs deposited Fig. 9b curve (19 nodes after t = 0) | rms 9.4e-4 (3-HP), 1.0e-3 (MA), 2.6e-5 (AAc), 2.1e-5 (OXA), 1.7e-4 (D); max 2.0e-3 |
| H2O2 model vs deposited Fig. 9a curve (23 nodes after t = 0) | rms at most 3.3e-4 per species; max 9.6e-4 |
| O2 model `x_OPA` vs digitized SI Fig. S12b inset (13 vertices) | max 3.0e-7 on a peak of 1.4e-5 |
| O2 model vs experiment, pooled R2 | 0.9945 on all 16 rows; 0.9957 without the 120 min row (SI Fig. S11a prints 0.997) |
| H2O2 model vs experiment, pooled R2 | 0.9452 on the 11 rows plotted in SI Fig. S11b, which prints 0.945; 0.9412 on all 13 rows |
| H2O2 model `eq_H2O2` vs deposited H2O2 trace | up to 0.38 equivalents above it (rms 0.22) |

The deposited curves have 20 (O2) and 24 (H2O2) nodes, i.e. 19 and 23 solver
steps: node counts in these notes include the node at t = 0, and a curve of n
nodes has n - 1 steps. The nodes are unevenly spaced, the pattern of an
adaptive solver's own steps, and Section 2.6 of the paper names `ode23s`. For
O2 the nodes are up to 2.0e-3 away from the converged solution, and the
straight segments drawn between them in the published figure up to 0.010. The
residual differences in the table have two sources:

- the solver error of the deposited curves. The repeat of the solver run
  (next section), with constants fitted to the nodes instead of the
  three-digit ones of Table 3, returns the deposited values to their five
  printed decimals (max 5.0e-6 for O2, 6.8e-6 for H2O2);
- the rounding of Table 3. With Table 3 as printed the repeated run comes
  within 4.0e-4 (O2) and 4.4e-4 (H2O2) of the nodes. For O2 the fitted
  constants agree with Table 3 within its rounding: each prints as Table 3
  prints it. For H2O2 they differ from Table 3 by up to 0.6 %, more than
  rounding allows for five of them (item 3 below).

## Emulation of the authors' solver run

The deposited model curves are not converged solutions: their nodes are
spaced like the steps of an adaptive solver, so an accurate solution can
agree with them only to their solver error (previous section; up to 2.7e-3
for Fig. 10). To compare more closely, `solver_emulation.py` repeats the
run:

- `ode23s` is a port of MATLAB's `ode23s`, the modified Rosenbrock (2,3)
  pair of Shampine and Reichelt (SIAM J. Sci. Comput. 18 (1997) 1-22), with
  the default options the posted script uses (`RelTol` 1e-3, `AbsTol` 1e-6,
  largest step a tenth of the time span, one output per step) and a
  finite-difference Jacobian. It was written from the paper and the
  documented behaviour of the MATLAB function. MATLAB itself was not run.
- `emulate_authors_run(model, t_end)` integrates the state vector of the
  posted script: the liquid mole fractions on the authors' scale (0.000339
  for the feed), then H2O, H2, O2 and CO2, from the feed to 1.1 times the
  last sampling time. Every rate of change is evaluated by the loaded model
  at the state the solver asks for, so the rate laws, constants and
  stoichiometry (and the rate rule of `yO2`) are those of the shipped
  Antimony file, and the plotted fractions are read from the model's
  assignment rules. States the model does not have follow the
  stoichiometric matrix of the posted script. They enter no rate and do not
  change a step (tested: the same steps to 1e-12 without them). For the
  H2O2 run, whose script was not posted, the states are the six liquid
  species and H2O2.
- `emulate_fig9a`, `emulate_fig9b` and `emulate_fig10abc` run it for the
  deposited curves.

Result, as measured by the package tests
(`test_pd_c_3hp_oxidation_solver_emulation.py`). "Values identical" counts
the values that equal the deposited ones after rounding as the file prints
them; it was counted with the package and is not pinned by a test, because
a value on a rounding boundary can fall either way on another machine.

| Deposited curve | Constant set | Steps (deposited) | Node times, max. difference | Values, rms / max. difference | Values identical |
|---|---|---|---|---|---|
| Fig. 9a, five species | `'fig9a_nodes'` | 23 (23) | 8.1e-6 min | 3.0e-6 / 6.8e-6 | 109 of 120 |
| Fig. 9b | `'fig9b_nodes'` | 19 (19) | 6.2e-6 min | 2.4e-6 / 5.0e-6 | 100 of 100 |
| H2O2 trace of SI Fig. S12a | `'figS12a_h2o2_nodes'` | 24 (23: the file ends one step early, at 1537 min) | 5.8e-6 min | 2.6e-6 / 4.4e-6 equivalents | 24 of 24 |
| Fig. 10a (25 C) | `'25C'` | 18 (18) | 0.0047 min | 2.1e-6 / 7.2e-6 | 91 of 95 |
| Fig. 10b (75 C) | `'75C'` | 19 (19) | 0.0054 min | 2.5e-6 / 6.3e-6 | 98 of 100 |
| Fig. 10c (100 C) | `'100C'` | 24 (24) | 0.0052 min | 2.3e-6 / 5.7e-6 | 113 of 125 |

The Fig. 9 file prints its node times to five decimals and the Fig. 10
files to two; all files print values to five decimals (small values to
three significant digits or six decimals). The differences are therefore at
the print rounding of the files (5e-6, or 0.005 min for Fig. 10). The tests
require the deposited number of steps, node times within 2e-5 min (Fig. 9)
or 0.01 min (Fig. 10) and values within 1e-5. For comparison, the same
models solved accurately differ from these emulated steps by 8.5e-4
(Fig. 9a), 2.0e-3 (Fig. 9b), 0.012 equivalents (H2O2 trace), 1.8e-4, 2.5e-3
and 2.7e-3 (Fig. 10a-c). That is the solver error of the deposited
curves, which the accurate-solver tests have to allow for; for the plotted
fractions the limit of 1e-5 is 18 to 270 times smaller.

The constant sets of the Fig. 9 nodes, in the units of Table 3 (`k_1` to
`k_6` in 1/(g_cat min), `k_7` in 1/min), with their ratio to Table 3 as
printed. All are **inferred** from the deposited nodes (fits that are not
shipped) and are selected with `parameter_set` in `load_o2_model` and
`load_h2o2_model`; the default of both loaders and of the Antimony files is
Table 3 as printed.

| | `'fig9b_nodes'` (O2) | `'fig9a_nodes'` (H2O2) | `'figS12a_h2o2_nodes'` (H2O2) |
|---|---|---|---|
| `k_1` | 7.1405267 (1.00007) | 5.3228596 (1.00242) | 5.5412148 (1.04354) |
| `k_2` | 138.68010 (0.99770) | 30.913113 (0.99399) | 32.725866 (1.05228) |
| `k_3` | 3.1608819e-4 (1.00028) | 7.3582251e-4 (1.00522) | 6.7028262e-4 (0.91569) |
| `k_4` | 0.058607778 (1.00013) | 1.1471899 (0.99756) | 1.1053986 (0.96122) |
| `k_5` | 0.011270915 (0.99743) | 0.0030247652 (0.99827) | 0.0023665310 (0.78103) |
| `k_6` | 0.0055212537 (1.00023) | 4.10e-9 (as printed) | 4.10e-9 (as printed) |
| `k_7` | - | 0.0064883315 (1.00129) | 0.0069365260 (1.07045) |

The O2 set is for `m_cat` = 0.2 g and `x_O2_ref` = 0.005085 (only the
products are determined). Each of its constants prints, to three digits, as
Table 3 prints it. Of the `'fig9a_nodes'` set only `k_4` does (item 3 of
the next section).

With **Table 3 as printed** the repeated run takes the deposited number of
steps and comes within 4.0e-4 of the Fig. 9b nodes (rms 1.4e-4; node times
within 0.041 min) and within 4.4e-4 of the Fig. 9a nodes (rms 2.0e-4; node
times within 0.47 min). This is what the rounding of Table 3 to three
digits (and, for H2O2, the deviations of item 3 below) costs.

Details of the solver that the deposited time grids fix (all tested):

- the tolerances: with a relative tolerance of 5e-4 or 2e-3 the node times
  of Fig. 9b are more than 20 min off, and with an absolute tolerance of
  5e-7 or 2e-6 the run takes 21 or 17 steps instead of 19;
- the first step, `1/h = 1.25*sqrt(0.5*max(|J f|/w))/RelTol^(1/3)`: it is
  the first deposited node of every curve but one (5.87666 min for Fig. 9b,
  12.547 min for Fig. 9a);
- the rule for a rejected step: at 25 C the first step (3.83 min) is
  rejected and retried with a tenth of it, and the first deposited node is
  at 0.38 min. It is the only rejected step of the six runs. No step of any
  run is accepted or rejected narrowly (the error estimate of every
  accepted step is below 0.6 of the tolerance, that of the rejected one
  1.4e4 times it), so rounding differences between machines cannot change
  the sequence of steps;
- the last step, which ends at 1.1 times the last sampling time.

Three rules of the solver are tested by these curves in part or not at all,
and the package tests check them on other problems:

- the rule for a rejected step, on `y' = cos(10 t)`. The rejected step of
  the 25 C run is a first retry, at the lower limit of a tenth; a retry by
  another factor and a second rejection do not occur in the six runs;
- the limit of five on the growth of a step from one to the next, on
  `y' = -100 y`. In the six runs an error estimate small enough for it
  occurs only in the last step and in steps that already have the largest
  size, a tenth of the run, so that the limit does not show. The largest
  growth of a step in the six runs is 4.17-fold (25 C, from the fifth to
  the sixth step);
- the stretching of a step that would end within 10 % of its length before
  the end of the run, on `y' = -y`. The last step of the 100 C run is
  stretched in this way (181.17 min, where the largest step size is
  171.6 min); the last step of each of the other five runs is cut short at
  the end of the run instead.

The 13 vertices of the OPA curve in the SI Fig. S12b inset (digitized)
agree with the first 13 nodes of the emulated Fig. 9b run, t = 0 and its
first 12 steps (`x_OPA` within 1.9e-7 on a peak of 1.4e-5; tested), so the
inset is consistent with the same solver run.

What this establishes, and what it does not:

- It shows that the shipped rate laws, with the constant sets named in the
  table, return the deposited nodes to their printed digits. It does not
  show that these are the only constants or the only rate laws that do.
  The fits that produced the constant sets are not shipped.
- The tests measure how tightly the nodes pin each constant on its own:
  the relative change that makes the run miss a limit at least twofold
  (the smallest power of ten).

  | | `k_1` | `k_2` | `k_3` | `k_4` | `k_5` | `k_6` | `k_7` |
  |---|---|---|---|---|---|---|---|
  | `'fig9b_nodes'` | 1e-6 | 1e-6 | 1e-3 | 1e-3 | 1e-6 | 1e-4 | - |
  | `'fig9a_nodes'` | 1e-6 | 1e-6 | 1e-3 | 1e-4 | 1e-5 | not detected | 1e-6 |
  | `'figS12a_h2o2_nodes'` | 1e-6 | 1e-6 | 1e-3 | 1e-4 | 1e-4 | not detected | 1e-6 |

  ("Not detected": a change of 10 % leaves the run within the limits;
  `k_6` of the H2O2 model is not identifiable.) The sets are therefore
  given to eight significant digits. The table for the reversible variant
  is in its section. These are changes of one constant at a time; several
  constants changed together can compensate each other in part.
- Each negative control of the tests changes one inferred choice and keeps
  the constants: `x_O2_ref` (section "`x_O2_ref`"), the O2 state, `k_7`
  and the catalyst mass, the normalization, and for the reversible variant
  the water stoichiometry and the reverse of reaction 1 (its section). A
  control of this kind shows that the shipped constants need the shipped
  choice. Whether another choice could do as well with other constants
  takes a refit; where one was made, it is not shipped, and the text says
  so.

The controls for the two Table 3 models (limit 1e-5; node times 2e-5 min):

| Change | Effect on the repeated run |
|---|---|
| `x_O2_ref` = 0.0051 (Fig. 9b) | misses the nodes by 9.4e-4, node times by 1.8 min; see "`x_O2_ref`" |
| O2 fraction held at 1 in the rates, as in the posted script (Fig. 9b) | 7.7e-5, node times 0.15 min |
| `k_7` multiplied by `m_cat` (Fig. 9a) | 16 steps instead of 23; 0.39 |
| OPA included in the sum of the normalization | same steps; 0.038 (Fig. 9b), 0.094 (Fig. 9a), 0.029, 0.098 and 0.157 (Fig. 10a-c) |
| H2O2 trace run with the `'fig9a_nodes'` set or with Table 3 | 23 steps instead of 24, node times more than 50 min off |

On the O2 fraction: the posted script integrates an O2 state and uses
`yO2 = 1.0` in its rates. The control shows that the constants of the node
set need the depleting state. In a refit of the six constants with the
fraction held at 1, on the deposited time grid (not shipped), the largest
difference from the nodes stayed at 1.5e-5, three times that of the model
with the depleting state (5e-6). The evidence for the coupling is therefore
this factor of three in a fit, not the control alone.

`figures.plot_solver_emulation` draws the six deposited curves with the
emulated steps on their nodes and, below each, the differences of the
emulated steps and of the accurate solution from the nodes.

## Paper versus deposit: discrepancies that affect these models

1. **Deposit file `10-1_kinetic_model_H2O2.csv`.** The block labelled `exp`
   holds the model nodes and the block labelled `model` holds the experimental
   points. The experimental t = 0 row is printed as zeros (it is the feed,
   1, 0, 0, 0, 0). The t = 0 cell of the H2O2 trace reads 1 (it is 15). The
   readers in `reference.py` correct all three; the shipped CSV is verbatim.
2. **H2O2 trace of SI Fig. S12a.** The deposited trace and the published
   inset agree with each other and not with the H2O2 of the Table 3 model,
   which lies up to 0.38 equivalents (rms 0.22) above them. The trace is
   reproduced as the H2O2 of another run of the same model: with the
   constant set `'figS12a_h2o2_nodes'` (from a fit to the trace that is not
   shipped) the repeated solver run returns its 24 nodes to their printed
   digits (previous section). That set has `k_7` = 6.94e-3 1/min instead of
   6.48e-3, a `k_5` 22 % lower and `k_1` to `k_4` 4 to 8 % off Table 3; the
   species curves it gives differ from the deposited Fig. 9a curve by up to
   0.05. So the trace is consistent with a parameter set other than the one
   behind the Fig. 9a curve and Table 3; how the authors computed it is not
   stated. The species curves in the same figure are the Fig. 9a curve,
   which is reproduced with constants within 0.6 % of Table 3 (the
   `'fig9a_nodes'` set; item 3).
3. **Table 3, H2O2 column, versus the Fig. 9a curve.** The constants with
   which the repeated solver run returns the deposited curve
   (`'fig9a_nodes'`; previous section; the fit itself is not shipped)
   differ from the printed ones by 0.1 to 0.6 %: +0.24 % (`k_1`), -0.60 %
   (`k_2`), +0.52 % (`k_3`), -0.24 % (`k_4`), -0.17 % (`k_5`) and +0.13 %
   (`k_7`). That is more than rounding to three digits allows for `k_1`,
   `k_2`, `k_3` and `k_7` (at most 0.09, 0.16, 0.07 and 0.08 %), marginally
   more for `k_5` (0.173 against 0.165 %), and within it only for `k_4`
   (0.43 %). Of the six, only `k_4` prints as Table 3 prints it (tested).
   The effect on the curve is 4.4e-4 (previous section). `k_6` (4.10e-9) is
   not identifiable: any value below about 1e-6 1/min gives the same curves.
4. **"H2O2 is fully consumed after 1000 min"** (Section 3.4). In the Table 3
   model H2O2 falls below 1 % of its charge at 671 min, and 0.016 equivalents
   remain at 1000 min.
5. **OPA.** The text gives the predicted OPA peak with O2 as 8.4e-4 mM near
   25 min and calls it negligible. The plotted quantity of SI Fig. S12b is a
   mole fraction that peaks at 1.4e-5, and 1.4e-5 times 60 mM is 8.4e-4 mM,
   which suggests that the printed value is that mole fraction multiplied by
   the feed concentration. As a concentration, the OPA peak of the model is
   2.5 mM, 4.2 % of the initial 3-HP, at 20.5 min; 25 min is a vertex of the
   plotted polyline. The H2O2 model predicts a larger undetected OPA pool:
   11.5 % of the initial 3-HP at 84 min and 7.8 % at 1440 min.
6. **Parity R2 (SI Fig. S11).** The text states "parity plots yielding
   R2 > 0.95 (Figure S11a and b)".
   - H2O2: SI Fig. S11b prints R2 = 0.945, which disagrees with the text.
     It plots 11 sampling times: all but t = 20 and 1205 min, whose 3-HP
     points would lie at (0.914, 0.962) and (0.402, 0.240) and are not in
     the figure. The H2O2 model reproduces the printed value on those 11
     rows (0.9452); on all 13 rows it gives 0.9412.
   - O2: SI Fig. S11a prints 0.997. The O2 model gives 0.9945 on all 16
     rows and 0.9957 without the 120 min row (whose 3-HP point is not in
     S11a either, as far as the image can be read). The origin of 0.997 is
     not identified: leaving out any one or two rows gives at most 0.9965,
     and the squared correlation coefficient is 0.9947 (0.9959 without the
     120 min row). Nothing is tested against it.
7. **The 120 min O2 sample** is plotted in Fig. 9b and listed in the
   deposit. A least-squares refit of the six constants to the experiment
   (not shipped) returns the published constants within 0.05 % when that
   row is left out, and not otherwise.
8. **Initial-rate ratio.** The paper's first Eq. 6 (Section 3.4; see the
   next section on the numbering) is
   `r0(H2O2)/r0(O2) = k_f(H2O2)*x_3HP*x_H2O2/(k_f(O2)*x_3HP*P_O2)` for
   reaction 1, and the text gives 0.204 "under the previously optimized
   conditions (3 bar O2 and 15 equimolar H2O2)", i.e. O2 about five times
   faster. This is not reproduced. In the shipped models the initial rates
   of reaction 1 are 3.24e-4 M/min (H2O2) and 4.36e-4 M/min (O2): both
   oxidant factors equal 0.005085 at t = 0, so the ratio is that of the
   Table 3 constants, 5.31/7.14 = 0.744. With the Table 3 constants, that
   equation gives 0.204 for `x_H2O2/P_O2` = 0.274. The quantities stated in
   the paper give other rate ratios: the H2O2 concentration of 0.9 M over
   3 bar (0.3 in place of `x_H2O2/P_O2`) gives 0.223, and the mole fraction
   of the models, 0.005085, over 3 bar or over `P/f0` = 2.96 gives 0.0013.
   How 0.204 was obtained is not known.

   The text goes on to the case of equal oxidant mole fractions, with the O2
   pressure from Henry's law: "Under the hypothetical condition where
   xH2O2 = xO2, the equivalent O2 pressure should be approximately 1000 bar,
   implying that r0O2 would be ~1300 times faster than r0H2O2". This is not
   reproduced either. At t = 0 the two shipped models have the same oxidant
   factor, 0.005085, which is that condition if `x_O2_ref` is read as a mole
   fraction, and their O2 rate is only 7.14/5.31 = 1.34 times their H2O2 rate
   (tested). In the paper's Eq. 6, `k_f(O2)` multiplies a pressure and not a
   mole fraction; scaled linearly from 3 to 1000 bar, the 0.204 would become
   a factor of 1000/(3*0.204) = 1630 in favour of O2 (our arithmetic; the
   paper gives ~1300 without a derivation). In the O2 model the Table 3
   constants multiply the factor `x_O2_ref`, which has the size of a mole
   fraction, because that is what reproduces the deposited curve (section
   "`x_O2_ref`").
9. **Solver.** Section 2.6 names `ode23s` and `ode15s`; the deposited curves
   are reproduced step for step by an emulation of `ode23s` at default
   tolerances (previous section; MATLAB itself was not run).
10. **Maximum MA yield at 75 C.** The text gives maximum MA yields of 52.6 %
    (75 C) and 56.8 % (100 C). The largest measured MA fractions in the
    deposit are 0.436 (75 C, 180 min) and 0.568 (100 C, 40 min); the 75 C
    value of the text is not in the deposited data. (At 25 C the text's
    2.3 % after 31.5 h is the deposited 0.023 at 1890 min.)
11. **Time axis of Fig. 10a-c.** The published panels end at 1440 min. The
    samples at 1740 and 1890 min (25 C) and 1560 min (100 C) are in the
    deposit but not in the figure. They are counted in the pooled R2 of the
    comparison with the O2 model below, whose figure
    (`plot_fig10_table3_o2`) ends at 1440 min too and does not draw them.
12. **Activation energies.** Two of the three activation energies of the
    text are not the slopes of the lines of Fig. 10d, from which the text
    says they were obtained (section "Arrhenius plot (Fig. 10d)").

## Temperature dependence of the O2 model (NSKinetics addition)

`k_i(T) = k_i_ref*exp(-(Ea_i/R_gas)*(1/T - 1/T_ref))` for all six steps, with
`T_ref` = 323.15 K, so that the default model is exactly the published 50 C
model. This is the Arrhenius relation of the paper.

**Equation numbering.** The paper numbers two displayed equations (6). The
first, in Section 3.4 (page 7), is the initial-rate ratio of item 8 above.
The second, printed at the end of Section 3.4 (page 9, before the
Conclusions), is the Arrhenius relation
`k(T) = k_ref*exp[-(Ea/R)*(1/T - 1/T_ref)]`; the text's "Arrhenius
relationship (Eq. 6)" refers to this one. (Fig. 6 of the paper is unrelated:
XPS profiles and a reusability test.) Wherever this package says "the
Arrhenius relation of the paper", the second Eq. 6 is meant.

The activation energies are those of the authors' heat-map script, which was
written for their simplified five-reaction model (3-HP -> MA in one step);
the "MA production" value is applied here to both R1 and R2. Three of them
are quoted in the text, which says that they were obtained from the
Arrhenius plot of Fig. 10d; the lines of that plot give other values for two
of the three (section "Arrhenius plot (Fig. 10d)").

None of the paper's own temperature results is reproduced by this
six-reaction model. The Fig. 11 heat maps are those of the simplified
reversible model of the deposited script
(`pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt`, next section),
and the points of the Arrhenius plot (Fig. 10d) are taken to be constants of
that model too. The model curves at 25, 75 and 100 C (Fig. 10a-c) are
reproduced by a reversible six-reaction variant with a separate constant set
for each curve (`pd_c_3hp_oxidation_o2_reversible_antimony.txt`; section
"Reversible O2 variant" below). Predictions of this file away from 50 C are
therefore an extrapolation. They are compared with the published heat maps
in the section after next and with the experiments at 25, 75 and 100 C in
the section "Temperature-extended Table 3 O2 model versus the 25, 75 and
100 C experiments"; they reproduce neither.

## Simplified reversible O2 model (the model of the heat-map script)

### Purpose

The text presents Fig. 11 as the O2 model extended with the Arrhenius
relation ("the O2 model was extended ... The heatmap of MA selectivity was
illustrated based on the Arrhenius relationship (Eq. 6), with referencing
kinetic parameter (kref) at Tref of 50 C"). The caption says that "Fig. 9b
data was used as a reference kref of the Arrhenius equation". The caption
matches the deposited script that the deposit's README lists as the "code
for generating heatmaps" (`APCATB_D_25_08770_heatmap.m`): the script fits a
model of its own to the 50 C O2 data of Fig. 9b and extrapolates that model.
The model is not the Table 3 model: it has five reactions, no OPA and a
reverse term in every step, and it is not written out in the paper or the
SI. It is shipped as
`pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt`
(`load_o2_simplified_reversible_model`) because, run as the script runs it,
it reproduces the published heat maps, and the Table 3 model does not. The
50 C point of the Arrhenius plot (Fig. 10d) supplies three of its constants.

### Structure (transcribed from the script)

Script states: liquid "mole fractions" x of 3-HP, MA, AAc, OXA, D, and gas
"fractions" yH2, yCO2 (it also integrates an H2O and an O2 state that no
rate reads; they are left out here). With `x_H2O = 1 - sum(x)`,
`yO2 = 1 - yH2 - yCO2` and `PF = 3.5/1.01325`:

| | Rate, 1/min | Stoichiometry |
|---|---|---|
| r1 | `k_1*(x_HP3*yO2*PF - x_MA*x_H2O/K_6)` | 3-HP -> MA |
| r2 | `k_2*(x_MA - x_AAc*yCO2*PF/K_7)` | MA -> AAc, +1 yCO2 |
| r3 | `k_3*(x_MA*yO2*PF - x_OXA*yH2*yCO2*PF^2/K_8)` | MA -> OXA, +1 yH2, +1 yCO2 |
| r4 | `k_4*(x_HP3 - x_D/K_9)` | 3-HP -> D |
| r5 | `k_5*(x_MA - x_D/K_10)` | MA -> D |

The numbering is the script's, not that of Table 3. In the Antimony file the
liquid species are concentrations in M as in the other models,
`x_i = x_per_M*c_i`, the rates above are the assignment rules `r_1`..`r_5`,
each reaction runs at `r_j/x_per_M` (M/min), and `yH2' = r_3`,
`yCO2' = r_2 + r_3`. This is the script's system exactly. The package tests
compare the rates and derivatives of the model with the script's equations,
written out again in the test module. In addition, the trajectories were
compared with an independent scipy implementation of the script's equations
(not shipped; same constants; DOP853, LSODA and Radau at rtol 1e-11 or
tighter) in both conventions below, for the three constant sets at 30, 50,
84, 95 and 100 C: at the loader's CVODE tolerances (rtol 1e-9, atol 1e-10)
the liquid fractions agree to 2e-8 of the feed and the gas fractions to
2e-9. For scale, rounding `k_2` to six digits (5.19998e-5) moves the
trajectories by 1e-7 of the feed.

Temperature: `p(T) = p_ref*exp(-(Ea/R_gas)*(1/T - 1/T_ref))` for all ten
parameters, `T_ref` = 323.15 K, with `Ea` = 59781.2, 91594.6, 10000,
51788.73 and 59549.8 J/mol for steps 1 to 5.

### Three properties of the script that are in the published figure

Each was tested against the contour lines digitized from Fig. 11b by removing
it (constants of the `'fig11b'` set; the rms offset is 0.27 px with all three
in place, 1 px = 3.07 min or 0.148 C). The offsets without refitting were
measured with the package. Offsets above about 1 px are not distances: they
are the linearized estimate of `contour_offsets` (value difference over
gradient), clipped at 30 px per sample (`max_offset=30`), and are
indicative only. The values "after refitting" are from fits of the constants
to the contour lines (not shipped). The package tests repeat neither.

1. **The Arrhenius factor of step i multiplies `k_i` and also `K_(i+5)`**
   (the script's `Ea` vector is the five values repeated), so `k_i/K_(i+5)`
   does not change with temperature. Scaling only the rate constants gives
   2.9 px (1.3 px after refitting four rate constants).
2. **The heat maps start from unit mole fractions.** The script has two
   initial conditions, selected here by `x_per_M` (the loader's
   `initial_condition`):
   - `'heatmap'` (default; script function `intg2`): x of 3-HP = 1,
     `x_per_M = 1/C0`. The five liquid fractions then always add up to 1, so
     `x_H2O = 0`: the reverse term of r1 vanishes, `K_6` has no effect, and
     `yH2 = x_OXA`, `yCO2 = x_AAc + x_OXA`, so that
     `yO2 = 1 - x_AAc - 2*x_OXA` falls as by-products form (to 0.92 at
     50 C and to 0.83 at 100 C, the lowest value on the map).
   - `'fitting'` (script functions `intg` and `plotfinal`, used to fit the
     constants to the 50 C data): x of 3-HP = 0.000339,
     `x_per_M = 0.000339/0.06`. Then `x_H2O = 0.999661`, r1 has its reverse
     term, and `yO2` stays above 0.9999.

   The same constants give different curves in the two conventions: at 50 C
   the MA fraction peaks at 0.571 (heat map) against 0.549 (fitting), and the
   curves differ by up to 0.03. The 50 C row of Fig. 11 is therefore not the
   50 C curve the same script fits and plots. With the fitting convention
   the heat maps miss the published contours by 14 px with the clip (91 px
   without it; either way far off) and by 9 px after refitting.
3. **The time axis of Fig. 11 is the reaction time plus 10 min.** The script
   places the feed at its first plotted time, 10 min. Model time here is the
   reaction time; `heat_maps.py` adds the offset (`FIG11_TIME_OFFSET`).
   Starting at plotted time 0 gives 3.8 px (6.5 px without the clip; 3.2 px
   after refitting).

Two further points: the script uses a pressure of 3.5 bar (the captions give
3 bar of O2), and what Fig. 11 calls "MA selectivity" and "MA yield" are the
water-free liquid fraction of MA, `X_MA`, and `(1 - X_HP3)*X_MA` (model
quantities `fig11a_selectivity` and `fig11b_yield`).

### Constants at 50 C

The fitted constants are not in the deposit. The script only prints the
vector it starts its fit from, to four decimals (`'k0_printed'`:
0.0021, 0.0001, 0.0000, 0.0023, 0.0010, 4.7840, 0.5918, 3.5116, 81.6870,
52.0472), and that vector does not reproduce the figure (4.6 and 7.1 px on
panels a and b). The two sets shipped were obtained as follows. Rate
constants in 1/min; they are lumped (no catalyst mass).

| Constant | `'fig11b'` (default) | `'fig11a'` | Source |
|---|---|---|---|
| `k_1` | 2.139999e-3 | same | exp(-6.14695): 50 C row of the deposited `11-4_Arrhenius.csv`, "ln k1" |
| `k_2` | 5.199984e-5 | same | exp(-9.86427): same row, "ln k2" |
| `k_3` | 1.435e-5 | 1.4912e-5 | **inferred** from the digitized contours of the panel |
| `k_4` | 2.3040e-3 | 2.2397e-3 | **inferred**, same fit |
| `k_5` | 1.062501e-3 | 1.1110e-3 | `'fig11b'`: exp(-6.84713), same row, the column labelled "ln k3"; `'fig11a'`: **inferred** |
| `K_6`, `K_7`, `K_8`, `K_9` | 4.7840, 0.5918, 3.5116, 81.6870 | same | as printed in the script |
| `K_10` | 52.0472 | 23.28 | `'fig11b'`: as printed in the script; `'fig11a'`: **inferred** |

The inferred values are least-squares fits of the model's contours to the
digitized contour lines, the other constants held fixed (not shipped). How
well they are determined (spread between independent fits of the same
contours): `'fig11b'`: `k_3` and `k_4` about 1 %. `'fig11a'`: `k_3` about
10 %, `k_4` 2 %, `k_5` 1 %, `K_10` 15 %. `K_6` has no effect on the maps,
and a 10 % change of `K_8`, `K_7`, `K_9` or `K_10` moves them by at most
8e-8, 4e-5, 9e-4 and 9e-4 (measured with the package). The authors' own
vectors cannot be recovered beyond this: the script obtains them from a
particle swarm (3 iterations) followed by five least-squares iterations.

That the column "ln k3" of `11-4_Arrhenius.csv` is this model's `k_5` is
supported by the figure: with all five rate constants fitted to the contours
of Fig. 11b (not shipped; the same optimum from five starting points, rms
offset 0.254 px), `k_5` came out as 1.0614e-3, 0.10 % below the value in the
table (0.30 % below with `k_1` and `k_2` held at their values in the 50 C
row). In the fit of all five, `k_1` came out 0.03 % above its value in that
row and `k_2` 2.2 % below. The contours determine `k_2` least well of the
three: moved by 2.2 % either way with the other four refitted, `k_2` raises
the rms offset from 0.254 to 0.255 px, `k_5` to 0.38 px and `k_1` to 0.9 px.
The value of `k_5` in the table, 1.0625e-3, would print as 0.0011, not as
the 0.0010 of the script's starting vector, so that vector and the table are
not from the same run.

**The two panels need different sets.** Since yield <= selectivity
everywhere, the late side of a yield contour must lie earlier than that of
the selectivity contour of the same level: with one set of constants
(`'fig11b'`) it lies 0.005 to 0.09 % earlier for the levels 0.1 and 0.3 above
62 C (measured with the package, row by row). In the published figure the
late 0.1 and 0.3 contours of panel b lie later than those of panel a above
62 C, by 2.2 % in the median and up to 4.4 % (a median of 1.2 to 2.4 px, up
to 4.6 px; panel a interpolated to the temperatures of the panel b samples).
No single run can produce both panels. With the `'fig11b'` set panel a is
off by 1.07 px, with the `'fig11a'` set panel b by 1.45 px. Why the panels
differ (separate fits, or other settings) is not known; the `'fig11a'` set
is one of several that fit panel a equally well.

### Agreement with the published Fig. 11

There is no table of the heat maps in the deposit. The reference is the
published figure: the contour lines 0.1, 0.3 and 0.5 of the two main panels,
digitized (`reference_data/`). Offsets are normal to the lines, in pixels of
the published raster (1 px = 3.07 min or 0.148 C); the digitization itself
is good to about 0.3 px. Measured with the package (CVODE, relative
tolerance 1e-9):

| Comparison | Panel a (`'fig11a'` set) | Panel b (`'fig11b'` set) |
|---|---|---|
| rms offset of the contour samples | 0.17 px | 0.27 px |
| largest mean offset of one contour branch | 0.26 px | 0.22 px |
| largest cell of the map | 0.5906 at 95 C, 26.76 min | 0.5483 at 84 C, 49.10 min |
| maximum of the 50 C row | 0.5705 at 255.8 min | 0.5331 at 322.8 min |
| red dashed line minus the digitized one, along time: mean / rms | +0.60 / 0.93 px | +0.15 / 1.48 px |

The red dashed line joins the row maxima (`heat_map_ridge`), so it is
defined to one grid column, 5.585 min or 1.82 px; a line that sat one column
off would be 1.82 px off in the mean. The digitized line has 118 (a) and 226
(b) samples. The tests require the mean offset within half a column and the
rms within one column. The largest single offsets, 2.1 px (a) and 6.5 px
(b), are at 31 to 32 C, where the row maxima are flat and the line advances
several columns per row. In panel b the agreement holds over the whole line
and not for every part of it: the mean offset is -1.3 px between 45 and
60 C and +1.1 px between 60 and 80 C (0.7 and 0.6 of a column).

Read from the colours of the published maps (not shipped): the maps agree
with the colour field to an rms of 0.004 (panel a) and 0.005 (panel b) where
the value is above 0.3, and the 50 C rows peak at 0.567 and 0.532 (+-0.01).

The script integrates each of the 299 grid intervals with `ode23s` at
default tolerances. Repeating those solver runs with the `ode23s` port of
`solver_emulation.py` and the script's equations changes the maps by at
most 1.1e-4 (panel a) and 1.4e-4 (panel b), so the solver does not matter
at the level of the figure. The difference grows with temperature and is
largest in the 100 C row (3e-5 at 84 C, below 1e-6 at 50 C). The package
tests repeat the 100 C row of each panel; the other rows were run outside
the tests (not shipped).

### Paper versus deposit for Fig. 11

1. **Model.** The text describes the heat maps as the O2 model extended with
   the Arrhenius relation. The caption says only that the Fig. 9b data gave
   `kref`, which is what the script does: it fits its own simplified
   reversible model to those data. A reader of the text alone would compute
   the maps with the Table 3 model; the next section gives the difference.
2. **Arrhenius scaling of the equilibrium constants**, with the activation
   energies of the rate constants, is not mentioned in the paper. Two of the
   five activation energies are not in the paper either (51.79 kJ/mol for
   3-HP -> D, and 10 kJ/mol for MA -> OXA, a round value).
3. **Time axis.** "27 min" and "49 min" in the text are the grid columns
   26.76 and 49.10 min of the plotted time, i.e. 16.8 and 39.1 min of
   reaction.
4. **Maximum selectivity, "59.1 % (95 C, 27 min)".** Reproduced: 0.5906 in
   that cell, which is the largest cell of the map.
5. **Maximum yield, "55.2 % (84 C, 49 min)".** The cell is reproduced as the
   largest of the map, but its value is 0.548. None of the constant sets found
   that match the contours gives 0.552, and the colour of the published map
   along its crest reads 0.549 (+-0.006) (both not shipped). The crest is
   flat: 13 cells between 76 and 90 C are within 0.001 of the maximum, so the
   cell that comes out largest is sensitive to small changes; the tests check
   the value in the cell, not that it is the largest.
6. **The two "observed" maxima, 56.9 % and 50.5 %.** The text sets the
   59.1 % against "the observed maximum of 56.9 % at 50 C (purple star on
   plot)" and the 55.2 % against an "experimentally observed maximum" of
   50.5 %. The abstract places both at one condition ("Malonic acid
   selectivity of 56.9 % and yield of 50.5 % were achieved at 3 bar oxygen,
   equimolar NaOH:3-HP ratio, and 50 C"), which is that of the Fig. 9b
   experiment. The paper does not say how either number was computed.
   Compared with the deposited points of that experiment (tested):
   - **50.5 % is matched.** The quantity of panel b, `(1 - X_HP3)*X_MA`, is
     largest at 360 min: (1 - 0.0647)*0.5401 = 0.505 (next largest 0.484).
     `X_MA` alone, which is how the yield of the paper's Eq. 3 (moles of MA
     per mole of 3-HP fed) would be read, is 0.517, 0.548, 0.540 and 0.493
     at 180, 240, 360 and 450 min, and 0.505 at no sampling time.
   - **56.9 % is not identified.** The quantity of panel a is `X_MA`, which
     is at most 0.548 (240 min). MA per 3-HP converted, `X_MA/(1 - X_HP3)`,
     is 0.621, 0.577 and 0.513 at 240, 360 and 450 min and higher before:
     0.569 at no sampling time. The 50 C row of the reproduced selectivity
     map peaks at 0.5705 at 255.8 min, next to the purple star (254 min),
     but that is a model value, not an observed one. The 8 h screening
     experiment at the same nominal condition has an MA selectivity of
     56.4 % (Fig. 3a), which is not 56.9 either.

   So the match of the 50.5 % is consistent with the reading of the "yield"
   given above, but it stands alone and is an inference: the same reading
   does not give the 56.9 %. The paper's use of "yield" is not uniform
   either. In its Fig. 10 passage the "MA yield" of an experiment is `X_MA`
   itself (2.3 % at 25 C and 56.8 % at 100 C are the deposited `X_MA` of
   0.023 and 0.568; the 52.6 % given for 75 C is not in the deposited
   75 C data, whose largest experimental `X_MA` is 0.436, see item 10 of
   "Paper versus deposit: discrepancies that affect these models"), and in
   the screening figures it is Eq. 3 (in the deposited tables it equals
   conversion times MA selectivity).
7. **Stars.** The red star of panel a sits on the largest cell (27.5 min,
   94.7 C). The red star of panel b is drawn at 55.9 min, 7 min right of the
   cell named in the text. The purple stars are at 254 min (a) and 336 min
   (b) on the 50 C row, whose maxima are at 255.8 and 322.8 min.
8. **Insets.** The insets of the published panels are not exact
   magnifications of the main panels: at equal temperature their lines are
   7 to 10 min later (a) or up to 15 min off (b; 41 min for the late 0.3
   contour). `figures.py` draws the insets as zooms of the same matrix,
   with the contour labels on the same lines as published (those of the
   three steep early lines are staggered by up to 1 C so that they do not
   touch).

`plot_fig11` departs from the published figure in two more things, on
purpose: the largest cell of each map is annotated with its value, and the
colour-bar label of panel a reads "MA Selectivity" (the published panel
prints "MA Selectvity").

## Temperature-extended Table 3 O2 model versus the published heat maps

The paper's text describes Fig. 11 as an Arrhenius extension of the O2
model; the deposited script computes it with the simplified reversible
model. This section measures how far the first is from the second.
`table3_o2_heat_maps()` computes the two Fig. 11 quantities with
`load_o2_model(T)` on the grid of the figure, each row started from the feed
at time 0, with the time axis taken as the reaction time (as a reader of the
paper would). Compared with the same digitized contours:

| | Table 3 O2 model | Simplified model (previous section) |
|---|---|---|
| rms offset of the contour samples, panel a / b | 7.6 / 7.6 px | 0.17 / 0.27 px |
| largest cell, "selectivity" | 0.587 at 91 C, 21 min | 0.591 at 95 C, 27 min |
| largest cell, "yield" | 0.542 at 84 C, 38 min | 0.548 at 84 C, 49 min |
| 50 C row maximum, "selectivity" / "yield" | 0.566 / 0.527 | 0.571 / 0.533 |

The offsets of the Table 3 model are linearized estimates (value difference
over gradient) clipped at 30 px per sample, as returned by
`contour_offsets(..., max_offset=30)`; without the clip panel a is unchanged
(7.6 px) and panel b reads 8.8 px. They are indicative only. Measured
directly on the contour lines (`contour_shifts`): the late side of every
contour of the Table 3 model is earlier than published, by 8 to 11 % of the
time in the median per line (28 to 74 min; 9 to 260 min over the whole
lines), which is the same as 1.3 to 1.8 C lower at a given time. The early
sides differ less, by a median of 1 to 6 min in either direction.

If the time axis of the Table 3 maps is also shifted by the 10 min of the
published figure, the rms offsets are 6.8 px (a) and 6.1 px (b), and the
largest cells are 0.587 at 90 C, 32 min and 0.542 at 81 C, 55 min. Either
way the Table 3 model gives maps of the same shape as the published ones,
with maxima 0.003 to 0.006 lower, but it does not reproduce their contour
lines: the difference is 20 times the reading uncertainty of the figure
(0.3 px) or more.

Part of the difference is in the structure and not only in the constants:
the Table 3 model has the OPA intermediate and no reverse terms, its O2
factor stays above 0.999 where the simplified model's falls to 0.83, and its
Arrhenius factors are those of the simplified model's steps assigned to six
steps.

## Reversible O2 variant (reproduces the curves of Fig. 10a-c)

### What it is

The deposit holds, for the O2 experiments at 25, 75 and 100 C, the model
curves drawn in Fig. 10a-c (`11-1_25C_Kinetic.csv`, `11-2_75C_Kinetic.csv`,
`11-3_100C_Kinetic.csv`). Neither posted script produces them as posted:
`APCATB_D_25_08770_kinetic_opt.m` is the irreversible six-reaction model of
Table 3, and `APCATB_D_25_08770_heatmap.m` is the simplified reversible model
without OPA. The README of the deposit says that the posted codes are "an
example" and that "details were adjusted for each time (H2O2, O2, parameters
to fit, etc)"; how the curves of Fig. 10a-c were computed is not stated.
The curves are reproduced by the reaction network of the Table 3 O2 model
with reactions 1, 2, 5 and 6 made reversible, with a separate constant set
for each curve. This variant is shipped as
`pd_c_3hp_oxidation_o2_reversible_antimony.txt`
(`load_o2_reversible_model(parameter_set)`, with the sets named `'25C'`,
`'75C'` and `'100C'`; `read_fig10_o2`, `emulate_fig10abc` and
`figures.draw_fig10abc` take the temperature of the experiment in C or the
same label). **Its structure and all of its constants are inferred from the
deposited curves; nothing of it is printed in the paper, the SI or the
deposit.** It is a reproduction aid.

Rates in M/min (the model is written in concentrations, as the other
models; `C0` = 0.06 M):

| | Rate | Stoichiometry |
|---|---|---|
| R1 | `k_1*[HP3] - k_1r*[OPA]*[W]/C0` | 3-HP -> OPA + W |
| R2 | `k_2*[OPA] - k_2r*[MA]` | OPA -> MA |
| R3 | `k_3*[MA]` | MA -> AAc |
| R4 | `k_4*[MA]` | MA -> OXA + W |
| R5 | `k_5*[HP3] - k_5r*[D]` | 3-HP -> D |
| R6 | `k_6*[MA] - k_6r*[D]` | MA -> D |

`W` is the water formed by the reactions, a state that starts at zero; it is
not the solvent. Every term is first order in each species it contains. There
is no oxidant factor, no catalyst mass and no temperature law: the constants
are effective ones in 1/min for the 0.2 g batch (`k_1r` in 1/min per unit of
`[W]/C0`, i.e. with OPA and W counted in mol per mol of initial 3-HP). The
plotted quantities are `X_i = c_i/([HP3]+[MA]+[AAc]+[OXA]+[D])` as in Fig. 9;
OPA is again left out of the sum. The deposited curves run from 0 to 1.1
times the last sampling time (2079, 1584 and 1716 min).

### Constants

In 1/min; all inferred (see the next subsection), given to six significant
digits.

| Constant | `'25C'` | `'75C'` (default) | `'100C'` |
|---|---|---|---|
| `k_1` (3-HP -> OPA + W) | 0.0465036 | 0.00918893 | 0.0414885 |
| `k_1r` (OPA + W -> 3-HP) | 51.0644 | 1.88278e-5 | 0.00474984 |
| `k_2` (OPA -> MA) | 0.0088933 | 0.0490036 | 0.122878 |
| `k_2r` (MA -> OPA) | 0.00136441 | 0.00193586 | 0.0103261 |
| `k_3` (MA -> AAc) | 2.43424e-4 | 4.9301e-4 | 0.00169028 |
| `k_4` (MA -> OXA + W) | 6.17632e-5 | 7.97777e-6 | 7.86978e-6 |
| `k_5` (3-HP -> D) | 1.01563e-6 | 0.00160383 | 9.32644e-7 |
| `k_5r` (D -> 3-HP) | 1.47217e-7 | 2.86632e-5 | 0 |
| `k_6` (MA -> D) | 0.00284329 | 0.00613222 | 0.0157177 |
| `k_6r` (D -> MA) | 3.12354e-4 | 2.56106e-5 | 1.05066e-4 |

How well the deposited curve determines each constant, in per cent of the
value. Each entry is the larger of two estimates from the fits of the next
subsection (not shipped; the package does not recompute them): the
linearized standard deviation of the least-squares fit, and the spread
between independent refits from different starting points. The next
subsection has the corresponding table measured by the package tests.

| Constant | `'25C'` | `'75C'` | `'100C'` |
|---|---|---|---|
| `k_1` | about 0.1 | < 0.01 | < 0.01 |
| `k_1r` | about 0.15 | 5 to 10 | 0.03 |
| `k_2` | 0.01 | < 0.01 | 0.01 |
| `k_2r` | 0.1 | 0.15 | 0.05 |
| `k_3` | 0.01 | < 0.01 | < 0.01 |
| `k_4` | 0.02 | 0.2 | 0.05 |
| `k_5` | 1 to 3 | < 0.01 | about 8 |
| `k_5r` | about 60 | 0.04 | not detectable, set to 0 |
| `k_6` | 0.05 | < 0.01 | < 0.01 |
| `k_6r` | 0.1 | 0.05 | 0.01 |

The weakly determined ones (`k_5r` and `k_5` at 25 C, `k_1r` at 75 C, `k_5`
at 100 C) move the curves by less than 1e-5 within their uncertainty
(measured with the package for `k_1r` at 75 C: a 10 % change moves the
curves by 7e-6, dropping the term by 7e-5). There is no set for 50 C: the
50 C curve (Fig. 9b) is the irreversible Table 3 model, and the loader
accepts only the three set names.

The three sets are independent fits and follow no temperature law. `k_1` is
0.047, 0.0092 and 0.041 1/min at 25, 75 and 100 C, `k_5` is 1.0e-6, 1.6e-3
and 9.3e-7 1/min, and the ratio `k_1/k_1r`, which would be an equilibrium
constant, is 9.1e-4, 488 and 8.7. The constants should therefore be read as
fit parameters of each curve, not as kinetic or thermodynamic quantities.
For that reason the loader takes the name of a set and not a temperature:
there is nothing to interpolate.

For comparison only, the forward constants in the units of Table 3
(`k/(m_cat*x_O2_ref)` for R1, R2 and R4, `k/m_cat` for the others, with
`m_cat` = 0.2 g and `x_O2_ref` = 0.005085), next to the Table 3 O2 values.
They are not interchangeable with Table 3: R1, R2, R5 and R6 also have
reverse terms here.

| | Table 3 (50 C) | `'25C'` | `'75C'` | `'100C'` |
|---|---|---|---|---|
| `k_1` | 7.14 | 45.7 | 9.04 | 40.8 |
| `k_2` | 139 | 8.74 | 48.2 | 121 |
| `k_3` | 3.16e-4 | 1.22e-3 | 2.47e-3 | 8.45e-3 |
| `k_4` | 5.86e-2 | 6.07e-2 | 7.84e-3 | 7.74e-3 |
| `k_5` | 1.13e-2 | 5.08e-6 | 8.02e-3 | 4.66e-6 |
| `k_6` | 5.52e-3 | 1.42e-2 | 3.07e-2 | 7.86e-2 |

### How the model was inferred and how firmly each feature is established

The deposited curves have 18, 19 and 24 unevenly spaced steps, like the
Fig. 9 curves. The constants above were obtained by fitting them through a
repeat of the authors' solver run to the deposited nodes, values and node
times (fits that are not shipped). What the package tests establish with
the shipped repeat of that run (`solver_emulation.py`; section "Emulation of
the authors' solver run"): with these constants the run, left to choose its
own steps from the feed onwards, takes the deposited number of steps and
returns the node times within 0.0054 min and the values within 7.2e-6 (rms
2.5e-6 or less); 91 of 95, 98 of 100 and 113 of 125 values are identical
after rounding as in the file. The residuals are the print rounding of the
files. MATLAB itself was not run.

How tightly the nodes pin each constant on its own, measured by the tests:
the relative change of one constant that makes the repeated run miss a
limit (0.01 min, 1e-5) at least twofold, as the smallest power of ten.

| Constant | `'25C'` | `'75C'` | `'100C'` |
|---|---|---|---|
| `k_1` | 1e-4 | 1e-4 | 1e-4 |
| `k_1r` | 1e-3 | 1 (a doubling) | 1e-2 |
| `k_2` | 1e-4 | 1e-3 | 1e-3 |
| `k_2r` | 1e-3 | 1e-2 | 1e-3 |
| `k_3` | 1e-2 | 1e-3 | 1e-3 |
| `k_4` | 1e-2 | 0.1 | 0.1 |
| `k_5` | 0.1 | 1e-3 | not detected |
| `k_5r` | not detected | 1e-2 | zero in the set |
| `k_6` | 1e-3 | 1e-4 | 1e-4 |
| `k_6r` | 1e-2 | 1e-2 | 1e-2 |

"Not detected": a change of 10 % leaves the run within the limits. One unit
of the last (sixth) printed digit, up or down, is not detectable for any of
the thirty constants, and neither is a change of 1e-5 of the value, the most
that such a unit can be (both tested). With the limits taken as they are, not
twofold, the smallest single change that is detected is about 1.5e-5 (`k_1`
at 75 and 100 C; then 2e-5 for `k_6` at 100 C and 3e-5 for `k_2` at 25 C and
`k_6` at 75 C; measured on a grid of changes, not asserted by a test). So six
digits are what the comparison needs, and the last digit of each constant is
not significant on its own. For comparison, the accurate-solver tests detect
a 10 % error in most constants only.

What the negative controls of the tests show. Each keeps the constants of
the set and changes one feature; the numbers are the largest difference
from the deposited nodes at 25 / 75 / 100 C (limit 1e-5):

- **Reverse of R1.** Without it (`k_1r` = 0): 0.94 / 7.6e-5 / 0.011. All
  three curves need the term.
- **The reverse of R1 is bilinear in OPA and W.** With `k_1r*[OPA]` in its
  place: 0.10 / 8.3e-5 / 0.0095.
- **W is formed by R1 and R4.** With water from R1 only: 9.1e-4 (node times
  1.6 min off) at 25 C; at 75 C the run is unchanged to the print rounding,
  and at 100 C the values stay within 8.2e-6 while the node times move by
  0.012 min, just beyond their limit of 0.01 min. So the 25 C curve, where
  the reverse of R1 matters most, shows where the water comes from; the
  100 C curve does so only marginally, and the 75 C curve not at all.
- **OPA is left out of the plotted sum.** With OPA in the sum: 0.029 /
  0.098 / 0.157 (the steps are the same; only the plotted quantity
  changes).

These controls do not refit the constants of the alternative. The
following comes from such refits, which are not shipped (rms of the best
fit of each alternative to the deposited curves, with an accurate solver,
at 25 / 75 / 100 C; this structure refitted the same way reaches 1.6e-5 /
2.5e-4 / 1.2e-4, the level of the deposited curves' solver error):

- **Not the Table 3 model with other constants.** The irreversible
  six-reaction model reaches 7.3e-3 / 2.7e-3 / 3.2e-3.
- **Not the simplified reversible model.** A model without OPA and with
  linear reverse terms reaches 1.1e-3 / 1.8e-2 / 2.2e-2.
- **Linear reverse of R1, refitted.** The fit is 1.1e-3 / 2.5e-4 / 1.2e-4.
  The 25 C curve excludes it: its best fit, reached from random starting
  points and from the shipped constants alike, is about 70 times worse than
  that of the bilinear form and no better than that of the model without
  OPA. The 75 and 100 C curves do not exclude it when solved accurately (a
  fit through the repeated solver run excluded it at 100 C as well).
- **Water from other reactions, refitted.** In the fit of the 25 C curve
  through the repeated solver run the misfit (chi-square) was 81 with water
  from R1 and R4, 15164 with water from R1 only and above 1e5 for the other
  combinations tried; the coefficient of R4 was 1.0 +- 0.05. R1 and R4 are
  the steps that form water in the O2 stoichiometry of SI Table S4 (the
  reaction list of Fig. 9b writes R4 without water).
- **Not identifiable:** how the authors' code factors the constants (oxidant
  fraction, pressure, equilibrium constants); only the effective constants
  are determined. Reverse terms of R3 and R4, or an additional linear reverse
  of R1, would not be detectable if their constants were small.

This is the only structure found that reproduces the three curves; that does
not prove that it is the one the authors used.

### Agreement with the deposit

Measured with the package (CVODE, relative tolerance 1e-9) at the node times
of the deposited curves; values are normalized concentrations.

| Comparison | `'25C'` | `'75C'` | `'100C'` |
|---|---|---|---|
| Model vs deposited curve, rms / max over all values | 3.8e-5 / 1.8e-4 | 7.1e-4 / 2.5e-3 | 6.9e-4 / 2.7e-3 |
| the same, max per species (3-HP, MA, AAc, OXA, D) | 7.7e-5, 1.8e-4, 7.9e-6, 4.9e-6, 9.4e-5 | 1.7e-3, 2.5e-3, 1.4e-4, 4.0e-6, 1.7e-3 | 1.7e-3, 2.7e-3, 1.9e-4, 1.0e-6, 1.7e-3 |
| Model vs the straight segments between the nodes, max | 7e-4 | 0.010 | 0.009 |
| Repeated solver run vs deposited curve, rms / max over all values | 2.1e-6 / 7.2e-6 | 2.5e-6 / 6.3e-6 | 2.3e-6 / 5.7e-6 |
| Model vs an accurate solution of the same equations, max over the plotted fractions | 1.5e-8 | 2.5e-8 | 2.4e-8 |
| Model vs experiment, pooled R2 / rms (all rows) | 0.9996 / 0.008 | 0.990 / 0.029 | 0.993 / 0.026 |

The differences of the first rows are the solver error of the deposited
curves, not a difference in the model: the published lines are not converged
solutions, and with the solver run repeated the differences are at the print
rounding of the files (fourth row). The published lines are drawn as
straight segments between the nodes, which at 75 and 100 C lie up to 0.010
from the solution (third row).

The accurate solution of the fifth row is computed in the tests
(`test_pd_c_3hp_oxidation_fig10.py`): the equations of the table above,
written out again in the test module on their own basis (mol per mol of
initial 3-HP) with the same six-digit constants, are integrated with scipy
(`solve_ivp`, DOP853, relative tolerance 1e-13; LSODA agrees to 1e-12),
independently of the Antimony file and of roadrunner. The row is therefore
the integration error of the package at the tolerances the loader sets.

The tests accept 1.5 to 2 times the solver error of the deposited curves
(limits in `test_pd_c_3hp_oxidation_fig10.py`). What these limits detect,
measured with the package at 25 C, where the solver error is smallest: water
formed by R1 only (rms 2.2e-4, max 9.0e-4), no reverse of R1 (max 0.97), and
a 10 % error in any constant except `k_5` and `k_5r`. At 75 and 100 C the
solver error of the deposited curves is larger and the limits are wider: a
10 % error in `k_1r`, `k_2r`, `k_5r` or `k_6r` (75 C) or in `k_1r`, `k_2r`
or `k_5` (100 C) passes them (measured with the package). The comparison
with the accurate solution (limit 1e-7) is the check that covers those: a
1 % error in any single constant of any of the three sets moves the
solution by at least 1.6e-7 (the smallest effects are `k_5r` at 25 C, `k_5`
at 100 C, 1.8e-7, and `k_1r` at 75 C, 6.6e-7; a 10 % error moves it by ten
times as much). That comparison guards the implementation, not the
inference. The comparison with the repeated solver run
(`test_pd_c_3hp_oxidation_solver_emulation.py`, previous subsection) is the
sharper check against the deposit.

### What the model implies

- **An undetected OPA pool.** OPA peaks at 2.9 % of the initial 3-HP after
  2.3 min at 25 C, 12.6 % after 42 min at 75 C and 20.1 % after 14.5 min at
  100 C (0.6, 0.07 and 0.04 % are left at the end of the curves). OPA was
  not detected in the experiments ("not detected with our analytical
  methods", Section 3.3) and is not in the plotted sum, so the plotted
  fractions add up to 1 while up to a fifth of the carbon of the model is
  in OPA.
- **Water of reaction.** `[W]/C0` ends at 0.133, 0.881 and 1.000.
- At 25 C the forward step of R1 is fast (`k_1` = 0.047 1/min) and is held
  back by its reverse: OPA and W reach about 3 % of the feed within minutes,
  and after that the net rate of R1 is small.

### Paper versus deposit for Fig. 10a-c

1. **Reversibility.** The text (Section 3.4) treats all reactions as
   irreversible and SI Table S4 lists equilibrium constants of 1e30 and
   more for R1 and R2. The curves of Fig. 10a-c are reproduced only with
   reverse terms in R1, R2, R5 and R6.
2. **Model and constants are not given.** The paper prints no constants for
   25, 75 or 100 C, and neither posted script produces the curves as posted.
   The Arrhenius plot (Fig. 10d) shows constants of a different model (next
   section).
3. **Units and orders.** The constants are per minute for the 0.2 g batch,
   not per gram of catalyst as in Table 3, and all terms are first order;
   the reaction orders of SI Table S4 are not used.
4. **R2.** The published panels print R2 = 0.99 for all three temperatures.
   The pooled R2 of this model against the deposited experimental points
   is 0.9996, 0.990 and 0.993 (`figures.py` prints the measured values).
5. **Time axis.** The published panels end at 1440 min. The ends of the
   model curves are in the deposit but not in the figure, like the late
   samples (item 11 of "Paper versus deposit: discrepancies that affect
   these models"; item 10 there is the maximum MA yield at 75 C). The
   redrawn panels (`draw_fig10abc`) end at 1440 min too: the late samples
   are not drawn but are counted in the pooled R2 of item 4.
6. **Solver.** The published lines carry a solver error of up to 2.7e-3 at
   75 and 100 C. (Section 2.6 names MATLAB R2024b, the README of the deposit
   R2024a.)

## Arrhenius plot (Fig. 10d)

### What is plotted

`11-4_Arrhenius.csv` holds the points and the lines of Fig. 10d. The points
are `ln k` at 25, 50, 75 and 100 C of three rate constants, labelled "k1",
"k2" and "k3" in the file and "3-HP -> MA", "MA -> AAc" and "MA -> Others"
in the figure. The paper does not say of which model, and the README of the
deposit names none where it describes the file: "Raw .csv data of Arrhenius
plot of O2 kinetic models" in its table of files, and "Arrhenius plot for
three reactions" in its item 11, which covers the four files of Fig. 10. The
README does describe the heat-map script as "MATLAB code for generating
heatmaps, expanded from 10. Arrhenius plot data" (the numbers in the
README's two script descriptions are one lower than in its list of the data,
where the kinetic model is item 10 and the Arrhenius plot item 11: the
posted O2 script is "used for 9. Kinetic model in xlsx dataset"). That is
consistent with the plotted constants following the rate law of that script,
but does not establish it. The package names the constants `k_1`, `k_2` and
`k_5`: the numbers these three reactions have in the authors' heat-map
script, whose activation energies for them are the ones the text quotes.
These are not the numbers of the O2 model, in which 3-HP -> MA is R1
followed by R2, MA -> AAc is R3 and MA -> D is R6. The abscissa is
`1/(T + 273.15)`; the file prints it to five decimals.

The points are taken here as `k_1` (3-HP -> MA), `k_2` (MA -> AAc) and `k_5`
(MA -> D) of the **simplified reversible model**, the model of that script,
fitted to the experiment of each temperature, for three reasons: the figure
labels name a one-step 3-HP -> MA reaction, which only that model has; with
the 50 C row as `k_1`, `k_2` and `k_5` and its other two rate constants
fitted, that model reproduces the contour lines of Fig. 11b, and with all five
rate constants fitted to those lines `k_1`, `k_2` and `k_5` come out 0.03 %
above, 2.2 % below and 0.10 % below the row (not shipped; of the three, `k_1`
is the one that multiplies the O2 activity of the script, and the lines
determine `k_2` least well; "Constants at 50 C" above); and refits of that
model to the experiments come close to the other points (below). Taking the
points as constants of that model means that they follow the rate law of its
script. The three reasons support that and do not establish it: the fits
behind the second and third are not shipped, and neither the paper nor the
deposit says it (previous paragraph). They are not constants of the model
that reproduces the curves of Fig. 10a-c: at 25 C, for example, the plotted
`k_1` is 7.75e-5 1/min, while the constant set of the Fig. 10a curve has a
forward constant of 0.047 1/min for 3-HP -> OPA with a fast reverse.

**The unit of k is inferred.** No source prints it: the file labels its
ordinate columns "exp" and "fitting" and gives a unit ("/K") for the
abscissa only, the README of the deposit gives none, and the axis of the
published figure reads "ln k". The constants are taken here to be in 1/min,
the time unit of the authors' two scripts. In that unit two of the three
50 C values lie within the confidence intervals of Table 3 as they stand,
and the third does with the O2 activity of the heat-map script (below);
read as 1/s or 1/h each value would be 60 times larger or smaller and none
would (tested). The slopes of the lines, and with them the activation
energies, do not depend on the unit. (In 1/min the 50 C values are also
three of the five rate constants with which the simplified reversible model
reproduces the contour lines of Fig. 11b, whose time axis is in min;
"Constants at 50 C" above.)

The constants whose logarithms are the deposited points (shortest decimals
that reproduce them within 5e-6), taken as 1/min:

| T (C) | `k_1` | `k_2` | `k_5` |
|---|---|---|---|
| 25 | 7.75e-5 | 6.64e-4 | 1.676e-3 |
| 50 | 2.14e-3 | 5.2e-5 | 1.0625e-3 |
| 75 | 2.43685e-3 | 4.44096e-4 | 5.4797e-3 |
| 100 | 1.00763e-2 | 1.59645e-3 | 1.41038e-2 |

The points cannot be recomputed: the script's fits are stochastic and
stopped after a few iterations. A refit of the simplified model to the
deposited experiments with an accurate solver (not shipped) gave an optimum
that, divided by the plotted value, was 1.000, 1.017, 0.998 (`k_1`, `k_2`,
`k_5`) at 100 C, 0.916, 1.085, 1.032 at 75 C and 0.996, 1.045, 0.917 at
50 C, and was ill-determined at 25 C. This supports the assignment of the
three columns and shows that the points are not converged optima. The points
are therefore taken from the file (`fig10d_points`).

**The 50 C points against the Table 3 constants.** The plotted 50 C values
are not numerically the Table 3 constants. In 1/min for the 0.2 g batch,
the first-order constants of the O2 model at 50 C (with `yO2` = 1) are
`k_1*m_cat*x_O2_ref` = 7.26e-3 for 3-HP -> OPA (R2, OPA -> MA, is 19.5
times faster), `k_3*m_cat` = 6.32e-5 for MA -> AAc and `k_6*m_cat` =
1.104e-3 for MA -> D. The rate laws behind the plotted constants are not
stated, so the two are compared on two bases (tested):

| | Plotted at 50 C | O2 model, 1/min | Table 3 interval, 1/min | Plotted over model, plain first order | The same with the O2 activity of the heat-map script |
|---|---|---|---|---|---|
| 3-HP -> MA (`k_1`; model R1) | 2.14e-3 | 7.261e-3 | 6.855e-3 to 7.668e-3 | 0.295, outside | 1.018, inside |
| MA -> AAc (`k_2`; model R3) | 5.2e-5 | 6.32e-5 | 3.52e-5 to 9.12e-5 | 0.823, inside | 0.823, inside |
| MA -> D (`k_5`; model R6) | 1.0625e-3 | 1.104e-3 | 1.044e-3 to 1.164e-3 | 0.962, inside | 0.962, inside |

- *Table 3 interval.* Table 3 prints the O2 constants as 7.14 (+- 0.4),
  3.16 (+- 1.4)e-4 and 5.52 (+- 0.3)e-3 per g and min and calls the +- values
  "the 95 % confidence interval". They are read here as half-widths in the
  units of the printed mantissa (the posted O2 script prints its intervals
  as a lower and an upper value of each constant), and scaled like the
  constants. Read as percentages they would contain none of the three
  plotted values on either basis.
- *Plain first order.* The plotted constant is taken as it stands. The O2
  model has 3.39, 1.22 and 1.04 times the plotted values.
- *O2 activity of the heat-map script.* In the deposited heat-map script,
  whose activation energies for these three reactions are the ones the text
  quotes, the constant of 3-HP -> MA multiplies an O2 activity `yO2*P/f0`
  and the other two constants multiply the MA mole fraction alone. The
  script starts with `yO2` = 1, pure O2, for which the activity is `P/f0` =
  3.5/1.01325 = 3.45 as posted. On that basis the plotted 2.14e-3
  corresponds to 7.39e-3 1/min. (These are the forward terms of `r_1`, `r_2`
  and `r_5` of the simplified reversible model.)

So the 50 C points of MA -> AAc and MA -> D agree with the O2 model within
the confidence intervals of Table 3 on either basis, and the 50 C point of
3-HP -> MA does if the plotted constant follows the rate law of the heat-map
script. That it follows that rate law is what the assignment to the
simplified reversible model above means; it is supported by the reasons
given for the assignment but not established (the README's "expanded from",
above, is consistent with it and no more).

### The lines

Each dashed line is the ordinary least-squares line of `ln k` on
`1/(T + 273.15)` through three of the four points
(`fig10d_regression_lines`, recomputed in the package):

| | Fitted to | Slope (K) | Intercept | Ea of the line | Ea in the text | Point left out, ln k above the line | Drawn from 100 C to |
|---|---|---|---|---|---|---|---|
| `k_1` | 25, 75, 100 C | -7208.69 | 14.7075 | 59.93 kJ/mol | 59.8 | 50 C, +1.45 | 25.0 C |
| `k_2` | 50, 75, 100 C | -8293.97 | 15.8974 | 68.96 kJ/mol | 91.6 | 25 C, +4.60 | 25.13 C |
| `k_5` | 50, 75, 100 C | -6265.49 | 12.6203 | 52.09 kJ/mol | 59.5 | 25 C, +2.00 | 23.50 C |

(Ea = -slope*R with R = 8.314 J/(mol K). The last column is the temperature
at which the recomputed line takes the `ln k` of the last row of the
deposited line; the 1/T printed in those rows, 0.00335, 0.00335 and 0.00337,
reads 25.36, 25.36 and 23.59 C.) The recomputed lines agree with the
deposited ones to 7e-6 in `ln k` at 100 C (the first row of each line) and
to 5e-7 at 25 C for `k_1`. At the abscissae as printed in the file the
difference is up to 0.04, which is the rounding of the printed 1/T; for the
same reason a regression of the deposited lines on the printed abscissae
gives 59.97, 69.32 and 52.38 kJ/mol. The fitted points lie within 0.21 of
their line. Which three points each line is fitted to is itself inferred:
only these subsets reproduce the deposited lines (any other subset of two to
four points misses the first row of the line by 0.013 or more; tested).

### Paper versus deposit for Fig. 10d

1. **Which model.** The caption calls the plot an "Arrhenius plot for three
   reactions", and the text does not name a model. The figure names a
   one-step 3-HP -> MA reaction, which the Table 3 model does not have (it
   has R1 followed by the faster R2). The 50 C points are not numerically
   the Table 3 constants, but they agree with them within the confidence
   intervals of Table 3 if the constant of 3-HP -> MA multiplies the O2
   activity of the heat-map script (above). The constants are taken to be
   those of the simplified reversible model ("What is plotted"); they are
   not those of the model that reproduces the curves in panels a to c.
2. **Each line leaves out one point** (table above); the figure shows all
   four points without saying so.
3. **Activation energies.** The text says that "Apparent activation energies
   (Ea) for key steps were obtained from Arrhenius plot in Fig. 10d" and
   gives 59.8 (MA production), 91.6 (MA to AAc) and 59.5 kJ/mol
   (decomposition). The drawn lines give 59.93, 68.96 and 52.09 kJ/mol. No
   subset of two, three or four of the plotted points gives 91.6 for `k_2`
   (nearest: 80.2 from 50 and 75 C) or 59.5 for `k_5` (nearest: 61.4 from 50
   and 75 C); for `k_1` the drawn line (59.93) is the nearest to 59.8. The
   text values are those of the heat-map script (59781.2, 91594.6 and
   59549.8 J/mol). Where they come from is not known; they are used as
   printed in the O2 model (`Ea_1` = `Ea_2`, `Ea_3` and `Ea_6`). The
   script's value for 3-HP -> D, 51788.73 J/mol (`Ea_5`), is not in the
   paper; it is close to what the 50 and 100 C points of the MA -> D
   constant give (51.85 kJ/mol). The simplified reversible model has the
   same four values under the script's numbers, in which `Ea_5` is the
   59549.8 J/mol of MA -> D and `Ea_4` the 51788.73 J/mol of 3-HP -> D
   (`Ea_1` and `Ea_2` are the other two); there, as in the script, they
   scale the 50 C constants (Fig. 11).
4. **The activation energies of the O2 model against the plotted points.**
   Lines with the text's activation energies through the 50 C points
   (`fig10d_text_lines`, drawn in `plot_fig10_table3_o2`) miss the other
   points: for `k_1` the 25, 75 and 100 C points lie 1.45, 1.47 and 1.43
   below the line (a factor of 4.2 to 4.3 in k), for `k_2` 5.41 above and
   0.30 and 1.14 below, for `k_5` 2.31 and 0.05 above and 0.38 below. Each
   of these lines is the plotted 50 C constant multiplied by the Arrhenius
   factor `k(T)/k_ref` that the O2 model applies to R1 and R2, R3 or R6
   (tested). How they compare with the model's own constants `k_j(T)`
   follows from the 50 C points (above), because the offset is the same at
   every temperature. For MA -> AAc and MA -> D the drawn lines run 0.20
   and 0.04 in `ln k` below the model's, inside the Table 3 intervals
   (-0.59 to +0.37 and -0.06 to +0.05 about the model's line). For
   3-HP -> MA the drawn line runs 1.22 below the model's if the plotted
   constant is taken as it stands, and 0.02 above it, inside the interval
   (-0.06 to +0.05), with the O2 activity of the heat-map script. So if
   the plotted constants follow the rate law of that script, which is
   supported by the reasons given for the assignment but not established
   (the README's "expanded from" is consistent with it and no more; section
   "What is plotted"), the solid lines of panel d are the O2 model's own
   `k_j(T)` to within the uncertainty of Table 3. In either case a constant
   factor does not change a slope, so the comparison of the activation
   energies holds for the model. These lines are also how the heat-map
   script extrapolates the three constants from 50 C (Fig. 11).

## Temperature-extended Table 3 O2 model versus the 25, 75 and 100 C experiments

The text describes the O2 model as extended with the Arrhenius relation.
`load_o2_model(T)` is that extension (an NSKinetics addition, see
"Temperature dependence of the O2 model"). The paper compares "modeled
data" with the experiments at other temperatures in Fig. 10a-c without
saying which model; those curves are reproduced only by separate constant
sets, not by this extension. This section compares the extension with
the deposited experimental points (all rows of each file, t = 0 included,
also the rows after the 1440 min at which the published panels end; item 11
of "Paper versus deposit: discrepancies that affect these models"), next to
those per-temperature sets (`load_o2_reversible_model`). Measured with the
package; the time at which `X_HP3` falls to 0.5 and the times of the maxima
are interpolated on a 0.1 min grid:

| | 25 C | 75 C | 100 C |
|---|---|---|---|
| Table 3 O2 model at T: pooled R2 / rms | 0.579 / 0.243 | 0.644 / 0.174 | 0.746 / 0.150 |
| Per-temperature set: pooled R2 / rms | 0.9996 / 0.008 | 0.990 / 0.029 | 0.993 / 0.026 |
| Time at which the plotted fraction `X_HP3` falls to 0.5, Table 3 at T / per-temperature set | 458 min / not reached (0.878 left at 1890 min) | 16.0 / 75.0 min | 4.2 / 22.0 min |
| Maximum MA fraction, Table 3 at T | 0.522 at 1479 min | 0.584 at 51.7 min | 0.587 at 12.8 min |
| Maximum MA fraction, per-temperature set | 0.028 at 374 min | 0.403 at 134.6 min | 0.549 at 42.2 min |
| Maximum MA fraction, experiment | 0.023 at 1890 min | 0.436 at 180 min | 0.568 at 40 min |

(The text gives a maximum MA yield of 52.6 % at 75 C, which is not in the
deposited data; item 10 of the same list.)

The temperature-extended Table 3 model does not describe these experiments.
It is too fast at all three temperatures:

- **25 C.** It predicts 0.051 of the 3-HP left at 1890 min and an MA
  fraction of 0.51; measured are 0.876 and 0.023. Its first-order constant
  of 3-HP loss (R1 plus R5), 1.6e-3 1/min, is about 22 times the measured
  average (7e-5 1/min, from 0.876 left at 1890 min).
- **75 C.** The predicted `X_HP3` falls to 0.5 after 16 min; measured, it
  is 0.53 at 60 min and 0.28 at 120 min. The predicted MA maximum is higher
  (0.58 against 0.44) and 3.5 times earlier.
- **100 C.** The height of the MA maximum is about right (0.59 against
  0.57), but it comes at 13 min instead of 40 min, and the final acetic
  acid fraction is 0.18 against a measured 0.09 to 0.11.

For 75 C this follows from the data: the measured 3-HP conversion at 75 C
is no faster than at 50 C (`X_HP3` falls to 0.5 between 60 and 90 min at
both temperatures), whereas the Arrhenius factor of the model from 50 to
75 C is 4.9 for R1 and R2. The paper's own Arrhenius plot shows the same for
all three temperatures: in Fig. 10d the 50 C point of 3-HP -> MA (`k_1`)
lies a factor of 4.3 above the line through the 25, 75 and 100 C points
(previous section). An extrapolation anchored at 50 C with the activation
energies of the text therefore overshoots at 25, 75 and 100 C.

That anchor accounts for part of the misfit only (measured with the package
and pinned by a test). With `k_1_ref` and `k_2_ref` lowered by the factor
4.3, the pooled R2 rises from 0.579, 0.644 and 0.746 to 0.818, 0.851 and
0.749: nothing is gained at 100 C. Most of the rest is the 3-HP -> D step
(R5). Its 50 C constant, `k_5*m_cat` = 2.3e-3 1/min, is extrapolated with an
activation energy that appears only in the script (51.79 kJ/mol): Fig. 10d
has no 3-HP -> D constant, and the value is close to the slope through the
50 and 100 C points of MA -> D (51.85 kJ/mol; item 3 of "Paper versus
deposit for Fig. 10d"). With `k_5_ref` lowered tenfold as well, the pooled
R2 is 0.957, 0.982 and 0.928. At 25 C, for example, the constant of 3-HP
loss of the model is 1.57e-3 1/min, of which 4.5e-4 is R5; dividing R1 by
4.3 still leaves 7.1e-4, ten times the measured 7.0e-5. The two changes are
illustrations, not fits. The constant sets that reproduce the curves at 25
and 100 C have about 1e-6 1/min for R5.

No value of `Ea_5` does what the tenfold lower constant does. With the 50 C
constant of R5 kept (and R1 and R2 lowered by 4.3), a low `Ea_5` fits 75 and
100 C and not 25 C, and a high one the reverse: the pooled R2 is 0.26, 0.99
and 0.93 for `Ea_5` = 0 and 0.96, 0.10 and 0.20 for 120 kJ/mol. From 0 to
150 kJ/mol the lowest of the three is 0.79 at best, between 47 and
48 kJ/mol, where the values at 25 and 100 C cross (tested in steps of
1 kJ/mol: 0.787, 0.879 and 0.794 at 47 kJ/mol, 0.794, 0.874 and 0.785 at
48). In the model this step is too fast on both sides of 50 C, which an
activation energy cannot mend.

The heat-map script extrapolates the simplified model in the same way
(previous section, item 4): its `k_1` at 75 and 100 C is 4.3 and 4.2 times
the value of the Fig. 10d points at those temperatures. How far the heat
maps of Fig. 11 are from the time courses at 75 and 100 C has not been
measured here; for the 8 h screening experiments see "Experiments of the
paper that the models were not built on".

`plot_fig10_table3_o2` in `figures.py` draws this comparison; the box in
each panel gives the pooled R2 of the Table 3 model and of the
per-temperature set, and panel d is the Arrhenius plot of the previous
section.

## Not modelled

Liquid volume and sampling withdrawals, gas-phase composition, pH, catalyst
deactivation, the loss of activity at higher O2 pressure and at higher H2O2
charge that the paper's screening experiments show (Figs. 2a and 3a), and any
dependence on O2 pressure other than through the fixed factor `x_O2_ref` (O2
model) or the factor `P/f0` of the script (simplified model). The Table 3
constants are per gram of catalyst for the paper's batch and are used with
`m_cat` = 0.2 g; the compartment volume is a nominal 1 L and cancels from all
concentrations (see "Operating conditions, units and the range of the data").

## Analyses that are not in the repository

The following results are quoted above from work done outside the package.
Its scripts are not shipped and the package tests do not repeat it, so
these numbers cannot be re-derived from the repository:

- **Fits through the repeated solver run.** The repeat of the authors'
  solver run is shipped and tested (`solver_emulation.py`). The
  least-squares fits made through such a repeat are not: they produced the
  constant sets of the Fig. 9 nodes (`'fig9a_nodes'`, `'figS12a_h2o2_nodes'`,
  `'fig9b_nodes'`), all thirty constants of the reversible O2 variant and
  their uncertainties, and the refits of the alternatives that were
  excluded, made through the repeated run or with an accurate solver (the
  irreversible model, the model without OPA, other water stoichiometries, a
  refitted linear reverse of R1, the O2 fraction held constant). The tests
  show that the shipped constants return the deposited nodes and how far
  each constant can be moved on its own; they do not show that no other
  constants or rate laws would return them.
- **Fits to the experiments:** the refit of the O2 constants with and
  without the 120 min row; the refits of the simplified model behind the
  assignment of the Fig. 10d points.
- **Fits to the digitized contour lines of Fig. 11:** `k_3` and `k_4` of the
  `'fig11b'` set and `k_3`, `k_4`, `k_5` and `K_10` of the `'fig11a'` set,
  their spreads, the offsets "after refitting", and the fits of the five
  rate constants to the contour lines of Fig. 11b quoted in "Constants at
  50 C" and for the assignment of the Fig. 10d points (all five free, and
  with one or two held).
- **Other solver runs.** The repeat of the heat-map script's solver runs
  for the rows of Fig. 11 other than 100 C (the tests repeat the 100 C
  row), and the independent scipy implementations quoted for the simplified
  model and for the Table S4 orders (every row of that table but the
  first). `solver_emulation.py` covers the deposited curves of Fig. 9 and
  Fig. 10a-c only. Its Jacobian is a plain forward difference; MATLAB's
  adapts its increments. Other increments, tried outside the package, moved
  the values by less than 1e-6. The margin of the node-time limits over
  rounding that the solver-emulation test module quotes (random relative
  errors of 1e-15 and 1e-14 added to the rates) was measured outside the
  tests too. MATLAB itself was never run.
- **Digitization.** The digitized curves that are shipped
  (`reference_data/`) were read from the published figures as described in
  `reference_data/README.md`; the colour field of Fig. 11 was read too and
  is not shipped.

What the package and its tests do establish: that the shipped models, with
the constants as documented, agree with the deposited curves, the digitized
lines and the experiments to the numbers given in the "Agreement" tables;
that, with the authors' solver run repeated, they return the deposited
model nodes of Fig. 9a, Fig. 9b, the H2O2 trace of SI Fig. S12a and
Fig. 10a-c to the digits the files print; what each inferred choice
(`x_O2_ref`, the O2 state, `k_7` without the catalyst mass, the
normalization without OPA, the water stoichiometry and the bilinear reverse
of R1 of the reversible variant) contributes to that with the shipped
constants; and that each model implements the rate laws written in these
notes.
