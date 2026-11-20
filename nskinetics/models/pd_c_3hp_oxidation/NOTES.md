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
was **fitted** here (one constant, `k_3` of the simplified reversible
model). They also list the places where the paper and the deposit disagree
and that matter for these models.

Two conventions are used throughout:

- How the authors produced a figure is known only where the paper or the
  deposit says so. Elsewhere these notes say what the numbers show ("is
  reproduced by", "agrees with") and not what the authors did.
- Several results quoted here come from analyses that are **not in the
  repository**: fits of constants to the deposited curves and to digitized
  contour lines, refits of alternative models, an emulation of the authors'
  random fitting procedure, and readings of the colour field of a published
  figure. Each such result is marked "(not shipped)". The last section
  lists them. Everything else can be recomputed with the package; the main
  results are pinned by its tests. In particular the repeat of the authors'
  solver run is shipped (`solver_emulation.py`): that the shipped models
  return the deposited model nodes to their printed digits is tested, and
  so is what each inferred choice contributes to that. The fit that
  defines one constant of the simplified reversible model is shipped as
  well (`fitting.py`), and the tests repeat it.

## The three models

The paper presents one kinetic model per oxidant, fitted at 50 C. Its
deposit holds a third model, that of the script with which the temperature
maps of Fig. 11 were computed. The package ships these three:

| Model (loader) | In the paper? | Compared with |
|---|---|---|
| H2O2 model (`load_h2o2_model`) | Table 3, H2O2 column | Fig. 9a (50 C), which it reproduces. Its predicted H2O2 is compared with SI Fig. S12a, whose H2O2 curve lies up to 0.38 equivalents lower and is matched only with other constants |
| O2 model (`load_o2_model`) | Table 3, O2 column, with one inferred factor | Fig. 9b and SI Fig. S12b (50 C), which it reproduces. It has no temperature parameter |
| Simplified reversible O2 model (`load_o2_simplified_reversible_model`; the name is this package's) | no; in the authors' deposited heat-map script | Fig. 11 (rms offset of the contour lines 1.05 px in panel a and 0.30 px in panel b); the 50 C experiment; the experiments at 25, 75 and 100 C, which it does not describe |

The two Table 3 models are the published ones, at 50 C. The only model of
the package with a temperature is the simplified reversible model, and its
temperature dependence is that of the script.

The constants that the script's run fitted for the simplified reversible
model are not in the deposit. The ten 50 C constants of the shipped model
are defined by one rule: **every constant the deposit gives is used as
given, the full-precision value where it gives one twice, and the one it
does not give is fitted to the 50 C data with the authors' fitting
set-up.** Nine are given: three rate constants from the 50 C row of the
deposited Arrhenius table, and `k_4` = 0.0023 and the five equilibrium
constants from the vector `k0` written in the script. `k_3`, which `k0`
writes as 0.0000, is fitted, by a function that is shipped ("Constants at
50 C" below). No constant of the package is fitted to a published figure.
Fig. 11 is therefore a comparison, not a fit: that the model comes close
to it shows that the constants are consistent with the published figure,
not that they are the authors'.

## Files

| File | Content |
|---|---|
| `pd_c_3hp_oxidation_h2o2_antimony.txt` (+ `_sbml.xml`) | H2O2 model: Table 3 H2O2 column, Fig. 9a |
| `pd_c_3hp_oxidation_o2_antimony.txt` (+ `_sbml.xml`) | O2 model: Table 3 O2 column, Fig. 9b, SI Fig. S12b |
| `pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt` (+ `_sbml.xml`) | the simplified reversible O2 model (the name is this package's; no OPA): the model of the authors' heat-map script (Fig. 11); not in the paper |
| `pd_c_3hp_oxidation.py` | loaders (`load_h2o2_model`, `load_o2_model`, `load_o2_simplified_reversible_model`), `simulate_at_times`, `TABLE_S4`, the constant sets of the simplified model, and the constant sets inferred from the deposited Fig. 9 nodes (`H2O2_NODE_PARAMETER_SETS`, `O2_NODE_PARAMETER_SETS`) |
| `fitting.py` | `fit_simplified_reversible_constants`: the least-squares fit to the 50 C data of the script that defines `k_3` of the simplified model |
| `solver_emulation.py` | a port of MATLAB `ode23s` and the repeat of the authors' solver run with the rates of a shipped model, for the comparison with the deposited model nodes of Fig. 9 (verification only) |
| `heat_maps.py` | the temperature-time maps of Fig. 11 and their comparison with digitized contour lines |
| `temperature_series.py` | the lines of the Arrhenius plot (Fig. 10d) and the measures used to compare a model with a table of experimental points |
| `reference.py`, `reference_data/` | deposited Fig. 9 and Fig. 10 tables, the 50 C data table of the heat-map script, digitized Fig. S12 inset curves and Fig. 11 lines, and their readers |
| `figures.py` | draws Fig. 9, SI Fig. S12 and Fig. 11 from the models, and the verification figures |

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
are the five figures of `figures.py`, written from the repository root by

    python -m nskinetics.models.pd_c_3hp_oxidation.figures \
        docs/source/_static/images/pd_c_3hp_oxidation --formats png --dpi 100

(the default resolution of the figure CLI is 600 dpi; the comment at the top
of the page gives the same command). They show the models with the deposited
and digitized data overlaid; no image of the paper is in the repository.

## Model basis

This section and the next six are mainly about the two Table 3 models.
The simplified reversible model has its own sections further down.

Conditions of both models: 60 mM 3-HP, NaOH/3-HP = 1:1, 0.2 g Pd/C, 50 C;
3 bar O2, or 15 equivalents of H2O2 under 3 bar Ar (Fig. 9 caption).

The authors' scripts work in liquid mole fractions and plot fractions of the
initial 3-HP. Here the liquid species are concentrations in M:
`c_i = C0 * n_i`, with `n_i` the moles of species i per initial mole of 3-HP
and `C0 = 0.06` M. The mapping is exact because every rate law is first order
in the liquid species. Time is in min. The two Table 3 models are at 50 C
and have no temperature; that of the simplified model is in K.

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
| all three | `C0` | M | 0.06 | initial 3-HP concentration (takes effect at the next `reset()`) |
| H2O2 | `eq_H2O2_0` | - | 15 | initial H2O2 in equivalents of `C0` (0.9 M; takes effect at the next `reset()`) |
| simplified | `T` | K | 323.15 | temperature. The model was fitted at 50 C only and does not describe the experiments at 25, 75 and 100 C (see "Limits of the temperature dependence") |
| simplified | `P` | bar | 3.5 | pressure of the heat-map script. Not for exploring pressure: the paper's Fig. 3a contradicts the model above 3 bar (see below) |

Rate constants: `k_1` to `k_6` of the Table 3 models are in 1/(g_cat min)
and `k_7` (H2O2) in 1/min. The rate constants of the simplified model are
lumped, in 1/min for the paper's batch; that model has no catalyst-mass
parameter.

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
`C0` or `eq_H2O2_0` is an extrapolation of the fitted rate laws. The
paper has experiments at other H2O2 charges, temperatures and O2 pressures;
the models were not built on them and do not describe them (next
subsection). The Table 3 models have no temperature parameter. That of the
simplified model extrapolates from 50 C, and the comparison with the time
courses of Fig. 10 is poor (section "Limits of the temperature
dependence"). For `m_cat` and `C0` there is
nothing to compare with: no experiment found in the paper or the deposit
varies the 3-HP concentration or the mass of the 10 % Pd/C (where captions
state them they are 60 mM and 0.2 g; SI Table S1 adjusts the amount of the
3 % and 5 % catalysts to the same metal-to-feed ratio). The simplified model
is fitted, by its script and here, to the 50 C experiment only (samples up
to 1680 min).

### Experiments of the paper that the models were not built on

The time courses behind the three models are those of Fig. 9 (50 C); the O2
one is at 3 bar of O2 (caption of Fig. 9), and the heat-map script uses
3.5 bar. The paper also has time courses with O2 at 25, 75 and 100 C
(Fig. 10a-c, 3 bar), which are compared with the simplified model in the
section "Limits of the temperature dependence", and it reports experiments
in which one condition is varied. Most are
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
| Temperature with O2, Fig. 3c | 25 C | 8.5 / 1.4 | simplified model at T: 53.5 / 36.3 |
| | 50 C | 72.8 / 41.1 | 98.9 / 49.1 (O2 model: 99.0 / 48.0) |
| | 75 C | 95.5 / 35.7 | 99.9 / 6.0 |
| | 100 C | 84.6 / 2.1 | 100.0 / 0.2 |
| H2O2 charge, Fig. 2a | 5 equivalents | 14.0 / 4.7 | H2O2 model, `eq_H2O2_0` set: 41.8 / 9.7 |
| | 10 | 38.1 / 20.3 | 54.6 / 24.3 |
| | 15 | 78.4 / 58.5 | 64.7 / 36.5 |
| | fourth bar (see below) | 56.4 / 28.4 | 72.6 / 45.5 (20 equivalents); 83.6 / 55.7 (30) |
| | fifth bar | 20.7 / 6.6 | 83.6 / 55.7 (30 equivalents); 94.2 / 58.2 (50) |

The simplified model is as its loader builds it (default constants,
heat-map start); with the script's fitting start the pressure column reads
85.5, 95.9, 97.9, 98.9 and 99.5 % conversion and 27.3, 44.5, 49.0, 51.1
and 49.8 % MA, and the temperature column 50.3, 96.7, 99.9 and 100.0 %
conversion and 33.0, 46.2, 5.8 and 0.1 % MA. The last two bars of Fig. 2a
are labelled 20 and 30 on the axis of the published figure and `1: 30` and
`1: 50` in the deposit; the text speaks of ratios "from 1:5 to 1:30". Both
readings are given.

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
  pressure**; leave it at the script's 3.5 bar. The O2 model has no
  pressure parameter. `x_O2_ref` of the O2 model
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
- **Temperature.** At 25 C the simplified model, the model of Fig. 11,
  converts 53 % in 8 h where 8.5 % was measured, as in the comparison with
  the 25 C time course. At 75 C it has 6 % MA left after 8 h where the
  screening experiment has a yield of 35.7 %.
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
  | O2 model started from 60 mM MA, 120 min | 13.7 % | no 100 C model |
  | O2 model, 1300 min | 79.7 % | no 100 C model |

  (Fig. 4b: 100 minus the deposited MA percentage; the 50 C series reads
  101.4 and 100.5 % MA at 1140 and 1230 min.) The published panels show the
  values of their own files, as far as the image can be read. The caption
  and both files give 1300 min; the text quotes the Fig. 4a values, 3.8 and
  89.8 %, as "after 24 h", and the 30.1 % of Fig. 4b. The model values are
  tested. The model is contradicted by the 50 C series (79.7 % against 1.2
  or 3.8 % at 1300 min; 13.7 against 0.8 % at 120 min). Its steps R3, R4
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
`eq_H2O2_0` and `T` are in the models because the rate laws, or the
authors' script, contain them, not because the response to them was
established.

What does and does not depend on `C0`, measured with the package:

- O2 model: the plotted fractions `X_i` hardly depend on `C0`, because
  every rate is first order in its liquid species. Halving or doubling `C0`
  changes them by less than 1e-4 (through the O2 depletion, which is
  proportional to `C0`).
- H2O2 model: the initial H2O2 concentration is `eq_H2O2_0*C0`, so changing
  `C0` alone changes the H2O2 concentration and with it the rates of R1, R2
  and R4 (`X_MA` at 600 min is 0.201, 0.397 and 0.564 for `C0` = 0.03, 0.06
  and 0.12 M at 15 equivalents).
- Simplified model: `x_per_M` holds the script's initial 3-HP mole fraction
  divided by `C0`; if `C0` is changed, `x_per_M` has to be changed with it.

`set_value` does not validate its input. Zero or negative `m_cat`, `C0` or
`eq_H2O2_0`, or a `T` set in degrees Celsius, give meaningless results
without an error. Only the loader arguments are checked (`T` of the
simplified model must lie within 250 to 450 K; unknown constant-set names
are refused).

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

| Term of Fig. S2 | H2O2 model | O2 model | Simplified reversible model |
|---|---|---|---|
| `r_j = k_jf*prod(a_i^n)`, irreversible | as Fig. S2 | as Fig. S2 | forward minus reverse term in every step |
| orders `n` (fitted, 0 to 3; Table S4) | all exactly 1 (inferred) | all exactly 1 (inferred) | all 1 |
| liquid activity `gamma*x` | `x_per_M*c`, `gamma = 1`: the authors' scale, 0.000339 for the 60 mM feed | the same | the script's mole fractions `x_per_M*c` (0.000339 scale in its fits, unit scale in its heat-map run) |
| oxidant activity | H2O2 is a liquid species: `x_per_M*[H2O2]`, 0.005085 at t = 0. As Fig. S2 | **departs**: Fig. S2 gives `y*P/f0` = 3/1.01325 = 2.96 for 3 bar of O2; the model has `x_O2_ref*yO2` = 0.005085*`yO2` | `yO2*P/f0` with `P` = 3.5 bar, `f0` = 1.01325 bar. As Fig. S2, with the script's pressure |
| `k_jf` | Table 3, times `m_cat` (not `k_7`) | Table 3, times `m_cat` | effective constants, 1/min |
| `Keq >> 1` | irreversible, as Fig. S2 | irreversible, as Fig. S2 | does not hold: "equilibrium constants" `K_6`..`K_10` of 0.6 to 82 |
| Eq. S1 | liquid species as concentrations, `c = C0*n` | the same; `yO2` loses one unit per unit of liquid mole fraction reacted in R1, R2, R4 | the script's states; `yH2' = r_3`, `yCO2' = r_2 + r_3` |

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
| `k_1`..`k_7` (H2O2), `k_1`..`k_6` (O2) | Table 3 as printed | printed |
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

The table is for the two Table 3 models. The numbers of the simplified
reversible model (its ten 50 C constants, five activation energies, `P`,
`f0`, `T_ref` and `R_gas`) are traced in its section ("Constants at 50 C"
and "Structure (transcribed from the script)").

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
agree with them only to their solver error (previous section). To compare
more closely, `solver_emulation.py` repeats the run:

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
- `emulate_fig9a` and `emulate_fig9b` run it for the deposited curves.

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

The files print their node times and values to five decimals (small values
to three significant digits). The differences are therefore at the print
rounding of the files (5e-6). The tests require the deposited number of
steps, node times within 2e-5 min and values within 1e-5. For comparison,
the same models solved accurately differ from these emulated steps by
8.5e-4 (Fig. 9a), 2.0e-3 (Fig. 9b) and 0.012 equivalents (H2O2 trace). That
is the solver error of the deposited curves, which the accurate-solver
tests have to allow for; for the plotted fractions the limit of 1e-5 is 85
to 200 times smaller.

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
  the first deposited node of each curve (5.87666 min for Fig. 9b,
  12.547 min for Fig. 9a);
- the last step, which ends at 1.1 times the last sampling time.

No step of the three runs is rejected, and none is accepted narrowly (the
error estimate of every step is below 0.6 of the tolerance), so rounding
differences between machines cannot change the sequence of steps.

Three rules of the solver are not tested by these curves, and the package
tests check them on other problems:

- the rule for a rejected step, on `y' = cos(10 t)`;
- the limit of five on the growth of a step from one to the next, on
  `y' = -100 y`. In the three runs an error estimate small enough for it
  occurs only in the last step and, in the run of the H2O2 trace, in the
  step before it, which already has the largest size, a tenth of the run,
  so that the limit does not show. The largest growth of a step in the
  three runs is 2.18-fold (Fig. 9b);
- the stretching of a step that would end within 10 % of its length before
  the end of the run, on `y' = -y`. The last step of each of the three runs
  is cut short at the end of the run instead.

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
  given to eight significant digits. These are changes of one constant at a
  time; several constants changed together can compensate each other in
  part.
- Each negative control of the tests changes one inferred choice and keeps
  the constants: `x_O2_ref` (section "`x_O2_ref`"), the O2 state, `k_7`
  and the catalyst mass, and the normalization. A control of this kind
  shows that the shipped constants need the shipped
  choice. Whether another choice could do as well with other constants
  takes a refit; where one was made, it is not shipped, and the text says
  so.

The controls for the two Table 3 models (limit 1e-5; node times 2e-5 min):

| Change | Effect on the repeated run |
|---|---|
| `x_O2_ref` = 0.0051 (Fig. 9b) | misses the nodes by 9.4e-4, node times by 1.8 min; see "`x_O2_ref`" |
| O2 fraction held at 1 in the rates, as in the posted script (Fig. 9b) | 7.7e-5, node times 0.15 min |
| `k_7` multiplied by `m_cat` (Fig. 9a) | 16 steps instead of 23; 0.39 |
| OPA included in the sum of the normalization | same steps; 0.038 (Fig. 9b), 0.094 (Fig. 9a) |
| H2O2 trace run with the `'fig9a_nodes'` set or with Table 3 | 23 steps instead of 24, node times more than 50 min off |

On the O2 fraction: the posted script integrates an O2 state and uses
`yO2 = 1.0` in its rates. The control shows that the constants of the node
set need the depleting state. In a refit of the six constants with the
fraction held at 1, on the deposited time grid (not shipped), the largest
difference from the nodes stayed at 1.5e-5, three times that of the model
with the depleting state (5e-6). The evidence for the coupling is therefore
this factor of three in a fit, not the control alone.

`figures.plot_solver_emulation` draws the three deposited curves with the
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
   tolerances (previous section; MATLAB itself was not run). Section 2.6
   names MATLAB R2024b, the README of the deposit R2024a.
10. **Maximum MA yield at 75 C.** The text gives maximum MA yields of 52.6 %
    (75 C) and 56.8 % (100 C). The largest measured MA fractions in the
    deposit are 0.436 (75 C, 180 min) and 0.568 (100 C, 40 min); the 75 C
    value of the text is not in the deposited data. (At 25 C the text's
    2.3 % after 31.5 h is the deposited 0.023 at 1890 min.)
11. **Time axis of Fig. 10a-c.** The published panels end at 1440 min. The
    samples at 1740 and 1890 min (25 C) and 1560 min (100 C) are in the
    deposit but not in the figure. They are counted in the pooled R2 of the
    comparison with the simplified model below ("Limits of the temperature
    dependence").
12. **Activation energies.** Two of the three activation energies of the
    text are not the slopes of the lines of Fig. 10d, from which the text
    says they were obtained (section "Arrhenius plot (Fig. 10d)").

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
(`load_o2_simplified_reversible_model`) because it is the model the
authors deposited for the temperature dependence, with its rate laws, its
data and its fitting set-up. The constants it fitted are not deposited;
"Constants at 50 C" says how the shipped ones are defined.

**Equation numbering.** The paper numbers two displayed equations (6). The
first, in Section 3.4 (page 7), is the initial-rate ratio of item 8 of
"Paper versus deposit: discrepancies that affect these models". The second,
printed at the end of Section 3.4 (page 9, before the Conclusions), is the
Arrhenius relation `k(T) = k_ref*exp[-(Ea/R)*(1/T - 1/T_ref)]`; the text's
"Arrhenius relationship (Eq. 6)" refers to this one. (Fig. 6 of the paper
is unrelated: XPS profiles and a reusability test.) Wherever this package
says "the Arrhenius relation of the paper", the second Eq. 6 is meant.

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
compared with an independent scipy implementation of the script's
equations (not shipped; same constants; LSODA and Radau at rtol 1e-12) in
both conventions below, for the two constant sets at 30, 50, 84, 95 and
100 C: at the loader's CVODE tolerances (rtol 1e-9, atol 1e-10) the liquid
fractions agree to 2e-8 of the feed and the gas fractions to 2e-9.

Temperature: `p(T) = p_ref*exp(-(Ea/R_gas)*(1/T - 1/T_ref))` for all ten
parameters, `T_ref` = 323.15 K and `R_gas` = 8.314 J/(mol K), with `Ea` =
59781.2, 91594.6, 10000, 51788.73 and 59549.8 J/mol for steps 1 to 5 (all
from the script). The first, second and fifth are printed to three digits
in Section 3.4 of the paper (59.8, 91.6 and 59.5 kJ/mol, for "MA
production", "MA conversion to AAc" and the "decomposition reaction"); the
text says they were obtained from the Arrhenius plot of Fig. 10d, whose
lines give 59.93, 68.96 and 52.09 kJ/mol (section "Arrhenius plot
(Fig. 10d)"). The third (MA -> OXA, a round value) and the fourth
(3-HP -> D) are in the script only.

### Three properties of the script that are in the published figure

Each was tested against the contour lines digitized from Fig. 11b by removing
it (default constants; the rms offset is 0.30 px with all three in place,
1 px = 3.07 min or 0.148 C). The offsets without refitting were measured
with the package; the tests pin those of the second and third. Offsets
above about 1 px are not distances: they are the linearized estimate of
`contour_offsets` (value difference over gradient), clipped at 30 px per
sample (`max_offset=30`), and are indicative only. The values "after
refitting" are from fits of the five rate constants to the contour lines of
Fig. 11b with the property removed (not shipped).

1. **The Arrhenius factor of step i multiplies `k_i` and also `K_(i+5)`**
   (the script's `Ea` vector is the five values repeated), so `k_i/K_(i+5)`
   does not change with temperature. Scaling only the rate constants gives
   2.8 px (1.3 px after refitting).
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
   Starting at plotted time 0 gives 3.8 px (6.6 px without the clip; 2.9 px
   after refitting).

Two further points: the script uses a pressure of 3.5 bar (the captions give
3 bar of O2), and what Fig. 11 calls "MA selectivity" and "MA yield" are the
water-free liquid fraction of MA, `X_MA`, and `(1 - X_HP3)*X_MA` (model
quantities `fig11a_selectivity` and `fig11b_yield`).

### Constants at 50 C

**What the deposit holds.** The constants the script's run fitted are not
in the deposit. The script writes one vector of ten numbers, to four
decimals, `k0 = [0.0021 0.0001 0.0000 0.0023 0.0010 4.7840 0.5918 3.5116
81.6870 52.0472]` (five rate constants, then five equilibrium constants;
`'k0_printed'`). In the script `k0` only sets the size and the bounds of
the global search: the particle swarm searches 0 to 1.5 times `k0`. For
`k_3` that interval is 0 to 0, so the swarm can only return `k_3` = 0;
`lsqcurvefit` then starts from the swarm's result with a lower bound of 0
and no upper bound and stops after at most five iterations. What `k_3` the
authors' run ended with therefore does not follow from the script text,
and whether `k0` holds the result of an earlier run is not stated. Run as
the heat-map model, `k0` is far from the figure (4.6 and 7.1 px on panels a
and b). The deposit also holds the Arrhenius table of Fig. 10d
(`11-4_Arrhenius.csv`), whose 50 C row gives the logarithms of three rate
constants, taken to be `k_1`, `k_2` and `k_5` of this model (section
"Arrhenius plot (Fig. 10d)").

**The rule.** Every constant the deposit gives is used as given. Where the
deposit gives a constant twice, the full-precision value is used: `k_1`,
`k_2` and `k_5` are exp() of the 50 C row of `11-4_Arrhenius.csv`, which
`k0` gives to four decimals as 0.0021, 0.0001 and 0.0010 (the last
disagrees with the table's 1.0625e-3, which would round to 0.0011).
`K_6` to `K_10` are those of `k0`, and so is `k_4` = 0.0023, its only
source (two significant digits). `k_3`, which `k0` writes as 0.0000 (below
5e-5), is the only constant fitted: to the 50 C data table of the script,
with the script's fitting convention. Rate constants in 1/min; they are
lumped (no catalyst mass).

| Constant | Value (`'deposit'`, the default) | Source |
|---|---|---|
| `k_1` | 2.139999e-3 | exp(-6.14695): 50 C row of the deposited `11-4_Arrhenius.csv`, "ln k1" (`k0`: 0.0021) |
| `k_2` | 5.199984e-5 | exp(-9.86427): same row, "ln k2" (`k0`: 0.0001) |
| `k_3` | 1.4215e-5 | **fitted** (below); standard error 4.4e-6 (31 %); `k0`: 0.0000 |
| `k_4` | 0.0023 | `k0`, its only source |
| `k_5` | 1.062501e-3 | exp(-6.84713): same row, the column labelled "ln k3" (`k0`: 0.0010) |
| `K_6`, `K_7`, `K_8`, `K_9`, `K_10` | 4.7840, 0.5918, 3.5116, 81.6870, 52.0472 | as written in the script (`k0`) |

In the model these 50 C constants are the parameters `k_1_ref` to
`K_10_ref`; `k_1` to `K_10` are their values at `T` (the Arrhenius
assignment rules). `fit_simplified_reversible_constants` accepts either
name and returns the `_ref` names.

**The fit.** `fit_simplified_reversible_constants` (`fitting.py`) fits
chosen constants of the model to the 50 C data by least squares, the others
held; by default `k_3` alone. Its set-up is that of the script's own fit:

- *Data:* the table of 16 samples (0 to 1680 min, five fractions each) that
  the script carries in its function `data()`, shipped as
  `reference_data/heatmap_script_50C_data.csv` (`read_heatmap_script_data`).
  It is the experimental block of the deposited `10-2_kinetic_model_O2.csv`
  (Fig. 9b) with one difference: the script gives the last two rows (1470
  and 1680 min) to four decimals where the file has five, of which they
  are the rounded values (largest difference 5e-5). The other 14 rows are
  equal (tested).
- *Initial condition:* the script's fitting convention, an initial 3-HP
  "mole fraction" of 0.000339 (previous subsection).
- *Residuals:* the five water-free liquid fractions of the model minus
  those of the table at the 16 sample times, all 80 values, unweighted (the
  five of the t = 0 row are zero).

With the other nine constants held at the values above, the fit gives
`k_3` = 1.42150e-5 1/min, shipped to five significant digits, with a sum
of squares of 0.033797 (pooled R2 0.9939, rms 0.0206). The same optimum is
reached from other starting points to 1e-5 of the value (7e-6 at most over
eleven starts from 1e-6 to 1e-3; three are tested), and an independent
scipy implementation of the script's equations (not shipped; Radau and
LSODA, scalar minimization and least squares) gives 1.42150e-5 as well.
The test suite checks that the shipped `k_3` is the rounded result of the
function. The result is at least 3.4e-5 (relative) from a rounding
boundary of the fifth digit; tolerances of the integration and the step of
the central-difference Jacobian move it by less than 1e-5, whereas scipy's
default two-point Jacobian stops the iteration early: 2.8e-4 lower
(1.4211e-5) from the default start, and from the starts 1e-6, 4.3e-6,
3e-5 and 1e-4 between 1.1e-3 lower and 3.6e-3 higher. A two-point Jacobian
with a step of 1e-6 comes closer but still crosses the rounding boundary
from the start 3e-5 (1.42145e-5, -3.8e-5). The standard error (from the
Jacobian at the optimum, 79 degrees of freedom; the script divides by
twice the number of residuals minus ten) says how well the data determine
`k_3` with the others held: to 31 % only. The product of step 3, OXA, is at
most 3.2 % of the liquid in the data, so they hold little information on
`k_3`.

**Alternative: `k_4` fitted as well** (`free=('k_3', 'k_4')`; tested, not
the shipped set). With the other eight held, the fit gives `k_3` =
1.42000e-5 and `k_4` = 2.32041e-3 1/min (standard errors 31 % and 3.7 %)
and a sum of squares of 0.0337728, 0.07 % below that of the shipped set.
The fitted `k_4` agrees with the 0.0023 of `k0` to its precision (it
rounds to it), and `k_3` changes by 0.1 %. With these two values the heat
maps lie 1.11 px (panel a) and 0.35 px (panel b) from the published
contours, with maxima of 0.592 and 0.548 (against 1.05 and 0.30 px and
0.593 and 0.548 with the shipped set).

**What is not the authors' procedure.** The script obtains its constants in
two stages, all ten constants free: a particle swarm (`particleswarm`, 100
particles, three iterations, bounds 0 to 1.5 times `k0`) and, from its
result, a least-squares refinement (`lsqcurvefit`, lower bound 0) with
options that set a maximum of five iterations. Every trajectory is integrated
with `ode23s` at default tolerances, restarted at each sample. The swarm is
random, and neither stage is run to convergence, so the authors' constants
cannot be recomputed from the deposit: another run of the same script gives
other constants. The fit of this package differs in four ways: it is
carried to convergence; its integration is accurate (CVODE, relative
tolerance 1e-12) and runs through all sample times at once; only the
constant without a deposited value is free; and it is deterministic.

**Why the given values are held, and what supports them.** The numbers
below were computed with the package unless marked "(not shipped)".

1. *A converged fit of all ten constants does not determine the equilibrium
   constants.* With all ten free (the same function), the five rate
   constants come out the same from every starting point tried (2.1310e-3,
   5.434e-5, 1.4636e-5, 2.3931e-3 and 9.747e-4 1/min; sum of squares
   0.03195) and so does `K_6` (4.44), but `K_7` and `K_8` end wherever the
   start leads (for example 24 to 6800 and 0.001 to 1.6e5 over ten starting
   points, and 2.1e3 and 2.7e5 from the function's default start; these are
   examples, not bounds), and `K_9` and `K_10` run to 1e5 and beyond, where
   the reverse terms of steps 4 and 5 vanish. The heat maps of these fits
   lie 2.7 to 2.8 px (panel a) and 1.8 to 1.9 px (panel b) from the
   published contours, with maxima of 0.596 and 0.554. The data alone
   therefore do not give the equilibrium constants, and the published
   contours are matched better by the equilibrium constants of `k0` with
   the rate constants that go with them (1.05 and 0.30 px, item 3) than by
   the converged fit. (Measured with the package; not pinned by a test.)
2. *With the equilibrium constants of `k0`, the data return the 50 C row of
   the Arrhenius table and the `k_4` of `k0`.* With `K_6` to `K_10` held as
   written and all five rate constants fitted, `k_1`, `k_2` and `k_5` come
   out as 2.1419e-3, 5.277e-5 and 1.0664e-3 1/min: 0.09, 1.5 and 0.37 %
   above the values of the row (2.1400e-3, 5.200e-5, 1.0625e-3), well
   inside their standard errors (2.5, 31 and 3.9 %). `k_4` comes out as
   2.3122e-3, 0.53 % above the 0.0023 of `k0` (standard error 5.4 %), and
   `k_3` as 1.4165e-5, 0.35 % below the shipped value; the sum of squares
   falls from 0.033797 to 0.033768 (tested). This is the support for
   reading that row as constants of this model: they are, within the
   accuracy of the data, what this model with these equilibrium constants
   gives when it is fitted to the script's table. It is an inference;
   neither the paper nor the deposit says so.
3. *The resulting heat maps are close to the published ones* (next
   subsection): rms offsets of 1.05 px (panel a) and 0.30 px (panel b) from
   the digitized contour lines, and maxima of 0.593 (95 C, 27 min) and
   0.548 (84 C, 49 min), where the text gives 59.1 % and 55.2 % for the
   same cells.
4. *Fits to the published contour lines agree with the data fit* (not
   shipped). With the other nine constants held, `k_3` fitted to the
   contour lines digitized from Fig. 11b comes out as 1.449e-5, 1.9 %
   above the value fitted to the data (rms offset 0.28 px in panel b,
   1.08 px in panel a); fitted to the lines of both panels together, as
   1.418e-5, 0.3 % below it. With `k_3` and `k_4` both fitted to the panel
   b lines (the other eight held), they come out as 1.433e-5 and
   2.305e-3 1/min, 0.8 % above the shipped `k_3` and 0.2 % above the `k_4`
   of `k0` (0.27 px in panel b, 1.07 px in panel a). With all five rate
   constants fitted to those lines, `k_1`, `k_2` and `k_5` come out 0.03 %
   above, 2.2 % below and 0.10 % below the 50 C row (rms offset 0.25 px).
   The figure was not used to define any shipped constant.
5. *The authors' procedure, emulated, scatters by more than these
   differences* (not shipped, and approximate: a plain global-best particle
   swarm and a bounded trust-region least squares stopped after five
   iterations stand in for MATLAB's `particleswarm` and `lsqcurvefit`,
   whose details differ, and the integration is accurate). Over 20 random
   seeds the sum of squares ended between 0.0320 and 0.0331, the rms
   offsets from the published contours between 1.7 and 3.9 px (panel a) and
   0.6 and 4.8 px (panel b), the largest "selectivity" between 0.592 and
   0.598 and the largest "yield" between 0.548 and 0.556; `K_9` and `K_10`
   ended between 137 and 9400 and between 48 and 3800. One run of the
   script is one draw from such a spread, which is why its constants cannot
   be recovered from the deposit.
6. *The solver of the script does not matter for the fit* (not shipped).
   With the trajectories integrated as the script integrates them (the
   `ode23s` port of `solver_emulation.py`, restarted at each sample, default
   tolerances), they differ from the accurate ones by at most 3.0e-5 at the
   samples, and the fitted `k_3` changes by +0.012 %. Fitting to the
   experimental block of `10-2_kinetic_model_O2.csv` instead of the
   script's table changes it by -0.015 % (tested).

How much the equilibrium constants matter for the heat maps: `K_6` has no
effect on them, and a 10 % change of `K_8`, `K_7`, `K_9` or `K_10` moves
them by at most 9e-8, 5e-5, 1.0e-3 and 1.1e-3 (measured with the package).

**The two published panels cannot come from one set of constants.** Since
yield <= selectivity everywhere, the late side of a yield contour must lie
earlier than that of the selectivity contour of the same level: in the
model it lies 0.005 to 0.09 % earlier for the levels 0.1 and 0.3 above
62 C (tested, row by row). In the published figure the late 0.1 and 0.3
contours of panel b lie later than those of panel a above 62 C, by 2.2 % in
the median and up to 4.4 % (a median of 1.2 to 2.4 px, up to 4.6 px; panel
a interpolated to the temperatures of the panel b samples; tested on the
digitized lines). No single run can produce both panels. The shipped set
is 0.30 px from panel b and 1.05 px from panel a. Fits to the contour
lines (not shipped) show what each panel would take: panel a alone is
matched to 0.17 px with `k_3`, `k_4`, `k_5` and `K_10` free (`K_10` then
falls to 23), which leaves panel b 1.45 px off; with `k_3` free and both
panels fitted together the offsets are 1.05 and 0.31 px, with `k_3` and
`k_4` free 0.98 and 0.37 px, and with all ten constants free still 0.60
and 0.48 px. Why the panels differ (separate runs of the script, or other
settings) is not known. The package ships one set, defined without the
figure, and reports the misfit of panel a
(`figures.plot_fig11_verification` draws it).

### Agreement with the published Fig. 11

There is no table of the heat maps in the deposit. The reference is the
published figure: the contour lines 0.1, 0.3 and 0.5 of the two main panels,
digitized (`reference_data/`). Offsets are normal to the lines, in pixels of
the published raster (1 px = 3.07 min or 0.148 C); the digitization itself
is good to about 0.3 px. Both panels are computed with the one default
constant set (`fig11_heat_maps`). Measured with the package (CVODE, relative
tolerance 1e-9) and pinned by its tests:

| Comparison | Panel a | Panel b |
|---|---|---|
| rms offset of the contour samples | 1.05 px | 0.30 px |
| largest rms offset of one contour branch | 1.28 px (0.5, late) | 0.39 px (0.5, late) |
| largest mean offset of one contour branch | 0.69 px (0.3, late) | 0.23 px (0.5, late) |
| largest offset of a single sample | 3.0 px (0.5, late, 31 C) | 1.6 px (0.5, late, 30 C) |
| largest cell of the map | 0.5931 at 95 C, 26.76 min | 0.5484 at 84 C, 49.10 min |
| maximum of the 50 C row | 0.5713 at 261.3 min | 0.5333 at 322.8 min |
| red dashed line minus the digitized one, along time: mean / rms | +0.25 / 0.89 px | +0.29 / 1.64 px |

Panel b is matched to about the reading uncertainty of the digitized lines.
Panel a is not: its early contours are within 0.56 px rms, but its late
0.1, 0.3 and 0.5 contours are 1.26, 1.03 and 1.28 px rms off, in a pattern
that changes sign along each line (the lower panels of the verification
figure). An offset is positive where the computed map is above the level
of the contour at the sample: on a late branch the published line then
lies earlier than the computed one, on an early branch later.
This is the misfit that one set of constants cannot remove (previous
subsection).

The red dashed line joins the row maxima (`heat_map_ridge`), so it is
defined to one grid column, 5.585 min or 1.82 px; a line that sat one column
off would be 1.82 px off in the mean. The digitized line has 118 (a) and 226
(b) samples. The tests require the mean offset within half a column and the
rms within one column. The largest single offsets, 1.8 px (a) and 7.5 px
(b), are at 30 to 31 C, where the row maxima are flat and the line advances
several columns per row. In panel b the agreement holds over the whole line
and not for every part of it: the mean offset is -1.3 px between 45 and
60 C and +1.2 px between 60 and 80 C (0.7 and 0.6 of a column).

Read from the colours of the published maps (not shipped): the maps agree
with the colour field to an rms of 0.005 in both panels where the value is
above 0.3, and the 50 C rows of the published maps peak at 0.567 and 0.532
(+-0.01).

The script integrates each of the 299 grid intervals with `ode23s` at
default tolerances. Repeating those solver runs with the `ode23s` port of
`solver_emulation.py` and the script's equations changes the maps by at
most 1.1e-4 (panel a) and 1.4e-4 (panel b), so the solver does not matter
at the level of the figure. The difference grows with temperature and is
largest in the 100 C row (3e-5 at 84 C, below 1e-6 at 50 C). The package
tests repeat the 100 C row; the other rows were run outside the tests (not
shipped).

### Paper versus deposit for Fig. 11

1. **Model.** The text describes the heat maps as the O2 model extended with
   the Arrhenius relation. The caption says only that the Fig. 9b data gave
   `kref`, which is what the script does: it fits its own simplified
   reversible model to those data. The Table 3 O2 model of this package is
   not extended to other temperatures.
2. **Arrhenius scaling of the equilibrium constants**, with the activation
   energies of the rate constants, is not mentioned in the paper. Two of the
   five activation energies are not in the paper either (51.79 kJ/mol for
   3-HP -> D, and 10 kJ/mol for MA -> OXA, a round value).
3. **Time axis.** "27 min" and "49 min" in the text are the grid columns
   26.76 and 49.10 min of the plotted time, i.e. 16.8 and 39.1 min of
   reaction.
4. **Maximum selectivity, "59.1 % (95 C, 27 min)".** The model has 0.593 in
   that cell, which is the largest cell of its map.
5. **Maximum yield, "55.2 % (84 C, 49 min)".** The cell is the largest of the
   model's map as well, with a value of 0.548. The colour of the published
   map along its crest reads 0.549 (+-0.006) (not shipped). The converged
   fit of all ten constants to the 50 C data has 0.554 in that cell, with
   contour lines 1.9 px from the published ones (item 1 of "Constants at
   50 C"). The crest is flat: 13 cells between 76 and 90 C are within 0.001
   of the maximum, so the cell that comes out largest is sensitive to small
   changes; the tests check the value in the cell, not that it is the
   largest.
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
     0.569 at no sampling time. The 50 C row of the model's selectivity map
     peaks at 0.571 at 261 min, near the purple star (254 min), but that is
     a model value, not an observed one. The 8 h screening experiment at
     the same nominal condition has an MA selectivity of 56.4 % (Fig. 3a),
     which is not 56.9 either.

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
   (b) on the 50 C row; the 50 C row of the model peaks at 261.3 and
   322.8 min.
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

## Experiments at 25, 75 and 100 C (Fig. 10a-c)

The deposit holds, for the O2 experiments at 25, 75 and 100 C, the
experimental points and the model curves drawn in Fig. 10a-c
(`11-1_25C_Kinetic.csv`, `11-2_75C_Kinetic.csv`, `11-3_100C_Kinetic.csv`).
The three files are shipped verbatim, and `read_fig10_o2` returns both
blocks. The package uses the experimental points (section "Limits of the
temperature dependence").

**The model curves of Fig. 10a-c are not reproduced by this package.** The
paper prints no constants for 25, 75 or 100 C and does not say which model
the curves are. Neither posted script produces them as posted:
`APCATB_D_25_08770_kinetic_opt.m` is the irreversible six-reaction model of
Table 3 (posted with the 75 C data as its example), and
`APCATB_D_25_08770_heatmap.m` is the simplified reversible model, with the
50 C data. The deposited curves themselves exclude the first: in the model
block of `11-2_75C_Kinetic.csv` 3-HP falls from 1 to 0.00265 by 707 min and
then stays at 0.00245 from 1024 to 1584 min, and in that of
`11-3_100C_Kinetic.csv` it stays at 5.20e-5 from 677 to 1716 min. In the
posted O2 script 3-HP is consumed by two irreversible first-order steps
(`r(1) = k(1)*x3HP*yO2` with `yO2 = 1.0`, and `r(5) = k(5)*x3HP`), so it
falls by the same factor in equal times and cannot level off at a value
above zero; a reverse term, as in the heat-map script, can. The README of the
deposit says that the posted codes are "an example" and that "details were
adjusted for each time (H2O2, O2, parameters to fit, etc)"; how the curves
of Fig. 10a-c were computed is not stated. The published panels print
R2 = 0.99 for all three temperatures. They are not curves of the
simplified reversible model with its shipped constants, whose pooled R2
against the same points is 0.56 to 0.76 (section "Limits of the temperature
dependence").

What the package states about these experiments comes from the deposited
points (tested): the largest measured MA fraction is 0.023 at 1890 min
(25 C), 0.436 at 180 min (75 C) and 0.568 at 40 min (100 C); the published
panels end at 1440 min, and the samples at 1740 and 1890 min (25 C) and
1560 min (100 C) are in the deposit but not in the figure (items 10 and 11
of "Paper versus deposit: discrepancies that affect these models").

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
fitted to the experiment of each temperature, for three reasons:

- the figure labels name a one-step 3-HP -> MA reaction, which of the
  models in the paper and the deposit only that one has;
- fitted to the 50 C data of the script with its equilibrium constants as
  printed, that model returns the 50 C row: `k_1`, `k_2` and `k_5` come out
  0.09, 1.5 and 0.37 % above it, inside their standard errors (tested;
  item 2 of "Constants at 50 C"). Fitted to the contour lines of Fig. 11b
  instead, they come out 0.03 % above, 2.2 % below and 0.10 % below the row
  (not shipped; item 4 there);
- refits of that model to the experiments at 75 and 100 C come close to
  the other points (below).

Taking the points as constants of that model means that they follow the
rate law of its script. The three reasons support that and do not establish
it: neither the paper nor the deposit says it (previous paragraph). This
assignment is what makes the 50 C row the source of three constants of the
shipped model.

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
what the simplified reversible model gives when it is fitted to the 50 C
data of the script, whose time axis is in min; "Constants at 50 C" above.)

The constants whose logarithms are the deposited points (shortest decimals
that reproduce them within 5e-6), taken as 1/min:

| T (C) | `k_1` | `k_2` | `k_5` |
|---|---|---|---|
| 25 | 7.75e-5 | 6.64e-4 | 1.676e-3 |
| 50 | 2.14e-3 | 5.2e-5 | 1.0625e-3 |
| 75 | 2.43685e-3 | 4.44096e-4 | 5.4797e-3 |
| 100 | 1.00763e-2 | 1.59645e-3 | 1.41038e-2 |

The points cannot be recomputed: the script's fits are stochastic and
stopped after a few iterations. A converged fit of all ten constants of the
simplified model to the deposited experiment of each temperature
(`fit_simplified_reversible_constants` with the experiment as data;
measured with the package, not pinned by a test) gave an optimum that,
divided by the plotted value, was 1.000, 1.02, 0.998 (`k_1`, `k_2`, `k_5`)
at 100 C, 0.916, 1.085, 1.032 at 75 C and 0.996, 1.045, 0.917 at 50 C, and
was ill-determined at 25 C (different starting points ended at different
optima). This supports the assignment of the three columns and shows that
the points are not converged optima. The points are therefore taken from
the file (`fig10d_points`).

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
   those of the simplified reversible model ("What is plotted"). Which
   model the curves in panels a to c are is not stated either (section
   "Experiments at 25, 75 and 100 C (Fig. 10a-c)").
2. **Each line leaves out one point** (table above); the figure shows all
   four points without saying so.
3. **Activation energies.** The text says that "Apparent activation energies
   (Ea) for key steps were obtained from Arrhenius plot in Fig. 10d" and
   gives 59.8 (MA production), 91.6 (MA to AAc) and 59.5 kJ/mol
   (decomposition). The drawn lines give 59.93, 68.96 and 52.09 kJ/mol. No
   subset of two, three or four of the plotted points of MA -> AAc (`k_2`)
   gives 91.6 (nearest: 80.2 from 50 and 75 C), and none of MA -> Others
   (`k_5`) gives 59.5 (nearest: 61.4 from 50 and 75 C); for `k_1` the drawn
   line (59.93) is the nearest to 59.8. The
   text values are those of the heat-map script (59781.2, 91594.6 and
   59549.8 J/mol). Where they come from is not known; they are used as the
   script has them in the simplified reversible model (`Ea_1`, `Ea_2` and
   `Ea_5`), where they scale the 50 C constants (Fig. 11). The script's
   value for 3-HP -> D, 51788.73 J/mol (`Ea_4`), is not in the paper; it is
   close to what the 50 and 100 C points of the MA -> D constant give
   (51.85 kJ/mol).
4. **The model's constants against the plotted points.** Lines with the
   text's activation energies through the 50 C points (`fig10d_text_lines`)
   are the constants `k_1(T)`, `k_2(T)` and `k_5(T)` of the simplified
   reversible model with its shipped constants, and how the heat-map script
   extrapolates the three constants from 50 C (tested). They miss the
   other points: for `k_1` the 25, 75 and 100 C points lie 1.45, 1.47 and
   1.43 below the line (a factor of 4.2 to 4.3 in k), for `k_2` 5.41 above
   and 0.30 and 1.14 below, for `k_5` 2.31 and 0.05 above and 0.38 below.
   So the model that is extrapolated in Fig. 11 does not pass through the
   25, 75 and 100 C points of Fig. 10d: its `k_1` is 4.2 to 4.3 times the
   plotted value at all three temperatures (next section).

## Limits of the temperature dependence

The simplified reversible model is the only model of the package with a
temperature. Its temperature dependence is the script's: the ten constants
fitted at 50 C, each multiplied by an Arrhenius factor. It was fitted to
the 50 C experiment only. This section says how far the extrapolation
holds, with the model as shipped (default constants). All numbers are
measured with the package and pinned by its tests.

**At 50 C** the model describes the experiment it was fitted to: pooled R2
0.9939 and rms difference 0.0206 against the 16 deposited rows of Fig. 9b in
the script's fitting convention (0.9921 and 0.0234 in its heat-map
convention; the Table 3 O2 model has 0.9945).

**At 25, 75 and 100 C** it does not describe the experiments of Fig. 10a-c.
Against all deposited experimental rows of each temperature (t = 0
included, also the rows after the 1440 min at which the published panels
end), with `pooled_r2` and `rms_difference` of `temperature_series.py`:

| | 25 C | 75 C | 100 C |
|---|---|---|---|
| pooled R2 / rms, fitting convention | 0.653 / 0.220 | 0.641 / 0.175 | 0.757 / 0.147 |
| pooled R2 / rms, heat-map convention | 0.561 / 0.248 | 0.637 / 0.176 | 0.758 / 0.146 |
| time at which `X_HP3` falls to 0.5, model | 475 min | 15.2 min | 3.9 min |
| the same, experiment | not reached (0.876 left at 1890 min) | between 60 and 90 min | between 20 and 40 min |
| maximum MA fraction, model | 0.427 at 1203 min | 0.586 at 51 min | 0.594 at 13 min |
| maximum MA fraction, experiment | 0.023 at 1890 min | 0.436 at 180 min | 0.568 at 40 min |

(Model rows three and five: fitting convention, interpolated on a 0.1 min
grid. The text gives a maximum MA yield of 52.6 % at 75 C, which is not in
the deposited data; item 10 of "Paper versus deposit: discrepancies that
affect these models".)

The model is too fast at all three temperatures:

- **25 C.** It has 0.16 of the 3-HP left at 1890 min and an MA fraction of
  0.40; measured are 0.876 and 0.023. Its forward first-order constant of
  3-HP loss (`k_1*P/f0 + k_4`), 1.6e-3 1/min, is about 23 times the
  measured average (7e-5 1/min, from 0.876 left at 1890 min).
- **75 C.** The predicted `X_HP3` falls to 0.5 after 15 min; measured, it
  is 0.53 at 60 min and 0.28 at 120 min. The predicted MA maximum is higher
  (0.59 against 0.44) and 3.5 times earlier.
- **100 C.** The height of the MA maximum is about right (0.59 against
  0.57), but it comes at 13 min instead of 40 min, and the final acetic
  acid fraction is 0.17 against a measured 0.09 to 0.11.

For 75 C this follows from the data: the measured 3-HP conversion at 75 C
is no faster than at 50 C (`X_HP3` falls to 0.5 between 60 and 90 min at
both temperatures), whereas the Arrhenius factor of step 1 from 50 to 75 C
is 4.9. The paper's own Arrhenius plot shows the same for all three
temperatures: in Fig. 10d the 50 C point of 3-HP -> MA (`k_1`) lies a
factor of 4.3 above the line through the 25, 75 and 100 C points (previous
section). An extrapolation anchored at 50 C with the activation energies
of the text therefore overshoots at 25, 75 and 100 C, and the model's `k_1`
is 4.3, 4.3 and 4.2 times the plotted value there (`k_2`: 0.0045, 1.35 and
3.14 times; `k_5`: 0.10, 0.95 and 1.47 times).

**The 8 h screening experiments** at 25, 75 and 100 C are not described
either (25 C: 53 % conversion against a measured 8.5 %; 75 C: 6 % MA left
against a yield of 35.7 %), and the model's response to the pressure `P` is
contradicted by the paper's pressure series. Both comparisons are in
"Experiments of the paper that the models were not built on".

**What this means for Fig. 11.** The heat maps are the model fitted at 50 C
and extrapolated to 30 to 100 C. Their 50 C row rests on data. Their other
rows are an extrapolation that the paper's own experiments at 25, 75 and
100 C do not support. This holds for the published maps as for the maps of
this package: the published ones are the same model, with constants that
give nearly the same maps. `T` is in the model because the script contains
it, not because the response to it was established.

## Not modelled

Liquid volume and sampling withdrawals, gas-phase composition, pH, catalyst
deactivation, the loss of activity at higher O2 pressure and at higher H2O2
charge that the paper's screening experiments show (Figs. 2a and 3a), any
dependence on O2 pressure other than through the fixed factor `x_O2_ref` (O2
model) or the factor `P/f0` of the script (simplified model), and any
temperature other than 50 C in the two Table 3 models. The Table 3
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
  `'fig9b_nodes'`), and the refit of the O2 model with the O2 fraction held
  constant. The tests show that the shipped constants return the deposited
  nodes and how far each constant can be moved on its own; they do not show
  that no other constants or rate laws would return them.
- **Fits to the experiments:** the refit of the O2 constants with and
  without the 120 min row, and the independent scipy fit of `k_3` quoted in
  "Constants at 50 C". (The fits of the simplified model to the 50 C data
  are made with the shipped `fit_simplified_reversible_constants`. Those of
  `k_3`, of `k_3` and `k_4`, and of the five rate constants are tested;
  those of all ten constants, at 50 C and at the other temperatures, were
  run with the same function and are not pinned by a test.)
- **Fits to the digitized contour lines of Fig. 11:** the fits of `k_3`,
  of `k_3` and `k_4`, of the five rate constants, of `k_3`, `k_4`, `k_5`
  and `K_10`, and of all ten constants to the contour lines of one or both
  panels, quoted
  in "Constants at 50 C" and for the assignment of the Fig. 10d points, and
  the offsets "after refitting" of the three properties of the script. No
  shipped constant comes from these fits.
- **The emulation of the authors' fitting procedure** (item 5 of "Constants
  at 50 C"): 20 runs of a particle swarm of three iterations followed by
  five least-squares iterations. It is approximate (MATLAB's `particleswarm`
  and `lsqcurvefit` were not run and differ in detail) and is quoted for
  the size of the scatter only.
- **Other solver runs.** The fit of `k_3` with the trajectories
  integrated by the `ode23s` port (item 6 there); the repeat of the
  heat-map script's solver runs for the rows of Fig. 11 other than 100 C
  (the tests repeat the 100 C row); and the independent scipy
  implementations quoted for the simplified model and for the Table S4
  orders (every row of that table but the first). `solver_emulation.py`
  covers the deposited curves of Fig. 9 only. Its Jacobian is a plain
  forward difference; MATLAB's adapts its increments. Other increments,
  tried outside the package, moved the values by less than 1e-6. The
  margin of the node-time limits over rounding that the solver-emulation
  test module quotes (random relative errors of 1e-15 and 1e-14 added to
  the rates) was measured outside the tests too. MATLAB itself was never
  run.
- **Digitization.** The digitized curves that are shipped
  (`reference_data/`) were read from the published figures as described in
  `reference_data/README.md`; the colour field of Fig. 11 was read too and
  is not shipped.

What the package and its tests do establish: that the two Table 3 models,
with the constants as documented, agree with the deposited curves and the
experiments to the numbers given in the "Agreement" tables; that, with the
authors' solver run repeated, they return the deposited model nodes of
Fig. 9a, Fig. 9b and the H2O2 trace of SI Fig. S12a to the digits the file
prints; what each inferred choice (`x_O2_ref`, the O2 state, `k_7` without
the catalyst mass, the normalization without OPA) contributes to that with
the shipped constants; that the simplified reversible model implements the
equations of the authors' script; that its constants are the values of the
deposit and, for `k_3`, the result of the shipped fit; how far
its heat maps are from the published contour lines; that it does not
describe the experiments at 25, 75 and 100 C; and that each model
implements the rate laws written in these notes.
