# Pd/C 3-HP oxidation models: provenance notes

Models of the selective oxidation of 3-hydroxypropionic acid (3-HP) to malonic
acid (MA) over 10 % Pd/C, after

> S. Y. Kim, J. Chavarrio, J. Install, G. W. Huber, "Selective oxidation of
> 3-hydroxypropionic acid to malonic acid over Pd/C: Mechanistic and kinetic
> study", Appl. Catal. B: Environ. Energy 386 (2026) 126403,
> https://doi.org/10.1016/j.apcatb.2026.126403

and the authors' data deposit (Dryad, https://doi.org/10.5061/dryad.jsxksn0q4,
CC0), referred to below as "the deposit".

These notes separate what is **printed** in the paper or its SI, what was
**inferred** from the deposit because the paper does not state it, and what
is an **NSKinetics addition**. They also list the places where the paper and
the deposit disagree and that matter for these models.

## Files

| File | Content |
|---|---|
| `pd_c_3hp_oxidation_h2o2_antimony.txt` (+ `_sbml.xml`) | H2O2 model: Table 3 H2O2 column, Fig. 9a, SI Fig. S12a |
| `pd_c_3hp_oxidation_o2_antimony.txt` (+ `_sbml.xml`) | O2 model: Table 3 O2 column, Fig. 9b, SI Fig. S12b, plus Arrhenius temperature dependence |
| `pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt` (+ `_sbml.xml`) | the authors' simplified reversible O2 model (no OPA), which generated Fig. 11; not in the paper |
| `pd_c_3hp_oxidation.py` | loaders (`load_h2o2_model`, `load_o2_model`, `load_o2_simplified_reversible_model`), `simulate_at_times`, `TABLE_S4`, the constant sets of the simplified model |
| `heat_maps.py` | the temperature-time maps of Fig. 11 and their comparison with digitized contour lines |
| `reference.py`, `reference_data/` | deposited Fig. 9 and Arrhenius tables, digitized Fig. S12 inset curves and Fig. 11 lines, and their readers |
| `figures.py` | regenerates Fig. 9, SI Fig. S12 and Fig. 11 |

The SBML files are generated from the Antimony files (the Antimony is the
source; `//` comments exist only there).

## Model basis

This section and the next four are about the two Table 3 models. The
simplified reversible model has its own section further down.

Conditions of both models: 60 mM 3-HP, NaOH/3-HP = 1:1, 0.2 g Pd/C, 50 C;
3 bar O2, or 15 equivalents of H2O2 under 3 bar Ar (Fig. 9 caption).

The authors' scripts work in liquid mole fractions and plot fractions of the
initial 3-HP. Here the liquid species are concentrations in M:
`c_i = C0 * n_i`, with `n_i` the moles of species i per initial mole of 3-HP
and `C0 = 0.06` M. The mapping is exact because every rate law is first order
in the liquid species. Time is in min.

States: `HP3`, `OPA` (3-oxopropanoic acid, the aldehyde intermediate; never
measured), `MA`, `AAc` (acetic acid), `OXA` (oxalic acid), `D` (lumped other
by-products), plus `H2O2` (H2O2 model) or the dimensionless O2 activity `yO2`
(O2 model). Co-products (H2O, H2, O2, CO2) are not tracked.

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

## Where each number comes from

| Quantity | Value | Status |
|---|---|---|
| `k_1`..`k_7` (H2O2), `k_1_ref`..`k_6_ref` (O2) | Table 3 as printed | printed |
| `m_cat` | 0.2 g | printed (Fig. 9 caption) |
| `C0`, `eq_H2O2_0` | 0.06 M, 15 | printed (Fig. 9 caption) |
| `k_7` not multiplied by `m_cat` | - | printed (Table 3 footnote b) |
| `x_per_M` | 0.000339/0.06 per M | inferred: the authors' scripts map the 60 mM feed to a mole fraction of 0.000339. This is not the mole fraction of 60 mM in water (about 1.08e-3); it is taken as used |
| `x_O2_ref` | 15*0.000339 = 0.005085 | **inferred**; not in the paper, the SI, the README or the posted script. See below |
| `yO2` as a depleting state, coefficient 1 for R1, R2, R4 | - | inferred from the deposited curve (the posted script sets yO2 = 1.0) |
| first-order, irreversible rate laws | - | inferred from the deposited curves; the paper says orders were fitted (Table S4) |
| H2O2 model structure | - | reconstructed: the authors posted only the O2 version of their script |
| `Ea_1` = `Ea_2`, `Ea_3`, `Ea_6` | 59781.2, 91594.6, 59549.8 J/mol | printed to three digits in Section 3.4 (59.8, 91.6, 59.5 kJ/mol); full values from the authors' heat-map script |
| `Ea_4` | 10000 J/mol | heat-map script only (a placeholder there); not in the paper |
| `Ea_5` | 51788.73 J/mol | heat-map script only; not in the paper |
| `T_ref`, `R_gas` | 323.15 K, 8.314 J/(mol K) | heat-map script; 50 C reference is printed |
| Arrhenius form on all six steps of the O2 model | - | **NSKinetics addition** (paper Eq. 6 form). See below |

### `x_O2_ref`

The Table 3 O2 entries for R1, R2 and R4 reproduce the deposited Fig. 9b model
curve only when multiplied by `m_cat` and by a further factor of 0.005085,
which equals the H2O2 mole fraction of the 15-equivalent H2O2 run
(15*0.000339). The effective constants recovered from the deposited curve,
divided by Table 3 times 0.2 g, give 0.0050854 (R1), 0.0050857 (R4) and
0.0050733 (R2); the rounding of the printed constants allows 0.00508 to
0.00509, and a rounded 0.0051 is excluded. How the factor entered the authors'
code is unknown. Only the products `k_i*m_cat*x_O2_ref` are identifiable.
The posted script `APCATB_D_25_08770_kinetic_opt.m`, run as posted (yO2 = 1.0
and its printed starting constants, no such factor), misses the deposited
curve by an rms of 0.42 to 0.76 for 3-HP, MA, OXA and D.

### Table S4 is not used

SI Table S4 lists reaction orders (H2O2: 0.99, 0.99, 1.49, 1.14, 1.02, 1.01
and 1 for R7; O2: 0.99, 0.95, 0.78, 0.87, 0.98, 1.01) and DFT equilibrium
constants. Neither enters the rate laws: the deposited curves are reproduced
only by first-order irreversible steps, and the posted script uses neither.
Applying the Table S4 orders moves the curves by up to 0.17 (O2) and 0.04
(H2O2). Table S4 also writes 0.5, 0.5 and 1.5 O2 for R1, R2 and R4, whereas
the O2 depletion in the deposited run is 1 for each; with `yO2` falling only
to 0.9995 this is immaterial (dropping the depletion changes the plotted
values by at most 7e-5). The values are kept as `TABLE_S4` in
`pd_c_3hp_oxidation.py` for documentation.

## Agreement with the deposit

Measured with these models (CVODE, relative tolerance 1e-9), evaluated at the
deposited node or sampling times; values are normalized concentrations.

| Comparison | Result |
|---|---|
| O2 model vs deposited Fig. 9b curve (19 nodes, t > 0) | rms 9.4e-4 (3-HP), 1.0e-3 (MA), 2.6e-5 (AAc), 2.1e-5 (OXA), 1.7e-4 (D); max 2.0e-3 |
| H2O2 model vs deposited Fig. 9a curve (24 nodes) | rms at most 3.3e-4 per species; max 9.6e-4 |
| O2 model `x_OPA` vs digitized SI Fig. S12b inset (13 vertices) | max 3.0e-7 on a peak of 1.4e-5 |
| O2 model vs experiment, pooled R2 (16 rows) | 0.9945 |
| H2O2 model vs experiment, pooled R2 | 0.9452 on 11 rows (t = 20 and 1205 min left out); 0.9412 on all 13 |
| H2O2 model `eq_H2O2` vs deposited H2O2 trace | up to 0.38 equivalents above it (rms 0.22) |

The residual differences from the deposited curves are the authors' solver
error and the rounding of Table 3, not a difference in the model: the
deposited curves are the raw steps of MATLAB `ode23s` at its default
tolerances (20 nodes for O2, 24 for H2O2). For O2 the nodes are up to 2.0e-3
away from the converged solution, and the straight segments drawn between
them in the published figure up to 0.010. When these models were
reconstructed, emulating the authors' solver with unrounded constants
reproduced the deposited values to their five printed decimals (max 5e-6 for
O2, 7e-6 for H2O2); the package tests do not repeat that emulation.

## Paper versus deposit: discrepancies that affect these models

1. **Deposit file `10-1_kinetic_model_H2O2.csv`.** The block labelled `exp`
   holds the model nodes and the block labelled `model` holds the experimental
   points. The experimental t = 0 row is printed as zeros (it is the feed,
   1, 0, 0, 0, 0). The t = 0 cell of the H2O2 trace reads 1 (it is 15). The
   readers in `reference.py` correct all three; the shipped CSV is verbatim.
2. **H2O2 trace of SI Fig. S12a.** The deposited trace (and the published
   inset) was generated with a different parameter set than the Fig. 9a curve
   and Table 3: recovered from the trace, its `k_7` is 6.94e-3 1/min instead
   of 6.48e-3 and its `k_5` is 22 % lower. The species curves in the same
   figure are the Table 3 ones. The Table 3 model therefore predicts up to
   0.38 equivalents more H2O2 than the published trace.
3. **Table 3, H2O2 column, versus the Fig. 9a curve.** The constants that
   generated the deposited curve differ from the printed ones by 0.1 to 0.6 %:
   +0.24 % (`k_1`), -0.60 % (`k_2`), +0.52 % (`k_3`), -0.24 % (`k_4`),
   -0.17 % (`k_5`) and +0.13 % (`k_7`). That is more than rounding to three
   digits allows for `k_1`, `k_2`, `k_3` and `k_7` (at most 0.09, 0.16, 0.07
   and 0.08 %), marginally more for `k_5` (0.173 against 0.165 %), and within
   it only for `k_4` (0.43 %). The effect on the curve is below 5e-4. `k_6`
   (4.10e-9) is not identifiable: any value below about 1e-6 1/min gives the
   same curves.
4. **"H2O2 is fully consumed after 1000 min"** (Section 3.4). In the Table 3
   model H2O2 falls below 1 % of its charge at 671 min, and 0.016 equivalents
   remain at 1000 min.
5. **OPA.** The text gives the predicted OPA peak with O2 as 8.4e-4 mM near
   25 min and calls it negligible. That number is the mole fraction on the
   authors' scale (1.4e-5) times 60 mM. In the model the OPA peak is 4.2 % of
   the initial 3-HP (2.5 mM) at 20.5 min; 25 min is a vertex of the plotted
   polyline. The H2O2 model predicts a larger unmeasured OPA pool: 11.5 % of
   the initial 3-HP at 84 min and 7.8 % at 1440 min.
6. **Parity R2.** The text states R2 > 0.95 for both oxidants and SI Fig. S11a
   prints 0.997 for O2. These models give 0.9945 for O2 (0.9957 without the
   120 min row) and 0.945 for H2O2. The source of 0.997 could not be
   identified, so nothing here is tested against it.
7. **The 120 min O2 sample** is plotted in Fig. 9b and listed in the deposit,
   but refitting shows it was left out of the fit that produced Table 3
   (leaving it out returns the published constants within 0.05 %).
8. **Initial-rate ratio.** The text reports r0(H2O2)/r0(O2) = 0.204 from
   Eq. 6. In these models both oxidant factors equal 0.005085 at t = 0, so
   the ratio of the modelled initial rates is 5.31/7.14 = 0.74.
9. **Solver.** Section 2.6 names `ode23s` and `ode15s`; the deposited curves
   come from `ode23s` at default tolerances.

## Temperature dependence of the O2 model (NSKinetics addition)

`k_i(T) = k_i_ref*exp(-(Ea_i/R_gas)*(1/T - 1/T_ref))` for all six steps, with
`T_ref` = 323.15 K, so that the default model is exactly the published 50 C
model. The activation energies are those of the authors' heat-map script,
which was written for their simplified five-reaction model (3-HP -> MA in one
step); the "MA production" value is applied here to both R1 and R2.

The paper's own temperature results (Fig. 10 curves at 25, 75 and 100 C and
the Fig. 11 heat maps) were computed with that simplified, reversible model
(`pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt`, next section),
not with this six-reaction model. Predictions of this file away from 50 C are
therefore an extrapolation. They are compared with the published heat maps in
the section after next; they have not been compared with the Fig. 10 data.

## Simplified reversible O2 model (the model behind Fig. 11)

### Purpose

The text presents Fig. 11 as the O2 model of Fig. 9b extended with the
Arrhenius equation. The deposited script that draws the figure
(`APCATB_D_25_08770_heatmap.m`) uses a different model: five reactions, no
OPA, a reverse term in every step. That model is not written out in the
paper or the SI. It is shipped as
`pd_c_3hp_oxidation_o2_simplified_reversible_antimony.txt`
(`load_o2_simplified_reversible_model`) because it, run as the script runs
it, reproduces the published heat maps, and the Table 3 model does not.
The 50 C point of the Arrhenius plot (Fig. 10d) supplies three of its
constants.

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
liquid species are concentrations in M as in the other two models,
`x_i = x_per_M*c_i`, the rates above are the assignment rules `r_1`..`r_5`,
each reaction runs at `r_j/x_per_M` (M/min), and `yH2' = r_3`,
`yCO2' = r_2 + r_3`. This is the script's system exactly. When the model
was built, its trajectories were compared with an independent scipy
implementation of the script's equations (same constants; DOP853, LSODA and
Radau at rtol 1e-11 or tighter) in both conventions below, for the three
constant sets at 30, 50, 84, 95 and 100 C. At the loader's CVODE tolerances
(rtol 1e-9, atol 1e-10) the liquid fractions agree to 2e-8 of the feed and
the gas fractions to 2e-9; with CVODE at rtol 1e-12, atol 1e-14 they agree
to 1e-10 (LSODA, Radau), so the remaining difference is solver tolerance.
For scale, rounding `k_2` to six digits (5.19998e-5) moves the trajectories
by 1e-7 of the feed. The package tests compare the rates and derivatives
instead.

Temperature: `p(T) = p_ref*exp(-(Ea/R_gas)*(1/T - 1/T_ref))` for all ten
parameters, `T_ref` = 323.15 K, with `Ea` = 59781.2, 91594.6, 10000,
51788.73 and 59549.8 J/mol for steps 1 to 5.

### Three properties of the script that are in the published figure

Each was tested against the contour lines digitized from Fig. 11b by removing
it (constants of the `'fig11b'` set; the rms offset is 0.27 px with all three
in place, 1 px = 3.07 min or 0.148 C). The offsets without refitting were
measured with this package; those "after refitting" are from the
reconstruction of this model. The package tests repeat neither.

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
   the heat maps miss the published contours by 14 px (9 px after
   refitting).
3. **The time axis of Fig. 11 is the reaction time plus 10 min.** The script
   places the feed at its first plotted time, 10 min. Model time here is the
   reaction time; `heat_maps.py` adds the offset (`FIG11_TIME_OFFSET`).
   Starting at plotted time 0 gives 3.8 px (3.2 px after refitting).

Two further points: the script uses a pressure of 3.5 bar (the captions give
3 bar of O2), and what Fig. 11 calls "MA selectivity" and "MA yield" are the
water-free liquid fraction of MA, `X_MA`, and `(1 - X_HP3)*X_MA` (model
quantities `fig11a_selectivity` and `fig11b_yield`).

### Constants at 50 C

The fitted constants are not in the deposit. The script only prints the
vector it starts its fit from, to four decimals (`'k0_printed'`:
0.0021, 0.0001, 0.0000, 0.0023, 0.0010, 4.7840, 0.5918, 3.5116, 81.6870,
52.0472), and that vector does not reproduce the figure (4.6 and 7.1 px on
panels a and b). The two sets shipped were reconstructed as follows. Rate
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

How well the inferred values are determined (spread between independent fits
of the same contours): `'fig11b'`: `k_3` and `k_4` about 1 %. `'fig11a'`:
`k_3` about 10 %, `k_4` 2 %, `k_5` 1 %, `K_10` 15 %. `K_6` has no effect on
the maps, and a 10 % change of `K_8`, `K_7`, `K_9` or `K_10` moves them by
at most 8e-8, 4e-5, 9e-4 and 9e-4. The authors' own vectors cannot be
recovered beyond this: they come from a particle swarm followed by five
least-squares iterations.

That the column "ln k3" of `11-4_Arrhenius.csv` is this model's `k_5` follows
from the figure: fitted freely to the contours of Fig. 11b (during the
reconstruction), `k_5` came out as 1.0613e-3 (-0.1 %). Its value, 1.0625e-3,
would print as 0.0011, not as the 0.0010 of the script's starting vector, so
that vector and the table are not from the same run.

**The two panels need different sets.** Since yield <= selectivity
everywhere, the late side of a yield contour must lie earlier than that of
the selectivity contour of the same level (one set of constants predicts
0.02 to 0.05 % earlier for the levels 0.1 and 0.3 above 62 C). In the
published figure the late 0.1 and 0.3 contours of panel b lie later than
those of panel a above 62 C, by 2.2 % in the median and up to 4.4 % (a
median of 1.2 to 2.4 px, up to 4.6 px). No single run can produce both
panels. With the `'fig11b'` set panel a is off by 1.07 px, with the
`'fig11a'` set panel b by 1.45 px. Why the panels differ (separate fits, or
other settings) is not known; the `'fig11a'` set is one of several that fit
panel a equally well.

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

Read from the colours of the published maps during the reconstruction (data
not shipped, not tested here): the maps agree with the colour field to an
rms of 0.004 (panel a) and 0.005 (panel b) where the value is above 0.3, and
the 50 C rows peak at 0.567 and 0.532 (+-0.01).

The script integrates each of the 299 grid intervals with `ode23s` at
default tolerances. An emulation of that differed from the converged
solution by at most 1.2e-4, so the solver does not matter at the level of
the figure.

### Paper versus deposit for Fig. 11

1. **Model.** The text and the caption describe the heat maps as the O2
   model of Fig. 9b with Arrhenius scaling. The script computes them with the
   simplified reversible model, with constants of its own. The next section
   gives the difference.
2. **Arrhenius scaling of the equilibrium constants**, with the activation
   energies of the rate constants, is not mentioned in the paper. Two of the
   five activation energies are not in the paper either (51.79 kJ/mol for
   3-HP -> D, and 10 kJ/mol for MA -> OXA, a placeholder).
3. **Time axis.** "27 min" and "49 min" in the text are the grid columns
   26.76 and 49.10 min of the plotted time, i.e. 16.8 and 39.1 min of
   reaction.
4. **Maximum selectivity, "59.1 % (95 C, 27 min)".** Reproduced: 0.5906 in
   that cell, which is the largest cell of the map.
5. **Maximum yield, "55.2 % (84 C, 49 min)".** The cell is reproduced as
   the largest of the map, but its value is 0.548. No constant set that
   matches the contours gives 0.552, and the colour of the published map
   along its crest reads 0.549 (+-0.006). The crest is flat: 13 cells
   between 76 and 90 C are within 0.001 of the maximum, so the cell that
   comes out largest is sensitive to small changes; the tests check the
   value in the cell, not that it is the largest.
6. **Stars.** The red star of panel a sits on the largest cell (27.5 min,
   94.7 C). The red star of panel b is drawn at 55.9 min, 7 min right of the
   cell named in the text. The purple stars are at 254 min (a) and 336 min
   (b) on the 50 C row, whose maxima are at 255.8 and 322.8 min.
7. **Insets.** The insets of the published panels are not exact
   magnifications of the main panels: at equal temperature their lines are
   7 to 10 min later (a) or up to 15 min off (b; 41 min for the late 0.3
   contour). `figures.py` draws the insets as zooms of the same matrix,
   with the contour labels on the same lines as published (those of the
   three steep early lines are staggered by up to 1 C so that they do not
   touch).

## Temperature-extended Table 3 O2 model versus the published heat maps

The paper's text describes Fig. 11 as an Arrhenius extension of the O2 model
of Fig. 9b; the deposited script computes it with the simplified reversible
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

The offsets above 1 px are linearized estimates (value difference over
gradient, capped at 30 px). Measured directly on the contour lines: the late
side of every contour of the Table 3 model is earlier than published, by 8
to 11 % of the time in the median per line (28 to 74 min; 9 to 260 min
over the whole lines), which is the same as 1.3 to 1.8 C lower at a given
time. The early sides differ less, by a median of 1 to 6 min in either
direction.

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

## Not modelled

Liquid volume and sampling withdrawals, gas-phase composition, pH, catalyst
deactivation, and any dependence on O2 pressure other than through the fixed
factor `x_O2_ref`. The Table 3 constants are per gram of catalyst for the
paper's batch and are used with `m_cat` = 0.2 g; the compartment volume is a
nominal 1 L and cancels from all concentrations.
