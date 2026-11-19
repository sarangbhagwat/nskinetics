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
| `pd_c_3hp_oxidation.py` | loaders (`load_h2o2_model`, `load_o2_model`), `simulate_at_times`, `TABLE_S4` |
| `reference.py`, `reference_data/` | deposited Fig. 9 tables, digitized Fig. S12 inset curves, and their readers |
| `figures.py` | regenerates Fig. 9 and SI Fig. S12 |

The SBML files are generated from the Antimony files (the Antimony is the
source; `//` comments exist only there).

## Model basis

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
the Fig. 11 heat maps) were computed with that simplified, reversible model,
not with this six-reaction model. Predictions of this file away from 50 C are
therefore an extrapolation that has not been compared with the paper's
temperature data.

## Not modelled

Liquid volume and sampling withdrawals, gas-phase composition, pH, catalyst
deactivation, and any dependence on O2 pressure other than through the fixed
factor `x_O2_ref`. The Table 3 constants are per gram of catalyst for the
paper's batch and are used with `m_cat` = 0.2 g; the compartment volume is a
nominal 1 L and cancels from all concentrations.
