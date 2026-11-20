# Reference data for the Pd/C 3-HP oxidation models

Data used to check and plot the models in this folder against S. Y. Kim,
J. E. Chavarrio, J. Install and G. W. Huber, Appl. Catal. B: Environ. Energy
386 (2026) 126403 (https://doi.org/10.1016/j.apcatb.2026.126403). Read the
files through `reference.py` in the parent folder, which returns tidy tables
and corrects the quirks listed below.

## Deposited by the authors (verbatim copies)

Source: S. Y. Kim, J. Chavarrio, J. Install, G. W. Huber, "Data and code from:
Selective oxidation of 3-hydroxypropionic acid to malonic acid over Pd/C:
Mechanistic and kinetic study", Dryad, https://doi.org/10.5061/dryad.jsxksn0q4.
Licence: CC0 1.0 (public domain dedication), under Dryad's policy for all
deposited data; the README of the deposit states no licence of its own. The
paper's data-availability statement does not cite the deposit, which is
identified by its title, its authors and the manuscript number
APCATB-D-25-08770 in its file names. The six files are byte-for-byte copies
(UTF-8 with byte-order mark, CRLF line ends); `.gitattributes` at the
repository root marks them (`10-*.csv` and `11-*.csv` of this folder)
`-text`, so that git never converts their line ends, and the package tests
check their SHA-256.

| File | Content |
|---|---|
| `10-1_kinetic_model_H2O2.csv` | Fig. 9a (H2O2, 50 C): experimental points, the authors' model curve, and their predicted H2O2 trace of SI Fig. S12a |
| `10-2_kinetic_model_O2.csv` | Fig. 9b (O2, 50 C): experimental points and the authors' model curve |
| `11-1_25C_Kinetic.csv` | Fig. 10a (O2, 25 C): experimental points and the authors' model curve |
| `11-2_75C_Kinetic.csv` | Fig. 10b (O2, 75 C): the same |
| `11-3_100C_Kinetic.csv` | Fig. 10c (O2, 100 C): the same |
| `11-4_Arrhenius.csv` | Fig. 10d (Arrhenius plot, O2): logarithms of three rate constants at 25, 50, 75 and 100 C, and the straight lines drawn through them. The paper does not say of which model the constants are (see `../NOTES.md`) |

In `10-1`, `10-2` and `11-1` to `11-3` all values are normalized
concentrations (each species over the sum of 3-HP, MA, AAc, OXA and D); the
H2O2 trace is in equivalents of the initial 3-HP. Time is in min.

The model curves of `10-1` and `10-2` are solver steps, not converged
solutions: they are reproduced step for step by the method of the solver the
paper names (MATLAB `ode23s`) at default tolerances, in a Python port; MATLAB
itself was not run. `solver_emulation.py` in the parent folder repeats those
runs with the shipped models, and the package tests check that it returns
every deposited node: the same number of steps (one fewer than the node
counts below, which include t = 0), node times within 2e-5 min (printed to
five decimals), and values within 1e-5, one unit of the last printed decimal
(see `../NOTES.md`, "Emulation of the authors' solver run"). The constants
used for the species curves of `10-1` and `10-2` are sets inferred from the
nodes, 0.01 to 0.6 % off Table 3 as printed; with Table 3 as printed those
values agree within 4.4e-4. The H2O2 trace of `10-1` needs a separate
inferred set, 4 to 22 % off Table 3 (`k_6` unchanged); with Table 3 as
printed that run takes 23 steps instead of 24.

Layout and known quirks:

- `10-2`: columns 1-6 are the experiment (16 rows, 0 to 1680 min) and columns
  7-12 the model curve (20 solver nodes, 0 to 1848 min). Labelled correctly.
- `10-1`: two caption lines come first. Columns 1-6, labelled `exp`, hold the
  **model** curve (24 solver nodes, 0 to 1584 min). Columns 9-14, labelled
  `model`, hold the **experimental** points (13 rows, 0 to 1440 min); their
  t = 0 row is printed as zeros and stands for the feed (1, 0, 0, 0, 0).
  Columns 7-8 hold the H2O2 trace on a time grid of its own (24 nodes, 0 to
  1537 min); its t = 0 cell reads 1 and stands for 15. The trace is not the
  H2O2 of the model curve in the same file: it lies up to 0.38 equivalents
  below the H2O2 of the Table 3 model and is matched only with other
  constants (see `../NOTES.md`). It is reproduced as a solver run of its
  own, of 24 steps to 1584 min, of which the file holds all nodes but the
  last.
- `11-1`, `11-2`, `11-3`: laid out like `10-2` and labelled correctly; no
  quirks. Columns 1-6 are the experiment (14, 12 and 11 rows; last sampling
  times 1890, 1440 and 1560 min) and columns 7-12 the model curve (19, 20 and
  25 nodes; it ends at 1.1 times the last sampling time: 2079, 1584
  and 1716 min). The model times are printed to two decimals and the model
  values to five, or to six decimals (shown with three digits) below 1e-3.
  These model curves are produced by neither of the two models of this
  package; only the experimental points are used (see `../NOTES.md`).
- `11-4`: two header rows. Columns 1-4 hold the points: 1/T in 1/K, printed
  to five decimals (0.00335, 0.00309, 0.00287, 0.00268, i.e. 25, 50, 75 and
  100 C), and the natural logarithms of three rate constants, labelled
  `ln k1`, `ln k2` and `ln k3`. The second header row gives a unit for 1/T
  only (`/K`); the constants are taken to be in 1/min, which is inferred
  (`../NOTES.md`). Columns 5-6, 7-8 and 9-10 hold the lines (10, 18 and 20
  rows). The paper does not say of which model the three constants are; the
  published figure names them 3-HP to MA, MA to AAc and MA to Others. The
  reader returns them as `k_1`, `k_2` and `k_5`, the numbers these reactions
  have in the authors' heat-map script, so the column labelled `k3` is read
  as `k_5`. These are not the numbers of the O2 model of this package (3-HP
  to MA is its R1 followed by R2, MA to AAc its R3 and MA to D its R6), and
  the 50 C row is not numerically the Table 3 constants (`../NOTES.md`
  compares the two).

## Digitized from SI Fig. S12 (not deposited data)

These two curves are not in the deposit. They were read off SI Fig. S12 of
the article (the embedded 1430 x 549 px raster), twice and independently,
and the two readings were merged: axes calibrated on the tick marks by
least squares, line centres taken per pixel column (per row where the line is
steep). The scripts that read the figures are not shipped.

| File | Content | Reading uncertainty |
|---|---|---|
| `figS12a_h2o2_inset_digitized.csv` | Predicted H2O2 of SI Fig. S12a, in equivalents of the initial 3-HP. `source_axes` tells whether a point was read on the 0-15 inset or on the main axes, where the same curve is drawn below 1. Not readable between 0.99 and about 1.2 (342 to 368 min) | inset: 0.02 equivalents, or 1.5 min on the steep part; main axes: 0.3 % of the value above 0.1 and 0.0005 below 0.05 |
| `figS12b_opa_inset_digitized.csv` | Predicted OPA of the SI Fig. S12b inset as liquid mole fraction (the inset axis is in units of 1e-5): the 13 vertices of the plotted polyline, with the 3-HP model line of the main axes at the same times | `OPA`: 0.006e-5 on the decay and 0.015e-5 from 6 to 26 min; times as in column `t_unc` |

Columns of `figS12a_h2o2_inset_digitized.csv`: `t_min`, `h2o2_norm`,
`source_axes`, `kind` (how the point was located), `flag` (reading caveats).
Columns of `figS12b_opa_inset_digitized.csv`: `t_min`, `t_unc`, `OPA`, `HP3`,
`ratio_OPA_over_3HP`.
