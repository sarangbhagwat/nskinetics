# Reference data for the Pd/C 3-HP oxidation models

Data used to check and plot the models in this folder against Kim, Chavarrio,
Install and Huber, Appl. Catal. B: Environ. Energy 386 (2026) 126403
(https://doi.org/10.1016/j.apcatb.2026.126403). Read the files through
`reference.py` in the parent folder, which returns tidy tables and corrects
the quirks listed below.

## Deposited by the authors (verbatim copies)

Source: S. Y. Kim, J. Chavarrio, J. Install, G. W. Huber, "Data and code from:
Selective oxidation of 3-hydroxypropionic acid to malonic acid over Pd/C:
Mechanistic and kinetic study", Dryad, https://doi.org/10.5061/dryad.jsxksn0q4.
Licence: CC0 1.0 (public domain dedication). The files are byte-for-byte
copies (UTF-8 with byte-order mark, CRLF line ends); `.gitattributes` at the
repository root marks the `*.csv` files of this folder `-text`, so that git
never converts their line ends.

| File | Content |
|---|---|
| `10-1_kinetic_model_H2O2.csv` | Fig. 9a (H2O2, 50 C): experimental points, the authors' model curve, and their predicted H2O2 trace of SI Fig. S12a |
| `10-2_kinetic_model_O2.csv` | Fig. 9b (O2, 50 C): experimental points and the authors' model curve |
| `11-4_Arrhenius.csv` | Fig. 10d (Arrhenius plot, O2): logarithms of three rate constants of the authors' simplified reversible model at 25, 50, 75 and 100 C, and the straight lines drawn through them |

In `10-1` and `10-2` all values are normalized concentrations (each species
over the sum of 3-HP, MA, AAc, OXA and D); the H2O2 trace is in equivalents
of the initial 3-HP. Time is in min.

Layout and known quirks:

- `10-2`: columns 1-6 are the experiment (16 rows, 0 to 1680 min) and columns
  7-12 the model curve (20 solver nodes, 0 to 1848 min). Labelled correctly.
- `10-1`: two caption lines come first. Columns 1-6, labelled `exp`, hold the
  **model** curve (24 solver nodes, 0 to 1584 min). Columns 9-14, labelled
  `model`, hold the **experimental** points (13 rows, 0 to 1440 min); their
  t = 0 row is printed as zeros and stands for the feed (1, 0, 0, 0, 0).
  Columns 7-8 hold the H2O2 trace on a time grid of its own (24 nodes, 0 to
  1537 min); its t = 0 cell reads 1 and stands for 15. The trace was generated
  with a different parameter set than the model curve and Table 3 of the
  paper (see `../NOTES.md`).
- `11-4`: two header rows. Columns 1-4 hold the points: 1/T in 1/K, printed
  to five decimals (0.00335, 0.00309, 0.00287, 0.00268, i.e. 25, 50, 75 and
  100 C), and the natural logarithms of three rate constants in 1/min,
  labelled `ln k1`, `ln k2` and `ln k3`. Columns 5-6, 7-8 and 9-10 hold the
  lines (10, 18 and 20 rows). The three constants are `k_1` (3-HP to MA),
  `k_2` (MA to AAc) and `k_5` (MA to D) of the simplified reversible model:
  the column labelled `k3` is that model's `k_5`. The 50 C row supplies
  `k_1`, `k_2` and `k_5` of the `'fig11b'` constant set of that model.

## Digitized from SI Fig. S12 (not deposited data)

These two curves are not in the deposit. They were read off SI Fig. S12 of
the article (the embedded 1430 x 549 px raster) by two independent
digitizations that were then reconciled: axes calibrated on the tick marks by
least squares, line centres taken per pixel column (per row where the line is
steep).

| File | Content | Reading uncertainty |
|---|---|---|
| `figS12a_h2o2_inset_digitized.csv` | Predicted H2O2 of SI Fig. S12a, in equivalents of the initial 3-HP. `source_axes` tells whether a point was read on the 0-15 inset or on the main axes, where the same curve is drawn below 1. Not readable between 0.99 and about 1.2 (342 to 368 min) | inset: 0.02 equivalents, or 1.5 min on the steep part; main axes: 0.3 % of the value above 0.1 and 0.0005 below 0.05 |
| `figS12b_opa_inset_digitized.csv` | Predicted OPA of the SI Fig. S12b inset as liquid mole fraction (the inset axis is in units of 1e-5): the 13 vertices of the plotted polyline, with the 3-HP model line of the main axes at the same times | `OPA`: 0.006e-5 on the decay and 0.015e-5 from 6 to 26 min; times as in column `t_unc` |

Columns of `figS12a_h2o2_inset_digitized.csv`: `t_min`, `h2o2_norm`,
`source_axes`, `kind` (how the point was located), `flag` (reading caveats).
Columns of `figS12b_opa_inset_digitized.csv`: `t_min`, `t_unc`, `OPA`, `HP3`,
`ratio_OPA_over_3HP`.

## Digitized from Fig. 11 (not deposited data)

The deposit contains the script that draws the heat maps of Fig. 11 but not
the heat-map matrices. The lines of the two main panels were therefore read
off the published figure (the embedded 1604 x 596 px raster) by two
independent digitizations that were then reconciled; they agree within 0.02
to 0.5 px on the contour lines. Axes were calibrated on the tick marks by
least squares (residuals below 0.07 px):

| Panel | Time axis | Temperature axis |
|---|---|---|
| a | t = 3.0723 col - 289.16 min | T = -0.14760 row + 105.244 C |
| b | t = 3.0763 col - 2809.3 min | T = -0.14801 row + 105.227 C |

(0-based pixel centres of the full image.) One pixel is 3.07 min or 0.148 C;
the comparison of a computed map with these lines uses 3.0723 min/px and
0.1476 C/px for both panels.

| File | Content |
|---|---|
| `fig11a_selectivity_digitized.csv` | Fig. 11a (MA "selectivity"): contour lines 0.1, 0.3 and 0.5, red dashed line, red and purple star |
| `fig11b_yield_digitized.csv` | Fig. 11b (MA "yield"): the same |

Columns: `series` (`contour`, `ridge` for the red dashed line, `red_star`,
`purple_star`), `level` (contours only), `branch` (contours only: `early` is
the part of a contour left of the red dashed line, `late` the part right of
it), `t_min` (time as plotted, which is the reaction time plus 10 min, see
`../NOTES.md`) and `T_C`.

What is and is not in these files:

- Every second sample of the digitization is kept (one sample per about
  2 px along a line), to keep the files small.
- Samples that the digitization flagged as uncertain are left out. They lie
  where lines run within 3 px of each other or within 4 to 5 px of the left
  frame (short times, mostly above 75 C in panel a and 70 C in panel b), and
  next to the red star. Their reading uncertainty is about 1 px (3 min).
- There are gaps where the published contour labels and the stars interrupt
  the lines.
- Reading uncertainty of the samples kept: about 0.3 px normal to a line
  (1 min where a line is steep, 0.05 C where it is flat); the steep early
  branches are quantized to pixel columns (0.5 px, 1.5 min). Stars: 0.3 px.
- The insets of the published panels are not included. They are not exact
  magnifications of the main panels: at equal temperature the inset lines of
  panel a lie 7 to 10 min later than those of the main panel, and in panel b
  the inset and main lines differ by up to 15 min (41 min for the late 0.3
  contour).
- The colour field of the maps was also digitized but is not shipped (several
  MB); numbers read from it are quoted in `../NOTES.md`.
