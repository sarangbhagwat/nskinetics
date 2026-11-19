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
copies (UTF-8 with byte-order mark, CRLF line ends).

| File | Content |
|---|---|
| `10-1_kinetic_model_H2O2.csv` | Fig. 9a (H2O2, 50 C): experimental points, the authors' model curve, and their predicted H2O2 trace of SI Fig. S12a |
| `10-2_kinetic_model_O2.csv` | Fig. 9b (O2, 50 C): experimental points and the authors' model curve |

All values are normalized concentrations (each species over the sum of 3-HP,
MA, AAc, OXA and D); the H2O2 trace is in equivalents of the initial 3-HP.
Time is in min.

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

## Digitized from the published SI figure (not deposited data)

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
