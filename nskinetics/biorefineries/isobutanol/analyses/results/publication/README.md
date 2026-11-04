# Manuscript figures: code and data

One folder per figure of Bhagwat et al. (NSKinetics + isobutanol biorefinery
manuscript). Each folder holds everything needed to redraw that figure:

| Folder   | Manuscript figure                           | Drawn by (snapshot in `code\`)                                                                              | Inputs (`data\`)                                        | Regenerated vs manuscript                                               |
| -------- | ------------------------------------------- | ----------------------------------------------------------------------------------------------------------- | ------------------------------------------------------- | ----------------------------------------------------------------------- |
| `Fig1`   | Fig. 1 framework overview                   | claude.ai Design artifact "Figure 1 – NSKinetics framework" (`Main.dc.html`), rendered with headless Chrome | none                                                    | same text and geometry; HTML text sits 1–5 px lower (browser rendering) |
| `Fig2`   | Fig. 2 kinetic-model network                | `models\...\conceptual_diagram.py` (network-only variant)                                                   | none                                                    | byte-identical                                                          |
| `Fig3`   | Fig. 3 block flow diagram                   | `plots\plot_block_flow_diagram.py`                                                                          | none                                                    | byte-identical                                                          |
| `Fig4`   | Fig. 4 Monte Carlo bivariate KDEs           | `plots\plot_uncertainty_MPSP_vs_TCI.py`                                                                     | 6,000-sim uncertainty workbook (2026.9.29-23.45)        | byte-identical                                                          |
| `Fig5`   | Fig. 5 IRR, glycolysis × ethanol inhibition | `plots\plot_glycolysis_inhib_ethanol_A.py`                                                                  | 40×40 sweep CSVs                                        | byte-identical                                                          |
| `Fig6`   | Fig. 6 optimization campaigns + proteome    | `plots\plot_kin_opt_parameter_sets.py --variants irr` (placed in PowerPoint)                                | 8 split_12d campaign CSVs                               | identical panel; PowerPoint callouts/heading added                      |
| `EDFig1` | ED Fig. 1 feeding strategy + trajectory     | `plots\plot_feeding_strategy_baseline_A.py` (cache stage only)                                              | feeding-strategy sweep CSVs + baseline trajectory cache | byte-identical                                                          |
| `EDFig2` | ED Fig. 2 Sobol'/Shapley                    | `analyses\sensitivity\analyze_sobol_split12d.py --replot`                                                   | Sobol' index table                                      | byte-identical                                                          |
| `EDFig3` | ED Fig. 3 TRY trajectories + proteome       | `plots\plot_kin_opt_parameter_sets.py --variants process` (placed in PowerPoint)                            | 8 campaign CSVs                                         | identical panel; heading swapped in PowerPoint                          |
| `EDFig4` | ED Fig. 4 k3 sweep                          | `plots\plot_axis_try_informed_sweep.py --param k_3` (placed in PowerPoint)                                  | k3 sweep CSV + anchor JSON                              | curves identical; legend cropped, label added                           |
| `EDFig5` | ED Fig. 5 TEA breakdowns                    | `plots\plot_tea_breakdowns_split12d.py` (stage 2)                                                           | breakdown JSON (2026.09.24-17.48)                       | byte-identical                                                          |
| `EDFig6` | ED Fig. 6 Spearman tornado                  | `plots\plot_uncertainty_spearman_tornado.py`                                                                | same 6,000-sim workbook                                 | byte-identical                                                          |
| `EDFig7` | ED Fig. 7 final parameters                  | `plots\plot_kin_opt_parameter_sets.py --variants parameters`                                                | 8 campaign CSVs + A/B workbooks                         | identical (PowerPoint export at 300 dpi)                                |
| `EDFig8` | ED Fig. 8 IRR, k13 × isobutanol inhibition  | `plots\plot_k13_inhib_isobutanol_opt_PI_TRY_informed.py`                                                    | 40×40 sweep CSVs (2026-10-26)                           | byte-identical                                                          |
| `EDFig9` | ED Fig. 9 TRY-informed optimum trajectory   | `plots\plot_trajectory_opt_PI_TRY_informed.py` (cache stage only)                                           | trajectory cache CSV/JSON                               | byte-identical                                                          |

Each folder contains:

- `make_figure.py`: regenerates the figure into `output\`, using only the folder's own `code\` and `data\`.
- `manuscript\`: the image as placed in the manuscript, plus the PowerPoint source for composites.
- `comparison\`: diffs and notes.

## Regenerating

Run a folder's script from any working directory:

```
python make_figure.py
```

These scripts don't run any simulation -- they plot saved results from the analyses used to generate the data for these figures, scripts which are in `nskinetics\biorefineries\isobutanol\analyses`.
