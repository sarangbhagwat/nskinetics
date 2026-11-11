# Manuscript figures: code and data

One folder per figure in Bhagwat et al., 2026. Each folder holds everything needed to redraw that figure:

- `make_figure.py`: regenerates the figure into `output\`, using only the folder's own `code\` and `data\`.
- `make_figure_frames.py` (four folders) and `make_figure_combined_frames.py` (`Fig6`): the figure as progressive frames for a presentation, see below.
- `make_figure_by_metric.py` (`EDFig2`): the sensitivity bars one metric per figure, see below.
- `manuscript\`: the image as placed in the manuscript, plus the PowerPoint source for composites.
- `comparison\`: diffs and notes.

## Regenerating

Run a folder's script from any working directory:

```
python make_figure.py
```

These scripts don't run any simulation -- they plot saved results from the analyses used to generate the data for these figures, scripts which are in `nskinetics\biorefineries\isobutanol\analyses`.

## Presentation frames

Four folders (`Fig6`, `EDFig3`, `EDFig4`, `EDFig8`) can also draw their figure as a sequence of progressive frames, for slides that reveal the figure a step at a time:

```
python make_figure_frames.py
```

This writes `output\frames\<name of make_figure.py's output>_frameNN_<frame name>.png`. The frames are cumulative and the last one is the complete figure. The figure is drawn once, by the same code as `make_figure.py`, and each frame only switches parts of it visible, so every frame has the same size and nothing moves between frames: placed at the same position on consecutive slides, they overlay exactly.

| Folder | Frames | Sequence |
| --- | --- | --- |
| `Fig6` | 13 | baseline; uninformed profitability campaign; process-level (TRY) campaigns; TRY trials as scouts; TRY-informed campaign, each followed by its panel-A annotation; then the proteome rows |
| `Fig6`, combined | 11 | `python make_figure_combined_frames.py`, frames named `kinBO_parameter_sets_combined_frameNN_<frame name>.png`: the optimization trajectories of `Fig6` (panel A, with its annotations and the campaign legend) and of `EDFig3` (panel A) side by side, without panel B: baseline; uninformed profitability campaign; TRY campaigns in their own metrics, the other campaigns' incumbents, then the TRY campaigns' profitability; TRY trials as scouts; TRY-informed campaign, with each annotation after what it describes |
| `EDFig3` | 10 | baseline; isobutanol campaigns; ethanol campaigns; other campaigns' incumbents; uninformed and TRY-informed profitability campaigns; then the proteome rows |
| `EDFig4` | 5 | x axis; IRR; profitability index; ethanol yield; growth derating |
| `EDFig8` | 7 | axes; IRR map; break-even line; profitability optimum; ethanol optima; isobutanol optima; cell density and TCI |

To change what each frame adds, edit the `FRAMES` table at the top of the folder's script in `code\frames\`; the groups it can name are listed there. `--dpi N` and `--formats png pdf svg` are passed on to that script.

The frames follow the manuscript versions of the figures where those differ from `make_figure.py`'s render: panel A of `Fig6` and `EDFig3` is titled "Optimization trajectories", `Fig6` carries the four panel-A annotations on a canvas extended at the top, and `EDFig4` has the label beside the growth-derating line in place of the legend.

## Sensitivity bars by metric

`EDFig2` can also draw its Sobol' / Shapley bars one metric per figure, with the decision variables grouped on the x axis instead of sorted by their effect on the profitability index: the six enzyme expression levels (glycolysis, ethanol production, isobutanol production), the three inhibition levels and the three substrate feeding variables.

```
python make_figure_by_metric.py
python make_figure_by_metric.py --metric "EtOH yield"
python make_figure_by_metric.py --list
```

Without `--metric` this writes the five metrics of the manuscript figure to `output\by_metric\EDFig2_bars_<metric>.png` (the metric's name with hyphens for spaces and symbols, e.g. `EtOH-yield`); `--metric` takes one or more of the metrics that `--list` prints, by name or title. The bars are the ones `make_figure.py` draws, from the same saved index table. The x order and the canvas are the same for every metric, so the figures overlay exactly on consecutive slides. To change the grouping, edit the `GROUPS` table at the top of `code\by_metric\plot_sobol_bars_by_metric.py`. `--dpi N` and `--formats png pdf svg` are passed on to that script.
