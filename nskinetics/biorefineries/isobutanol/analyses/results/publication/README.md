# Manuscript figures: code and data

One folder per figure in Bhagwat et al., 2026. Each folder holds everything needed to redraw that figure:

- `make_figure.py`: regenerates the figure into `output\`, using only the folder's own `code\` and `data\`.
- `manuscript\`: the image as placed in the manuscript, plus the PowerPoint source for composites.
- `comparison\`: diffs and notes.

## Regenerating

Run a folder's script from any working directory:

```
python make_figure.py
```

These scripts don't run any simulation -- they plot saved results from the analyses used to generate the data for these figures, scripts which are in `nskinetics\biorefineries\isobutanol\analyses`.
