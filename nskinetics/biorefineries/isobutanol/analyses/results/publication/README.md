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

## Build GIFs for presentations

Fig6, EDFig3, EDFig4 and EDFig8 also have a `make_gifs.py`, the animated version of `make_figure.py`. It draws the same figure with the same `code\` and `data\` and writes one GIF per step towards the full figure into `output\gifs\`, to show a little of the figure at a time:

```
python make_gifs.py
```

- `<figure>_stepNN_<what it adds>.gif`: one step. It starts on the last frame of the step before it, plays once and rests on its last frame, so the steps can go one per slide or click.
- `<figure>_stepNN_<what it adds>.png`: the finished state of that step, as a still (`<figure>_step00_start.png` is the state before the first step).
- `<figure>_all_steps.gif`: the whole build in a loop, for previewing.
- `<figure>_steps.txt`: what each step shows.

The last step ends on the figure `make_figure.py` draws, with these differences:

- Fig6 carries the four panel A annotations of the manuscript version, with headroom above panel A for them and the panel title ("Optimization trajectories") moved above them.
- Fig6 and EDFig3 drop the blank band below panel B.
- Where Arial is not installed, the GIFs use a metric-compatible face (e.g. Liberation Sans) so the layout holds, while `make_figure.py` falls back to DejaVu Sans.

`python make_gifs.py --help` lists the options: frame width (smaller frames make smaller files), frame rate, speed, looping. The figure-specific storyboard is `code\animate_figure.py`. `code\gif_engine.py` is the shared engine, the same file in each of the four folders; `code\kin_opt_animation.py`, the same file in Fig6 and EDFig3, holds what those two share (their panel B). Needs Pillow in addition to the plotting packages.

These scripts don't run any simulation -- they plot saved results from the analyses used to generate the data for these figures, scripts which are in `nskinetics\biorefineries\isobutanol\analyses`.
