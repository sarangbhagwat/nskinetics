# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Step-by-step build GIFs of a finished matplotlib figure, for presentations.

The figure is drawn once, by the unmodified plot script. Its artists are then
sorted into named groups, every group is hidden, and a storyboard of *steps*
reveals the groups again, one step per GIF:

    board = Storyboard(fig)
    board.group('yield', yield_line, yield_axis)
    board.group('irr', irr_line, irr_axis)
    board.step('ethanol-yield', 'Ethanol yield rises with Pdc capacity',
               wipe('yield', region=Region(ax)))
    board.step('irr', 'IRR peaks and then falls',
               wipe('irr', region=Region(ax)))
    board.render(out_dir, 'my_figure')

Each step writes one GIF that starts on the last frame of the step before it
and ends on its own finished state, so the GIFs chain across slides with no
jump. The last step ends on the complete figure; render() checks that this
final state is pixel-identical to the figure with nothing hidden.

How a frame is made: the figure is rendered once per *state* (the artists
visible before and after a phase). Every in-between frame is a per-pixel blend
of the two states, ``A*(1 - m) + B*m``, where the mask ``m`` is uniform for a
fade and a soft-edged moving front for a wipe. Since B is A plus the newly
shown artists, a wipe reads as the new artists being drawn along the axis
(trials accumulating, a curve being traced), and no frame ever needs a
matplotlib render of its own.

Files written by render() into ``out_dir``:

    <stem>_step00_start.png            the state before the first step
    <stem>_stepNN_<slug>.gif           one build step (plays once by default)
    <stem>_stepNN_<slug>.png           the finished state of that step
    <stem>_all_steps.gif               every step in sequence, looping
    <stem>_steps.txt                   step number -> file -> caption

GIF notes:

* By default a step GIF carries no loop instruction, so players that follow
  the file play it once and rest on its last frame. Players that loop
  regardless restart only after the last frame's hold, ``end_hold`` seconds
  (default 60).
* Every GIF of a figure shares one 256-colour palette, built from all the
  states, so a step's first frame equals the previous step's last frame
  exactly. The background and pure black are kept exact.
* The final state is compared with the complete figure as this script renders
  it. It matches make_figure.py's image wherever both use the same font.

Fonts: the plot scripts ask for Arial. Where Arial is not installed (most
Linux machines), use_arial_or_metric_twin() registers a metric-compatible
face (Liberation Sans, Arimo, ...) under the name Arial, so the layout tuned
for Arial holds; call it before building the figure.

Sim-safe: imports only the standard library, numpy, matplotlib and Pillow.
"""
import contextlib
import dataclasses
import io
import math
import os
import re
from typing import Optional, Sequence

import numpy as np
import matplotlib
from matplotlib import font_manager
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.transforms import Bbox
from PIL import Image

__all__ = ['Storyboard', 'Region', 'fade', 'wipe', 'cut', 'pause',
           'capture_figure', 'use_arial_or_metric_twin', 'add_cli_arguments',
           'render_kwargs']

ARIAL_TWINS = ('Liberation Sans', 'Arimo', 'Nimbus Sans', 'Nimbus Sans L',
               'FreeSans', 'Helvetica')


# -----------------------------------------------------------------------------
# fonts and figure capture
# -----------------------------------------------------------------------------
def use_arial_or_metric_twin():
    """Make the family name 'Arial' resolve on this machine.

    Returns 'Arial' when it is installed, the name of the metric-compatible
    face registered in its place otherwise, or None when neither exists (the
    plot scripts then fall back to DejaVu Sans, whose wider glyphs can make
    labels collide).
    """
    manager = font_manager.fontManager
    names = {entry.name for entry in manager.ttflist}
    if 'Arial' in names:
        return 'Arial'
    for twin in ARIAL_TWINS:
        entries = [e for e in manager.ttflist if e.name == twin]
        if entries:
            manager.ttflist.extend(
                dataclasses.replace(e, name='Arial') for e in entries)
            # forget lookups made before the alias existed
            cache = getattr(manager, '_findfont_cached', None)
            if hasattr(cache, 'cache_clear'):
                cache.cache_clear()
            return twin
    return None


def capture_figure(draw, *args, **kwargs):
    """Call a plot function that saves and closes its figure, and return the
    figure instead, still open and with nothing written to disk.

    ``draw(*args, **kwargs)`` runs with Figure.savefig and pyplot.close
    switched off; the figure it created last is returned. What the function
    prints is dropped (it reports files as written that are not). The plot
    script itself stays unmodified.
    """
    before = set(plt.get_fignums())
    savefig, close = Figure.savefig, plt.close
    Figure.savefig = lambda self, *a, **k: None
    plt.close = lambda *a, **k: None
    try:
        with open(os.devnull, 'w') as sink, contextlib.redirect_stdout(sink):
            draw(*args, **kwargs)
    finally:
        Figure.savefig, plt.close = savefig, close
    new = sorted(set(plt.get_fignums()) - before)
    if not new:
        raise RuntimeError('the plot function created no figure')
    for num in new[:-1]:
        plt.close(num)
    return plt.figure(new[-1])


# -----------------------------------------------------------------------------
# regions and phases
# -----------------------------------------------------------------------------
@dataclasses.dataclass
class Region:
    """A rectangle of the figure that a wipe sweeps across.

    axes     one Axes or several (their bounding boxes are joined, e.g. the
             two halves of a broken axis).
    x, y     optional (low, high) data limits on the first axes, to sweep
             only part of it (trials 1000-2000, one bar row).
    artists  alternatively (or additionally), artists whose drawn extents are
             joined into the rectangle.
    pad      margin in points added on every side, so markers and line caps
             that poke out of the axes are swept with the rest.
    """
    axes: object = None
    x: Optional[Sequence[float]] = None
    y: Optional[Sequence[float]] = None
    artists: Sequence = ()
    pad: float = 3.0

    def bbox(self, fig, renderer):
        """The rectangle in display pixels (origin bottom-left)."""
        boxes = []
        axes = self.axes
        if axes is not None:
            axes = list(axes) if isinstance(axes, (list, tuple)) else [axes]
            box = Bbox.union([a.get_window_extent(renderer) for a in axes])
            x0, y0, x1, y1 = box.extents
            first = axes[0]
            if self.x is not None:
                lo, hi = sorted(first.transData.transform(
                    [(self.x[0], 0), (self.x[1], 0)])[:, 0])
                x0, x1 = max(x0, lo), min(x1, hi)
            if self.y is not None:
                # the y limits keep their data order (an inverted axis is fine)
                lo, hi = sorted(first.transData.transform(
                    [(first.get_xlim()[0], self.y[0]),
                     (first.get_xlim()[0], self.y[1])])[:, 1])
                y0, y1 = max(y0, lo), min(y1, hi)
            boxes.append(Bbox.from_extents(x0, y0, x1, y1))
        for artist in self.artists:
            boxes.append(artist.get_window_extent(renderer))
        if not boxes:
            raise ValueError('a Region needs axes or artists')
        if axes is None and (self.x is not None or self.y is not None):
            raise ValueError('Region x / y limits need axes')
        pad = self.pad * fig.dpi / 72.0
        return Bbox.union(boxes).padded(pad)


@dataclasses.dataclass
class Phase:
    kind: str                       # 'fade' | 'wipe' | 'cut' | 'pause'
    groups: tuple = ()
    seconds: float = 0.0
    regions: tuple = ()
    direction: str = 'right'
    edge: float = 7.0               # soft front width, points
    rest: str = 'fade'              # outside the regions: 'fade' | 'end'


def fade(*groups, seconds=0.7):
    """Show the groups by fading them in together."""
    return Phase('fade', groups, seconds)


def wipe(*groups, region, seconds=2.0, direction='right', edge=7.0,
         rest='fade'):
    """Show the groups by sweeping a front across ``region`` (one Region or a
    list swept in step with each other).

    direction: 'right', 'left', 'up' or 'down', the way the front travels.
    Anything of these groups lying outside the regions (a legend entry, an
    axis label) fades in over the same time (rest='fade') or appears only
    when the sweep is done (rest='end').
    """
    regions = tuple(region) if isinstance(region, (list, tuple)) else (region,)
    if direction not in ('right', 'left', 'up', 'down'):
        raise ValueError(f'unknown wipe direction {direction!r}')
    if rest not in ('fade', 'end'):
        raise ValueError(f'unknown rest mode {rest!r}')
    return Phase('wipe', groups, seconds, regions, direction, edge, rest)


def cut(*groups):
    """Show the groups at once, with no transition."""
    return Phase('cut', groups, 0.0)


def pause(seconds):
    """Hold the current state."""
    return Phase('pause', (), seconds)


@dataclasses.dataclass
class Step:
    slug: str
    caption: str
    phases: tuple


def _ease(t):
    """Smoothstep: slow in, slow out."""
    return t * t * (3.0 - 2.0 * t)


def _flatten(items):
    for item in items:
        if item is None:
            continue
        if isinstance(item, (list, tuple, set)):
            yield from _flatten(item)
        elif hasattr(item, 'set_visible'):
            yield item
        elif hasattr(item, '__iter__'):   # matplotlib containers
            yield from _flatten(list(item))
        else:
            raise TypeError(f'not an artist: {item!r}')


# -----------------------------------------------------------------------------
# the storyboard
# -----------------------------------------------------------------------------
class Storyboard:
    """Groups of artists of one figure and the steps that reveal them.

    savefig_kwargs: the keyword arguments the plot script passes to
    fig.savefig besides the file name and dpi (``bbox_inches='tight'``,
    ``pad_inches``, ``facecolor`` ...), so the frames are framed like the
    static figure. A tight bounding box is measured once, on the complete
    figure, and reused for every state.
    """

    def __init__(self, fig, savefig_kwargs=None):
        self.fig = fig
        self.savefig_kwargs = dict(savefig_kwargs or {})
        self.groups = {}
        self.steps = []

    # -- declaration ----------------------------------------------------------
    def group(self, name, *artists):
        """Register artists under a group name (repeat calls add to it).
        Artists that are already invisible are left out, so they stay hidden;
        an artist can belong to one group only. Returns the group's artist
        list."""
        members = self.groups.setdefault(name, [])
        owner = {id(a): other for other, group in self.groups.items()
                 for a in group}
        for artist in _flatten(artists):
            if not artist.get_visible() or owner.get(id(artist)) == name:
                continue
            if id(artist) in owner:
                raise ValueError(f'{artist!r} is already in group '
                                 f'{owner[id(artist)]!r}, cannot join {name!r}')
            owner[id(artist)] = name
            members.append(artist)
        return members

    def require_grouped(self, *artists):
        """Raise if any of these (visible) artists is in no group. Call it
        with everything the plot script draws on an axes, so that an artist
        the storyboard does not know about (a new landmark, say) is reported
        instead of sitting on screen from the first frame."""
        grouped = {id(a) for group in self.groups.values() for a in group}
        loose = [a for a in _flatten(artists)
                 if a.get_visible() and id(a) not in grouped]
        if loose:
            raise ValueError(f'artists in no group: {loose}')

    def step(self, slug, caption, *phases):
        """Add a step: one GIF, made of one or more phases played in order."""
        if not phases:
            raise ValueError(f'step {slug!r} has no phases')
        self.steps.append(Step(slug, caption, tuple(phases)))

    # -- checks ---------------------------------------------------------------
    def _check(self):
        shown = {}
        for step in self.steps:
            for phase in step.phases:
                for name in phase.groups:
                    if name not in self.groups:
                        raise KeyError(f'step {step.slug!r} shows the unknown '
                                       f'group {name!r}')
                    if name in shown:
                        raise ValueError(
                            f'group {name!r} is shown twice (steps '
                            f'{shown[name]!r} and {step.slug!r})')
                    shown[name] = step.slug
        never = [name for name in self.groups if name not in shown]
        if never:
            raise ValueError(f'groups never shown: {never}')
        empty = [name for name, members in self.groups.items() if not members]
        if empty:
            raise ValueError(f'empty groups: {empty}')
        slugs = [step.slug for step in self.steps]
        if len(set(slugs)) != len(slugs):
            raise ValueError('step slugs must be unique')

    # -- rendering ------------------------------------------------------------
    def _set(self, names, visible):
        for name in names:
            for artist in self.groups[name]:
                artist.set_visible(visible)

    def _snapshot(self, dpi, crop):
        """The current state as an RGB uint8 array."""
        buffer = io.BytesIO()
        kwargs = dict(self.savefig_kwargs)
        kwargs.pop('bbox_inches', None)
        kwargs.pop('pad_inches', None)
        self.fig.savefig(buffer, format='png', dpi=dpi, bbox_inches=crop,
                         **kwargs)
        buffer.seek(0)
        image = Image.open(buffer).convert('RGBA')
        canvas = Image.new('RGBA', image.size, (255, 255, 255, 255))
        canvas.alpha_composite(image)
        return np.asarray(canvas.convert('RGB'))

    def _crop(self, renderer):
        """The saved area in inches: the whole canvas, or the tight box of
        the complete figure when the plot script saves with
        bbox_inches='tight'."""
        fig = self.fig
        requested = self.savefig_kwargs.get('bbox_inches')
        if requested == 'tight':
            pad = self.savefig_kwargs.get(
                'pad_inches', matplotlib.rcParams['savefig.pad_inches'])
            return fig.get_tightbbox(renderer).padded(pad)
        if isinstance(requested, Bbox):
            return requested
        return Bbox.from_bounds(0, 0, *fig.get_size_inches())

    def _mask(self, phase, t, shape, crop, renderer):
        """The blend mask of a phase at progress t in [0, 1]: float32 (H, W),
        0 = the state before, 1 = the state after."""
        height, width = shape
        if phase.kind == 'fade':
            return np.full(shape, _ease(t), np.float32)
        dpi = self.fig.dpi
        if phase.rest == 'fade':
            mask = np.full(shape, _ease(t), np.float32)
        else:
            mask = np.full(shape, 1.0 if t >= 1.0 else 0.0, np.float32)
        swept = np.full(shape, -1.0, np.float32)   # -1 = outside every region
        soft = max(1.0, phase.edge * dpi / 72.0)
        left_px = crop.x0 * dpi
        top_px = crop.y1 * dpi
        for region in phase.regions:
            box = region.bbox(self.fig, renderer)
            c0 = int(math.floor(box.x0 - left_px))
            c1 = int(math.ceil(box.x1 - left_px))
            r0 = int(math.floor(top_px - box.y1))
            r1 = int(math.ceil(top_px - box.y0))
            c0, c1 = max(c0, 0), min(c1, width)
            r0, r1 = max(r0, 0), min(r1, height)
            if c1 <= c0 or r1 <= r0:
                continue
            horizontal = phase.direction in ('right', 'left')
            n = (c1 - c0) if horizontal else (r1 - r0)
            # position along the sweep, 0 at the side the front starts from
            along = np.arange(n, dtype=np.float32) + 0.5
            if phase.direction in ('left', 'up'):
                along = along[::-1]
            # the front travels from -soft to n, so the region is empty at
            # t = 0 and complete at t = 1
            front = -soft + t * (n + soft)
            ramp = np.clip((front - along) / soft + 1.0, 0.0, 1.0)
            if t >= 1.0:
                ramp[:] = 1.0
            # rows run top to bottom; 'up' already reversed `along`
            ramp = ramp[None, :] if horizontal else ramp[:, None]
            # where regions overlap, the front that is further along wins
            swept[r0:r1, c0:c1] = np.maximum(swept[r0:r1, c0:c1], ramp)
        return np.where(swept >= 0.0, swept, mask)

    def render(self, out_dir, stem, width_px=1600, fps=20, speed=1.0,
               loop=False, lead_in=0.4, end_hold=None, combined=True,
               combined_pause=1.5, dither=False, stills=True, strict=True,
               verbose=True):
        """Write the step GIFs (and stills, the combined GIF and the step
        list) into out_dir. Returns the list of step GIF paths.

        width_px        frame width; the height follows the figure.
        fps, speed      frame rate (at most 50, the fastest a GIF can play),
                        and a factor on every phase duration (speed=2 halves
                        them).
        loop            False: play once and rest on the last frame; True:
                        loop forever.
        lead_in         seconds the first frame is held before motion starts.
        end_hold        seconds the last frame of a step GIF is held. Default:
                        60 for a play-once GIF (it matters only where the
                        player loops regardless), 2 with loop=True. A pause()
                        that ends a step is replaced by this hold.
        combined        also write <stem>_all_steps.gif (looping), holding
                        combined_pause seconds after each step.
        dither          Floyd-Steinberg dithering of the frames (smoother
                        gradients, larger files).
        strict          raise if the final state differs from the complete
                        figure.
        """
        if not 0 < fps <= 50:
            raise ValueError(f'fps must be in (0, 50], got {fps}')
        if speed <= 0:
            raise ValueError(f'speed must be positive, got {speed}')
        if width_px < 16:
            raise ValueError(f'width_px must be at least 16, got {width_px}')
        if end_hold is None:
            end_hold = 2.0 if loop else 60.0
        self._check()
        fig = self.fig
        # the quarter pixel keeps float rounding from costing a pixel column
        dpi = (width_px + 0.25) / fig.get_size_inches()[0]
        fig.set_dpi(dpi)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        crop = self._crop(renderer)
        requested = self.savefig_kwargs.get('bbox_inches')
        if requested is not None:
            # width_px is the width of the saved (cropped) frame; a tight box
            # moves slightly with the dpi, so settle it in two passes
            for _ in range(2 if requested == 'tight' else 1):
                dpi = (width_px + 0.25) / crop.width
                fig.set_dpi(dpi)
                fig.canvas.draw()
                renderer = fig.canvas.get_renderer()
                crop = self._crop(renderer)
            # what is left of the drift (a pixel or two) comes off the padding
            excess = crop.width - (width_px + 0.25) / dpi
            crop = Bbox.from_extents(crop.x0 + excess / 2, crop.y0,
                                     crop.x1 - excess / 2, crop.y1)
        reference = self._snapshot(dpi, crop)

        # every state, in order: all groups hidden, then phase by phase
        self._set(self.groups, False)
        states = [self._snapshot(dpi, crop)]
        plan = []                       # (step index, phase, state before, after)
        for i, step in enumerate(self.steps):
            for phase in step.phases:
                before = len(states) - 1
                if phase.groups:
                    self._set(phase.groups, True)
                    states.append(self._snapshot(dpi, crop))
                plan.append((i, phase, before, len(states) - 1))
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        if not np.array_equal(states[-1], reference):
            differing = float(np.mean(np.any(states[-1] != reference, axis=2)))
            message = (f'the final state differs from the complete figure in '
                       f'{100 * differing:.3f} % of the pixels')
            if strict:
                raise RuntimeError(message)
            print('  WARNING ' + message)
        shape = reference.shape[:2]

        quantize = _Quantizer(_shared_palette(states), dither)
        frame_s = 1.0 / fps
        written, lines = [], []
        everything = []                 # [frame, seconds] of the combined GIF

        os.makedirs(out_dir, exist_ok=True)
        _remove_previous(out_dir, stem)
        if stills:
            start = os.path.join(out_dir, f'{stem}_step00_start.png')
            Image.fromarray(states[0]).save(start)
        for i, step in enumerate(self.steps):
            # [frame, seconds]; the first frame is the state the step starts
            # from, held for the lead-in plus any pause that opens the step
            shots = [[quantize(states[plan_before(plan, i)]), 0.0]]
            for index, phase, before, after in plan:
                if index != i:
                    continue
                seconds = phase.seconds / speed
                if phase.kind == 'pause':
                    shots[-1][1] += seconds
                    continue
                a = states[before].astype(np.float32)
                b = states[after].astype(np.float32)
                n = 0 if phase.kind == 'cut' else max(
                    1, int(round(seconds * fps)))
                for k in range(1, n):
                    m = self._mask(phase, k / n, shape, crop, renderer)
                    blend = a + (b - a) * m[:, :, None]
                    shots.append([quantize(
                        np.clip(blend + 0.5, 0, 255).astype(np.uint8)),
                        frame_s])
                shots.append([quantize(states[after]), frame_s])
            opening_pause = shots[0][1] if len(shots) > 1 else 0.0
            name = f'{stem}_step{i + 1:02d}_{step.slug}'
            path = os.path.join(out_dir, name + '.gif')
            frames = [frame for frame, _ in shots]
            seconds = [lead_in + opening_pause] + [s for _, s in shots[1:]]
            seconds[-1] = end_hold
            count = _save_gif(path, frames, seconds, loop)
            written.append(path)
            if stills:
                Image.fromarray(states[plan_after(plan, i)]).save(
                    os.path.join(out_dir, name + '.png'))
            lines.append(f'{i + 1:02d}  {name}.gif  {step.caption}')
            if verbose:
                print(f'wrote {path} ({count} frames, '
                      f'{os.path.getsize(path) / 1e6:.1f} MB)')
            if combined:
                last = i == len(self.steps) - 1
                hold = combined_pause * (2.0 if last else 1.0)
                if not everything:
                    everything.append([frames[0], lead_in])
                # a later step's first frame repeats the frame already there,
                # which takes over the pause that opens (or is) the step
                everything[-1][1] += shots[0][1]
                everything.extend([frame, s] for frame, s in shots[1:])
                if len(shots) > 1:
                    everything[-1][1] = hold
                elif last:
                    everything[-1][1] = max(everything[-1][1], hold)
        if combined:
            path = os.path.join(out_dir, f'{stem}_all_steps.gif')
            count = _save_gif(path, [frame for frame, _ in everything],
                              [s for _, s in everything], True)
            if verbose:
                print(f'wrote {path} ({count} frames, '
                      f'{os.path.getsize(path) / 1e6:.1f} MB)')
        with open(os.path.join(out_dir, f'{stem}_steps.txt'), 'w',
                  encoding='utf-8') as handle:
            handle.write('\n'.join(lines) + '\n')
        return written


def _remove_previous(out_dir, stem):
    """Delete what an earlier render() of this stem wrote, so a changed
    storyboard leaves no outdated step files behind."""
    pattern = re.compile(re.escape(stem) + r'_(step\d\d_.*\.(gif|png)'
                         r'|all_steps\.gif|steps\.txt)$')
    for name in os.listdir(out_dir):
        if pattern.match(name):
            os.remove(os.path.join(out_dir, name))


def plan_before(plan, step_index):
    """Index of the state a step starts from."""
    return next(before for index, _, before, _ in plan if index == step_index)


def plan_after(plan, step_index):
    """Index of the state a step ends on."""
    return [after for index, _, _, after in plan if index == step_index][-1]


# -----------------------------------------------------------------------------
# palette and GIF writing
# -----------------------------------------------------------------------------
def _shared_palette(states, colours=256, sample_cap=1_500_000):
    """One palette, an (n, 3) uint8 array, for every frame of the figure.

    Built from the distinct colours of all states and of the blends between
    consecutive states (what fades and wipe fronts pass through), each weighted by the square root of its pixel
    count: thin coloured lines and markers keep accurate colours next to the
    large flat areas. The most common colour (the background) and pure black
    are then stored exactly.
    """
    samples = [state[::2, ::2].reshape(-1, 3) for state in states]
    for a, b in zip(states[:-1], states[1:]):
        a, b = a[::3, ::3].astype(np.uint16), b[::3, ::3].astype(np.uint16)
        for weight in (1, 2, 3):      # quarter, half, three-quarter blends
            blend = (a * (4 - weight) + b * weight) // 4
            samples.append(blend.astype(np.uint8).reshape(-1, 3))
    pixels = np.concatenate(samples).astype(np.uint32)
    packed = (pixels[:, 0] << 16) | (pixels[:, 1] << 8) | pixels[:, 2]
    unique, counts = np.unique(packed, return_counts=True)
    weights = np.sqrt(counts)
    repeats = np.maximum(
        1, np.round(weights * min(1.0, sample_cap / weights.sum()))
    ).astype(np.int64)
    expanded = np.repeat(unique, repeats)
    rgb = np.stack([(expanded >> 16) & 255, (expanded >> 8) & 255,
                    expanded & 255], axis=1).astype(np.uint8)
    side = int(math.ceil(math.sqrt(len(rgb))))
    tile = np.zeros((side * side, 3), np.uint8)
    tile[:len(rgb)] = rgb
    tile[len(rgb):] = rgb[-1]
    quantized = Image.fromarray(tile.reshape(side, side, 3)).quantize(
        colors=colours - 2, method=Image.Quantize.MEDIANCUT, kmeans=1)
    table = np.asarray(quantized.getpalette()[:3 * (colours - 2)],
                       np.uint8).reshape(-1, 3)
    background = unique[np.argmax(counts)]
    exact = np.array([[(background >> 16) & 255, (background >> 8) & 255,
                       background & 255], [0, 0, 0]], np.uint8)
    return np.concatenate([exact, table])


class _Quantizer:
    """Maps RGB frames onto the shared palette.

    Every colour goes to its nearest palette entry exactly (Pillow's own
    palette mapping works on a coarse colour grid and can move the background
    off white). The colour -> entry table is filled as colours are met, so
    later frames cost little. With dither, Pillow's Floyd-Steinberg mapping
    is used instead, with the background and black put back exactly (entries
    0 and 1 of the palette).
    """

    def __init__(self, table, dither=False):
        self.table = np.asarray(table, np.uint8)
        self.dither = dither
        self.lookup = np.full(1 << 24, -1, np.int16)
        self.palette = Image.new('P', (1, 1))
        self.palette.putpalette(self.table.reshape(-1).tolist())

    def __call__(self, rgb):
        if self.dither:
            entries = np.array(Image.fromarray(rgb).quantize(
                palette=self.palette, dither=Image.Dither.FLOYDSTEINBERG))
            # error diffusion must not tint the background or solid black
            for entry in (0, 1):
                entries[np.all(rgb == self.table[entry], axis=2)] = entry
            return self._image(entries)
        packed = ((rgb[:, :, 0].astype(np.uint32) << 16)
                  | (rgb[:, :, 1].astype(np.uint32) << 8) | rgb[:, :, 2])
        entries = self.lookup[packed]
        missing = entries < 0
        if missing.any():
            new = np.unique(packed[missing])
            colours = np.stack([(new >> 16) & 255, (new >> 8) & 255,
                                new & 255], axis=1).astype(np.int32)
            table = self.table.astype(np.int32)
            for start in range(0, len(new), 20000):
                block = colours[start:start + 20000]
                distance = ((block[:, None, :] - table[None, :, :]) ** 2
                            ).sum(axis=2)
                self.lookup[new[start:start + 20000]] = distance.argmin(axis=1)
            entries = self.lookup[packed]
        return self._image(entries)

    def _image(self, entries):
        image = Image.fromarray(entries.astype(np.uint8))
        image.putpalette(self.table.reshape(-1).tolist())   # 'L' becomes 'P'
        return image


def _milliseconds(seconds):
    """Frame durations in ms on the GIF's 10 ms grid. The rounding error is
    carried forward, so the total stays right (20 fps is exact; 30 fps
    alternates 30 and 40 ms). No frame is shorter than 20 ms, below which
    players fall back to a slow default."""
    out, elapsed, shown = [], 0.0, 0
    for s in seconds:
        elapsed += s
        ms = max(20, int(round(elapsed * 100)) * 10 - shown)
        out.append(ms)
        shown += ms
    return out


def _save_gif(path, frames, seconds, loop):
    """Write the frames; consecutive identical frames become one, shown for
    their summed time. Returns the number of frames stored."""
    kept, times = [], []
    for frame, s in zip(frames, seconds):
        data = frame.tobytes()
        if kept and data == kept[-1][1]:
            times[-1] += s
        else:
            kept.append((frame, data))
            times.append(s)
    images = [frame for frame, _ in kept]
    kwargs = dict(save_all=True, append_images=images[1:],
                  duration=_milliseconds(times), optimize=False, disposal=1)
    if loop:
        kwargs['loop'] = 0          # forever; omitted = play once
    images[0].save(path, **kwargs)
    return len(images)


# -----------------------------------------------------------------------------
# shared command line
# -----------------------------------------------------------------------------
def add_cli_arguments(parser):
    """The options every figure's GIF script takes."""
    parser.add_argument('--out-dir', default=None,
                        help='where the GIFs go')
    parser.add_argument('--width-px', type=int, default=1600,
                        help='frame width in pixels (default 1600)')
    parser.add_argument('--fps', type=float, default=20.0,
                        help='frames per second, at most 50 (default 20)')
    parser.add_argument('--speed', type=float, default=1.0,
                        help='factor on the animation speed (2 = twice as fast)')
    parser.add_argument('--loop', action='store_true',
                        help='loop each step GIF forever instead of playing '
                             'it once')
    parser.add_argument('--end-hold', type=float, default=None,
                        help='seconds the last frame of a step GIF is held '
                             '(default 60, for players that loop regardless; '
                             '2 with --loop)')
    parser.add_argument('--no-combined', dest='combined', action='store_false',
                        help='skip the all-steps GIF')
    parser.add_argument('--dither', action='store_true',
                        help='dither the frames (smoother gradients, larger '
                             'files)')
    return parser


def render_kwargs(args):
    """Storyboard.render keyword arguments from the parsed shared options."""
    return dict(width_px=args.width_px, fps=args.fps, speed=args.speed,
                loop=args.loop, end_hold=args.end_hold,
                combined=args.combined, dither=args.dither)
