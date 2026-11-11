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
  the file (PowerPoint, browsers) play it once and rest on its last frame.
  Players that loop regardless (Google Slides, Keynote) restart only after
  the last frame's hold, ``end_hold`` seconds (default 60).
* Every GIF of a figure shares one 256-colour palette, built from all the
  states, so a step's first frame equals the previous step's last frame
  exactly. The background and pure black are kept exact.

Fonts: the plot scripts ask for Arial. Where Arial is not installed (most
Linux machines), use_arial_or_metric_twin() registers a metric-compatible
face (Liberation Sans, Arimo, ...) under the name Arial, so the layout tuned
for Arial holds; call it before building the figure.

Sim-safe: imports only the standard library, numpy, matplotlib and Pillow.
"""
import dataclasses
import io
import math
import os
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
    switched off; the figure it created last is returned. The plot script
    itself stays unmodified.
    """
    before = set(plt.get_fignums())
    savefig, close = Figure.savefig, plt.close
    Figure.savefig = lambda self, *a, **k: None
    plt.close = lambda *a, **k: None
    try:
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
            if horizontal:
                mask[r0:r1, c0:c1] = ramp[None, :]
            else:
                # rows run top to bottom; 'up' already reversed `along`
                mask[r0:r1, c0:c1] = ramp[:, None]
        return mask

    def render(self, out_dir, stem, width_px=1600, fps=20, speed=1.0,
               loop=False, lead_in=0.4, end_hold=60.0, combined=True,
               combined_pause=1.5, dither=False, stills=True, strict=True,
               verbose=True):
        """Write the step GIFs (and stills, the combined GIF and the step
        list) into out_dir. Returns the list of step GIF paths.

        width_px        frame width; the height follows the figure.
        fps, speed      frame rate, and a factor on every phase duration
                        (speed=2 halves them).
        loop            False: play once and rest on the last frame; True:
                        loop forever.
        lead_in         seconds the first frame is held before motion starts.
        end_hold        seconds the last frame of a step GIF is held (matters
                        only where the player loops).
        combined        also write <stem>_all_steps.gif (looping), holding
                        combined_pause seconds after each step.
        dither          Floyd-Steinberg dithering of the frames (smoother
                        gradients, larger files).
        strict          raise if the final state differs from the complete
                        figure.
        """
        self._check()
        fig = self.fig
        os.makedirs(out_dir, exist_ok=True)
        dpi = float(width_px) / fig.get_size_inches()[0]
        fig.set_dpi(dpi)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        crop = self._crop(renderer)
        if self.savefig_kwargs.get('bbox_inches') == 'tight':
            # width_px is the width of the saved (cropped) frame
            dpi *= fig.get_size_inches()[0] / crop.width
            fig.set_dpi(dpi)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            crop = self._crop(renderer)
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
        written, lines, everything, everything_ms = [], [], [], []

        if stills:
            start = os.path.join(out_dir, f'{stem}_step00_start.png')
            Image.fromarray(states[0]).save(start)
        for i, step in enumerate(self.steps):
            frames = [quantize(states[plan_before(plan, i)])]
            durations = [lead_in]
            for index, phase, before, after in plan:
                if index != i:
                    continue
                seconds = phase.seconds / speed
                if phase.kind == 'pause':
                    durations[-1] += seconds
                    continue
                a = states[before].astype(np.float32)
                b = states[after].astype(np.float32)
                n = 0 if phase.kind == 'cut' else max(
                    1, int(round(seconds * fps)))
                for k in range(1, n):
                    m = self._mask(phase, k / n, shape, crop, renderer)
                    blend = a + (b - a) * m[:, :, None]
                    frames.append(quantize(
                        np.clip(blend + 0.5, 0, 255).astype(np.uint8)))
                    durations.append(frame_s)
                frames.append(quantize(states[after]))
                durations.append(frame_s)
            name = f'{stem}_step{i + 1:02d}_{step.slug}'
            path = os.path.join(out_dir, name + '.gif')
            step_ms = _milliseconds(durations)
            _save_gif(path, frames, step_ms[:-1] + [int(end_hold * 1000)],
                      loop)
            written.append(path)
            if stills:
                Image.fromarray(states[plan_after(plan, i)]).save(
                    os.path.join(out_dir, name + '.png'))
            lines.append(f'{i + 1:02d}  {name}.gif  {step.caption}')
            if verbose:
                print(f'wrote {path} ({len(frames)} frames, '
                      f'{os.path.getsize(path) / 1e6:.1f} MB)')
            if combined:
                last = len(self.steps) - 1
                hold = combined_pause * (2.0 if i == last else 1.0)
                # the first frame of a later step repeats the frame before it
                skip = 0 if i == 0 else 1
                everything.extend(frames[skip:])
                everything_ms.extend(step_ms[skip:-1] + [int(hold * 1000)])
        if combined:
            path = os.path.join(out_dir, f'{stem}_all_steps.gif')
            _save_gif(path, everything, everything_ms, True)
            if verbose:
                print(f'wrote {path} ({len(everything)} frames, '
                      f'{os.path.getsize(path) / 1e6:.1f} MB)')
        with open(os.path.join(out_dir, f'{stem}_steps.txt'), 'w',
                  encoding='utf-8') as handle:
            handle.write('\n'.join(lines) + '\n')
        return written


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

    Built from the distinct colours of all states and of the half-way blends
    between consecutive states, each weighted by the square root of its pixel
    count: thin coloured lines and markers keep accurate colours next to the
    large flat areas. The most common colour (the background) and pure black
    are then stored exactly.
    """
    samples = [state[::2, ::2].reshape(-1, 3) for state in states]
    for a, b in zip(states[:-1], states[1:]):
        half = ((a[::3, ::3].astype(np.uint16) + b[::3, ::3]) // 2)
        samples.append(half.astype(np.uint8).reshape(-1, 3))
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
    is used instead.
    """

    def __init__(self, table, dither=False):
        self.table = np.asarray(table, np.uint8)
        self.dither = dither
        self.lookup = np.full(1 << 24, -1, np.int16)
        self.palette = Image.new('P', (1, 1))
        self.palette.putpalette(self.table.reshape(-1).tolist())

    def __call__(self, rgb):
        if self.dither:
            return Image.fromarray(rgb).quantize(
                palette=self.palette, dither=Image.Dither.FLOYDSTEINBERG)
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
        image = Image.fromarray(entries.astype(np.uint8))
        image.putpalette(self.table.reshape(-1).tolist())   # 'L' becomes 'P'
        return image


def _milliseconds(seconds):
    """Frame durations in ms, on the GIF's 10 ms grid."""
    return [max(20, int(round(s * 100)) * 10) for s in seconds]


def _save_gif(path, frames, durations_ms, loop):
    kwargs = dict(save_all=True, append_images=frames[1:],
                  duration=list(durations_ms), optimize=False, disposal=1)
    if loop:
        kwargs['loop'] = 0          # forever; omitted = play once
    frames[0].save(path, **kwargs)


# -----------------------------------------------------------------------------
# shared command line
# -----------------------------------------------------------------------------
def add_cli_arguments(parser):
    """The options every figure's GIF script takes."""
    parser.add_argument('--out-dir', default=None,
                        help='where the GIFs go')
    parser.add_argument('--width-px', type=int, default=1600,
                        help='frame width in pixels (default 1600)')
    parser.add_argument('--fps', type=float, default=20.0)
    parser.add_argument('--speed', type=float, default=1.0,
                        help='factor on the animation speed (2 = twice as fast)')
    parser.add_argument('--loop', action='store_true',
                        help='loop each step GIF forever instead of playing '
                             'it once')
    parser.add_argument('--end-hold', type=float, default=60.0,
                        help='seconds the last frame of a step GIF is held, '
                             'for players that loop regardless (default 60)')
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
