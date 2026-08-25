# Quickstart demo — source

Source of truth for the interactive **quickstart demo player** that is embedded
on the docs landing page (`docs/source/index.rst`) and linked, via a poster,
from the project `README.rst`.

Nothing here is shipped or imported at runtime; these are authoring tools. The
*built* artifacts they produce live elsewhere in the tree (see below) and are
what gets committed and served.

## Files

| File | Role |
| --- | --- |
| `quickstart_demo_template.html` | The player itself (HTML/CSS/JS) with image **placeholders** (`__LOGO_LIGHT__`, `__FIG_KINETICS__`, …). Edit the demo here. |
| `build_demo.py` | Fills the placeholders and writes the built HTML (two modes). |
| `make_poster.py` | Renders the README poster PNG from real repo figures. |
| `make_loop_gif.py` | Renders the animated landing-page loop GIFs (light + dark) and their reduced-motion stills. |

## Regenerate

Run with the project interpreter (the `IBO_2026` env — matplotlib is needed for
the poster and the loop GIFs; the loop GIFs also need Pillow):

```bash
# docs copy: images referenced by relative path (small; embedded via <iframe>)
python docs/_demo_src/build_demo.py
#   -> docs/source/_static/quickstart_demo.html

# README poster (a *.png, so gitignored — commit it with `git add -f`)
python docs/_demo_src/make_poster.py
#   -> docs/source/_static/images/examples/quickstart_demo_poster.png

# landing-page loop GIFs + their frame-0 reduced-motion stills
# (GIFs commit normally; the stills are *.png, so gitignored — `git add -f` them)
python docs/_demo_src/make_loop_gif.py
#   -> docs/source/_static/images/demo/loop_light.gif
#   -> docs/source/_static/images/demo/loop_dark.gif
#   -> docs/source/_static/images/demo/loop_light_still.png
#   -> docs/source/_static/images/demo/loop_dark_still.png
```

The stills are what the landing page serves under
`prefers-reduced-motion: reduce` (the `<picture>` blocks in
`docs/source/index.rst`), so regenerate and re-commit them together with the
GIFs.

### Self-contained build (the shareable Artifact)

For a standalone page with every image inlined as a base64 `data:` URI (no
external files to load):

```bash
python docs/_demo_src/build_demo.py --self-contained quickstart_demo.selfcontained.html
```

This is the version published as a Claude Artifact; it is intentionally **not**
committed (the docs copy is the in-repo one).

## Notes

- The demo reuses the quickstart figures already under
  `_static/images/examples/` (`tutorial_01_quickstart_*.png`) and the logos
  under `_static/images/logo/` — no duplicate assets.
- Every number and figure shown in the demo is real quickstart output; keep it
  in sync with `docs/source/tutorial/01_quickstart.rst`, `index.rst`, and the
  `README.rst` quickstart if those change.
- Non-ASCII glyphs are entity-encoded by `build_demo.py` so the page renders
  correctly no matter what charset the host serves.
