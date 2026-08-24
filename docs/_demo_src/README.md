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

## Regenerate

Run with the project interpreter (the `IBO_2026` env — matplotlib is only needed
for the poster):

```bash
# docs copy: images referenced by relative path (small; embedded via <iframe>)
python docs/_demo_src/build_demo.py
#   -> docs/source/_static/quickstart_demo.html

# README poster (a *.png, so gitignored — commit it with `git add -f`)
python docs/_demo_src/make_poster.py
#   -> docs/source/_static/images/examples/quickstart_demo_poster.png
```

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
