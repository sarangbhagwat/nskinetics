# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Build the NSKinetics quickstart interactive demo from its template.

The template ``quickstart_demo_template.html`` (in this directory) is the single
source of truth for the demo player. It carries image placeholders
(``__LOGO_LIGHT__`` etc.) that this script fills in two different ways:

* ``--docs`` (default) writes ``docs/source/_static/quickstart_demo.html`` with
  every image referenced by a *relative* path. The quickstart figures already
  ship under ``_static/images/``, so the committed file stays tiny (~29 KB) and
  embeds via ``<iframe>`` on the docs landing page / Read the Docs.

* ``--self-contained PATH`` writes a standalone HTML with every image inlined as
  a base64 ``data:`` URI. That is the version published as the shareable
  Artifact (nothing external to load).

Non-ASCII glyphs in the template (``τ``, ``≈``, ``→`` …) are
entity-encoded on the way out so the page renders identically regardless of the
charset the host serves it with.

Usage
-----
::

    python docs/_demo_src/build_demo.py                 # -> docs/source/_static/quickstart_demo.html
    python docs/_demo_src/build_demo.py --self-contained out.html

Run ``docs/_demo_src/make_poster.py`` to regenerate the README poster.
"""
import argparse
import base64
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent      # docs/_demo_src
ROOT = HERE.parents[1]                               # repo root
STATIC = ROOT / "docs" / "source" / "_static"
TEMPLATE = HERE / "quickstart_demo_template.html"
DOCS_OUT = STATIC / "quickstart_demo.html"

# placeholder token -> image path relative to the _static directory
IMAGES = {
    "__LOGO_LIGHT__":    "images/logo/logo_nskinetics_light.png",
    "__LOGO_DARK__":     "images/logo/logo_nskinetics_dark.png",
    "__FIG_FLOWSHEET__": "images/examples/tutorial_01_quickstart_flowsheet.png",
    "__FIG_KINETICS__":  "images/examples/tutorial_01_quickstart_kinetics.png",
    "__FIG_RETUNED__":   "images/examples/tutorial_01_quickstart_kinetics_retuned.png",
}

# non-ASCII -> HTML entity, so the page is charset-independent
_ENTITIES = {
    "—": "&mdash;", "–": "&ndash;", "τ": "&tau;", "≈": "&asymp;",
    "×": "&times;", "³": "&sup3;", "·": "&middot;", "→": "&rarr;",
    "←": "&larr;", "’": "&rsquo;", "“": "&ldquo;", "”": "&rdquo;",
    "≤": "&le;", "≥": "&ge;", "…": "&hellip;",
}


def _entity_encode(html):
    for ch, ent in _ENTITIES.items():
        html = html.replace(ch, ent)
    leftover = sorted({hex(ord(c)) for c in html if ord(c) > 127})
    if leftover:
        raise ValueError(f"un-encoded non-ASCII characters remain: {leftover}")
    return html


def _fill(replace):
    """Load the template and apply ``replace(token) -> str`` for each image."""
    html = TEMPLATE.read_text(encoding="utf-8")
    for token, rel in IMAGES.items():
        src = STATIC / rel
        if not src.exists():
            raise FileNotFoundError(f"referenced image is missing: {src}")
        if token not in html:
            raise ValueError(f"template is missing placeholder: {token}")
        html = html.replace(token, replace(token, rel, src))
    return _entity_encode(html)


def build_docs(out=DOCS_OUT):
    """Write the docs copy (relative image paths)."""
    html = _fill(lambda token, rel, src: rel)
    out.write_text(html, encoding="utf-8")
    # every remaining image ref must resolve on disk
    for rel in re.findall(r'"(images/[^"]+\.png)"', html):
        if not (STATIC / rel).exists():
            raise FileNotFoundError(f"unresolved image ref in output: {rel}")
    return out


def build_self_contained(out):
    """Write a standalone copy with images inlined as base64 data URIs."""
    def as_data_uri(token, rel, src):
        return "data:image/png;base64," + base64.b64encode(src.read_bytes()).decode()

    out = pathlib.Path(out)
    out.write_text(_fill(as_data_uri), encoding="utf-8")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-contained", metavar="PATH",
                    help="write a standalone base64-inlined HTML to PATH "
                         "(the Artifact build) instead of the docs copy")
    args = ap.parse_args(argv)
    if args.self_contained:
        out = build_self_contained(args.self_contained)
    else:
        out = build_docs()
    print(f"wrote {out} ({out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
