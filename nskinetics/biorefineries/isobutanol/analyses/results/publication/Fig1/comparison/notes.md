Fig1 comparison (2026-11-04)
- The manuscript image1.png == Figures/Figure1_framework_overview.png (md5 cb33a78ea98cc7453634e91edb639c6a).
- The source is the claude.ai Design artifact eceea534-5348-4131-ad8c-d64b91c638d5 v1793486932-ef0d. Its
  project/Main.dc.html is saved verbatim as code/Main.dc.html (sha256 matches the artifact listing).
- Regenerated output/Fig1.png: headless Chrome render, 1440x846 window, device scale factor 3 = 4320x2538.
- Whole-image diff vs the manuscript: mean abs 2.56, max 243, 3.36 % of pixels differ, 2.54 % by more than 32
  (diff_regen_vs_manuscript.png; black marks the >32 pixels).
- Best-alignment offsets (regen relative to the manuscript, device px), each with the mean abs residual after alignment:

  | Region | Offset | Residual |
  |---|---|---|
  | SVG text, icons, inset curve, seed arrow | dx=0, dy=0 | 0.2-1.7 |
  | HTML box titles, italic examples, inset captions, Inform div | dy=+1 | 0.9-2.3 |
  | HTML 19 px stage headings (top and bottom) | dy=+5 | 0.6-1.0 |
  | HTML 16 px inset title | dy=+4 | 1.0 |

- The step-1 dashed border differs in dash phase and length. The numbered circles' borders antialias differently.
- Crops (manuscript on top, regen below): crop_heading_*, crop_box1_*, crop_inset_*.
- Verdict: MATCH-with-allowed-diffs. The geometry and text content are identical, and the glyphs are the same Arimo.
  Some HTML text lines sit 1-5 device px lower, which is attributed to the Design runtime CSS or browser layout, not to content.
- Earlier version: final vs OLD_Figure1_framework_overview.png (2026-10-02). 2.08 % of pixels differ, all in text
  glyphs (diff_final_vs_OLD.png); the step-5 example text was reworded.
