# EDFig4 comparison

- output/k_3_try_informed_sweep_etoh_yield.png (1716x1146, 300 dpi) is byte-identical
  (md5 d2e524b90f6be95a372a404f2b697858) to
  publication/_archive_pre_reorg/Objective-landscape/k_3_try_informed_sweep_etoh_yield.png
  (formerly publication/Objective-landscape/), to manuscript/raw_panel_k_3_..._etoh_yield.png
  and to the Bioindustrial-Park copy of it (mtime 2026-09-25 05:42).
- The manuscript image (docx image10 = Figures\Extended_Data_Fig_2.png, md5 60b5dc50..., 1888x1260,
  330 dpi) is that render placed in PowerPoint and re-exported. Measured from the red/blue/purple
  spine columns and the axes top/bottom rows: scale 1.100 (= 330/300), offset about 0 px. The
  bottom legend "Growth derated by enzyme burden (from 4.85)" is cropped away, and a rotated text
  label "growth derated at higher k3 values" is added next to the dotted line.
- Pixel diff (regen LANCZOS-upscaled to 1888x1260 vs manuscript, RGB composited on white):
  - whole image: max 255, mean 1.60, 0.69 % px > 32 (all in the legend row and the added label);
  - masked (rows >= 1100 and the label box x 1135-1245, y 520-900 excluded): max 20, mean 0.28,
    0.000 % px > 32, i.e. resampling noise only. These numbers depend on the resampling filter:
    an independent reviewer's masked diff (same 1.1x scale, 0 offset, different filter) gave
    max 136, mean 1.19, 1.3 % px > 32 and 0.28 % > 64, all on antialiased edges, with zoomed
    crops of the PI dip and the right-end curves visually identical;
  - curve masks (PI red, ethanol-yield blue, IRR purple): IoU 0.988-0.993; 100 % of the
    manuscript's curve pixels lie within 2 px of the regenerated curve.
- Files: absdiff_x3_manuscript_vs_regen_upscaled1.1.png (inverted |diff| x3) and
  regen_upscaled_1888x1260.png.
- Verdict: MATCH-with-allowed-diffs (legend removed, annotation added in PowerPoint).
