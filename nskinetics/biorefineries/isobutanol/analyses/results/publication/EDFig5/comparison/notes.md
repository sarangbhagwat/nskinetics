# EDFig5 comparison

- md5(output/EDFig5.png) = f620c56afef55b7beb091aa60737e2a6 (3450x4650, 300 dpi)
- md5(manuscript/EDFig5_manuscript.png) (= docx image11
  = Figures\tea_breakdowns_split12d_figure_2026.10.02-12.09.png
  = Bioindustrial-Park results\tea_breakdowns_split12d_figure_2026.10.02-12.09.png)
  = f620c56afef55b7beb091aa60737e2a6
- Verdict: MATCH (byte-identical).
- Note: the older stage-1 JSON tea_breakdowns_split12d_2026.09.24-16.52.json has no product
  revenue, and the current script refuses it. The 17.48 JSON is the right input.
- Re-checked 2026-11-04 after the polish edits to make_figure.py (PYTHONDONTWRITEBYTECODE,
  SOURCE_DATE_EPOCH): output/EDFig5.png md5 is still f620c56afef55b7beb091aa60737e2a6, and the
  nine breakdown CSVs are unchanged. The PDF md5 varied per run before (it carried the render-time
  CreationDate). It is now pinned to D:20261002190900Z.
- None of these files lives in the old theme folders (now
  publication/_archive_pre_reorg/<theme>/). The manuscript-side copy is in
  C:\Users\saran\Documents\Academia\Manuscripts_Postdoc\01_NSKinetics\Figures and the BP results root.
