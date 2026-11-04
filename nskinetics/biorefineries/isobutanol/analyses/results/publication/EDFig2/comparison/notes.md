# EDFig2 comparison

- md5(output/EDFig2_bars_stack.png) = 5444a50eaa61dcb16e4120b6369e0316
- md5(manuscript/EDFig2_manuscript.png) (= docx word/media/image8.png) = 5444a50eaa61dcb16e4120b6369e0316
- Also identical to the Bioindustrial-Park render
  `...\Bioindustrial-Park\biorefineries\isobutanol\analyses\results\kin_sobol_..._custom9f1df84b_aA_burden_seed20260922_sobol_bars_stack.png`
  (mtime 2026-09-30 22:57, = BP commit 532dc606 time).
- Verdict: MATCH (byte-identical, 4440x6595 px, 600 dpi). No pixel diff needed.
- Q2 in the panels come from the index table: PI 0.857, EtOH titer 0.891, EtOH yield 0.848,
  IBO titer 0.898, IBO yield 0.854 (all HGB surrogates, reliable).
  (two decimals as drawn: 0.86, 0.89, 0.85, 0.90, 0.85).
- Only the PNG is compared; the output PDFs carry a creation timestamp and differ in md5 per run.
