// NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
// Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
// MIT open-source license -- see https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
//
// Bundle entry for build_amcharts.mjs. amCharts5 ships ES modules (for
// bundlers), not browser UMD/script bundles, so a classic <script> tag on
// node_modules/@amcharts/amcharts5/index.js cannot expose the `am5` global.
// esbuild bundles this entry into an IIFE (vendor/am5.bundle.js) that
// attaches the amCharts core and hierarchy namespaces to `window`, which is
// what treemap.html loads.
import * as am5 from "@amcharts/amcharts5";
import * as am5hierarchy from "@amcharts/amcharts5/hierarchy";

window.am5 = am5;
window.am5hierarchy = am5hierarchy;
