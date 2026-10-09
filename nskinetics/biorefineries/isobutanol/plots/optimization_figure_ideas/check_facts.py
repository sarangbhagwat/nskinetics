# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Offline test of the optimization-figure folder (sim-safe; no load(), no
simulation): recompute every fact the figures quote from the campaign CSVs
(`_common.compute_facts()`), print each computed-vs-expected value, assert
them against `_common.EXPECTED`, confirm the data scope (only the eight
campaigns of 2026-09-23/24 are registered), scan the figure scripts
(main, S2) for hard-coded annotation numbers (spec 6.A.7), check the
palette (`_style.cvd_check()`) and confirm no simulation package was
imported. Exit 0 = all checks passed,
1 = any mismatch.

When the campaign / Sobol' CSVs are absent (only
analyses/results/publication ships with the repository) the facts are
SKIPPED with a message naming the missing inputs and the analyses that write
them; the data-independent checks (data scope, literal scan, palette,
sim-safety of loading pk / ko / eb) still run, and the exit code is 0 if they
pass.

Usage::

    python check_facts.py [--failures-only] [--no-literal-scan] [--no-palette]
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common as C    # noqa: E402


def _fmt(v):
    return C._fmt_val(v)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--failures-only', action='store_true',
                    help='print only the facts that differ from EXPECTED')
    ap.add_argument('--no-literal-scan', action='store_true',
                    help='skip the hard-coded-literal scan of the figure '
                         'scripts')
    ap.add_argument('--no-palette', action='store_true',
                    help='skip the colour-vision / grayscale palette check')
    args = ap.parse_args(argv)

    # a fresh process (script) gets the strict process-wide sim-safety check;
    # a shared one (e.g. pytest) checks only what this run imports
    before = None if not C._PRELOADED else frozenset(sys.modules)
    t0 = time.time()
    failures = []
    missing = C.missing_inputs()
    if missing:
        print(f'facts: SKIPPED -- {len(missing)} input file(s) missing:')
        for path, produced_by in missing:
            print(f'  {path}\n      (produced by {produced_by})')
        print(f'  run those analyses first, or set {C._paths.RESULTS_DIR_ENV} '
              'to a results folder that has them')
        with C.sim_safe_block('loading pk / ko / eb'):
            C._pk()
    else:
        with C.sim_safe_block('compute_facts()'):
            facts = C.compute_facts()
        rows = C.fact_table(facts)
        n_bad = 0
        print(f'{"":4s} {"fact":60s} {"computed":>14s} {"expected":>14s}  tol')
        for r in rows:
            if not r.ok:
                n_bad += 1
            if args.failures_only and r.ok:
                continue
            flag = 'ok  ' if r.ok else 'FAIL'
            note = f'  [{r.note}]' if r.note else ''
            print(f'{flag} {r.path:60s} {_fmt(r.computed):>14s} '
                  f'{_fmt(r.expected):>14s}  {r.tol:.2g}{note}')
        print(f'\nfacts: {len(rows) - n_bad} of {len(rows)} match EXPECTED '
              f'({time.time() - t0:.1f} s)')
        if n_bad:
            failures.append(f'{n_bad} fact mismatch(es)')

    # data scope: only this work's eight campaigns (run 2026-09-23/24; the
    # seed-350 uninformed + six scouts and the _rl15c111dc relay)
    allowed = ('_rs350_burden', '_rl15c111dc_burden')
    stray = [k for k, c in C.CAMPAIGNS.items()
             if not c.stem.endswith(allowed)]
    print(f'\ndata scope: {len(C.CAMPAIGNS)} campaigns '
          f'({", ".join(C.CAMPAIGNS)}) + the relay manifest; '
          f'{len(stray)} outside it')
    if len(C.CAMPAIGNS) != 8 or stray:
        failures.append(f'data scope: {len(C.CAMPAIGNS)} campaigns, '
                        f'stray {stray}')

    if not args.no_literal_scan:
        hits = C.literal_scan()
        present = [f for f in C.FIGURE_SCRIPTS
                   if os.path.exists(os.path.join(C.FIGDIR, f))]
        print(f'\nliteral scan of {present or "no figure scripts yet"}: '
              f'{len(hits)} hit(s) of {list(C.FORBIDDEN_LITERALS)}')
        for f, line, snip in hits:
            print(f'  {f}:{line}: {snip!r}')
        if hits:
            failures.append(f'{len(hits)} hard-coded annotation literal(s)')

    if not args.no_palette:
        import _style as S                 # matplotlib (Agg) only
        ok, report = S.cvd_check(raise_on_fail=False, verbose=True)
        if not ok:
            failures.append('palette check')

    try:
        C.assert_sim_safe(since=before)
        print('\nsim-safety: no forbidden module imported '
              f'({", ".join(C.FORBIDDEN_MODULES[:5])}, the isobutanol '
              'build modules)' + ('' if before is None else
                                  ' by this run (shared process)'))
    except RuntimeError as e:
        failures.append(str(e))

    if failures:
        print('\nCHECK FACTS FAILED: ' + '; '.join(failures))
        return 1
    if missing:
        print('\nCHECKS PASSED (facts SKIPPED: campaign / Sobol\' CSVs absent)')
        return 0
    print('\nALL CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
