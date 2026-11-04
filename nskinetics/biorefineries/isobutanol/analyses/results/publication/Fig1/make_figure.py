# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Regenerate manuscript Figure 1 (conceptual framework overview).

The source is the claude.ai Design artifact "Figure 1 - NSKinetics framework"
(slug eceea534-5348-4131-ad8c-d64b91c638d5, version 1793486932-ef0d). Its board
file project/Main.dc.html is saved verbatim as code/Main.dc.html.

This script:
1. derives code/Fig1_standalone.html from code/Main.dc.html by stripping the
   Design-tool wrapper (support.js <script>, <x-dc>, <helmet>, and the
   text/x-dc component script), keeping the Arimo Google Fonts link and the
   body style;
2. renders it with a local headless Chrome (or Edge) at a 1440 x 846 CSS-px
   viewport with device scale factor 3, i.e. 4320 x 2538 px, to output/Fig1.png.

Sim-safe: no Python package beyond the standard library is imported. The
Arimo web font is fetched by the browser from fonts.googleapis.com; offline,
the page falls back to Arial (metric-compatible), which changes glyph shapes
slightly. Works from any cwd.

For uniformity with the other figure folders, the child environment gets
NSKINETICS_IBO_RESULTS_DIR = <Fig1>/data (empty: Figure 1 has no data inputs),
MPLBACKEND=Agg and PYTHONDONTWRITEBYTECODE=1; the browser ignores all three.
"""
import os
import re
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True   # leave no __pycache__ in the figure folder

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, 'code')
OUT = os.path.join(HERE, 'output')
DATA = os.path.join(HERE, 'data')
SRC = os.path.join(CODE, 'Main.dc.html')
STANDALONE = os.path.join(CODE, 'Fig1_standalone.html')
PNG = os.path.join(OUT, 'Fig1.png')
W, H, SCALE = 1440, 846, 3

BROWSERS = (
    r'C:\Program Files\Google\Chrome\Application\chrome.exe',
    r'C:\Program Files (x86)\Google\Chrome\Application\chrome.exe',
    r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe',
    r'C:\Program Files\Microsoft\Edge\Application\msedge.exe',
)


def build_standalone():
    with open(SRC, encoding='utf-8') as f:
        html = f.read()
    head = re.search(r'<helmet>(.*?)</helmet>', html, re.S).group(1).strip()
    body = re.search(r'</helmet>(.*?)</x-dc>', html, re.S).group(1).strip()
    page = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<title>Figure 1 - Framework overview (standalone)</title>\n'
            f'{head}\n</head>\n<body>\n{body}\n</body>\n</html>\n')
    with open(STANDALONE, 'w', encoding='utf-8', newline='\n') as f:
        f.write(page)
    return STANDALONE


def find_browser():
    for b in BROWSERS:
        if os.path.exists(b):
            return b
    raise FileNotFoundError('no local Chrome/Edge found; tried: ' + ', '.join(BROWSERS))


def render(html_path):
    os.makedirs(OUT, exist_ok=True)
    profile = os.path.join(OUT, '_chrome_profile')
    url = 'file:///' + html_path.replace('\\', '/')
    cmd = [find_browser(), '--headless=new', '--disable-gpu', '--hide-scrollbars',
           '--no-first-run', '--no-default-browser-check',
           f'--user-data-dir={profile}',
           f'--force-device-scale-factor={SCALE}', f'--window-size={W},{H}',
           '--virtual-time-budget=10000', '--run-all-compositor-stages-before-draw',
           f'--screenshot={PNG}', url]
    env = dict(os.environ, NSKINETICS_IBO_RESULTS_DIR=DATA, MPLBACKEND='Agg',
               PYTHONDONTWRITEBYTECODE='1')
    print('browser:', cmd[0])
    subprocess.run(cmd, check=True, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(20):          # the browser may hold the profile briefly
        try:
            shutil.rmtree(profile)
            break
        except OSError:
            time.sleep(0.5)
    if not os.path.exists(PNG):
        sys.exit('render failed: no ' + PNG)
    print('->', PNG)


def main():
    render(build_standalone())


if __name__ == '__main__':
    main()
