# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
from setuptools import setup, find_packages

setup(
    name='nskinetics',

    packages=find_packages(),
    license='MIT',
    version='0.6.0',
    description='Simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena',
    long_description=open('README.rst', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    author='Sarang S. Bhagwat',
    install_requires=['IPython>=7.9.0',
                      # Floors below align with thermosteam 0.53.5 (pulled in
                      # by biosteam 2.53.11); the stack is validated on
                      # numpy 2.5.2 / scipy 1.18.0 / pandas 3.0.5 /
                      # matplotlib 3.11.1 / Python 3.14.7.
                      'numpy>=2.3.5,<3.0.0',
                      'numba>=0.62.1,<1.0.0',
                      'scipy>=1.15.2,<2.0.0',
                      'matplotlib>=3.10.1,<4.0.0',
                      'xlsxwriter>=3.2.5,<4.0.0',
                      'scikit-learn',
                      'pandas>=2.2.2,<4.0.0',
                      'tellurium',
                      'biosteam>=2.53.11,<2.55.0',
                      # 'python-libsbml>=5.20.5,<6.0.0',
		     ],
    extras_require={
        # Optional dependencies of the isobutanol biorefinery case study
        # (nskinetics.biorefineries.isobutanol); `import nskinetics` never
        # needs them. Install with: pip install "nskinetics[isobutanol]".
        # `biorefineries` (Bioindustrial-Park: the corn / cellulosic / ethanol
        # modules the flowsheet reuses) and `hensmith` (heat-exchanger-network
        # synthesis; requires biosteam>=2.54.0) are both on PyPI, but the
        # locally validated stack imports the Bioindustrial-Park
        # (biorefineries 2.35.1) and hensmith (0.1.2) git clones from
        # PYTHONPATH, not these PyPI releases.
        'isobutanol': [
            'biorefineries',
            'hensmith',
            'contourplots>=0.4.0',
            'optuna>=4.9',
            'SALib>=1.5',
            'openpyxl',  # pandas.read_excel engine for the parameter workbooks
        ],
    },
    include_package_data=True,
    package_data={
        # Shipped model files loaded by name at import/run time (e.g. the
        # isobutanol model module loads its antimony .txt on import, and the
        # tutorial loads the shipped SBML .xml). These must be in the wheel so
        # installed (non-editable) copies can find them.
        'nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo': ['*.txt', '*.xml'],
        # Reference antimony used by the events tests, plus the isobutanol
        # tests' fixture and the script-style checks their wrappers run in
        # fresh processes (scripts/ has no __init__.py, so it ships as data).
        'nskinetics.tests': ['data/*.txt',
                             'isobutanol/fixtures/*.csv',
                             'isobutanol/scripts/*.py'],
        # Isobutanol biorefinery case study: the parameter-distribution
        # workbooks read at run time (the 'parameter-distributions_*' glob
        # never matches Excel '~$' lock files), small JSON/CSV inputs, the
        # analysis scripts (plain folders without __init__.py, so they ship
        # as data) and the plotting extras. analyses/results/ (including the
        # ~87 MB local-only publication set) is never shipped: see
        # exclude_package_data below and the prunes in MANIFEST.in.
        'nskinetics.biorefineries.isobutanol': [
            'data/parameter_distributions/parameter-distributions_*.xlsx',
            'data/**/*.json',
            'analyses/**/*.py',
            'analyses/**/*.json',
            'analyses/**/*.md',
            'plots/optimization_figure_ideas/*.md',
            'plots/optimization_figure_ideas/*.py',
            'plots/voronoi/package.json',
            'plots/voronoi/package-lock.json',
            'plots/voronoi/*.mjs',
            'plots/voronoi/*.js',
            'plots/voronoi/treemap.html',
        ],
    },
    # Belt and braces: whatever the globs above or MANIFEST.in match, never
    # ship analysis outputs, node tooling, Excel lock files or bytecode.
    # (These are fnmatch patterns, so '*' also spans folders.) plots/ is a
    # package of its own, so it gets its own key.
    exclude_package_data={
        'nskinetics.biorefineries.isobutanol': [
            'analyses/results/*',
            'analyses/*/results/*',
            'plots/results/*',
            'plots/voronoi/node_modules/*',
            'plots/voronoi/vendor/*',
            '*~$*',
            '*__pycache__*',
        ],
        'nskinetics.biorefineries.isobutanol.plots': [
            'results/*',
            'voronoi/node_modules/*',
            'voronoi/vendor/*',
            '*~$*',
            '*__pycache__*',
        ],
    },
    python_requires='>=3.12',  # floor set by biosteam 2.53.11 / thermosteam 0.53.5
    platforms=['Windows', 'Mac', 'Linux'],
    author_email='sarangbhagwat.developer@gmail.com',
    url='https://github.com/sarangbhagwat/nskinetics',
    download_url='https://github.com/sarangbhagwat/nskinetics.git',
    classifiers=['License :: OSI Approved :: MIT License',
                 'Development Status :: 3 - Alpha',
                 'Environment :: Console',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Mathematics',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Intended Audience :: Developers',
                 'Intended Audience :: Education',
                 'Intended Audience :: Manufacturing',
                 'Intended Audience :: Science/Research',
                 'Natural Language :: English',
                 'Operating System :: MacOS',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Operating System :: POSIX :: BSD',
                 'Operating System :: POSIX :: Linux',
                 'Operating System :: Unix',
                 'Programming Language :: Python :: 3.12',
                 'Programming Language :: Python :: 3.13',
                 'Programming Language :: Python :: 3.14',
                 'Programming Language :: Python :: Implementation :: CPython',
                 'Topic :: Education'],
    keywords=['reaction kinetics', 'biocatalysis', 'biomanufacturing', 'bioprocess engineering', 'mass and energy balance', 'process simulation', 'biorefinery', 'biofuel', 'bioproducts'],
)
