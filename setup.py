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
    version='0.5.0',
    description='Simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena',
    long_description=open('README.rst', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    author='Sarang S. Bhagwat',
    install_requires=['IPython>=7.9.0',
                      'numpy>=1.26.4,<2.0.0', 
                      'numba>=0.60.0,<1.0.0',
                      'scipy>=1.13.1,<2.0.0',
                      'matplotlib>=3.5.2,<4.0.0',
                      'xlsxwriter>=3.2.5,<4.0.0',
                      'scikit-learn',
                      'pandas>=2.2.2,<3.0.0',
                      'tellurium',
                      'biosteam==2.47.0',
                      # 'python-libsbml>=5.20.5,<6.0.0',
		     ],
    # extras_require={
    #     'dev': [
    #     ]
    # },
    include_package_data=True,
    package_data={
        # Shipped model files loaded by name at import/run time (e.g. the
        # isobutanol model module loads its antimony .txt on import, and the
        # tutorial loads the shipped SBML .xml). These must be in the wheel so
        # installed (non-editable) copies can find them.
        'nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo': ['*.txt', '*.xml'],
        # Reference antimony used by the events tests.
        'nskinetics.tests': ['data/*.txt'],
    },
    python_requires='>=3.9',
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
                 'Programming Language :: Python :: 3.9',
                 'Programming Language :: Python :: 3.10',
                 'Programming Language :: Python :: Implementation :: CPython',
                 'Topic :: Education'],
    keywords=['reaction kinetics', 'biocatalysis', 'biomanufacturing', 'bioprocess engineering', 'mass and energy balance', 'process simulation', 'biorefinery', 'biofuel', 'bioproducts'],
)
