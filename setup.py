# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
from setuptools import setup

setup(
    name='nskinetics',
    packages=['nskinetics'],
    license='MIT',
    version='0.1.0',
    description='Simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena',
    long_description=open('README.rst', encoding='utf-8').read(),
    author='Sarang S. Bhagwat',
    install_requires=['IPython>=7.9.0',
                      'numpy>=1.26.4,<2.0.0', 
                      'numba>=0.60.0,<1.0.0',
                      'scipy>=1.13.1,<2.0.0',
                      'matplotlib>=3.5.2,<4.0.0'],
    # extras_require={ 
    #     'dev': [
    #     ]
    # }, 
    # package_data={
    #     'nskinetics': []
    # },
    # exclude_package_data={
    # },
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
