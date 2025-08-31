# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import os
import nskinetics

nskinetics_filepath = nskinetics.__file__

# def test_RxnSys_from_SBML(filename='michaelis_menten_valid.xml',
#                           plot=True, 
#                           show_progress=True,
#                           show_output=True,
#                           show_warnings=True):
#     delimiter = None
#     for i in ['/', '\\']:
#         if i in nskinetics_filepath:
#             delimiter=i
    
#     SBML_file_location = nskinetics_filepath.replace('__init__.py', '')
#     SBML_file_location += 'tests'
#     # os.chdir(nskinetics_filepath+SBML_file_location)
#     rxn_sys = nskinetics.ReactionSystem.from_SBML(filepath=SBML_file_location+delimiter+filename)
    
#     return rxn_sys
