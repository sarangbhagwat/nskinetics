# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import nskinetics as nsk
import numpy as np
from warnings import filterwarnings

def test_simple_ESP_log_transform(plot=False, 
                    show_progress=False,
                    show_output=False,
                    show_warnings=False):
    
    if not show_warnings: filterwarnings("ignore")
    
    # Create a SpeciesSystem object
    sp_sys = nsk.SpeciesSystem('sp_sys', 
                           ['E', 'S', 'ES', 'P'], # enzyme, substrate, enzyme-substrate complex, product
                           concentrations=[1e-4, 1e-4, 0., 0.])
    
    # Describe reactions by writing chemical equations and kinetic parameter info
    reactions = [
                'E + S <-> ES; kf = 12.0, kb = 10.0', # kf = kon, kb = koff
                'ES -> E + P; kf = 32.0' # kf = kcat (enzyme turnover number)
                ]
    
    # Generate a ReactionSystem from strings
    rxn_sys = nsk.ReactionSystem(ID='rxn_sys', 
                                     reactions=reactions,
                                     species_system=sp_sys)
    
    # Simulate the ReactionSystem
    rxn_sys.solve(t_span=[0, 2*24*3600], # I want to simulate the system over 2 days
                     sp_conc_for_events={'S':1e-6}, # In addition to a full simulation,
                     log_transform_concs=True)                              # I want to know the time at which [S] drops to 1e-6
    
    # Plot results
    if plot:
        rxn_sys.plot_solution()
        
        rxn_sys.plot_solution(sps_to_include=['ES'])
    
    # Tests
    
    # assert np.allclose(rxn_sys._solution['t_events'], 
    #                    np.array([42219.44616989]), 
    #                    rtol=1e-5, atol=1e-8)
    
    # assert np.allclose(rxn_sys._solution['y_events'], 
    #                    np.array([[9.99998909e-05],
    #                             [1.00000000e-06],
    #                             [1.09091871e-10],
    #                             [9.89998909e-05]]), 
    #                    rtol=1e-5, atol=1e-8)
    
    # assert np.allclose(rxn_sys._solution['y'][:,50],
    #                    np.array([9.99958946e-05, 
    #                              3.76342167e-05, 
    #                              4.10542307e-09, 
    #                              6.23616779e-05]),
    #                    rtol=1e-5, atol=1e-8)
    
    return rxn_sys