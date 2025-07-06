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

def test_simple_ESP(plot=False, 
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
                     )                              # I want to know the time at which [S] drops to 1e-6
    
    # Plot results
    if plot:
        rxn_sys.plot_solution() 
        
        rxn_sys.plot_solution(sps_to_include=['ES'])
    
    # Tests
    print(rxn_sys._solution['t_events'], rxn_sys._solution['y_events'])
    assert np.allclose(rxn_sys._solution['t_events'], 
                       np.array([42218.33254329]), 
                       rtol=1e-3, atol=1e-3)
    
    assert np.allclose(rxn_sys._solution['y_events'], 
                       np.array([[9.99998660e-05],
                                 [1.00000000e-06],
                                 [1.09091845e-10],
                                 [9.86354523e-05]]), 
                       rtol=1e-3, atol=1e-8)

    return rxn_sys


def test_simple_ESP_inhib(plot=False, 
                          show_progress=False,
                          show_output=False,
                          show_warnings=False):
    
    if not show_warnings: filterwarnings("ignore")
    
    # Create a SpeciesSystem object
    sp_sys = nsk.SpeciesSystem('sp_sys', 
                           ['E', 'S', 'ES', 'P',
                            'I_CI', 'EI_CI', 'Q',
                            'I_MBI', 'EI_MBI_unstable', 'EI_MBI_stable'], # mechanism-based_inhibitor, unstable enzyme-MBI complex, stable enzyme-MBI complex 
                           concentrations=[1e-4, 1e-4, 0, 0,
                                           5e-5, 0, 0,
                                           3e-5, 0, 0])
    
    # Describe reactions by writing chemical equations and kinetic parameter info
    reactions = [
                'E + S <-> ES; kf = 12, kb = 10.0',
                'ES -> E + P; kf = 32.0',
                'E + I_CI <-> EI_CI; kf=12, kb=10.0',
                'EI_CI -> E + Q; kf=32',
                'E + I_MBI <-> EI_MBI_unstable; kf=12.0, kb=10',
                'EI_MBI_unstable -> EI_MBI_stable; kf = 32'
                ]
    
    # Generate a ReactionSystem from strings
    rxn_sys = nsk.ReactionSystem(ID='rxn_sys', 
                                     reactions=reactions,
                                     species_system=sp_sys)
    
    # Simulate the ReactionSystem
    rxn_sys.solve(t_span=[0, 2*24*3600],
                  sp_conc_for_events={'S':1e-6})
    
    # Plot results
    if plot:
        rxn_sys.plot_solution()
    
    # Tests
    
    assert np.allclose(rxn_sys._solution['t_events'], 
                       np.array([55693.31551888]), 
                       rtol=1e-3, atol=1e-3)
    
    assert np.allclose(rxn_sys._solution['y_events'], 
                       np.array([[7.03279223e-05],
                                 [1.00000000e-06],
                                 [7.67219072e-11],
                                 [9.86283203e-05],
                                 [5.00000000e-07],
                                 [3.83609536e-11],
                                 [4.93141601e-05],
                                 [3.00000000e-07],
                                 [2.30165722e-11],
                                 [2.95884961e-05]]), 
                       rtol=1e-3, atol=1e-8)

    return rxn_sys


def test_simple_ESP_fed_batch(plot=False, 
                              show_progress=False,
                              show_output=False,
                              show_warnings=False):
    
    if not show_warnings: filterwarnings("ignore")
    
    # Create a SpeciesSystem object
    sp_sys = nsk.SpeciesSystem('sp_sys', 
                           ['E', 'S', 'ES', 'P',],
                           concentrations=[1e-4, 1e-4, 0, 0,])
    
    # Describe reactions by writing chemical equations and kinetic parameter info
    reactions = [
                'E + S <-> ES; kf = 12, kb = 10.0',
                'ES -> E + P; kf = 32.0',
                ]
    
    # Generate a ReactionSystem from strings
    rxn_sys = nsk.ReactionSystem(ID='rxn_sys', 
                                     reactions=reactions,
                                     species_system=sp_sys)
    
    
    # Describe forced concentration spikes for any species 
    # (e.g., from feeding substrate in a fed-batch regime)
    spikes = {20000: 'Target; S; 1e-4', # at t=40000, add enough S to achieve [S]=1e-4
              50000: 'Target; S; 1e-4', # at t=50000, add enough S to to achieve [S]=1e-4
              80000: 'Target; S; 1e-4', # at t=80000, add enough S to achieve [S]=1e-4
              100000: 'Change; S; 2e-4',# at t=100000, add enough S to increase [S] by 2e-4
              }
    
    # Simulate the ReactionSystem
    rxn_sys.solve(t_span=[0, 2*24*3600],
                  sp_conc_for_events={'S':1e-6},
                  spikes=spikes)
    
    # Plot results
    if plot:
        rxn_sys.plot_solution()
    
    # Tests
    
    assert np.allclose(rxn_sys._solution['t_events'], 
                       np.array([149077.3233837]), 
                       rtol=1e-3, atol=1e-3)
    
    assert np.allclose(rxn_sys._solution['y_events'], 
                       np.array([[9.99996948e-05],
                                 [1.00000000e-06],
                                 [1.09091658e-10],
                                 [5.78368281e-04]]), 
                       rtol=1e-3, atol=1e-8)

    return rxn_sys


def test_simple_ESP_inhib_fed_batch(plot=False, 
                                    show_progress=False,
                                    show_output=False,
                                    show_warnings=False):
    
    if not show_warnings: filterwarnings("ignore")
    
    # Create a SpeciesSystem object
    sp_sys = nsk.SpeciesSystem('sp_sys', 
                           ['E', 'S', 'ES', 'P',
                            'I_CI', 'EI_CI', 'Q',
                            'I_MBI', 'EI_MBI_unstable', 'EI_MBI_stable'], # mechanism-based_inhibitor, unstable enzyme-MBI complex, stable enzyme-MBI complex 
                           concentrations=[1e-4, 1e-4, 0, 0,
                                           5e-5, 0, 0,
                                           3e-5, 0, 0])
    
    # Describe reactions by writing chemical equations and kinetic parameter info
    reactions = [
                'E + S <-> ES; kf = 12, kb = 10.0',
                'ES -> E + P; kf = 32.0',
                'E + I_CI <-> EI_CI; kf=12, kb=10.0',
                'EI_CI -> E + Q; kf=32',
                'E + I_MBI <-> EI_MBI_unstable; kf=12.0, kb=10',
                'EI_MBI_unstable -> EI_MBI_stable; kf = 32'
                ]
    
    
    # Generate a ReactionSystem from strings
    rxn_sys = nsk.ReactionSystem(ID='rxn_sys', 
                                     reactions=reactions,
                                     species_system=sp_sys)
    
    # Describe forced concentration spikes for any species 
    # (e.g., from feeding substrate in a fed-batch regime)
    spikes = {20000: 'Target; S; 1e-4', # at t=40000, add enough S to achieve [S]=1e-4
              50000: 'Target; S; 1e-4', # at t=50000, add enough S to increase [S] by 1e-4
              80000: 'Target; S; 1e-4',
              100000: 'Change; S; 2e-4',
              }
    
    # Simulate the ReactionSystem
    rxn_sys.solve(t_span=[0, 2*24*3600],
                  sp_conc_for_events={'S':1e-6},
                  spikes=spikes)
    
    # Plot results
    if plot:
        rxn_sys.plot_solution()
    
    # Tests
    print(rxn_sys._solution['t_events'], rxn_sys._solution['y_events'])
    
    assert np.allclose(rxn_sys._solution['t_events'], 
                       np.array([170730.89673069]), 
                       rtol=1e-3, atol=1e-3)
    
    assert np.allclose(rxn_sys._solution['y_events'], 
                       np.array([[7.00060582e-05],
                                 [1.00000000e-06],
                                 [7.63707755e-11],
                                 [5.59863488e-04],
                                 [7.62382049e-11],
                                 [5.82237083e-15],
                                 [4.98832953e-05],
                                 [4.57429229e-11],
                                 [3.49342250e-15],
                                 [2.99299772e-05]]), 
                       rtol=1e-3, atol=1e-8)
    
    return rxn_sys
