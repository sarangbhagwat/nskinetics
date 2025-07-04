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
    
    assert np.allclose(rxn_sys._solution['t_events'], 
                       np.array([42219.44616989]), 
                       rtol=1e-5, atol=1e-8)
    
    assert np.allclose(rxn_sys._solution['y_events'], 
                       np.array([[9.99998909e-05],
                                [1.00000000e-06],
                                [1.09091871e-10],
                                [9.89998909e-05]]), 
                       rtol=1e-5, atol=1e-8)
    
    assert np.allclose(rxn_sys._solution['y'][:,50],
                       np.array([9.99958946e-05, 
                                 3.76342167e-05, 
                                 4.10542307e-09, 
                                 6.23616779e-05]),
                       rtol=1e-5, atol=1e-8)


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
                       np.array([55698.02533593]), 
                       rtol=1e-5, atol=1e-8)
    
    assert np.allclose(rxn_sys._solution['y_events'], 
                       np.array([[7.02998850e-05],
                              [1.00000000e-06],
                              [7.66929745e-11],
                              [9.89999233e-05],
                              [5.00000000e-07],
                              [3.83464834e-11],
                              [4.94999617e-05],
                              [3.00000000e-07],
                              [2.30078923e-11],
                              [2.96999770e-05]]), 
                       rtol=1e-5, atol=1e-8)
    
    assert np.allclose(rxn_sys._solution['y'][:,50],
                       np.array([8.56115744e-05, 5.20628936e-05, 4.86241697e-09, 4.79322440e-05,
                              2.60314468e-05, 2.43120849e-09, 2.39661220e-05, 1.56188681e-05,
                              1.45872509e-09, 1.43796732e-05]),
                       rtol=1e-5, atol=1e-8)


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
                       np.array([149078.391953]), 
                       rtol=1e-5, atol=1e-8)
    
    assert np.allclose(rxn_sys._solution['y_events'], 
                       np.array([[9.99998909e-05],
                              [1.00000000e-06],
                              [1.09094016e-10],
                              [5.80130395e-04]]), 
                       rtol=1e-5, atol=1e-8)
    
    assert np.allclose(rxn_sys._solution['y'][:,150],
                       np.array([9.99890949e-05, 
                                 9.99809209e-05, 
                                 1.09051014e-08, 
                                 1.84930371e-04]),
                       rtol=1e-5, atol=1e-8)


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
    
    assert np.allclose(rxn_sys._solution['t_events'], 
                       np.array([170738.8596293]), 
                       rtol=1e-5, atol=1e-8)
    
    assert np.allclose(rxn_sys._solution['y_events'], 
                       np.array([[6.99999694e-05],
                              [1.00000000e-06],
                              [7.63641321e-11],
                              [5.63156018e-04],
                              [7.62706609e-11],
                              [5.82434282e-15],
                              [4.99999237e-05],
                              [4.57623965e-11],
                              [3.49460569e-15],
                              [2.99999542e-05]]), 
                       rtol=1e-5, atol=1e-8)
    
    assert np.allclose(rxn_sys._solution['y'][:,150],
                       np.array([7.04568689e-05, 9.99897217e-05, 7.65819642e-09, 1.74218180e-04,
                              7.74311035e-07, 5.95157329e-11, 4.92256294e-05, 4.64586621e-07,
                              3.57094397e-11, 2.95353777e-05]),
                       rtol=1e-5, atol=1e-8)
