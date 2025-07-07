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

def test_gui_simple_ESP(plot=False, 
                    show_progress=False,
                    show_output=False,
                    show_warnings=True,
                    close_immediately=True):
    
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
    
    # Plot results
    if plot:
        rxn_sys.plot_solution() 
        
        rxn_sys.plot_solution(sps_to_include=['ES'])
    
    # GUI
    rxn_sys.GUI(close_immediately=True)
    

def test_gui_simple_ESP_inhib(plot=False, 
                          show_progress=False,
                          show_output=False,
                          show_warnings=False,
                          close_immediately=True):
    
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
    
    
    # GUI
    rxn_sys.GUI(close_immediately=True)


def test_gui_simple_ESP_fed_batch(plot=False, 
                              show_progress=False,
                              show_output=False,
                              show_warnings=False,
                              close_immediately=True):
    
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
    
    
    # GUI
    rxn_sys.GUI(close_immediately=True)


def test_gui_simple_ESP_inhib_fed_batch(plot=False, 
                                    show_progress=False,
                                    show_output=False,
                                    show_warnings=False,
                                    close_immediately=True):
    
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
    
    
    # GUI
    rxn_sys.GUI(close_immediately=True)
