# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 12:54:07 2025

@author: saran
"""

import nskinetics as nsk

# Create a SpeciesSystem object
sp_sys = nsk.SpeciesSystem('sp_sys', 
                       ['E', 'S', 'ES', 'P'], # enzyme, substrate, enzyme-substrate complex, product
                       concentrations=[1e-4, 1e-4, 0., 0.])

# Describe reactions by writing chemical equations and 
# kinetic parameter info
reactions = [
            'E + S <-> ES; kf = 12.0, kb = 10.0', # kf = kon, kb = koff
            'ES -> E + P; kf = 32.0' # kf = kcat (enzyme turnover number)
            ]

# Generate a ReactionSystem from strings
ESP_rxn_sys = nsk.ReactionSystem(ID='ESP_rxn_sys', 
                                 reactions=reactions,
                                 species_system=sp_sys)

# Simulate the ReactionSystem
ESP_rxn_sys.solve(t_span=[0, 2*24*3600], # I want to simulate the system over 2 days
                 sp_conc_for_events={'S':1e-6}, # In addition to a full simulation,
                 )                              # I want to know the time at which [S] drops to 1e-6

# Plot results
ESP_rxn_sys.plot_solution() 

