# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import nskinetics as nsk
import numpy as np

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
rxn_sys.plot_solution() 

rxn_sys.plot_solution(sps_to_include=['ES'])

#%% Design of experiments
param_keys = rxn_sys.reaction_kinetic_param_keys
candidate_initials = {'E': np.linspace(1e-4, 0.1, 30), 'S': np.linspace(1e-4, 1., 30), 'ES': [0.0], 'P': [0.0]}
spike_options = None
t_eval = np.linspace(0, 30, 50)
output_idx = [1, 3]  # Measuring S and P

best_expts = nsk.doe.fim.design_experiments(
    rxn_sys,
    param_keys,
    candidate_initials,
    t_eval,
    spike_options=spike_options,
    output_idx=output_idx,
    epsilon=1e-4,
    top_n=3
)

for i, expt in enumerate(best_expts):
    print(f"--- Experiment {i+1} ---")
    print("Initial concentrations:", expt['y0'])
    print("Spikes:", expt['spikes'])
    print("D-optimality:", expt['score'])
    