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

np.random.seed(3221)

#%% Define batches to run

batch_concs = [
               np.array([1e-4, 1e-4, 0., 0.]),
               np.array([1e-6, 1e-2, 0., 0.]),
               np.array([1e-5, 5e-2, 0., 0.]),
               # np.array([1e-3, 5e-2, 0., 0.]),
               ]
batch_t_spans = [
                 [0, 2*24*3600],
                 [0, 2*24*3600],
                 [0, 2*24*3600],
                 # [0, 10000],
                 ]

n_batches = len(batch_concs)
batch_indices = range(n_batches)

#%% Create a SpeciesSystem object with defined kinetic parameters
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

#%% Simulate all batches and save results
for i, concentrations, t_span in zip(batch_indices, 
                                     batch_concs, 
                                     batch_t_spans):
    sp_sys.concentrations = concentrations
    # Simulate the ReactionSystem
    rxn_sys.solve(t_span=t_span,
                  sp_conc_for_events={'S':1e-6},
                  filename=f'solution{i}',
                  save_events=False)                              
    
    # Plot results
    rxn_sys.plot_solution()
    rxn_sys.plot_solution(sps_to_include=['ES'])

#%% Fit kinetic parameters to saved results from simulating with original set of kinetic parameters

filterwarnings("ignore")
rxn_sys.fit_reaction_kinetic_parameters_to_data(data=[f'solution{i}.xlsx'
                                                      for i in batch_indices[:3]],
                                                p0=1*np.ones(len(rxn_sys.reaction_kinetic_params)),
                                                use_only=[
                                                          # 'E', 
                                                          'S', 
                                                          'P'
                                                          ],
                                                options={'disp':True},
                                                # plot_during_fit=True,
                                                show_progress=True,
                                                )
filterwarnings("default")

#%% Simulate batches with inverse-modeled rxn_sys
for i, concentrations, t_span in zip(range(len(batch_concs)), 
                                     batch_concs, 
                                     batch_t_spans):
    sp_sys.concentrations = concentrations
    # Simulate the ReactionSystem
    rxn_sys.solve(t_span=t_span,
                  sp_conc_for_events={'S':1e-6},
                  # filename=f'solution{i}',
                  save_events=False)                              
    
    # Plot results
    rxn_sys.plot_solution()
    rxn_sys.plot_solution(sps_to_include=['ES'])