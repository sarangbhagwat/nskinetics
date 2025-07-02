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
rxn_sys.solve(t_span=[0, 2*24*3600],
              sp_conc_for_events={'S':1e-6})                              

# Plot results
rxn_sys.plot_solution() 

rxn_sys.plot_solution(sps_to_include=['ES'])

# Fit to results from old set of kinetic parameters

filterwarnings("ignore")
rxn_sys.fit_reaction_kinetic_parameters_to_data(data=rxn_sys._solution_dfs[0],
                                                p0=np.ones(len(rxn_sys.reaction_kinetic_params)),
                                                use_only=['S', 'P'],)
filterwarnings("default")

# Simulate the ReactionSystem

sp_sys.concentrations = [1e-4, 1e-4, 0., 0.]
rxn_sys.solve(t_span=[0, 2*24*3600], # I want to simulate the system over 2 days
              sp_conc_for_events={'S':1e-6}, # In addition to a full simulation, I want to know the time at which [S] drops to 1e-6
              )  

# Plot new results
rxn_sys.plot_solution() 

rxn_sys.plot_solution(sps_to_include=['ES'])
