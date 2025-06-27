# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import nskinetics as nsk

# Create a SpeciesSystem object
sp_sys = nsk.SpeciesSystem('sp_sys', 
                       ['E', 'S', 'ES', 'P',
                        'I_CI', 'EI_CI', 'Q'], # competitive_inhibitor, enzyme-competitive_inhibitor complex, byproduct
                       concentrations=[1e-4, 1e-4, 0, 0,
                                       5e-5, 0, 0])

# Describe reactions by writing chemical equations and kinetic parameter info
reactions = [
            'E + S <-> ES; kf = 12, kb = 10.0',
            'ES -> E + P; kf = 32.0',
            'E + I_CI <-> EI_CI; kf=12, kb=10.0',
            'EI_CI -> E + Q; kf=32'
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
