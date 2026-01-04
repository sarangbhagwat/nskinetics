# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import nskinetics as nsk
import biosteam as bst
import tellurium as te
import simplesbml

__all__ = ('te_r', 'reset_kinetic_reaction_system',)

#%% Import SBML
# filename = 'BIOMD0000000245_url.xml' # Lei et al., 2001
# filepath = f'{filename}'
# model = simplesbml.loadSBMLFile(filepath)

sbml_url = 'https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1003250000/3/BIOMD0000000245_url.xml' # Lei et al., 2001
r = te.loadSBMLModel(sbml_url)
model = simplesbml.loadSBMLStr(r.getSBML())

#%% Update model

# no updates here

#%% Load updated model SBML in Tellurium
r = te.loadSBMLModel(model.toSBML())

#%% Update parameter values through Tellurium model
r.D = 0 # dilution rate # note r.S_f (inlet glucose conc) does not matter when D=0 

#%% Create NSKinetics Reaction System
te_r = nsk.TelluriumReactionSystem(r)
te_r._units['time'] = 'h'
te_r._units['conc'] = 'g/L'
    
def reset_kinetic_reaction_system(r):
    r.reset()
    r._te.n_glu_spikes = 0

#%% Set initial concs
    
r.s_glu = 100 # initial glucose conc
r.x = 2 # initial biomass conc

#%% Export to antimony

r.exportToAntimony('s_cerevisiae_ferm_antimony.txt', current=True)

#%%  simulations

simulate = False
if simulate:
    r.reset()
    r.s_glu = 100 # initial glucose conc
    r.x = 2 # initial biomass conc
    print(r.s_glu, r.x, r.X_a, r.s_EtOH, r.s_acetate, r.s_acetald)
    
    r.simulate(0, 300, 300, 
               # ['time', 'X_a', 'X_AcDH', 
                # 'a',
                # ],
               )
    
    print(r.s_glu, r.x, r.X_a, r.s_EtOH, r.s_acetate, r.s_acetald)
    r.plot()
    
    bst.settings.set_thermo(['Water', 'Carbon', 'Glucose', 'Ethanol', 'Sucrose'])
    
    feed = bst.Stream('feed', Water=100, Carbon=1, Glucose=1)
    
    map_chemicals_nsk_to_bst = {'s_glu': 'Glucose',
                                'x': 'Carbon',
                                's_EtOH': 'Ethanol',}

        
    R302 = nsk.units.NSKFermentation('R302', 
                                     ins=feed, 
                                     kinetic_reaction_system=te_r,
                                     map_chemicals_nsk_to_bst=map_chemicals_nsk_to_bst,
                                     n_simulation_steps=None,
                                     f_reset_kinetic_reaction_system=reset_kinetic_reaction_system,
                                     tau=3*24)
    
    R302.simulate()
    r.plot()

