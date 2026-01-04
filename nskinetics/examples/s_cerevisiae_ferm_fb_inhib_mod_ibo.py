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

# te_r_ethanol = nsk.examples.saccharomyces_cerevisiae_fermentation.te_r

__all__ = ('te_r', 'reset_kinetic_reaction_system',)

#%%

nsk_filepath = nsk.__file__.replace('\\__init__.py', '')
nsk_examples_filepath = nsk_filepath + '\\examples\\'

#%%
antimony_filename = 's_cerevisiae_ferm_fb_inhib_mod_ibo_antimony.txt'
r = te.loadAntimonyModel(f'{nsk_examples_filepath}\\{antimony_filename}')

#%% Create NSKinetics Reaction System
te_r = nsk.TelluriumReactionSystem(r)
te_r._units['time'] = 'h'
te_r._units['conc'] = 'g/L'
    
def reset_kinetic_reaction_system(r):
    r.reset()
    r_te = r._te
    r_te.n_glu_spikes = 0
    r_te.last_vol_glu_feed_added = 0.
    r_te.tot_vol_glu_feed_added = 0.
    r_te.env = 1.

#%% Document references for IBO pathway kinetic parameters

# r.k_6 = 2.9 * 88.06 * 1e-6 * 60 # 2.9 uM/min # Zhao 2021 https://doi.org/10.3390/foods10051013
# r.K_6 = 31.8 * 88.06 * 1e-3 # 31.8 mM # Zhao 2021 https://doi.org/10.3390/foods10051013
#%%
simulate = False
if simulate:
    reset_kinetic_reaction_system(te_r)
    r.s_glu = 100 # initial glucose conc
    r.x = 1 # initial biomass conc
    print(r.s_glu, r.x, r.X_a, r.s_EtOH, r.s_acetate, r.s_acetald, r.s_AL, r.s_DHI, r.s_KIV, r.s_IBO)
    
    r.simulate(0, 300, 2000,
               # ['time', 's_glu', '[s_glu]'],
               # ['time', 'X_a', 'X_AcDH', 
               #  'a',
               #  ],
               )
    
    print(r.s_glu, r.x, r.X_a, r.s_EtOH, r.s_acetate, r.s_acetald, r.s_AL, r.s_DHI, r.s_KIV, r.s_IBO)
    r.plot()
