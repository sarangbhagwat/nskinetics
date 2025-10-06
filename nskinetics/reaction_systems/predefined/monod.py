# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

from ..reaction_system import RxnSys
from ...reactions import Rxn
from math import exp
import numpy as np

__all__ = ('MonodCellGrowthReactionSystem', 'Monod',)

#%%
class MonodCellGrowthReactionSystem(RxnSys):
    def __init__(self, 
                 C, S,
                 umax, KS,
                 species_system,
                 stoich=1.0,
                 k_inhibs=None,
                 inhibitors=None,
                 ID=None,
                 freeze_species_index_params=True):
        if inhibitors is None:
            inhibitors = []
        if k_inhibs is None:
            k_inhibs = []
        if ID is None:
            ID = f'MCG_{C}_{S}'
            if not inhibitors==[]:
                ID += '_inhib'
                for i in inhibitors:
                    ID += f'_{i}'
        self.C = C
        self.S = S
        self.inhibitors = inhibitors
        
        S_index = species_system.index(S)
        inhib_indices = [species_system.index(i) for i in inhibitors]
        
        # inhib_params = {}
        # for k, inhib in zip(k_inhibs, inhibitors):
        #     inhib_index = species_system.index(inhib)
        #     inhib_params[f'k_{inhib}'] = k
        #     inhib_params[f'index_{inhib}'] = inhib_index
        # rate_params={'umax':umax, 'KS':KS, 'S_index':S_index,}
        # rate_params.update(inhib_params)
        
        def rate_f(species_concs_vector, rxn_stoichs, rl_exps, reactant_indices, product_indices, 
                   umax, KS, k_inhibs, S_index, inhib_indices, # rate_params
                   ):
            S_conc = species_concs_vector[S_index]
            u = umax * S_conc /(KS + S_conc)
            u *= exp(sum([-k*species_concs_vector[i] for k, i in zip(k_inhibs, inhib_indices)]))
            return u
        
        self.reaction = reaction = Rxn.from_equation(ID=f'r_{ID}', 
                               chem_equation=f'{C} + {S} <-> {1.0 + stoich} {C}', 
                               rate_f=rate_f,
                               rate_params={'umax':umax, 'KS':KS, 
                                            'k_inhibs':np.array(k_inhibs),
                                            'S_index':S_index, 'inhib_indices':inhib_indices},
                               species_system=species_system, 
                               is_multicompartment=False)
        RxnSys.__init__(self, ID=ID, reactions=[reaction], species_system=species_system)
    
        self.freeze_species_index_params = freeze_species_index_params
        if freeze_species_index_params:
            reaction._freeze_params.add('S_index')
            for i in range(len(inhib_indices)):
                reaction._freeze_params.add(f'inhib_indices_{i}')
            
    @property
    def kcat(self):
        return self.reaction.rate_params['kcat']
    @kcat.setter
    def kcat(self, val):
        self.reaction.rate_params['kcat'] = float(val)
        
    @property
    def KM(self):
        return self.reaction.rate_params['KM']
    @KM.setter
    def KM(self, val):
        self.reaction.rate_params['KM'] = float(val)
        
Monod = MonodCellGrowthReactionSystem
