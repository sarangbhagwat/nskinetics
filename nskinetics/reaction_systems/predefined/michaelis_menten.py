# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

from ..reaction_system import RxnSys
from ...reactions import Rxn

__all__ = ('MichaelisMentenReactionSystem', 'MichaelisMenten',
           'MichaelisMentenCellularReactionSystem', 'MichaelisMentenCellular')

# def get_ESP_reactions(E, S, ES, P, kon, koff, kcat, species_system):
#     return [Rxn.from_equation(ID=f'r_{E}_{S}_binding', 
#                            chem_equation=f'{E} + {S} <-> {ES}; kf={kon}, kb={koff}', 
#                            species_system=species_system, 
#                            is_multicompartment=False),
#             Rxn.from_equation(ID=f'r_{ES}_{P}_catalysis', 
#                                    chem_equation=f'{ES} -> {E} + {P}; kf={kcat}, kb=0.0', 
#                                    species_system=species_system, 
#                                    is_multicompartment=False),
#             ]

#%%
class MichaelisMentenReactionSystem(RxnSys):
    def __init__(self, E, S, P, 
                 kcat, KM,
                 species_system,
                 stoich=1.0,
                 ID=None):
        if ID is None:
            ID = f'MM_{E}_{S}_{P}'
        self.E = E
        self.S = S
        self.P = P
        E_index = species_system.index(E)
        S_index = species_system.index(S)
        def rate_f(species_concs_vector, rxn_stoichs, rl_exps, reactant_indices, product_indices, 
                   kcat, KM, E_index, S_index, # rate_params
                   ):
            S_conc = species_concs_vector[S_index]
            return kcat * species_concs_vector[E_index] * S_conc /(KM + S_conc)
        
        self.reaction = reaction = Rxn.from_equation(ID=f'r_{ID}', 
                               chem_equation=f'{E} + {S} <-> {E} + {stoich} {P}', 
                               rate_f=rate_f,
                               rate_params={'kcat':kcat, 'KM':KM, 'E_index':E_index, 'S_index':S_index},
                               species_system=species_system, 
                               is_multicompartment=False)
        RxnSys.__init__(self, ID=ID, reactions=[reaction], species_system=species_system)
    
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
        
MichaelisMenten = MichaelisMentenReactionSystem

#%%
class MichaelisMentenCellularReactionSystem(RxnSys):
    def __init__(self, C, C_E_frac, S, P, 
                 kcat, KM,
                 species_system,
                 stoich=1.0,
                 ID=None):
        if ID is None:
            ID = f'MMC_{C}_{S}_{P}'
        self.C = C
        self.S = S
        self.P = P
        C_index = species_system.index(C)
        S_index = species_system.index(S)
        def rate_f(species_concs_vector, rxn_stoichs, rl_exps, reactant_indices, product_indices, 
                   kcat, KM, C_index, C_E_frac, S_index, # rate_params
                   ):
            S_conc = species_concs_vector[S_index]
            return kcat * species_concs_vector[C_index] * C_E_frac * S_conc /(KM + S_conc)
        
        self.reaction = reaction = Rxn.from_equation(ID=f'r_{ID}', 
                               chem_equation=f'{C} + {S} <-> {C} + {stoich} {P}', 
                               rate_f=rate_f,
                               rate_params={'kcat':kcat, 'KM':KM, 
                                            'C_index':C_index, 'C_E_frac':C_E_frac,
                                            'S_index':S_index},
                               species_system=species_system, 
                               is_multicompartment=False)
        
        RxnSys.__init__(self, ID=ID, reactions=[reaction], species_system=species_system)
    
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
        
    @property
    def C_E_frac(self):
        return self.reaction.rate_params['C_E_frac']
    @C_E_frac.setter
    def C_E_frac(self, val):
        self.reaction.rate_params['C_E_frac'] = float(val)
        
MichaelisMentenCellular = MichaelisMentenCellularReactionSystem
