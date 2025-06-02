# -*- coding: utf-8 -*-
"""
Created on Thu May 29 17:32:54 2025

@author: sarangbhagwat
"""
from .nonsteady import SpeciesSystem
import numpy as np
from numba import njit
from matplotlib import pyplot as plt

__all__ = ('MichaelisMenten',)

#%% Michaelis Menten (STEADY-state only)
# (Do NOT use in non-steady state systems)

#!!!
@njit(cache=True)
def steady_MM_dconcs_dt(kcat, KM, species_concs_vector, 
                        e_index, s_index,
                        rxn_stoichs):
    # for a Michaelis-Menten RxnSys (STEADY-state),
    # returns array of temporal rates of change 
    # in the concentrations of given species
    # (- denotes decrease, + denotes increase)
    s_conc = species_concs_vector[s_index]
    change = kcat*species_concs_vector[e_index]*species_concs_vector[s_index]/(s_conc + KM)
    temp_vector = np.copy(species_concs_vector)
    temp_vector += change*rxn_stoichs
    # breakpoint()
    if np.any(temp_vector<0):
        tv_stoich_adj = temp_vector/np.abs(rxn_stoichs)
        tv_stoich_adj[np.where(np.isinf(tv_stoich_adj))] = 0.
        tv_stoich_adj[np.where(np.isnan(tv_stoich_adj))] = 0.
        limiting_reactant_index = np.where(tv_stoich_adj==np.min(tv_stoich_adj))[0][0]
        # print(limiting_reactant_index)
        change = species_concs_vector[limiting_reactant_index]/np.abs(rxn_stoichs[limiting_reactant_index]) # limiting reactant conc. adjusted by stoichiometry
    return change * rxn_stoichs

class MichaelisMenten():
    """
    Abstract class for a STEADY-STATE system of 
    reactions involving a single enzyme, 
    substrate, and product.
    
    Parameters
    ----------
    ID : str
        ID.
    enzyme : str
        Enzyme.
    substrate : str
        Substrate.
    product : str
        Product.
    kcat : float
        Rate constant for the dissociation of
        enzyme-substrate complex to form enzyme
        and product.
    KM: float.
        Michaelis constant. Equals (koff + kcat)/kon.
    species_system : SpeciesSystem
        A SpeciesSystem object containing all species
        involved in this system of reactions.
    """
    def __init__(self, 
                 ID, 
                 enzyme, 
                 substrate, 
                 product,
                 kcat, KM,
                 species_system):
        self.ID = ID
        self.enzyme = enzyme
        self.substrate = substrate
        self.product = product
        self.kcat = kcat
        self.KM = KM
        all_sps_IDs = [i.ID for i in species_system.all_sps]
        self._enzyme_index = all_sps_IDs.index(enzyme)
        self._substrate_index = all_sps_IDs.index(substrate)
        # self.product_index = all_sps.index(product)
        self.species_system = species_system
        
        stoichiometry = []
        
        reactants = {}
        if isinstance(substrate, dict):
            reactants.update(substrate)
        else:
            reactants = [substrate]
        
        products = {}
        if isinstance(product, dict):
            products.update(product)
        else:
            products = [product]
 
        if isinstance(reactants, dict) and isinstance(products, dict):
            reactant_keys = reactants.keys()
            product_keys = products.keys()
            for i in species_system.all_sps:
                if i.ID in reactant_keys: 
                    stoichiometry.append(reactants[i])
                elif i.ID in product_keys: 
                    stoichiometry.append(products[i])
                else:
                    stoichiometry.append(0.)
        elif isinstance(reactants, list) and isinstance(products, list):
            for i in species_system.all_sps:
                if i.ID in reactants: 
                    stoichiometry.append(-1)
                elif i.ID in products: 
                    stoichiometry.append(1)
                else:
                    stoichiometry.append(0.)
        # else:
        #     raise TypeError('\nGiven reactants and products must both be either lists or dicts.\n')
        
        self.stoichiometry = np.array(stoichiometry)
        
    def get_dconcs_dt(self):
        return steady_MM_dconcs_dt(kcat=self.kcat, 
                             KM=self.KM,
                             e_index=self._enzyme_index, 
                             s_index=self._substrate_index,
                             species_concs_vector=self.species_system.concentrations, 
                             rxn_stoichs=self.stoichiometry)
    
SteadyStateCompetitiveInhibition = MichaelisMenten
