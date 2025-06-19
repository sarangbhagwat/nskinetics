# -*- coding: utf-8 -*-
"""
Created on Thu May 29 17:39:02 2025

@author: sarangbhagwat
"""

import numpy as np
from numba import njit

__all__ = ('Reaction', 'Rxn', 'IrreversibleReaction', 'IrrevRxn',
           'ReversibleReaction', 'RevRxn',)

#%% Abstract chemical equation class
# class ChemicalEquation():
#     def __init__(self, ID, )
#%% Abstract reaction class
class AbstractReaction():
    """ 
    Abstract class for a chemical reaction.
    Parameters
    ----------
    ID : str, optional
        ID.
    reactants : dict, optional
        Key: reactant name as string; Value: stoichiometry as float or integer.
    products : dict, optional
        Key: product name as string; Value: stoichiometry as float or integer.
    chem_equation: Object, optional
        Must have a parameter 'stoichiometry' containing an array (of length 
        equal to species_system.concentrations) of stoichiometry float values.
    stoichiometry: array, optional
        Array (of length equal to species_system.concentrations) of 
        stoichiometry float values.
    species_system: SpeciesSystem object
        System of chemical species including at least those involved in this
        reaction.
    """
    def __init__(self, ID, 
                 species_system,
                 reactants=None, products=None,
                 chem_equation=None,
                 stoichiometry=None,
                 ):
        self.ID = ID
        self.reactants = reactants
        self.products = products
        self.chem_equation = chem_equation
        self.species_system = species_system
        
        if stoichiometry is None:
            if chem_equation is None:
                stoichiometry = []
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
                else:
                    raise TypeError('\nGiven reactants and products must both be either lists or dicts.\n')
                
                self.stoichiometry = np.array(stoichiometry)
            
            else:
                self.stoichiometry = chem_equation.stoichiometry
        else:
            self.stoichiometry = stoichiometry
            
#%% One-way reaction

@njit(cache=True)
def oneway_dconcs_dt(k, species_concs_vector, rxn_stoichs):
    # -----------------------------------------
    # for a one-way reaction,
    # returns array of temporal rates of change 
    # in the concentrations of given species
    # (- denotes decrease, + denotes increase)
    # -----------------------------------------
    # 1-stoichiometry-equivalent change for reactant concs
    change = k*np.prod(species_concs_vector[np.where(rxn_stoichs<0.)])
    temp_vector = np.copy(species_concs_vector)
    temp_vector += change*rxn_stoichs
    # breakpoint()
    if np.any(temp_vector<0):
        tv_stoich_adj = temp_vector/np.abs(rxn_stoichs)
        tv_stoich_adj[np.where(np.isinf(tv_stoich_adj))] = 0.
        tv_stoich_adj[np.where(np.isnan(tv_stoich_adj))] = 0.
        limiting_reactant_index = np.where(tv_stoich_adj==np.min(tv_stoich_adj))[0][0]
        # print(limiting_reactant_index)
        # try: 
        change = species_concs_vector[limiting_reactant_index]/np.abs(rxn_stoichs[limiting_reactant_index]) # limiting reactant conc. adjusted by stoichiometry
        # except:
        #     breakpoint()
    return change * rxn_stoichs


class Reaction(AbstractReaction):
    """ 
    Class for a one-way reaction.
    Using this class to make reversible reactions
    is not recommended as this will slow computation.
    Parameters
    ----------
    ID : str, optional
        ID.
    reactants : dict, optional
        Key: reactant name as string; Value: stoichiometry as float or integer.
    products : dict, optional
        Key: product name as string; Value: stoichiometry as float or integer.
    chem_equation: Object, optional
        Must have a parameter 'stoichiometry' containing an array (of length 
        equal to species_system.concentrations) of stoichiometry float values.
    stoichiometry: array, optional
        Array (of length equal to species_system.concentrations) of 
        stoichiometry float values.
    k : float
        rate constant.
    species_system: SpeciesSystem object
        System of chemical species including at least those involved in this
        reaction.
    """
    def __init__(self, ID, 
                 k, species_system,
                 reactants=None, products=None,
                 chem_equation=None,
                 stoichiometry=None,
                 ):
        AbstractReaction.__init__(self, ID, 
                     species_system,
                     reactants=reactants, products=products,
                     chem_equation=chem_equation,
                     stoichiometry=stoichiometry,)
        self.k = k
        
    def get_dconcs_dt(self):
        return oneway_dconcs_dt(k=self.k, 
                           species_concs_vector=self.species_system.concentrations, 
                           rxn_stoichs=self.stoichiometry)

Rxn = IrreversibleReaction = IrrevRxn = Reaction

#%% Reversible reaction
@njit(cache=True)
def rev_dconcs_dt(kf, kb, species_concs_vector, rxn_stoichs):
    # -----------------------------------------
    # for a one-way reaction,
    # returns array of temporal rates of change 
    # in the concentrations of given species
    # (- denotes decrease, + denotes increase)
    # -----------------------------------------
    # 1-stoichiometry-equivalent change for reactant concs
    change = kf*np.prod(species_concs_vector[np.where(rxn_stoichs<0.)]) - kb*np.prod(species_concs_vector[np.where(rxn_stoichs>0.)])
    # if change is too great, cap it to
    # the limiting reactant conc
    temp_vector = np.copy(species_concs_vector)
    temp_vector += change*rxn_stoichs
    # breakpoint()
    if np.any(temp_vector<0):
        # breakpoint()
        tv_stoich_adj = temp_vector/np.abs(rxn_stoichs)
        tv_stoich_adj[np.where(np.isinf(tv_stoich_adj))] = 0.
        tv_stoich_adj[np.where(np.isnan(tv_stoich_adj))] = 0.
        limiting_reactant_index = np.where(tv_stoich_adj==np.min(tv_stoich_adj))[0][0]
        # print(limiting_reactant_index)
        # try: 
        change = species_concs_vector[limiting_reactant_index]/np.abs(rxn_stoichs[limiting_reactant_index]) # limiting reactant conc. adjusted by stoichiometry
        # except:
        #     breakpoint()
    return  change * rxn_stoichs


class ReversibleReaction(AbstractReaction):
    """
    Class for a reversible reaction. 
    Using this class for irreversible reactions by setting kb=0
    is not recommended as this will slow computation (use Rxn instead).
    
    Parameters
    ----------
    ID : str, optional
        ID.
    reactants : dict, optional
        Key: reactant name as string; Value: stoichiometry as float or integer.
    products : dict, optional
        Key: product name as string; Value: stoichiometry as float or integer.
    chem_equation: Object, optional
        Must have a parameter 'stoichiometry' containing an array (of length 
        equal to species_system.concentrations) of stoichiometry float values.
    stoichiometry: array, optional
        Array (of length equal to species_system.concentrations) of 
        stoichiometry float values.
    kf : float
        rate constant for forward reaction.
    kb : float
        rate constant for backward reaction.
    species_system: SpeciesSystem object
        System of chemical species including at least those involved in this
        reaction.
    """
    def __init__(self, ID, 
                 kf, kb, species_system,
                 reactants=None, products=None,
                 chem_equation=None,
                 stoichiometry=None,
                 ):
        AbstractReaction.__init__(self, ID, 
                     species_system,
                     reactants=reactants, products=products,
                     chem_equation=chem_equation,
                     stoichiometry=stoichiometry,)
        self.kf = kf
        self.kb = kb
        
    def get_dconcs_dt(self):
        return rev_dconcs_dt(kf=self.kf, 
                             kb=self.kb,
                             species_concs_vector=self.species_system.concentrations, 
                             rxn_stoichs=self.stoichiometry)
        
RevRxn = ReversibleReaction

