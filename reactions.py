# -*- coding: utf-8 -*-
"""
Created on Thu May 29 17:39:02 2025

@author: sarangbhagwat
"""

import numpy as np
from numba import njit
from warnings import warn

__all__ = ('Reaction', 'Rxn', 'IrreversibleReaction', 'IrrevRxn',
           'ReversibleReaction', 'RevRxn',
           'ChemicalEquation')

#%% Abstract chemical equation class
class ChemicalEquation():
    def __init__(self, ID, 
                 species_system, 
                 stoichiometry=None, 
                 paired_obj=None, # can pass an object that has a parameter 'stoichiometry'
                 ):
        self.ID = ID
        self.species_system = species_system
        self._stoichiometry = stoichiometry if stoichiometry is not None\
            else paired_obj.stoichiometry
        self.paired_obj = paired_obj
        
    def __str__(self):
        lhs = ''
        rhs = ''
        for chem, stoich in zip(self.species_system.all_sps, self.stoichiometry):
            
            if stoich<0: 
                if not lhs=='': lhs+= ' + '
                lhs+= str(-stoich) + ' ' + chem.ID
            elif stoich>0: 
                if not rhs=='': rhs+= ' + '
                rhs+= str(stoich) + ' ' + chem.ID
        return 'ChemicalEquation(' + lhs + ' -> ' + rhs + ')'
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def stoichiometry(self):
        return self._stoichiometry
    @stoichiometry.setter
    def stoichiometry(self, new_stoichiometry):
        self._stoichiometry = new_stoichiometry
        if self.paired_obj is not None:
            if not self.paired_obj.stoichiometry==new_stoichiometry:
                warn(f'Replaced {self.ID} stoichiometry with {new_stoichiometry}, but this does not match the paired_obj stoichiometry: {self.paired_obj.stoichiometry}.',
                     RuntimeWarning)

#%% Abstract reaction class
class AbstractReaction():
    """ 
    Abstract class to define stoichiometries for a chemical reaction.
    In addition to required arguments, initialization must include 
    exactly ONE of the following sets of optional arguments:
        (a) reactants and products; OR
        (b) chem_equation; OR
        (c) stoichiometry
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
        
        if self.chem_equation is None:
            self.chem_equation = ChemicalEquation(ID=ID+'_eqn', 
                                                  species_system=species_system,
                                                  stoichiometry=self.stoichiometry)
        
#%% Reaction class

@njit(cache=True)
def dconcs_dt(kf, kb, species_concs_vector, rxn_stoichs, rxn_exps,
              reactant_indices, product_indices):
    # -----------------------------------------
    # for a one-way reaction,
    # returns array of temporal rates of change 
    # in the concentrations of given species
    # (- denotes decrease, + denotes increase)
    # -----------------------------------------
    
    # 'change' is the 1-stoichiometry-equivalent change for reactant concs
    
    # change = 0
    # if np.all(rxn_exps==1.):
    # change = kf*np.prod(species_concs_vector[reactant_indices]) -\
    #     kb*np.prod(species_concs_vector[product_indices])
    # else:
    change = kf*np.prod(np.power(species_concs_vector[reactant_indices], rxn_exps[reactant_indices])) -\
        kb*np.prod(np.power(species_concs_vector[product_indices], rxn_exps[product_indices]))
    
    # if change is too great, cap it to
    # the limiting reactant conc
    temp_vector = np.copy(species_concs_vector)
    # breakpoint()
    temp_vector += change*rxn_stoichs
    if np.any(temp_vector<0):
        tv_stoich_adj = temp_vector/np.abs(rxn_stoichs)
        tv_stoich_adj[np.where(np.isinf(tv_stoich_adj))] = 0.
        tv_stoich_adj[np.where(np.isnan(tv_stoich_adj))] = 0.
        limiting_reactant_index = np.where(tv_stoich_adj==np.min(tv_stoich_adj))[0][0]
        change = species_concs_vector[limiting_reactant_index]/np.abs(rxn_stoichs[limiting_reactant_index]) # limiting reactant conc. adjusted by stoichiometry
    return change * rxn_stoichs

class Reaction(AbstractReaction):
    """ 
    Class for a chemical reaction with defined kinetics.
    Can describe a one-way/irreversible reaction 
    (using kf; kb set to zero by default) OR a 
    reversible reaction (using kf and kb).
    
    In addition to the above, initialization must include 
    exactly ONE of the following sets of optional arguments:
        (a) reactants and products; OR
        (b) chem_equation; OR
        (c) stoichiometry
        
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
    exponents: array, optional
        Array (of length equal to species_system.concentrations) of 
        rate law exponent float values.
    kf : float
        forward reaction rate constant.
    kb : float, optional
        backward reaction rate constant.
    species_system: SpeciesSystem object
        System of chemical species including at least those involved in this
        reaction.
    """
    def __init__(self, ID, 
                 species_system,
                 kf,
                 kb=0.,
                 reactants=None, products=None,
                 chem_equation=None,
                 stoichiometry=None,
                 exponents=None,
                 get_exponents_from_stoich=False,
                 ):
        AbstractReaction.__init__(self, ID, 
                     species_system,
                     reactants=reactants, products=products,
                     chem_equation=chem_equation,
                     stoichiometry=stoichiometry,)
        stoich=self.stoichiometry
        self.kf = kf
        self.kb = kb
        if exponents is None:
            if get_exponents_from_stoich:
                self.exponents = np.abs(len(stoich))
            else:
                self.exponents = np.ones(len(stoich))
        else:
            self.exponents = exponents
        self.reactant_indices = np.where(stoich<0)
        self.product_indices = np.where(stoich>0)
        
    def get_dconcs_dt(self):
        return dconcs_dt(kf=self.kf, 
                         kb=self.kb,
                         species_concs_vector=self.species_system.concentrations, 
                         rxn_stoichs=self.stoichiometry,
                         rxn_exps=self.exponents,
                         reactant_indices=self.reactant_indices,
                         product_indices=self.product_indices)
    
    def __str__(self):
        lhs = ''
        rhs = ''
        arrow = '->' if self.kb==0. else '<->'
        for chem, stoich in zip(self.species_system.all_sps, self.stoichiometry):
            
            if stoich<0: 
                if not lhs=='': lhs+= ' + '
                lhs+= str(-stoich) + ' ' + chem.ID
            elif stoich>0: 
                if not rhs=='': rhs+= ' + '
                rhs+= str(stoich) + ' ' + chem.ID
        return 'Reaction(' + lhs + ' ' + arrow + ' ' + rhs + ')'
    
    def __repr__(self):
        return self.__str__()
    
Rxn = IrreversibleReaction = ReversibleReaction = IrrevRxn = RevRxn = Reaction

