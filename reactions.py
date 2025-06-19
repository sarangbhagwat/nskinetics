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

def read_equation_str(equation_str, all_sps):
    conjugations = ['+', '->', '<->']
    stoichiometry = []
    
    all_sp_IDs = [i.ID for i in all_sps]
    param_info_start = ';' # after this character (if present), there will be information on kinetic parameters
    param_info_junk = [param_info_start, '=', ',', '.',]
        
    split_str, param_info = None, None
    if param_info_start in equation_str:
        split_str = equation_str[:equation_str.index(param_info_start)].split(' ')
        param_info = equation_str[equation_str.index(param_info_start):].replace('=', ' ')
        for i in param_info_junk:
            param_info = param_info.replace(i, ' ')
        param_info = param_info.split(' ')
    else:
        split_str = equation_str.split(' ')
    
    arrow_ind = split_str.index('->') if '->' in split_str else split_str.index('<->')
    
    # ----- #
    # first, handle instances of stoichiometry number immediately next to species name
    # e.g., 2A rather than 2 A; make sure 2A isn't a chemical
    replace={}
    last_str_was_float = False
    
    for str_ in split_str:
        if str_ not in all_sp_IDs + conjugations: # is it a sp_ID or conjugation
            resolved=False
            try:
                if last_str_was_float: 
                    raise ValueError(f'Equation string "{equation_str}" has at least numbers one after the other that with neither being part of a registered chemical name in the species system {all_sp_IDs}.')
                float(str_) # is it a number by itself (stoichiometry)
                last_str_was_float = True
                resolved = True
            except Exception as e:
                use_index_upto=1
                for i in range(1, len(str_)):
                    try:
                        float(str[:i])
                        use_index_upto+=1
                    except:
                        if i>1:
                            stoich_part = str_[:use_index_upto]
                            sp_ID_part = str_[use_index_upto:]
                            replace[str_] = stoich_part, sp_ID_part
                            resolved=True
                            if last_str_was_float:
                                raise e
                            break
                        else:
                            RuntimeError(f'Error reading equation string "{equation_str}"')
                if not use_index_upto==len(str_)-1:
                    last_str_was_float = False
            if not resolved:
                raise ValueError(f'Equation string "{equation_str}" contains element "{str_}" that was not identified as a species, stoichiometry number, conjugation')
    for k, v in replace.items():
        ind = split_str.index(k)
        split_str[ind] = v[0]
        split_str.insert(ind+1, v[1])
    # ----- #
    
    # ----- #
    # Get stoichiometry array
    for sp_ID in all_sp_IDs:
        if sp_ID in split_str:
            try:
                stoichiometry.append(float(split_str[split_str.index(sp_ID)-1]))
            except:
                stoichiometry.append(1.)
            if split_str.index(sp_ID)<arrow_ind:
                stoichiometry[-1] *= -1
        else:
            stoichiometry.append(0.)
    # ----- #
    
    # ----- #
    # Get kinetic parameters, if any
    
    kf_kb = [None, None]
    if param_info is not None:
        param_info_clean = [i for i in param_info if not i in param_info_junk + ['']]
        kf_chars = ['kf',]
        kb_chars = ['kb',]
        for i in range(len(param_info_clean)-1):
            if param_info_clean[i] in kf_chars:
                kf_kb[0] = float(param_info_clean[i+1])
            elif param_info_clean[i] in kb_chars:
                kf_kb[1] = float(param_info_clean[i+1])
    
    return np.array(stoichiometry), kf_kb[0], kf_kb[1]

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
    
    def from_string(ID, equation_str, species_system, paired_obj=None):
        stoichiometry, kf, kb = read_equation_str(equation_str, 
                                                   species_system.all_sps)
        return ChemicalEquation(ID=ID, 
                                species_system=species_system,
                                stoichiometry=stoichiometry,
                                paired_obj=paired_obj),\
               kf, kb
        
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
        if kf is None:
            kf = 0.
        if kb is None:
            kb = 0.
            
        self.kf = kf
        self.kb = kb
        self.get_exponents_from_stoich = get_exponents_from_stoich
        
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
                if not np.abs(stoich)==1.: lhs+= str(-stoich) + ' '
                lhs+= chem.ID
            elif stoich>0: 
                if not rhs=='': rhs+= ' + '
                if not np.abs(stoich)==1.: rhs+= str(stoich) + ' '
                rhs+= chem.ID
        return 'Reaction(' + lhs + ' ' + arrow + ' ' + rhs + ')'
    
    def __repr__(self):
        return self.__str__()
    
    def from_equation(ID, chem_equation, species_system, 
                      kf=None, # overrides any parameter info in the chem_equation string
                      kb=None, # overrides any parameter info in the chem_equation string
                      exponents=None, get_exponents_from_stoich=None):
        kf_, kb_ = None, None
        if isinstance(chem_equation, str):
            chem_equation, kf_, kb_ = ChemicalEquation.from_string(ID=ID+'_eqn', 
                                                         equation_str=chem_equation,
                                                         species_system=species_system)
        if kf is None:
            kf = kf_
        if kb is None:
            kb = kb_
        
        return Reaction(ID=ID,
                        species_system=species_system,
                        chem_equation=chem_equation,
                        kf=kf,
                        kb=kb,
                        exponents=exponents,
                        get_exponents_from_stoich=get_exponents_from_stoich,)
        
Rxn = IrreversibleReaction = ReversibleReaction = IrrevRxn = RevRxn = Reaction

