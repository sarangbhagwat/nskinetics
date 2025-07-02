# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np
from numba import njit
from warnings import warn
from .utils import read_equation_str
from ..species.enzyme import enzyme_complex_joiner

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
        return f'{self.ID}: ChemicalEquation(' + lhs + ' -> ' + rhs + ')'
    
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
                                                   species_system)
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
def dconcs_dt(kf, kb, species_concs_vector, rxn_stoichs, rl_exps,
              reactant_indices, product_indices):
    # Compute forward rate
    forward = 1.0
    for i in reactant_indices:
        forward *= species_concs_vector[i] ** rl_exps[i]
    forward *= kf

    # Compute backward rate
    backward = 1.0
    for i in product_indices:
        backward *= species_concs_vector[i] ** rl_exps[i]
    backward *= kb

    change = forward - backward

    # Calculate temp vector
    temp_vector = np.empty_like(species_concs_vector)
    for i in range(species_concs_vector.size):
        temp_vector[i] = species_concs_vector[i] + change * rxn_stoichs[i]
            
    # Check for negative concentrations
    limiting_found = False
    min_ratio = 1e20
    limiting_index = -1
    for i in range(temp_vector.size):
        if temp_vector[i] < 0.0:
            stoich = rxn_stoichs[i]
            if stoich != 0.0:
                ratio = species_concs_vector[i] / abs(stoich)
                if ratio < min_ratio:
                    min_ratio = ratio
                    limiting_index = i
                    limiting_found = True

    if limiting_found:
        change = min_ratio

    # Compute final change in concentrations
    result = np.empty_like(species_concs_vector)
    for i in range(species_concs_vector.size):
        result[i] = change * rxn_stoichs[i]

    return result

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
                 freeze_kf=False,
                 freeze_kb=False,
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
            
        self._kf = kf
        self._kb = kb
        
        self._freeze_kf = freeze_kf
        self._freeze_kb = freeze_kb
        
        self.get_exponents_from_stoich = get_exponents_from_stoich
        
        if exponents is None:
            if get_exponents_from_stoich:
                self.exponents = np.abs(len(stoich))
            else:
                self.exponents = np.ones(len(stoich))
        else:
            self.exponents = exponents
        self.reactant_indices = np.where(stoich<0)[0]
        self.product_indices = np.where(stoich>0)[0]
        self._load_full_string()
    
    @property
    def kf(self):
        return self._kf
    @kf.setter
    def kf(self, new_kf):
        if not self._freeze_kf:
            self._kf = new_kf
        else:
            warn(f'kf for Reaction {self.ID} was not changed to {new_kf} as _freeze_kf was True.\n',
                 RuntimeWarning)
    
    @property
    def kb(self):
        return self._kb
    @kb.setter
    def kb(self, new_kb):
        if not self._freeze_kb:
            self._kb = new_kb
        else:
            warn(f'kb for Reaction {self.ID} was not changed to {new_kb} as _freeze_kb was True.\n',
                 RuntimeWarning)
            
    def get_dconcs_dt(self):
        kf, kb = self.kf, self.kb
        if kf==kb==0:
            return 0.
        return dconcs_dt(kf=kf, 
                         kb=kb,
                         species_concs_vector=self.species_system._concentrations, 
                         rxn_stoichs=self.stoichiometry,
                         rl_exps=self.exponents,
                         reactant_indices=self.reactant_indices,
                         product_indices=self.product_indices)
    
    def _load_full_string(self):
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
        param_info = f'kf={self.kf}'
        if arrow=='<->':
            param_info += ', ' + f'kb={self.kb}'
        self._lhs_string = lhs
        self._rhs_string = rhs
        self._param_info_string = param_info
        self._full_string = f'{self.ID}: Reaction(' + lhs + ' ' + arrow + ' ' + rhs + '; ' + param_info + ')'
    
    def __str__(self):
        self._load_full_string()
        return self._full_string
    
    def __repr__(self):
        return self.__str__()
    
    def from_equation(ID, chem_equation, species_system, 
                      kf=None, # overrides any parameter info in the chem_equation string
                      kb=None, # overrides any parameter info in the chem_equation string
                      exponents=None, get_exponents_from_stoich=None):
        kf_, kb_ = None, None
        freeze_kb = False
        if isinstance(chem_equation, str):
            freeze_kb = '->' in chem_equation and not '<->' in chem_equation
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
                        freeze_kf=False,
                        freeze_kb=freeze_kb,
                        get_exponents_from_stoich=get_exponents_from_stoich,)
        
Rxn = IrreversibleReaction = ReversibleReaction = IrrevRxn = RevRxn = Reaction
