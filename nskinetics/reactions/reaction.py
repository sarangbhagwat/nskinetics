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
    """
    Represents a chemical equation with defined stoichiometry.
    
    This class stores and formats the stoichiometry of a reaction,
    and is used to represent the equation symbolically and in string form.
    It can optionally be paired with an external object that also holds stoichiometric information.
    
    Parameters
    ----------
    ID : str
        Unique identifier for the chemical equation.
    species_system : SpeciesSystem
        System containing all chemical species.
    stoichiometry : array-like, optional
        Stoichiometric coefficients of all species (negative for reactants, positive for products).
    paired_obj : object, optional
        An object containing its own 'stoichiometry' attribute that this equation can reference 
        instead of passing 'stoichiometry' as an initialization argument.
    
    """
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
                if not np.abs(stoich)==1.: lhs+= str(-stoich) + ' '
                lhs+= chem.ID
            elif stoich>0: 
                if not rhs=='': rhs+= ' + '
                if not np.abs(stoich)==1.: rhs+= str(stoich) + ' '
                rhs+= chem.ID
                
        return f'{self.ID}: ChemicalEquation(' + lhs + ' -> ' + rhs + ')'
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def stoichiometry(self):
        """
        The stoichiometric coefficients for this chemical equation.
        
        Returns
        -------
        np.ndarray
            Array of stoichiometric values (negative for reactants, positive for products).
        
        """
        return self._stoichiometry
    @stoichiometry.setter
    def stoichiometry(self, new_stoichiometry):
        self._stoichiometry = new_stoichiometry
        if self.paired_obj is not None:
            if not self.paired_obj.stoichiometry==new_stoichiometry:
                warn(f'Replaced {self.ID} stoichiometry with {new_stoichiometry}, but this does not match the paired_obj stoichiometry: {self.paired_obj.stoichiometry}.',
                     RuntimeWarning)
    
    def from_string(ID, equation_str, species_system, paired_obj=None):
        """
        Create a ChemicalEquation object from a string representation of a reaction.
        
        Parameters
        ----------
        ID : str
            Identifier for the chemical equation.
        equation_str : str
            String representing the chemical equation, e.g., 
            "A + B -> C",
            "A + B <-> C",
            "A + B -> C, kf=20",
            "A + B <-> C, kf=20, kb=40.5".
        species_system : SpeciesSystem
            System containing all species involved in the reaction.
        paired_obj : object, optional
            Object that may contain a stoichiometry attribute.

        Returns
        -------
        ChemicalEquation
            Instantiated ChemicalEquation object.
        float
            Forward rate constant (kf) parsed from the string, if present.
        float
            Backward rate constant (kb) parsed from the string, if present.
        """
        
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
    Abstract base class to define the stoichiometry of a chemical reaction.
    
    This class does not define any kinetic behavior (i.e., no rate constants),
    but sets up the stoichiometry from either a ChemicalEquation, lists/dicts of
    reactants and products, or a raw stoichiometry array.
    
    One and only one of the following must be provided:
        - reactants and products (as dict or list)
        - chem_equation
        - stoichiometry
        
    Parameters
    ----------
    ID : str
        Identifier for the reaction.
    species_system : SpeciesSystem
        System containing the chemical species.
    reactants : list or dict, optional
        Reactant species. If dict, keys are species IDs and values are stoichiometric coefficients.
    products : list or dict, optional
        Product species. Same format as `reactants`.
    chem_equation : ChemicalEquation, optional
        Equation object defining reaction stoichiometry.
    stoichiometry : array-like, optional
        Array of stoichiometric coefficients.
        
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
def dconcs_dt_v0_1(kf, kb, species_concs_vector, rxn_stoichs, rl_exps,
              reactant_indices, product_indices):
    """
    (Old, unused.)
    
    Computation of species concentration changes using 
    forward and backward reaction rate constants and rate law exponents,
    with logic to identify limiting reactants and avoid negative concentrations.
    
    Uses numpy functions rather than loops and built-in Python operators.

    Parameters
    ----------
    kf : float
        Forward rate constant.
    kb : float
        Backward rate constant.
    species_concs_vector : ndarray
        Current species concentrations.
    rxn_stoichs : ndarray
        Stoichiometric coefficients.
    rl_exps : ndarray
        Rate law exponents.
    reactant_indices : ndarray
        Indices of reactant species.
    product_indices : ndarray
        Indices of product species.

    Returns
    -------
    ndarray
        Time derivatives of species concentrations.
    """
    # -----------------------------------------
    # for a one-way reaction,
    # returns array of temporal rates of change 
    # in the concentrations of given species
    # (- denotes decrease, + denotes increase)
    # -----------------------------------------
    
    # 'change' is the 1-stoichiometry-equivalent change for reactant concs
    
    # change = 0
    # if np.all(rl_exps==1.):
    # change = kf*np.prod(species_concs_vector[reactant_indices]) -\
    #     kb*np.prod(species_concs_vector[product_indices])
    # else:
    change = kf*np.prod(np.power(species_concs_vector[reactant_indices], rl_exps[reactant_indices])) -\
        kb*np.prod(np.power(species_concs_vector[product_indices], rl_exps[product_indices]))
    
    # if change is too great, cap it to the limiting reactant conc
    temp_vector = species_concs_vector + change*rxn_stoichs
    if np.any(temp_vector<0):
        tv_stoich_adj = temp_vector/np.abs(rxn_stoichs)
        tv_stoich_adj[np.where(np.isinf(tv_stoich_adj))] = 0.
        tv_stoich_adj[np.where(np.isnan(tv_stoich_adj))] = 0.
        limiting_reactant_index = np.where(tv_stoich_adj==np.min(tv_stoich_adj))[0][0]
        change = species_concs_vector[limiting_reactant_index]/np.abs(rxn_stoichs[limiting_reactant_index]) # limiting reactant conc. adjusted by stoichiometry
    return change * rxn_stoichs


@njit(cache=True)
def dconcs_dt_v0_2(kf, kb, species_concs_vector, rxn_stoichs, rl_exps,
              reactant_indices, product_indices):
    """
    Computation of species concentration changes using 
    forward and backward reaction rate constants and rate law exponents,
    with logic to identify limiting reactants and avoid negative concentrations.
    
    Parameters
    ----------
    kf : float
        Forward rate constant.
    kb : float
        Backward rate constant.
    species_concs_vector : ndarray
        Current species concentrations.
    rxn_stoichs : ndarray
        Stoichiometric coefficients.
    rl_exps : ndarray
        Rate law exponents.
    reactant_indices : ndarray
        Indices of reactants in species_concs_vector.
    product_indices : ndarray
        Indices of products in species_concs_vector.

    Returns
    -------
    ndarray
        Rate of change of concentrations of each species.
    """
    
    # Check if any reactant is exhausted:
    for ind in reactant_indices:
        if species_concs_vector[ind] <= 1e-20:
            return np.zeros(species_concs_vector.shape, dtype=species_concs_vector.dtype)
    
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

@njit(cache=True)
def dconcs_dt_v0_3(
    kf, kb, species_concs_vector, rxn_stoichs, rl_exps,
    reactant_indices, product_indices
    ):
    n_species = species_concs_vector.size
    
    # Early exit if any reactant is exhausted
    for idx in reactant_indices:
        if species_concs_vector[idx] <= 1e-20:
            return np.zeros(n_species, dtype=species_concs_vector.dtype)
        
    # Compute forward rate
    forward = kf
    for idx in reactant_indices:
        conc = species_concs_vector[idx]
        exp  = rl_exps[idx]
        forward *= conc ** exp
        
    # Compute backward rate
    backward = kb
    for idx in product_indices:
        conc = species_concs_vector[idx]
        exp  = rl_exps[idx]
        backward *= conc ** exp
        
    change = forward - backward
    
    # Evaluate limiting reactant (to avoid negative concentrations)
    min_ratio = 1e20
    limiting_found = False
    
    for i in range(n_species):
        stoich = rxn_stoichs[i]
        projected = species_concs_vector[i] + change * stoich
        if projected < 0.0:
            if stoich != 0.0:
                ratio = species_concs_vector[i] / abs(stoich)
                if ratio < min_ratio:
                    min_ratio = ratio
                    limiting_found = True
                    
    if limiting_found:
        change = min_ratio
        
    # Compute final concentration rate of change
    result = np.empty(n_species, dtype=species_concs_vector.dtype)
    for i in range(n_species):
        result[i] = change * rxn_stoichs[i]
        
    return result

class Reaction(AbstractReaction):
    """
    A complete chemical reaction object that includes kinetic parameters and rate law information.

    This class extends `AbstractReaction` by adding support for reaction rate constants,
    rate law exponents, reversible/irreversible behavior, and dynamic concentration updates.
    Includes logic to compute the rate of change in concentrations.
    
    One and only one of the following must be provided to define stoichiometry:
        - reactants and products (as dict or list)
        - chem_equation
        - stoichiometry
        
    Parameters
    ----------
    ID : str
        Identifier for the reaction.
    species_system : SpeciesSystem
        System containing the chemical species.
    kf : float
        Forward reaction rate constant.
    kb : float, optional
        Backward reaction rate constant (defaults to 0 for irreversible reactions).
    reactants : list or dict, optional
        Reactant species. If dict, keys are species IDs and values are stoichiometric coefficients.
    products : list or dict, optional
        Product species. Same format as `reactants`.
    chem_equation : ChemicalEquation, optional
        Equation object defining reaction stoichiometry.
    stoichiometry : array-like, optional
        Stoichiometric coefficients for all species.
    exponents : array-like, optional
        Rate law exponents for each species (same length as stoichiometry).
    get_exponents_from_stoich : bool, optional
        If True, use the absolute value of stoichiometry as the rate law exponents.
    freeze_kf : bool, optional
        If True, prevents modification of the forward rate constant, kf.
    freeze_kb : bool, optional
        If True, prevents modification of the backward rate constant, kb.
            
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
        """
        The forward reaction rate constant, kf.
    
        Returns
        -------
        float
            Forward reaction rate constant, kf.
            
        """
        return self._kf
    @kf.setter
    def kf(self, new_kf):
        """
        Set the forward reaction rate constant, kf.
        
        If `freeze_kf` is True, the value is not updated and a warning is issued.
        
        Parameters
        ----------
        new_kf : float
            New value for the forward reaction rate constant, kf.
            
        """
        if not self._freeze_kf:
            self._kf = new_kf
        else:
            warn(f'kf for Reaction {self.ID} was not changed to {new_kf} as _freeze_kf was True.\n',
                 RuntimeWarning)
    
    @property
    def kb(self):
        """
        The backward reaction rate constant, kf.
    
        Returns
        -------
        float
            Forward reaction rate constant, kb.
            
        """
        return self._kb
    @kb.setter
    def kb(self, new_kb):
        """
        Set the backward reaction rate constant, kb.
        
        If `freeze_kb` is True, the value is not updated and a warning is issued.
        
        Parameters
        ----------
        new_kb : float
            New value for the backward reaction rate constant, kb.
            
        """
        if not self._freeze_kb:
            self._kb = new_kb
        else:
            warn(f'kb for Reaction {self.ID} was not changed to {new_kb} as _freeze_kb was True.\n',
                 RuntimeWarning)
            
    def get_dconcs_dt(self):
        """
        Calculate the net rate of change in species concentrations for this reaction,
        using the current kinetic parameters and species concentrations.
        If self.species_system.log_transformed is True, returns rate 
        of change of np.log(concentrations) instead.
        
        Returns
        -------
        ndarray
            Change in concentration; or 
            zeros if kf and kb are both zero; or
            change in log(concentrations) if self.log_transformed is True.
        
        """
        kf, kb = self.kf, self.kb
        if kf==kb==0:
            return np.zeros(shape=self.species_system.concentrations.shape)
        
        return dconcs_dt_v0_3(kf=kf, 
                             kb=kb,
                             species_concs_vector=self.species_system._concentrations, 
                             rxn_stoichs=self.stoichiometry,
                             rl_exps=self.exponents,
                             reactant_indices=self.reactant_indices,
                             product_indices=self.product_indices)
    
    def _load_full_string(self):
        """
        Internal method to format and cache string representations of 
        the reaction's left-hand side, right-hand side, and kinetic parameters.
        
        """
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
        """
        Create a Reaction object from a ChemicalEquation object or string.
        
        Parameters
        ----------
        ID : str
            Identifier for the Reaction.
        chem_equation : str or ChemicalEquation
            Reaction as a string or ChemicalEquation object.
            String represents the chemical equation and may or may not
            include kf and kb (if not included, must pass at least kf 
                               -- kb defaults to zero -- or
                               both kf and kb as arguments). 
            E.g., 
            "A + B -> C",
            "A + B <-> C",
            "A + B -> C, kf=20",
            "A + B <-> C, kf=20, kb=40.5".
        species_system : SpeciesSystem
            System containing the involved species.
        kf : float, optional
            Forward rate constant (overrides any value parsed from str chem_equation).
        kb : float, optional
            Backward rate constant (overrides any value parsed from str chem_equation).
        exponents : ndarray, optional
            Rate law exponents.
        get_exponents_from_stoich : bool, optional
            Whether to use absolute stoichiometry as rate law exponents.
            
        Returns
        -------
        Reaction
            Instantiated Reaction object.
            
        """
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
