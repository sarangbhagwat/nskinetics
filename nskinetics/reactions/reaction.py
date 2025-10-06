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
    
    @classmethod
    def from_string(cls, ID, equation_str, species_system, paired_obj=None):
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
            "A + B -> C; kf=20",
            "A + B <-> C; kf=20, kb=40.5",
            "A + B <-> C; kf=20, kb = 10.5",
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
        str
            Rate expression (rate) parsed from the string, if present.
            
        """
        
        stoichiometry, kf, kb = read_equation_str(equation_str, 
                                                   species_system)
        return cls(ID=ID, 
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

@njit(cache=True)
def dconcs_dt_custom_rate_f(
    change,
    species_concs_vector, rxn_stoichs, rl_exps,
    reactant_indices, product_indices
    ):
    
    n_species = species_concs_vector.size
    
    # Early exit if any reactant is exhausted
    for idx in reactant_indices:
        if species_concs_vector[idx] <= 1e-20:
            return np.zeros(n_species, dtype=species_concs_vector.dtype)
        
    # # Compute forward rate
    # forward = kf
    # for idx in reactant_indices:
    #     conc = species_concs_vector[idx]
    #     exp  = rl_exps[idx]
    #     forward *= conc ** exp
        
    # # Compute backward rate
    # backward = kb
    # for idx in product_indices:
    #     conc = species_concs_vector[idx]
    #     exp  = rl_exps[idx]
    #     backward *= conc ** exp
        
    # change = forward - backward
    
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
    rate_f : function, optional
        Custom rate function (overrides any value parsed from str chem_equation).
        Must accept species_concs_vector, rxn_stoichs, rl_exps, reactant_indices, 
        and product_indices as arguments (see Reaction.get_dconcs_dt for example use)
        and additional keyword arguments for parameters from Reaction.rate_params.
    rate_params: dict, optional
        Dictionary of parameters passed as keyword arguments to Reaction.rate_f.
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
                 rate_f=None,
                 rate_params=None,
                 reactants=None, products=None,
                 chem_equation=None,
                 stoichiometry=None,
                 exponents=None,
                 is_transport=False, # note all reactants must have the same compartment and all products must have the same compartment
                 # also, a transport reaction should only dictate how concentrations change
                 is_multicompartment=False, # all reactants must have the same compartment; products may all have mutually different compartments
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
            
        self.rate_f = rate_f
        
        if rate_params is None:
            rate_params = {}
            rate_params['kf'] = float(kf)
            rate_params['kb'] = float(kb)
        
        self.rate_params = rate_params
        
        self._freeze_params = set()
        
        # self._freeze_kf = freeze_kf
        # self._freeze_kb = freeze_kb
        
        self.get_exponents_from_stoich = get_exponents_from_stoich
        
        if exponents is None:
            if get_exponents_from_stoich:
                self.exponents = np.abs(len(stoich))
            else:
                self.exponents = np.ones(len(stoich))
        else:
            self.exponents = exponents
        self.reactant_indices = reactant_indices = np.where(stoich<0)[0]
        self.product_indices = product_indices = np.where(stoich>0)[0]
        self._load_full_string()
        
        all_sps = species_system.all_sps
        try:
            self.source = all_sps[reactant_indices[0]].compartment
        except:
            breakpoint()
            
        self.is_transport = is_transport
        self.is_multicompartment = is_multicompartment
        
        if not (is_transport or is_multicompartment):
            self.destination = self.source
        elif is_transport:
            self.destination = all_sps[product_indices[0]].compartment
        
    @property
    def kf(self):
        """
        The forward reaction rate constant, kf.
    
        Returns
        -------
        float
            Forward reaction rate constant, kf.
            
        """
        return self.rate_params['kf']
    @kf.setter
    def kf(self, new_kf):
        """
        Set the forward reaction rate constant, kf.
        
        If 'kf' is `_freeze_params`, the value is not updated and a warning is issued.
        
        Parameters
        ----------
        new_kf : float
            New value for the forward reaction rate constant, kf.
            
        """
        if not 'kf' in self._freeze_params:
            self.rate_params['kf'] = new_kf
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
        return self.rate_params['kb']
    @kb.setter
    def kb(self, new_kb):
        """
        Set the backward reaction rate constant, kb.
        
        If 'kB' is `_freeze_params`, the value is not updated and a warning is issued.
        
        Parameters
        ----------
        new_kb : float
            New value for the backward reaction rate constant, kb.
            
        """
        if not 'kb' in self._freeze_params:
            self.rate_params['kb'] = new_kb
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
        rate_f = self.rate_f
        sp_sys = self.species_system
        species_concs_vector = sp_sys._concentrations
        stoichiometry = self.stoichiometry
        exponents = self.exponents
        reactant_indices = self.reactant_indices
        product_indices = self.product_indices
        rate_params = self.rate_params
        
        dconcs_dt_vector = None
        
        if rate_f is not None:
            change = rate_f(species_concs_vector=species_concs_vector, 
            rxn_stoichs=stoichiometry,
            rl_exps=exponents,
            reactant_indices=reactant_indices,
            product_indices=product_indices,
            **rate_params)
            
            dconcs_dt_vector =  dconcs_dt_custom_rate_f(change=change, 
            species_concs_vector=species_concs_vector, 
            rxn_stoichs=stoichiometry,
            rl_exps=exponents,
            reactant_indices=reactant_indices,
            product_indices=product_indices,)
        
        else:
            kf, kb = self.kf, self.kb
            if not kf==kb==0:
                dconcs_dt_vector = dconcs_dt_v0_3(kf=kf, 
                                     kb=kb,
                                     species_concs_vector=species_concs_vector, 
                                     rxn_stoichs=stoichiometry,
                                     rl_exps=exponents,
                                     reactant_indices=reactant_indices,
                                     product_indices=product_indices)
            
            else:
                dconcs_dt_vector = np.zeros(shape=species_concs_vector.shape)
            
            if self.is_transport: # multiply by volume ratio for conservation of mass 
                get_comp_vol = sp_sys.get_compartment_volume
                vol_ratio = get_comp_vol(self.source)/get_comp_vol(self.destination)
                for i in product_indices:
                    dconcs_dt_vector[i] *= vol_ratio
            elif self.is_multicompartment: # multiply each species conc change by volume ratio for conservation of mass
                get_comp_vol = sp_sys.get_compartment_volume
                src_vol = get_comp_vol(self.source)
                all_sps = sp_sys.all_sps
                for i in product_indices:
                    dconcs_dt_vector[i] *= src_vol/get_comp_vol(all_sps[i].compartment)
                    
        return dconcs_dt_vector
    
    def _load_full_string(self):
        """
        Internal method to format and cache string representations of 
        the reaction's left-hand side, right-hand side, and kinetic parameters.
        
        """
        lhs = ''
        rhs = ''
        rate_f = self.rate_f
        arrow = ''
        if rate_f is not None:
            arrow = '<->'
        elif rate_f is None and self.kb==0.:
            arrow = '->' 
        else:
            arrow = '<->'
            
        for chem, stoich in zip(self.species_system.all_sps, self.stoichiometry):
            
            if stoich<0: 
                if not lhs=='': lhs+= ' + '
                if not np.abs(stoich)==1.: lhs+= str(-stoich) + ' '
                lhs+= chem.ID
            elif stoich>0: 
                if not rhs=='': rhs+= ' + '
                if not np.abs(stoich)==1.: rhs+= str(stoich) + ' '
                rhs+= chem.ID
        # param_info = f'kf={self.kf}'
        param_info = ', '.join([f'{k}={v}' for k, v in self.rate_params.items()])
        
        if not arrow=='<->':
            kb = self.kb
            param_info = param_info.replace(', ' + f'kb={kb}', '')
            param_info = param_info.replace(f'kb={kb}', '')
            
        self._lhs_string = lhs
        self._arrow_string = arrow
        self._rhs_string = rhs
        self._param_info_string = param_info
        full_string = f'{self.ID}: Reaction(' + lhs + ' ' + arrow + ' ' + rhs + '; '
        if rate_f is not None:
            full_string += 'custome rate_f; '
        full_string+= param_info + ')'
        self._full_string = full_string
        
    def get_equation_str(self):
        return self._lhs_string + ' ' + self._arrow_string + ' ' + self._rhs_string
        
    def __str__(self):
        self._load_full_string()
        return self._full_string
    
    def __repr__(self):
        return self.__str__()
    
    @classmethod
    def from_equation(cls, ID, chem_equation, species_system, 
                      kf=None, # overrides any parameter info in the chem_equation string
                      kb=None, # overrides any parameter info in the chem_equation string
                      rate_f=None, # overrides any parameter info in the chem_equation string
                      rate_params=None, # overrides any parameter info in the chem_equation string
                      exponents=None, get_exponents_from_stoich=False,
                      is_transport=False, is_multicompartment=False):
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
        rate_f : function, optional
            Custom rate function (overrides any value parsed from str chem_equation).
            Must accept species_concs_vector, rxn_stoichs, rl_exps, reactant_indices, 
            and product_indices as arguments (see Reaction.get_dconcs_dt for example use)
            and additional keyword arguments for parameters from Reaction.rate_params.
        rate_params: dict, optional
            Dictionary of parameters passed as keyword arguments to Reaction.rate_f.
            Each parameter value can be a float or 1-d array.
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
            
        return cls(ID=ID,
                species_system=species_system,
                chem_equation=chem_equation,
                kf=kf,
                kb=kb,
                rate_f=rate_f,
                rate_params=rate_params,
                exponents=exponents,
                freeze_kf=False,
                freeze_kb=freeze_kb,
                get_exponents_from_stoich=get_exponents_from_stoich,
                is_transport=is_transport,
                is_multicompartment=is_multicompartment)
        
Rxn = IrreversibleReaction = ReversibleReaction = IrrevRxn = RevRxn = Reaction
