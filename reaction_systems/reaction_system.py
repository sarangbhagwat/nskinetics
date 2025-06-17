# -*- coding: utf-8 -*-
"""
Created on Thu May 29 17:39:02 2025

@author: sarangbhagwat
"""

import numpy as np
from scipy.integrate import solve_ivp

__all__ = ('ReactionSystem', 'RxnSys')

#%% Utility functions

def codify(statement):
    statement = replace_apostrophes(statement)
    statement = replace_newline(statement)
    return statement

def replace_newline(statement):
    statement = statement.replace('\n', ';')
    return statement

def replace_apostrophes(statement):
    statement = statement.replace('’', "'").replace('‘', "'").replace('“', '"').replace('”', '"')
    return statement

def create_function(code, namespace):
    def wrapper_fn(statement):
        def f(t, concs):
            namespace['t'] = t
            namespace['concs'] = concs
            exec(codify(statement), namespace)
            return namespace['y']
        return f
    function = wrapper_fn(code)
    return function

#%% Reaction system

class ReactionSystem():
    """
    Abstract class for a system of reactions.
    
    Parameters
    ----------
    ID : str
        ID.
    reactions : list
        List of Rxn, ReversibleRxn, or RxnSystem objects.
    species_system : SpeciesSystem
        A SpeciesSystem object containing all species
        involved in this system of reactions.
    """
    def __init__(self, ID, reactions, species_system):
        self.ID = ID
        self.reactions = reactions
        self.species_system = species_system
    
    def get_dconcs_dt(self):
        reactions = self.reactions
        # species_concs_vector = self.species_system.all_sps
        # breakpoint()
        return np.sum([r.get_dconcs_dt() for r in reactions], axis=0)
    
    def solve(self, 
              t_span,
              t_eval=None,
              method='LSODA',
              atol=None, rtol=1e-6, 
              events=None,
              sp_conc_for_events=None, # dict or None
              dense_output=False,
              y0=None):
        get_dconcs_dt = self.get_dconcs_dt
        sp_sys = self.species_system
        concentrations = sp_sys.concentrations
        if atol is None:
            atol = 1e-6*max(concentrations)
        if y0 is None:
            y0 = concentrations
        if sp_conc_for_events is not None:
            if events is None:
                events = []
            code = 'y = concs[index] - S'
            for sp, conc in sp_conc_for_events.items():
                index = sp_sys.index_from_ID(sp) if isinstance(sp, str) else sp_sys.index(sp)
                events.append(create_function(code=code, 
                                              namespace={'S': conc,
                                                         'index': index,
                                                         'y': None}
                                              ))
        def ode_system_RHS(t, concs):
            concs[np.where(concs<0)] = 0. # not strictly necessary with a low enough atol
            sp_sys.concentrations = concs
            return get_dconcs_dt()
        
        return solve_ivp(ode_system_RHS, 
                         t_span=t_span, 
                         y0=y0,
                         t_eval=t_eval,
                         atol=atol, # recommended: <= 1e-6*max(sp_sys.concentrations)
                         rtol=atol, # recommended: 1e-6
                         # the solver keeps the local error estimates less than atol + rtol * abs(y)
                         events=events,
                         method=method,
                         dense_output=dense_output)
    
RxnSys = ReactionSystem
