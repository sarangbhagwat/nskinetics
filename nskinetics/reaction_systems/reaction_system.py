# -*- coding: utf-8 -*-
"""
Created on Thu May 29 17:39:02 2025

@author: sarangbhagwat
"""

import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from ..reactions import Reaction

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
        List of Reaction, ReactionSystem, or str objects.
        If str, must include chemical equation and kinetic 
        parameter info.
    species_system : SpeciesSystem
        A SpeciesSystem object containing all species
        involved in this system of reactions.
    """
    def __init__(self, ID, reactions, species_system):
        self.ID = ID
        _reactions = []
        i = 0
        for r in reactions:
            if isinstance(r, str):
                _reactions.append(Reaction.from_equation(ID=ID+f'_r{i}', 
                                                         chem_equation=r, 
                                                         species_system=species_system))
            elif isinstance(r, Reaction) or isinstance(r, ReactionSystem):
                _reactions.append(r)
            i += 1
        self.reactions = _reactions
        
        self.species_system = species_system
        self._solution = None # stored solution from the most recent 'solve' call
        
    def get_dconcs_dt(self):
        reactions = self.reactions
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
        
        self._solution = sol = solve_ivp(ode_system_RHS, 
                         t_span=t_span, 
                         y0=y0,
                         t_eval=t_eval,
                         atol=atol, # recommended: <= 1e-6*max(sp_sys.concentrations)
                         rtol=atol, # recommended: 1e-6
                         # the solver keeps the local error estimates less than atol + rtol * abs(y)
                         events=events,
                         method=method,
                         dense_output=dense_output)
        return sol
    
    def plot_solution(self, show_events=True, sps_to_include=None):
        if sps_to_include is None:
            sps_to_include = [i.ID for i in self.species_system.all_sps]
        
        sol = self._solution
        t, y = sol.t, sol.y
        t_events, y_events = sol.t_events, sol.y_events
        all_sps = self.species_system.all_sps
        
        fig, ax = plt.subplots()
        
        for i, sp in zip(range(len(all_sps)), all_sps):
            if sp in sps_to_include or sp.ID in sps_to_include:
                ax.plot(t, y[i, :], label=sp.ID,
                        linestyle='solid',
                        linewidth=1.)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Concentration [mol/L]')
        plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
        
        if show_events:
            ylim = ax.get_ylim()
            ax.vlines(t_events, ylim[0], ylim[1], 
                      linestyles='dashed', linewidth=0.5)
            
        plt.show()
    
    def __str__(self):
        rxns = self.reactions
        str_ = f'{self.ID}: ReactionSystem(\n'
        for r in rxns:
            str_ += '    ' + r.__str__() + '\n'
        str_ += '    ' + ')'
        return str_
    
    def __repr__(self):
        return self.__str__()
    
RxnSys = ReactionSystem
