# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from ..reactions import Reaction
from ..utils import create_function

__all__ = ('ReactionSystem', 'RxnSys')

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
        self._reactions = _reactions
        
        self._species_system = species_system
        self._solution = None # stored solution from the most recent 'solve' call
    
    @property
    def reactions(self):
        return self._reactions
    
    @property
    def species_system(self):
        return self._species_system
    
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
    
    def add_reaction(self, reaction):
        r = reaction
        _reactions = self._reactions
        if isinstance(r, str):
            _reactions.append(Reaction.from_equation(ID=self.ID+f'_r{len(_reactions)}', 
                                                     chem_equation=r, 
                                                     species_system=self.species_system))
        elif isinstance(r, Reaction) or isinstance(r, ReactionSystem):
            _reactions.append(r)
    
    def change_reaction(self, index, 
                        new_equation_string=None,
                        new_kf=None, new_kb=None,
                        ):
        _reactions = self._reactions
        if new_equation_string is None and new_kf is None and new_kb is None:
            raise ValueError('Either new_equation_string or new_kf and new_kf must be provided.')
        elif new_equation_string is not None:
            _reactions[index] = Reaction.from_equation(chem_equation=new_equation_string,
                                                       species_system=self.species_system)
        else:
            r = _reactions[index]
            if new_kf is not None: r.kf = new_kf
            if new_kb is not None: r.kf = new_kb
            
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
