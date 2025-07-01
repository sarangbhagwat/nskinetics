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
from scipy.interpolate import interp1d
from ..reactions import Reaction
from ..utils import create_function, is_number, is_array_of_numbers, is_list_of_strings

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
        self._C_at_t_is_updated = False
        self._C_at_t_f_all = None
        self._C_at_t_fs_indiv_sps = None
        
    @property
    def reactions(self):
        return self._reactions
    
    @property
    def species_system(self):
        return self._species_system
    
    def get_dconcs_dt(self):
        reactions = self.reactions
        return np.sum([r.get_dconcs_dt() for r in reactions], axis=0)
    
    def _solve_single_phase(self, 
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
        
        sol = solve_ivp(ode_system_RHS, 
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
    
    def solve(self, 
              t_span,
              sp_conc_for_events=None, # dict or None
              events=None,
              spikes=None,
              method='LSODA',
              t_eval=None,
              atol=None, rtol=1e-6, 
              dense_output=False,
              y0=None,
              dt_spike=1e-6, # long dt_spike (e.g., slow feeding) not supported, only spikes in concentrations
              ):
        """
        Get concentration vs. time data.
        
        Arguments
        ---------
        t_span: list or tuple
            Two-item iterable consisting of the 
            desired start time and end time.
        sp_conc_for_events: dict, optional
            Keys are Species objects or Species ID strings.
            Values are concentrations. When a Species concentrations
            achieves the corresponding value, the time at which this
            event occurs is stored in ReactionSystem._solution['t_events'].
        spikes : dict, optional
            A dictionary with time stamps (t) as keys.
            Values for spikes can assume one of the following forms:
                (A) string or list of strings each of the form 
                    f'Change; {species_ID}; {conc}' 
                    (e.g., 'Change; Substrate; 0.001'), where conc is 
                    the desired change (spike) in the concentration 
                    of the corresponding Species at time t;
                (B) string or list of strings each of the form 
                    f'Target; {species_ID}; {conc}' 
                    (e.g., 'Target; Substrate; 0.001'), where conc is 
                    the concentration of the corresponding Species 
                    to be targeted via a change (spike) at time t;
                (C) string or list of strings each of the form 
                    f'Target (negative change allowed); {species_ID}; {conc}' 
                    (e.g., 'Target (negative change allowed); Substrate; 0.001'), 
                    where conc is the concentration of the corresponding Species 
                    to be targeted via a change (spike) at time t, with
                    negative values for the change (i.e., dips) permitted;
                (D) array of floats equal to the change (spike) in 
                    concentrations of all Species in the species_system
                    (0 if none) at time t; or
                (E) function that accepts the time and 
                    current concentrations array as arguments, and 
                    returns the desired change (spike) in concentrations
                    of all Species in the species_system (0 if none) at time t.
        
        """
        self._spikes = spikes
        self._spikes_list = sl = self._get_spikes_list_from_dict(spikes)
        tmin, tmax = t_span
        # dt_spike=max(1e-6, dt_spike)
        concentrations = self.species_system.concentrations
        if atol is None:
            atol = 1e-6*max(concentrations)
        if y0 is None:
            y0 = concentrations
            
        sols = []
        fs = [i() if callable(i) else i for i in sl]
        _solve_single_phase = self._solve_single_phase
        
        # if there are feed spikes
        if fs is not None and not fs==[]:
            tmin_curr = tmin
            y0_curr = y0
            dconcs_curr = np.zeros(len(concentrations))
            # simulate phases tmin -> feed spike time 1 -> feed spike time 2 ...
            for t, dconcs in self._spikes_list:
                assert tmin_curr < t
                y0_curr += dconcs_curr(t=t, concs=y0_curr)\
                    if callable(dconcs_curr)\
                    else np.array(dconcs_curr)
                
                sols.append(_solve_single_phase((tmin_curr, t),
                                                 t_eval=t_eval,
                                                 method=method,
                                                 atol=atol, rtol=rtol, 
                                                 events=events,
                                                 sp_conc_for_events=sp_conc_for_events, # dict or None
                                                 dense_output=dense_output,
                                                 y0=y0_curr))
                
                tmin_curr = t+dt_spike
                y0_curr = sols[-1].y[:, -1]
                dconcs_curr = dconcs
            
            # simulate last phase (last feed spike -> tmax)
            y0_curr += dconcs_curr(t=t, concs=y0_curr)\
                if callable(dconcs_curr)\
                else np.array(dconcs_curr)
            sols.append(_solve_single_phase((tmin_curr, tmax),
                                             t_eval=t_eval,
                                             method=method,
                                             atol=atol, rtol=rtol, 
                                             events=events,
                                             sp_conc_for_events=sp_conc_for_events, # dict or None
                                             dense_output=dense_output,
                                             y0=y0_curr))
            
        # if no feed spikes     
        else:
            # simulate single phase
            sols.append(_solve_single_phase((tmin, tmax),
                                             t_eval=t_eval,
                                             method=method,
                                             atol=atol, rtol=rtol, 
                                             events=events,
                                             sp_conc_for_events=sp_conc_for_events, # dict or None
                                             dense_output=dense_output,
                                             y0=y0))
            
        y_final = np.concatenate([sol.y.transpose() for sol in sols])
        t_final = list(sols[0].t)
        t_events = list(sols[0].t_events[0])
        y_events = list(sols[0].y_events[0])
        
        for sol in sols[1:]:
            # for arr1, arr2 in zip([t_final, y_final],
            #                       [sol.t, sol.y]):
            #     t_final = np.concatenate((t_final, sol.t.transpose()))
            #     y_final = np.concatenate((y_final, sol.y.transpose()))
            t_final += list(sol.t)
            t_events += list(sol.t_events[0])
            y_events += list(sol.y_events[0])
        
        if events is None:
            events = []
        
        solution = {'t': np.array(t_final).transpose(),
                    'y': np.array(y_final).transpose(),
                    't_events': np.array(t_events),
                    'y_events': np.array(y_events),
                    'sol': sols,
                    'events': events+['['+k+']' + ' = ' + str(v) 
                                      for k,v in sp_conc_for_events.items()]
                    }
        
        self._solution = solution
        self._C_at_t_is_updated = False
        return solution
    
    def C_at_t(self, t,  species=None):
        # Note on speed-up for C_at_t
        # Creating separate interp1d objects for each species concentration
        # results in faster C_at_t calls with species specified.
        # Creating a single interp1d object for all species concentrations
        # results in faster C_at_t calls with species not specified.
        # We create both here (if not self._C_at_t_is_updated; note this step
        # is slower as a result) and reference the faster object 
        # for the given call.
        
        species_system = self.species_system
        all_sps = species_system.all_sps
        index_f = self.species_system.index
        
        if not self._C_at_t_is_updated:
            _solution = self._solution
            _t, _y = _solution['t'], _solution['y']
            self._C_at_t_f_all = interp1d(_t, _y)
            self._C_at_t_fs_indiv_sps = [interp1d(_t, _y[index_f(sp), :]) 
                                         for sp in all_sps]
            self._C_at_t_is_updated = True
        
        if species is not None:
            ind = self.species_system.index(species)
            return self._C_at_t_fs_indiv_sps[ind](t)
        else:
            return self._C_at_t_f_all(t)
        
    def plot_solution(self, show_events=True, sps_to_include=None):
        if sps_to_include is None:
            sps_to_include = [i.ID for i in self.species_system.all_sps]
        
        sol = self._solution
        t, y = sol['t'], sol['y']
        t_events, y_events = sol['t_events'], sol['y_events']
        events = sol['events']
        all_sps = self.species_system.all_sps
        
        fig, ax = plt.subplots()
        
        for i, sp in zip(range(len(all_sps)), all_sps):
            if sp in sps_to_include or sp.ID in sps_to_include:
                ax.plot(t, y[i, :], label=sp.ID,
                        linestyle='solid',
                        linewidth=1.)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Concentration [mol/L]')
        
        if show_events:
            ylim = ax.get_ylim()
            for t, e in zip(t_events, events):
                label = str(e)
                if callable(e):
                    label = label.split(' at ')[0].remove('<')
                    label += ' = 0'
                label = 't | ' + label
                ax.vlines(t, ylim[0], ylim[1], 
                          linestyles='dashed', linewidth=0.5,
                          color='blue',
                          # label=label,
                          )
                xlim = ax.get_xlim()
                ax.annotate(label, 
                            xy=((t + 0.04*(xlim[1]-xlim[0])), 
                                (ylim[1]-ylim[0])/2), 
                            verticalalignment='center', 
                            horizontalalignment='right' , 
                            rotation = -270,
                            color='blue',
                            )
        plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
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
    
    def _get_spikes_list_from_dict(self, spikes_dct):
        sd = spikes_dct
        if sd is None:
            return []
        sl = []
        for t, dconcs in sd.items():
            assert is_number(t)
            if is_array_of_numbers(dconcs):
                pass
            elif isinstance(dconcs, str):
                dconcs = self._get_dconcs_from_str(dconcs)
            elif isinstance(dconcs, list):
                if is_array_of_numbers(np.array(dconcs)):
                    pass
                elif is_list_of_strings(dconcs):
                    dconcs_f_list = [self._get_dconcs_from_str(i) for i in dconcs]
                    code = 'y = np.sum([i() for i in dconcs_f_list])'
                    dconcs = create_function(code=code,
                                             namespace={'np': np,
                                                        'dconcs_f_list':dconcs_f_list})
            # finally, append to list
            sl.append((t, dconcs))
        return sl
    
    def _get_dconcs_from_str(self, s):
        sp_sys = self.species_system
        dconcs = None
        split_s = [i.replace(' ', '')
                   for i in s.split(';')]
        action = split_s[0].lower()
        action, species, value = split_s
        action = action.lower()
        value = float(value)
        
        sp_ind = sp_sys.index(species)
        
        if action=='change':
            dconcs = np.zeros((len(sp_sys.concentrations)))
            dconcs[sp_ind] = value
        elif action=='target':
            code = 'y = np.zeros(len(concs)); y[sp_ind] = max(0, target_value-concs[sp_ind])'
            dconcs = create_function(code=code, 
                                     namespace={'sp_ind':sp_ind,
                                                'target_value':value,
                                                'np': np})
        elif action=='target (negative change allowed)':
            code = 'y = np.zeros(len(concs)); y[sp_ind] = target_value-concs[sp_ind]'
            dconcs = create_function(code=code, 
                                     namespace={'sp_ind':sp_ind,
                                                'target_value':value,
                                                'np': np})
        
        return dconcs
    
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
