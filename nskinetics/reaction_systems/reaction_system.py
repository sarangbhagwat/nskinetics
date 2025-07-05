# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np
import pandas as pd
import time

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from typing import Union, Tuple, List        
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoLocator, AutoMinorLocator


from ..reactions import Reaction
from ..utils import create_function, is_number, is_array_of_numbers,\
                    is_list_of_strings, fit_multiple_dependent_variables

__all__ = ('ReactionSystem', 'RxnSys')

np_array = np.array

#%% Reaction system

class ReactionSystem():
    """
    Represents a system of one or more chemical reactions involving 
    dynamic changes in species concentrations over time.

    This class allows the construction, simulation, and analysis of 
    reaction systems governed by mass action kinetics or other 
    user-defined rate laws. It supports features such as species spikes, 
    event-triggered integration, parameter fitting, and nesting of 
    reaction systems.
    
    Parameters
    ----------
    ID : str
        Unique identifier for the reaction system.
        
    reactions : list of str, Reaction, or ReactionSystem
        A list defining the reactions in the system. Each item can be:
        - A `Reaction` object
        - Another `ReactionSystem` object (for hierarchical systems)
        - A string describing a chemical equation with kinetic information
        
    species_system : SpeciesSystem
        The species system containing all species involved in the reaction system.
        
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
        
        self._exclude_frozen_params = True
        self._interp1d_fill_value = None
        self._timeout_solve_ivp = None
        
    @property
    def reactions(self):
        """
        List of reactions and nested reaction systems in this reaction system.
        
        Returns
        -------
        list
            List of Reaction or ReactionSystem objects.
        
        """
        return self._reactions
    
    @property
    def species_system(self):
        """
        The SpeciesSystem associated with this reaction system.
        
        Returns
        -------
        SpeciesSystem
            Object containing all species in the system.
            
        """
        return self._species_system
    
    def get_dconcs_dt(self):
        """
        Compute the total rate of change of concentrations for all species
        in the reaction system by summing the contributions from all reactions.
        
        Returns
        -------
        numpy.ndarray
            Array representing the net rate of change of concentrations for each species.
        
        """
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
        """
        Solve a single phase of the reaction system.
        
        Parameters
        ----------
        t_span : tuple
            Time interval (start, end) for integration.
        t_eval : array_like, optional
            Time points at which to store the computed solution.
        method : str, default 'LSODA'
            Integration method to use.
        atol : float, optional
            Absolute tolerance.
        rtol : float, optional
            Relative tolerance.
        events : list of callables, optional
            List of event functions.
        sp_conc_for_events : dict, optional
            Dictionary of species and target concentrations to trigger events.
        dense_output : bool, default False
            Whether to compute a continuous solution.
        y0 : array_like, optional
            Initial concentrations of species.
            
        Returns
        -------
        scipy.integrate.OdeResult
            Object containing time points, solution, and event information.
            
        """
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
        
        # print('Solving ...')
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
        # print('Solved.')
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
              remove_negative_concs=True, # can safely do this as negatives don't affect calculation with ode_system_RHS
              filename=None,
              save_events=True,
              ):
        """
        Simulate the reaction system and get concentration vs. time data.
        
        Parameters
        ----------
        t_span: list or tuple
            Two-item iterable consisting of the 
            desired start time and end time.
        sp_conc_for_events: dict, optional
            Keys are Species objects or Species ID strings.
            Values are concentrations. When a Species concentrations
            achieves the corresponding value, the time at which this
            event occurs is stored in ReactionSystem._solution['t_events']
            and the concentrations array at that time is stored in
            ReactionSystem._solution['y_events']. Event indices are stored
            in ReactionSystem._solution['events'].
        events: callable, list or tuple, optional
            Function or list of functions that accepts arguments t (time)
            and concentrations (ndarray of species concentrations ordered the
            same as ReactionSystem.species_system.all_sps). Each function
            represents an event (triggered when the function returns zero),
            When the function returns zero, the time at which this event
            event occurs is stored in ReactionSystem._solution['t_events']
            and the concentrations array at that time is stored in
            ReactionSystem._solution['y_events']. Event indices are stored
            in ReactionSystem._solution['events'].
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
        
        filename: str, optional
            Filename or path including filename to save results. Results will 
            be saved to an Excel file in the following format:
            Sheet: Main
                t     [Species ID 1]     [Species ID 2]     ...
                -     ------------       --------------        
                .     .                  .                  ...
                .     .                  .                  ...
                .     .                  .                  ...
            Sheet: Events
                event     t_event    [Species ID 1] at t_event    [Species ID 2] at t_event      ...
                -----     -------    -------------------------    -------------------------
                .         .          .                            .                              ...
                .         .          .                            .                              ...
                .         .          .                            .                              ...
        
        See Also
        --------
        `scipy.integrate.solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_
        
        Examples
        --------
        >>> import nskinetics as nsk
        >>> sp_sys = nsk.SpeciesSystem('sp_sys', 
        >>>                        ['E', 'S', 'ES', 'P'], # enzyme, substrate, enzyme-substrate complex, product
        >>>                        concentrations=[1e-4, 1e-4, 0., 0.])
        >>> reactions = [
        >>>             'E + S <-> ES; kf = 12.0, kb = 10.0', # kf = kon, kb = koff
        >>>             'ES -> E + P; kf = 32.0' # kf = kcat (enzyme turnover number)
        >>>             ]
        >>> rxn_sys = nsk.ReactionSystem(ID='rxn_sys', 
        >>>                                  reactions=reactions,
        >>>                                  species_system=sp_sys)
        >>> rxn_sys.solve(t_span=[0, 2*24*3600], # I want to simulate the system over 2 days
        >>>                  sp_conc_for_events={'S':1e-6}, # In addition to a full simulation,
        >>>                  )                              # I want to know the time at which [S] drops to 1e-6
        >>> np.allclose(rxn_sys._solution['t_events'], np.array([42219.44616989]), rtol=1e-5, atol=1e-8)
        True
        >>> np.allclose(rxn_sys._solution['y_events'], 
        >>>        np.array([[9.99998909e-05],
        >>>                 [1.00000000e-06],
        >>>                 [1.09091871e-10],
        >>>                 [9.89998909e-05]]), 
        >>>                 rtol=1e-5, atol=1e-8)
        True
        
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
        
        events_non_sp_conc = events if events is not None else []
        
        timeout = self._timeout_solve_ivp
        
        if timeout is not None:
            start_time = time.time()
            def timeout_event(t, y):
                """
                Callback function to check for timeout.
                Returns 0 to trigger termination if timeout exceeded.
                """
                elapsed_time = time.time() - start_time
                # print(start_time, elapsed_time, timeout)
                if elapsed_time > timeout:
                    return 0  # Event occurs when elapsed_time - timeout > 0
                else:
                    return 1
            # Set the event to terminate integration
            timeout_event.terminal = True
                
            events_non_sp_conc.append(timeout_event)
        
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
                    else np_array(dconcs_curr)
                
                sols.append(_solve_single_phase((tmin_curr, t),
                                                 t_eval=t_eval,
                                                 method=method,
                                                 atol=atol, rtol=rtol, 
                                                 events=events_non_sp_conc,
                                                 sp_conc_for_events=sp_conc_for_events, # dict or None
                                                 dense_output=dense_output,
                                                 y0=y0_curr))
                
                tmin_curr = t+dt_spike
                y0_curr = sols[-1].y[:, -1]
                dconcs_curr = dconcs
            
            # simulate last phase (last feed spike -> tmax)
            y0_curr += dconcs_curr(t=t, concs=y0_curr)\
                if callable(dconcs_curr)\
                else np_array(dconcs_curr)
            sols.append(_solve_single_phase((tmin_curr, tmax),
                                             t_eval=t_eval,
                                             method=method,
                                             atol=atol, rtol=rtol, 
                                             events=events_non_sp_conc,
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
                                             events=events_non_sp_conc,
                                             sp_conc_for_events=sp_conc_for_events, # dict or None
                                             dense_output=dense_output,
                                             y0=y0))
        
        y_final = np.concatenate([sol.y.transpose() for sol in sols])
        t_final = list(sols[0].t)
        
        t_events = None
        y_events = None
        if sols[0].t_events is not None: 
            t_events = list(sols[0].t_events[0])
            y_events = list(sols[0].y_events[0])
        else:
            t_events = []
            y_events = []
        
        for sol in sols[1:]:
            # for arr1, arr2 in zip([t_final, y_final],
            #                       [sol.t, sol.y]):
            #     t_final = np.concatenate((t_final, sol.t.transpose()))
            #     y_final = np.concatenate((y_final, sol.y.transpose()))
            t_final += list(sol.t)
            if sol.t_events is not None: 
                t_events += list(sol.t_events[0])
                y_events += list(sol.y_events[0])
            else:
                pass
        
        if sp_conc_for_events is None:
            sp_conc_for_events = {}
        
        if remove_negative_concs:
            y_final[np.where(y_final<0)] = 0.
        
        events = events if events is not None else []
        
        solution = {'t': np_array(t_final).transpose(),
                    'y': np_array(y_final).transpose(),
                    'events': events+['['+k+']' + ' = ' + str(v) 
                                      for k,v in sp_conc_for_events.items()],
                    't_events': np_array(t_events).transpose(),
                    'y_events': np_array(y_events).transpose(),
                    'sol': sols,
                    }
        
        self._solution = solution
        self._C_at_t_is_updated = False # this generates new interp1d objects the next time C_at_t is called
        self.save_solution(filename=filename, save_events=save_events) # if filename is None, saves only to RxnSys._solution_dfs and does not save file
        
        return solution
    
    def save_solution(self, filename, save_events=True):
        """
        Save the latest simulation solution to an Excel file and cache DataFrames.
        
        Parameters
        ----------
        filename : str
            Path to the output Excel file. If None, data is not saved to disk.
        save_events : bool, default True
            Whether to include event data in the output.
            
        """
        solution = self._solution
        df_dict = {'t': solution['t']}
        all_sp_IDs = self.species_system.all_sp_IDs
        y = solution['y']
        
        for ind, sp_ID in zip(range(len(all_sp_IDs)), 
                              all_sp_IDs):
            df_dict[sp_ID] = y[ind, :]
        
        df_main = pd.DataFrame.from_dict(df_dict)
        df_events = None
        if save_events and not solution['t_events'].shape == (0,):
            df_dict_events ={'event': solution['events'],
                            't_event': solution['t_events'],
                            }
            y_event = solution['y_events']
            for ind, sp_ID in zip(range(len(all_sp_IDs)), 
                                  all_sp_IDs):
                df_dict_events[sp_ID] = y_event[ind, :]
            
            df_events = pd.DataFrame.from_dict(df_dict_events)
                
        if filename is not None:
            if not '.xlsx' in filename:
                filename += '.xlsx'
            writer = pd.ExcelWriter(filename, engine = 'xlsxwriter')
            df_main.to_excel(writer, sheet_name='Main', index=False)
            if df_events is not None:
                df_events.to_excel(writer, sheet_name='Events', index=False)
            writer.close()
        
        self._solution_dfs = (df_main, df_events)

    def plot_solution(self, show_events=True, sps_to_include=None,
                  x_ticks=None, y_ticks=None):
        """
        Plot the concentrations of selected species over time.
        
        Parameters
        ----------
        show_events : bool, default True
            Whether to annotate event times on the plot.
        sps_to_include : list of str or Species, optional
            Species to include in the plot. If None, includes all.
        x_ticks : list or array-like, optional
            Tick values for the x-axis. If None, selected automatically.
        y_ticks : list or array-like, optional
            Tick values for the y-axis. If None, selected automatically.
        """
        if sps_to_include is None:
            sps_to_include = [i.ID for i in self.species_system.all_sps]
        
        sol = self._solution
        t, y = sol['t'], sol['y']
        t_events, y_events = sol['t_events'], sol['y_events']
        events = sol['events']
        all_sps = self.species_system.all_sps
        
        fig, ax = plt.subplots()
    
        y_max = 0  # for auto y-axis limit
    
        for i, sp in zip(range(len(all_sps)), all_sps):
            if sp in sps_to_include or sp.ID in sps_to_include:
                ax.plot(t, y[i, :], label=sp.ID,
                        linestyle='solid',
                        linewidth=1.)
                y_max = max(y_max, np.max(y[i, :]))
    
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Concentration [mol/L]')
    
        # === Automatic ticks and limits ===
        def auto_ticks(data_min, data_max, n_ticks=5, round_base=None):
            range_span = data_max - data_min
            raw_step = range_span / (n_ticks - 1)
            if round_base is None:
                exponent = np.floor(np.log10(raw_step))
                base = 10**exponent
                multiples = np.array([1, 2, 5, 10])
                step = multiples[np.searchsorted(multiples * base, raw_step, side='right')]
                step = step * base
            else:
                step = round_base
    
            tick_min = step * np.floor(data_min / step)
            tick_max = step * np.ceil(data_max / step)
            ticks = np.arange(tick_min, tick_max + 0.5 * step, step)
            return ticks, (tick_min, tick_max)
    
        # X-axis ticks and limits
        if x_ticks is None:
            x_ticks, xlim = auto_ticks(t.min(), t.max())
            ax.set_xticks(x_ticks)
            # ax.set_xlim(xlim)
        else:
            ax.set_xticks(x_ticks)
            # ax.set_xlim(min(x_ticks), max(x_ticks))
    
        # Y-axis ticks and limits
        if y_ticks is None:
            y_ticks, ylim = auto_ticks(0, y_max)
            ax.set_yticks(y_ticks)
            # ax.set_ylim(ylim)
        else:
            ax.set_yticks(y_ticks)
            # ax.set_ylim(min(y_ticks), max(y_ticks))
            
        # Tick directions
        ax.tick_params(axis='x', direction='inout', which='both', bottom=True, top=False)
        ax.tick_params(axis='x', direction='in', which='both', top=True)
        ax.tick_params(axis='y', direction='inout', which='both', left=True, right=False)
        ax.tick_params(axis='y', direction='in', which='both', right=True)
    
        # Minor ticks
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    
        if show_events:
            ylim = ax.get_ylim()
            for t_ev, e in zip(t_events, events):
                label = str(e)
                if callable(e):
                    label = label.split(' at ')[0].replace('<', '') + ' = 0'
                label = 't | ' + label
                ax.vlines(t_ev, ylim[0], ylim[1], 
                          linestyles='dashed', linewidth=0.5,
                          color='blue')
                xlim = ax.get_xlim()
                ax.annotate(label, 
                            xy=((t_ev + 0.04*(xlim[1]-xlim[0])), 
                                (ylim[1]-ylim[0])/2), 
                            verticalalignment='center', 
                            horizontalalignment='right', 
                            rotation=-270,
                            color='blue')
        
        plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
        plt.tight_layout()
        plt.show()
    
    def C_at_t(self, t,  species=None):
        """
        Get interpolated species concentrations at a given time.
        
        Parameters
        ----------
        t : float or list of floats
            Time or times at which to evaluate concentrations.
        species : Species, str or None
            Species object or ID of species to evaluate. 
            If None, returns all concentrations.

        Returns
        -------
        float or np.ndarray
            Concentration(s) at time t.
            
        """
        # Note on speed-up for C_at_t
        # Creating separate interp1d objects for each species concentration
        # results in faster C_at_t calls with species specified.
        # Creating a single interp1d object for all species concentrations
        # results in faster C_at_t calls with species not specified.
        # We create both here (if not self._C_at_t_is_updated; note this step
        # is slower as a result) and reference the faster object 
        # for the given call.
        
        if isinstance(t, list) or isinstance(t, np.ndarray):
            return np.array([self.C_at_t(ti) for ti in t])
        
        species_system = self.species_system
        
        if not self._C_at_t_is_updated:
            self._update_C_at_t()
        
        if species is not None:
            ind = species_system.index(species)
            return self._C_at_t_fs_indiv_sps[ind](t)
        else:
            return self._C_at_t_f_all(t)
    
    def _update_C_at_t(self):
        """
        Update the internal interpolation functions for species concentrations.
        This enables fast lookup of concentrations at arbitrary times.
        Not used for solve, but used for fit_reaction_kinetic_parameters_to_data.
        
        """
        species_system = self.species_system
        all_sps = species_system.all_sps
        index_f = species_system.index
        _solution = self._solution
        _t, _y = _solution['t'], _solution['y']
        interp1d_fill_value = self._interp1d_fill_value
        
        self._C_at_t_f_all = interp1d(_t, _y, fill_value=interp1d_fill_value)
        self._C_at_t_fs_indiv_sps = [interp1d(_t, _y[index_f(sp), :], fill_value=interp1d_fill_value) 
                                     for sp in all_sps]
        self._C_at_t_is_updated = True
        
    def _get_all_reactions_flattened(self):
        """
        Recursively flatten all reactions into a single list, including 
        nested ReactionSystem objects.
        
        Returns
        -------
        list
            List of all Reaction objects in the system.
        """
        reactions_flattened = []
        for r in self.reactions:
            if isinstance(r, Reaction):
                reactions_flattened.append(r)
            elif isinstance(r, ReactionSystem):
                reactions_flattened.extend(r._get_all_reactions_flattened())
        return reactions_flattened

    @property
    def reactions_flattened(self):
        """
        All reactions in the system, flattened to a single list.
        
        Returns
        -------
        list of Reaction
            All individual reactions, including those from nested ReactionSystem
            objects.
        
        """
        return self._get_all_reactions_flattened()
    
    def _get_reaction_kinetic_params(self):
        """
        Retrieve the kinetic parameters (kf and kb) for all reactions in the system.
        Excludes frozen parameters if self._exclude_frozen_params is True.
        
        Returns
        -------
        numpy.ndarray
            Flat array of kinetic parameters (kf and kb) for each reaction.
            Excludes frozen parameters if self._exclude_frozen_params is True.
            
        """
        rf = self.reactions_flattened
        param_vector = []
        param_keys = []
        
        if self._exclude_frozen_params:
            for r in rf:
                if not r._freeze_kf:
                    param_vector.append(r.kf)
                    param_keys.append((r.ID+': '+r.chem_equation.__str__(), 'kf'))
                if not r._freeze_kb:
                    param_vector.append(r.kb)
                    param_keys.append((r.ID+': '+r.chem_equation.__str__(), 'kb'))
        else:
            for r in rf:
                param_vector.extend([r.kf, r.kb])
                param_keys.append((r.ID+': '+r.chem_equation.__str__(), 'kf'))
                param_keys.append((r.ID+': '+r.chem_equation.__str__(), 'kb'))
        
        self._reaction_kinetic_param_keys = param_keys
        return np_array(param_vector)
    
    @property
    def reaction_kinetic_params(self):
        """
        Current kinetic parameter values for all reactions.
        Excludes frozen parameters if self._exclude_frozen_params is True.
        
        Returns
        -------
        numpy.ndarray
            Array of kinetic parameters (kf and kb) in the system.
            Excludes frozen parameters if self._exclude_frozen_params is True.
            
        """
        return self._get_reaction_kinetic_params()
    
    @property
    def reaction_kinetic_param_keys(self):
        """
        Keys corresponding to the kinetic parameters.
        Excludes frozen parameters if self._exclude_frozen_params is True.
        
        Returns
        -------
        list of str
            Parameter keys in the format 'reaction_ID.kf' or 'reaction_ID.kb'.
            Excludes frozen parameters if self._exclude_frozen_params is True.
            
        """
        self._get_reaction_kinetic_params()
        return self._reaction_kinetic_param_keys
    
    def set_reaction_kinetic_params(self, param_vector):
        """
        Set kinetic parameters (kf, kb) for all reactions from a flat vector.
        
        Parameters
        ----------
        param_vector : array_like
            Flat list or array of parameters.
            
        """
        
        rf = self.reactions_flattened
        # expected_length = 2 * len(rf)
        # if len(param_vector) != expected_length:
        #     raise ValueError(f"Expected vector of length {expected_length}, got {len(param_vector)}.\n")
        # for i, r in enumerate(rf):
        #     r.kf = param_vector[2*i]
        #     r.kb = param_vector[2*i + 1]
        
        curr_param_ind = 0
        for r in rf:
            if not r._freeze_kf:
                r.kf = param_vector[curr_param_ind]
                curr_param_ind += 1
            if not r._freeze_kb:
                r.kb = param_vector[curr_param_ind]
                curr_param_ind += 1
                
    def _extract_t_spIDs_y(self,
                           data: Union[str, dict, pd.DataFrame]) -> Tuple[pd.Series, pd.DataFrame, List[str]]:
        """
        Extracts 't' values, species names, and species data from the input.
        
        Parameters:
        -----------
            data: A pandas DataFrame, a dictionary, or a path to a .csv or .xlsx file.
        
        Returns:
        --------
            t: List of time values
            species_IDs: List of species IDs in column names (excluding 't')
            _y: List of lists, each corresponding to a species column
        
        """
        df = None
        # Load data into a pandas DataFrame
        if isinstance(data, str):
            if data.endswith('.csv'):
                df = pd.read_csv(data)
            elif data.endswith('.xlsx'):
                df = pd.read_excel(data)
            else:
                raise ValueError(f"Unsupported file type: {data}")
        elif isinstance(data, dict):
            df = pd.DataFrame(data)
        elif isinstance(data, pd.DataFrame):
            # df = data.copy()
            df = data
        else:
            raise TypeError("Input data must be a DataFrame, dict, or path to a .csv or .xlsx file.")
            
        if 't' not in df.columns:
            raise ValueError("The input data must contain a 't' column.")
            
        t = np_array(df['t'].tolist())
        
        species_IDs = [col for col in df.columns if col != 't']
        all_sp_IDs = self.species_system.all_sp_IDs
        for sp_ID in species_IDs:
            if not sp_ID in all_sp_IDs:
                raise ValueError(f"A Species ID '{sp_ID}' obtained from input data was not found in all_sp_IDs:\n{all_sp_IDs}\n")
            
        _y = np_array([df[name].tolist() for name in species_IDs])
        
        return t, species_IDs, _y

    def fit_reaction_kinetic_parameters_to_data(self,
                                                data,
                                                p0=None,
                                                all_species_tracked=False,
                                                use_only=None,
                                                normalize=True,
                                                n_de_runs=5,
                                                n_minimize_runs=2,
                                                minimize_kwargs=None,
                                                differential_evolution_kwargs=None,
                                                show_output=True,
                                                show_progress=False,
                                                plot_during_fit=False,
                                                call_before_each_solve=None,
                                                timeout_solve_ivp=0.5,
                                                **kwargs):
        """
        Fit reaction kinetic parameters to experimental time-series data.
        (i.e., inverse modeling).
        
        This method optimizes the reaction kinetic parameters of the system to best 
        fit provided concentration data over time using nonlinear least squares optimization.
        The function assumes a shared set of parameters governing all species.
        
        Parameters
        ----------
        data : pandas.DataFrame, dict, str, or list
            A pandas DataFrame, a dictionary, a path to a .csv or .xlsx file,
            or a list containing any combinations of those items.
            Each item represents data containing time ('t') and concentrations 
            of one or more species as columns. Time should be in the column labeled 't'; 
            all other columns are interpreted as species concentrations.
            If data is a list of data items, each data item must have the same format
            (i.e., the same columns and order of columns).
            
        p0 : array_like, optional
            Initial guess for the kinetic parameters. If not provided, the current 
            `reaction_kinetic_params` will be used if they are finite and not NaN.
            
        all_species_tracked : bool, default False
            If True, assumes all species are tracked and uses the full concentration vector
            when computing predictions. Otherwise, uses only individual species.
            
        show_output : bool, default True
            If True, prints fit summary including R² score and success status.
            
        use_only : list of str or Species, optional
            List of Species object or species IDs to fit against. If not provided,
            all species in the input data will be used. It is recommended to include
            at least one species involved in each reaction in `reactions`.
            
        method : str, default 'Powell'
            Optimization method to use for parameter fitting. Passed to 
            `scipy.optimize.minimize`.
        
        n_minimize_runs : int, optional
            Number of local optimization (minimization) runs to perform, each with a different
            starting point (first uses `p0`, others use DE output).
        
        n_de_runs : int, optional
            Number of differential evolution runs per minimization run to generate good starting
            guesses; the best DE result is used to initialize the local minimization.
        
        minimize_kwargs : dict
            Additional keyword arguments passed to `scipy.optimize.minimize`.
        
        differential_evolution_kwargs : dict
            Additional keyword arguments passed to `scipy.optimize.differential_evolution_kwargs`.
        
        timeout_solve_ivp: float, int, or None, optional
            Enforce timeout of `scipy.integrate.solve_ivp` when it exceeds this value.
            Creates and passes a timeout_function as an event to `scipy.integrate.solve_ivp`.
            Defaults to 0.5 (seconds).
            
        Returns
        -------
        None
            Updates the reaction kinetic parameters of the system in place and stores 
            the fit result in `self._fitsol` as a tuple: (best-fit parameters, R² score, success flag).
            
        Notes
        -----
        - The fit is performed by normalizing each species' concentration to its maximum
          observed value to ensure scale invariance.
        - The objective function maximizes the mean R² score across selected species.
        - This method relies on `solve()` and `C_at_t()` methods to simulate the system 
          under trial parameters and look up interpolated results, 
          and on `fit_multiple_dependent_variables()` to perform optimization.
          
        See Also
        --------
        `scipy.optimize.minimize <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`_
        `scipy.optimize.differential_evolution <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html>`_
        
        """
        _orig_timeout_solve_ivp = self._timeout_solve_ivp # reset to original value after fit, usually None
        self._timeout_solve_ivp = timeout_solve_ivp
        
        sp_sys = self.species_system
        all_sp_IDs = sp_sys.all_sp_IDs
        
        if call_before_each_solve is None:
            call_before_each_solve = []
            
        use_only_inds = [sp_sys.index(sp) for sp in use_only]
        
        if not (isinstance(data, list) or isinstance(data, tuple)):
            data = [data]
        
        dataset = [self._extract_t_spIDs_y(di)
                              for di in data]
            
        t_dataset, sp_IDs_dataset, y_dataset = [], [], []
        
        for d in dataset:
            t_dataset.append(d[0])
            sp_IDs_dataset.append(d[1])
            y_dataset.append(d[2])
            
        data_sp_IDs = sp_IDs_dataset[0]
        sp_IDs_to_use = use_only if use_only is not None else data_sp_IDs
        sp_inds = [sp_sys.index(sp) for sp in sp_IDs_to_use]
        
        structured_xdata = []
        
        
        for t_, sp_IDs, y_ in zip(t_dataset, sp_IDs_dataset, y_dataset):
            
            tdata = t_
            
            y_maxes = np.array([np.max(y_[ind, :]) for ind in sp_inds])
            y_maxes = np_array([y_maxes for i in range(y_.shape[1])]).transpose()
            
            # y_normalized = y_
            # if use_only:
            #     use_only_inds = [sp_sys.index(sp) for sp in use_only]
            #     y_normalized = y_normalized[use_only_inds, :]
            # y_normalized /= y_maxes
            
            y0 = np.zeros(len(all_sp_IDs))
            
            for spID in data_sp_IDs:
                # get initial concentrations of all species in the dataset,
                # even if using only those in use_only for actual fitting
                y0[sp_sys.index(spID)] = y_[data_sp_IDs.index(spID), 0]
            
            structured_xdata.append((tdata, y_maxes, y0))
        
        if p0 is None:
            rkp = self.reaction_kinetic_params
            if not np.any(np.isinf(rkp)) or np.any(np.isnan(rkp)):
                p0 = rkp
        
        
        
        # y_normalized = [yi_/sdata[1] for (yi_, sdata) in zip(y_dataset, structured_xdata)]
        # y_normalized = y_normalized.flatten()
        y_to_use_0 = y_dataset[0]
        if use_only:
            y_to_use_0 = y_to_use_0[use_only_inds, :]
        
        y_dataset_normalized = y_to_use_0/structured_xdata[0][1] if normalize else y_to_use_0
        for (yi_, sdata) in zip(y_dataset[1:], structured_xdata[1:]):
            y_to_use = yi_
            if use_only:
                y_to_use = y_to_use[use_only_inds, :]
            to_concat = y_to_use/sdata[1] if normalize else y_to_use
            y_dataset_normalized = np.concatenate((y_dataset_normalized.transpose(), to_concat.transpose()))
            y_dataset_normalized = y_dataset_normalized.transpose()

        set_rxn_kp = self.set_reaction_kinetic_params
        solve = self.solve
        _update_C_at_t = self._update_C_at_t
        plot_solution = self.plot_solution
        
        # zeros_like_y_dataset_normalized = np.zeros(shape=y_dataset_normalized.shape)
        # nans_like_y_dataset_normalized = np.full(y_dataset_normalized.shape, np.nan)
        
        def f_single_xdata(xdata, new_rxn_kp):
            tdata, y_maxes, y0 = xdata
            t_span = np.min(tdata), np.max(tdata)
            set_rxn_kp(new_rxn_kp)
            for c in call_before_each_solve:
                c(new_rxn_kp)
            solve(t_span=t_span, y0=y0, save_events=False)
            for sol in self._solution['sol']:
                if (not sol.success) or (sol.status==1):
                    return np.full((len(y_maxes), len(tdata)), 0.)
            if plot_during_fit: plot_solution()
            _update_C_at_t()
            if not all_species_tracked:
                ypred = np_array([self._C_at_t_fs_indiv_sps[ind](tdata)
                        for ind in sp_inds])
                if normalize: ypred/= y_maxes
                return ypred
            else:
                ypred = self._C_at_t_f_all(tdata)
                if normalize: ypred/= y_maxes
                return ypred
            
        def f(xdataset, new_rxn_kp):
            ypred_normalized_concat = f_single_xdata(xdataset[0], new_rxn_kp)
            for xdata in xdataset[1:]:
                ypred_normalized_concat = np.concatenate((ypred_normalized_concat.transpose(), f_single_xdata(xdata, new_rxn_kp).transpose()))
                ypred_normalized_concat = ypred_normalized_concat.transpose()
            return ypred_normalized_concat
        
        if minimize_kwargs is None:
            minimize_kwargs = {'method': 'Powell',
                               'bounds':[(0., None) for i in p0],
                               }
        if differential_evolution_kwargs is None:
            differential_evolution_kwargs = {'strategy': 'best1bin',
                                             'bounds': [(0, 1e8) for i in p0]}
        elif 'bounds' not in differential_evolution_kwargs:
            differential_evolution_kwargs['bounds'] = [(0, 1e8) for i in p0]
        fitsol = fit_multiple_dependent_variables(f=f,
                                                   xdata=structured_xdata,
                                                   ydata=y_dataset_normalized,
                                                   p0=p0, # 
                                                   minimize_kwargs=minimize_kwargs,
                                                   differential_evolution_kwargs=differential_evolution_kwargs,
                                                   n_de_runs=n_de_runs,
                                                   n_minimize_runs=n_minimize_runs,
                                                   show_progress=show_progress,
                                                   # options={'maxiter':5},
                                                   )
        
        set_rxn_kp(fitsol[0])
        self._fitsol = fitsol
        
        self._timeout_solve_ivp = _orig_timeout_solve_ivp # reset to original, usually None
        
        if show_output: 
            print('\n')
            print('Fit results')
            print('-----------\n')
            print(f'R^2={fitsol[1]}\n')
            print(self.__str__())
            print('\n\n')
            
    def add_reaction(self, reaction):
        """
        Add a reaction to the system.
        
        Parameters
        ----------
        reaction : Reaction, ReactionSystem, or str
            The reaction object or string to add.
            
        """
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
        """
        Modify a reaction at a given index.
        
        Parameters
        ----------
        index : int
            Index of the reaction to modify.
        new_equation_string : str, optional
            New chemical equation.
        new_kf : float, optional
            New forward rate constant.
        new_kb : float, optional
            New backward rate constant.
            
        """
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
        """
        Convert dictionary of concentration spikes into a list of time-change pairs.
        
        Parameters
        ----------
        spikes_dct : dict
            Dictionary mapping time to spike definition(s).
            
        Returns
        -------
        list
            List of (time, change) tuples.
        
        """
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
                if is_array_of_numbers(np_array(dconcs)):
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
        """
        Parse a string description of concentration change and return a change function or array.
        
        Parameters
        ----------
        s : str
            String of the form 'Change; Species ID; value' 
            or 'Target; Species ID; value'
            or 'Target (negative change allowed); Species ID; value'.
            
        Returns
        -------
        array or function
            Array of concentration changes or function returning such an array.
        
        """
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
        """
        Return a string representation of the reaction system.
        
        Returns
        -------
        str
            String listing the reaction system and its reactions.
            
        """
        rxns = self.reactions
        str_ = f'{self.ID}: ReactionSystem(\n'
        for r in rxns:
            str_ += '    ' + r.__str__() + '\n'
        str_ += '    ' + ')'
        return str_
    
    def __repr__(self):
        """
        Return the formal string representation of the object.
        
        Returns
        -------
        str
            Same as __str__.
            
        """
        return self.__str__()
    
RxnSys = ReactionSystem
