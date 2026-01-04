# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
References
----------
.. [1] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden
    National. Renewable Energy Laboratory Golden, Colorado. P. Schoen,
    J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon. Harris Group
    Inc. Seattle, Washington and Atlanta, Georgia. Process Design and Economics
    for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid
    Pretreatment and Enzymatic Hydrolysis of Corn Stover. May 2011. Technical
    Report NREL/TP-5100-47764
"""

import numpy as np
from biosteam.units import BatchBioreactor
from thermosteam import Reaction, ParallelReaction
from ..reaction_systems.tellurium_based.tellurium_sbml import TelluriumReactionSystem
from ..reaction_systems.reaction_system import ReactionSystem
from ..utils import get_index_nearest_element_from_sorted_array

__all__ = ('NSKFermentation',)


class NSKFermentation(BatchBioreactor):
    """
    Create a unit that models large-scale batch fermentation
    for the production of 1st-generation ethanol using Saccharomyces cerevisiae.
    Only sucrose and glucose are taken into account for conversion.
    Conversion is based on reaction time, `tau`. Cleaning and unloading time,
    `tau_0`, fraction of working volume, `V_wf`, and number of reactors,
    `N_reactors`, are attributes that can be changed. Cost of a reactor
    is based on the NREL batch fermentation tank cost assuming volumetric
    scaling with a 6/10th exponent [1]_. 
    
    Parameters
    ----------
    ins : 
        Inlet fluids to be mixed into the fermentor.
    outs : 
        * [0] Vent
        * [1] Effluent
    tau : float
        Reaction time [h].
    N : int, optional
        Number of batch reactors
    V : float, optional
        Target volume of reactors [m^3].
    T=305.15 : float
        Temperature of reactor [K].
    P=101325 : float
        Operating pressure of reactor [Pa].
    Nmin=2 : int
        Minimum number of fermentors.
    Nmax=36: int
        Maximum number of fermentors.  
    kinetic_reaction_system: nskinetics.ReactionSystem or nskinetics.TelluriumReactionSystem object.
        Reaction system used for kinetic simulation.
    map_chemicals_nsk_to_bst: dict.
        Dictionary with keys as nskinetics species ID strings and values as biosteam chemical ID strings.
    Notes
    -----
    Either N or V must be given.
    
    Examples
    --------
    
    
    
    """
    line = 'NSKFermentation'
    _ins_size_is_fixed = False
    # Ins size is not fixed; however, the first three ins streams must be as follows:
    # 0. Initial feed
    # 1. Seed culture
    # 2. Spike feed (optional; for fed-batch only)
    
    autoselect_N = True
    
    def _init(self, 
              tau, 
              kinetic_reaction_system, 
              map_chemicals_nsk_to_bst={}, 
              n_simulation_steps=1000,
              f_reset_kinetic_reaction_system=None,
              N=None, V=None, T=305.15, P=101325., Nmin=2, Nmax=36,
              sugar_IDs=('Sucrose', 'Glucose', 'Xylose'),
              tau_max=24.*7.,
              tau_update_policy=('max', 's_IBO'),
              perform_hydrolysis=True):
        
        BatchBioreactor._init(self, tau=tau, N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax)
        self._load_components()
        
        chemicals = self.chemicals
        if perform_hydrolysis:
            self.hydrolysis_reaction = Reaction('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00, chemicals)
        
        self.kinetic_reaction_system = kinetic_reaction_system
        if isinstance(kinetic_reaction_system, TelluriumReactionSystem):
            self.simulate_kinetics = self._nsk_te_simulate_kinetics
        elif isinstance(kinetic_reaction_system, ReactionSystem):
            self.simulate_kinetics = self._nsk_simulate
        
        self.map_chemicals_nsk_to_bst = map_chemicals_nsk_to_bst
        self.n_simulation_steps = n_simulation_steps
        self.f_reset_kinetic_reaction_system = f_reset_kinetic_reaction_system if f_reset_kinetic_reaction_system is not None else lambda model: model.reset()
        
        self.sugar_IDs = sugar_IDs
        self.perform_hydrolysis = perform_hydrolysis
        
        self.tau_max = tau_max
        self.tau_update_policy = tau_update_policy
        
    def _nsk_simulate_kinetics(self, feed, tau, feed_spike_condition=None, plot=False): 
        # !!!
        self.tau = tau
        effluent = feed.copy()
        return effluent
    
    def _nsk_te_simulate_kinetics(self, feed, tau, feed_spike_condition=None, plot=False):
        # self.tau = tau
        map_chemicals_nsk_to_bst = self.map_chemicals_nsk_to_bst
        kinetic_reaction_system = self.kinetic_reaction_system
        self.f_reset_kinetic_reaction_system(kinetic_reaction_system)
        te_r = kinetic_reaction_system._te
        chems_nsk = list(map_chemicals_nsk_to_bst.keys())
        # get unit conversion factors and unit-based material indexers
        
        time_units = kinetic_reaction_system._units['time']
        if time_units.lower() in ('min', 'm'):
            time_conv_factor = 60.0
        elif time_units.lower() in ('sec', 's'):
            time_conv_factor = 3600.0
        elif time_units.lower() in ('hr', 'h'):
            time_conv_factor = 1.0
            
        conc_units = kinetic_reaction_system._units['conc']
        if conc_units in ('M', 'mol/L', 'kg/m3', 'kg/m^3'):
            material_indexer = 'imol'
            volume_attribute = "ivol['Water']"
        elif conc_units in ('g/L', 'kg/m3', 'kg/m^3'):
            material_indexer = 'imass'
            volume_attribute = "ivol['Water']"
            
        self._nsk_initial_concentration = initial_concentrations = {}
        for c_nsk, c_bst in map_chemicals_nsk_to_bst.items():
            exec(f'te_r.{c_nsk.replace("[", "").replace("]", "")} = feed.{material_indexer}[c_bst]/feed.{volume_attribute}')
            exec(f'initial_concentrations[c_nsk] = te_r.{c_nsk.replace("[", "").replace("]", "")}')
        
        initial_conc_sugars = te_r.s_glu
        
        results_col_names = ['time']+chems_nsk
        try:
            self.results = results = te_r.simulate(0, self.tau_max*time_conv_factor, self.n_simulation_steps,
                                    results_col_names,
                                    )
        except Exception as e:
            print(str(e))
            breakpoint()
            
        tau_index = -1
        tau_update_policy = self.tau_update_policy
        
        if tau_update_policy is None:
            tau_index = get_index_nearest_element_from_sorted_array(results[:, results_col_names.index('time')])
            
        elif tau_update_policy[0]=='max':
            chem_to_max = tau_update_policy[1]
            index_chem_to_max = results_col_names.index(chem_to_max)
            results = np.array(results)
            index_tau_with_max_chem = np.where(results[:, index_chem_to_max] == results[:, index_chem_to_max].max())
            tau_index = index_tau_with_max_chem
            
        self.results_specific_tau = results_specific_tau = results[tau_index][0]
        
        effluent = feed.copy()
        for c_nsk, c_bst in map_chemicals_nsk_to_bst.items():
            # exec(f'effluent.{material_indexer}[c_bst] = te_r.{c_nsk.replace("[", "").replace("]", "")} * effluent.{volume_attribute}')
            # breakpoint()
            exec(f'effluent.{material_indexer}[c_bst] = results_specific_tau[results_col_names.index(c_nsk)] * effluent.{volume_attribute}')

        self.amt_sugars_spiked = amt_sugars_spiked = te_r.conc_glu_feed_spike * te_r.tot_vol_glu_feed_added
        self.amt_sugars_initial = amt_sugars_initial = initial_conc_sugars*(te_r.env - te_r.tot_vol_glu_feed_added)
        
        effluent.F_vol *= te_r.env/(te_r.env - te_r.tot_vol_glu_feed_added)
        
        # effluent.F_vol *= (amt_sugars_initial + amt_sugars_spiked)/amt_sugars_initial
        
        self.tau = results_specific_tau[results_col_names.index('time')]
        
        return effluent
    
    def _run(self):
        vent, effluent = self.outs
        ins = self.ins
        spike_feed = ins[2]
        initial_feed_seed_others = (i for i in ins if not i==spike_feed) # exclude spike feed initially
        effluent.mix_from(initial_feed_seed_others)
        
        if self.perform_hydrolysis:
            self.hydrolysis_reaction.force_reaction(effluent)
            self.hydrolysis_reaction.force_reaction(spike_feed)
        
        minimal_feed = effluent.copy()
        
        for i in minimal_feed.chemicals:
            if not i.ID in list(self.map_chemicals_nsk_to_bst.values()) + ['Water',]:
                minimal_feed.imol[i.ID] = 0.0
        minimal_effluent = self.simulate_kinetics(feed=minimal_feed, tau=self._tau)
        
        te_r = self.kinetic_reaction_system._te
        
        # effluent.F_vol *= (te_r.env)/(te_r.env-te_r.tot_vol_glu_feed_added)
        # effluent.F_vol *= (effluent.F_vol + spike_feed.F_vol)/effluent.F_vol
        
        # effluent.mix_from((effluent, spike_feed))
        # copy non-nskinetics chemicals to minimal_effluent from effluent and spike_feed
        for i in minimal_effluent.chemicals:
            if not i.ID in list(self.map_chemicals_nsk_to_bst.values()) + ['Water',]:
                minimal_effluent.imol[i.ID] += effluent.imol[i.ID]
                minimal_effluent.imol[i.ID] += spike_feed.imol[i.ID]
        
        effluent.copy_like(minimal_effluent)
        effluent.imol['NH3'] = 0. # NH3 in ins must be based on final Yeast mass
        
        effluent.empty_negative_flows()
        vent.empty()
        vent.receive_vent(effluent, energy_balance=False)

    # @property
    # def Hnet(self):
    #     return 0