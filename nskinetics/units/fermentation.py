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
    
    autoselect_N = True
    
    def _init(self, 
              tau, 
              kinetic_reaction_system, 
              map_chemicals_nsk_to_bst={}, 
              n_simulation_steps=1000,
              f_reset_kinetic_reaction_system=None,
              N=None, V=None, T=305.15, P=101325., Nmin=2, Nmax=36,
              efficiency=None):
        
        BatchBioreactor._init(self, tau=tau, N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax)
        self._load_components()
        
        chemicals = self.chemicals
        self.hydrolysis_reaction = Reaction('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00, chemicals)
        
        self.kinetic_reaction_system = kinetic_reaction_system
        if isinstance(kinetic_reaction_system, TelluriumReactionSystem):
            self.simulate_kinetics = self._nsk_te_simulate_kinetics
        elif isinstance(kinetic_reaction_system, ReactionSystem):
            self.simulate_kinetics = self._nsk_simulate
        
        self.map_chemicals_nsk_to_bst = map_chemicals_nsk_to_bst
        self.n_simulation_steps = n_simulation_steps
        self.f_reset_kinetic_reaction_system = f_reset_kinetic_reaction_system if f_reset_kinetic_reaction_system is not None else lambda model: model.reset()
        
    def _nsk_simulate_kinetics(self, feed, tau, feed_spike_condition=None, plot=False): 
        # !!!
        self.tau = tau
        effluent = feed.copy()
        return effluent
    
    def _nsk_te_simulate_kinetics(self, feed, tau, feed_spike_condition=None, plot=False):
        self.tau = tau
        map_chemicals_nsk_to_bst = self.map_chemicals_nsk_to_bst
        kinetic_reaction_system = self.kinetic_reaction_system
        self.f_reset_kinetic_reaction_system(kinetic_reaction_system)
        te_r = kinetic_reaction_system._te
        
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
            volume_indexer = 'F_vol'
        elif conc_units in ('g/L', 'kg/m3', 'kg/m^3'):
            material_indexer = 'imass'
            volume_indexer = 'F_vol'
            
        self._nsk_initial_concentration = initial_concentrations = {}
        for c_nsk, c_bst in map_chemicals_nsk_to_bst.items():
            exec(f'te_r.{c_nsk} = feed.{material_indexer}[c_bst]/feed.{volume_indexer}')
            exec(f'initial_concentrations[c_nsk] = te_r.{c_nsk}')
        
        
        te_r.simulate(0, tau*time_conv_factor, self.n_simulation_steps)
        
        effluent = feed.copy()
        initially_zero = []
        initially_nonzero = []
        # breakpoint()
        for c_nsk, c_bst in map_chemicals_nsk_to_bst.items():
            if initial_concentrations[c_nsk] > 0.0:
                exec(f'effluent.{material_indexer}[c_bst] *= te_r.{c_nsk}/initial_concentrations[c_nsk]')
                initially_nonzero.append((c_nsk, c_bst))
            else:
                initially_zero.append((c_nsk, c_bst))
        
        for c_nsk, c_bst in initially_zero:
            exec(f'material_factor = feed.{material_indexer}[initially_nonzero[0][1]]/initial_concentrations[initially_nonzero[0][0]]; effluent.{material_indexer}[c_bst] = material_factor * te_r.{c_nsk}')
                
        return effluent
    
    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins)
        self.hydrolysis_reaction.force_reaction(effluent)
        effluent.copy_like(self.simulate_kinetics(feed=effluent, tau=self._tau))
        effluent.empty_negative_flows()
        vent.empty()
        vent.receive_vent(effluent, energy_balance=False)
