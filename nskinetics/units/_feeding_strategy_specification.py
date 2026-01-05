# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import nskinetics as nsk
import biosteam as bst
import thermosteam as tmo
from flexsolve import IQ_interpolation

__all__ = ('FedBatchStrategySpecification',)

#%%
class FedBatchStrategySpecification:
    """
    Specification and orchestration of a fed-batch feeding strategy based on
    sugar concentration control for fermentation processes.
    
    This class defines target and threshold sugar concentrations and translates
    them into operating conditions for upstream unit operations (evaporators,
    mixers, splitter) that supply both a base feed and a spike feed to a
    fermentation reactor. The strategy supports dynamic reconfiguration of
    feed concentrations, spike magnitude, and reactor residence time, and is
    designed to work with multiple sugar species (e.g., glucose, sucrose,
    xylose).
    
    The specification logic includes:
    - Setting desired sugar concentrations for continuous feed and spike feed
      via evaporator volume and mixer dilution adjustments.
    - Defining a threshold sugar concentration that triggers spike feeding.
    - Updating maximum reactor residence time (tau_max) and propagating its effects
      upstream through splitter and feed unit simulations.
    - Coordinated simulation of all upstream units to ensure consistency
      between specified targets and achievable process conditions.
    
    Parameters
    ----------
    target_conc_sugars : float
        Desired target concentration of total sugars in the continuous feed
        to the fermentation reactor.
    threshold_conc_sugars : float
        Sugar concentration threshold in the reactor environment that triggers
        a spike feeding event.
    conc_sugars_feed_spike : float
        Target sugar concentration of the spike feed stream used during
        threshold-triggered feeding.
    tau_max : float
        Maximum residence time of the fermentation reactor.
    fermentation_reactor : object
        Fermentation reactor unit, expected to expose a kinetic reaction system
        and support simulation with spike feeding logic.
    splitter : object
        Splitter unit used to divide flow between base feed and spike feed
        pathways.
    feed_evaporator : object
        Evaporator used to concentrate sugars in the continuous feed stream.
    feed_mixer : object
        Mixer used to dilute or adjust the continuous feed stream composition.
    spike_evaporator : object
        Evaporator used to concentrate sugars in the spike feed stream.
    spike_mixer : object
        Mixer used to dilute or adjust the spike feed stream composition.
    sugar_IDs : list of str, optional
        Identifiers of sugar species used to compute total sugar concentration
        (default: ['Glucose', 'Sucrose', 'Xylose']).
    """
    
    def __init__(
        self,
        # kinetic_reaction_system: nsk.TelluriumReactionSystem,
        target_conc_sugars: float,
        threshold_conc_sugars: float,
        conc_sugars_feed_spike: float,
        tau_max: float,
        fermentation_reactor,
        splitter,
        feed_evaporator,
        feed_mixer,
        spike_evaporator,
        spike_mixer,
        feed_units_sequential=None,
        spike_units_sequential=None,
        sugar_IDs=['Glucose', 'Sucrose', 'Xylose']
        ):
        # self.kinetic_reaction_system = kinetic_reaction_system
        self.target_conc_sugars = target_conc_sugars
        self.threshold_conc_sugars = threshold_conc_sugars
        self.conc_sugars_feed_spike = conc_sugars_feed_spike
        self.sugar_IDs = sugar_IDs
        self.tau_max = tau_max
        
        self.fermentation_reactor = fermentation_reactor
        self.splitter = splitter
        self.feed_evaporator = feed_evaporator
        self.feed_mixer = feed_mixer
        self.feed_units_sequential = feed_units_sequential
        self.spike_evaporator = spike_evaporator
        self.spike_mixer = spike_mixer
        self.spike_units_sequential = spike_units_sequential
        
        self._validate_parameters()
        
    def _validate_parameters(self) -> None:
        """
        Basic validation for initialization parameters.
        """
        if self.target_conc_sugars <= 0:
            raise ValueError("target_conc_sugars must be positive.")
            
        if self.threshold_conc_sugars < 0:
            raise ValueError("threshold_conc_sugars must be positive.")
            
        if self.conc_sugars_feed_spike <= 0:
            raise ValueError("conc_sugars_feed_spike must be positive.")
            
        if self.threshold_conc_sugars > self.target_conc_sugars:
            raise ValueError("threshold_conc_sugars should not exceed target_conc_sugars.")
        
        if self.target_conc_sugars > self.conc_sugars_feed_spike:
            raise ValueError("target_conc_sugars should not exceed conc_sugars_feed_spike.")
        
    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"target_conc_sugars={self.target_conc_sugars}, "
            f"threshold_conc_sugars={self.threshold_conc_sugars}, "
            f"conc_sugars_feed_spike={self.conc_sugars_feed_spike})"
        )
    
    def load_specifications(self,
                            target_conc_sugars, 
                            conc_sugars_feed_spike,
                            threshold_conc_sugars,
                            tau_max,
                            evaporator_V_ub=0.8, evaporator_V_lb=0.0,
                            mixer_dil_lb=0., mixer_dil_ub=100_000
                            ):
        self.load_desired_concs_sugars(target_conc_sugars=target_conc_sugars, 
            conc_sugars_feed_spike=conc_sugars_feed_spike,
            evaporator_V_ub=evaporator_V_ub, evaporator_V_lb=evaporator_V_lb,
            mixer_dil_lb=mixer_dil_lb, mixer_dil_ub=mixer_dil_ub)
        
        self.load_threshold_conc_sugars_and_tau_max(threshold_conc_sugars=threshold_conc_sugars,
                                                tau_max=tau_max)
        
    def load_desired_concs_sugars(self, 
                                target_conc_sugars, 
                                conc_sugars_feed_spike,
                                evaporator_V_ub, evaporator_V_lb,
                                mixer_dil_lb, mixer_dil_ub):
        # clear_units([V301, K301])
        self.target_conc_sugars = target_conc_sugars
        self.conc_sugars_feed_spike = conc_sugars_feed_spike
        
        te_r = self.fermentation_reactor.kinetic_reaction_system._te
        te_r.conc_glu_feed_spike = conc_sugars_feed_spike
        te_r.target_conc_glu_spike = target_conc_sugars
        
        feed_evaporator = self.feed_evaporator
        feed_mixer = self.feed_mixer
        spike_evaporator = self.spike_evaporator
        spike_mixer = self.spike_mixer
        
        fermentation_reactor = self.fermentation_reactor
        
        self._estimate_set_V_and_dil(desired_conc_sugars=target_conc_sugars,
                               evaporator=feed_evaporator,
                               mixer=feed_mixer,
                               evaporator_V_lb=evaporator_V_lb, evaporator_V_ub=evaporator_V_ub,
                               mixer_dil_lb=mixer_dil_lb, mixer_dil_ub=mixer_dil_ub,
                               f_simulate_units_sequential=self._simulate_feed_units,
                               terminal_stream=fermentation_reactor.ins[0])
        
        self._estimate_set_V_and_dil(desired_conc_sugars=conc_sugars_feed_spike,
                               evaporator=spike_evaporator,
                               mixer=spike_mixer,
                               evaporator_V_lb=evaporator_V_lb, evaporator_V_ub=evaporator_V_ub,
                               mixer_dil_lb=mixer_dil_lb, mixer_dil_ub=mixer_dil_ub,
                               f_simulate_units_sequential=self._simulate_spike_units,
                               terminal_stream=fermentation_reactor.ins[2])
    
    def load_threshold_conc_sugars_and_tau_max(self,
                                           threshold_conc_sugars,
                                           tau_max):
        self.threshold_conc_sugars = threshold_conc_sugars
        fermentation_reactor = self.fermentation_reactor
        te_r = fermentation_reactor.kinetic_reaction_system._te
        
        fermentation_reactor.tau_max = tau_max
        te_r.threshold_conc_glu_spike = threshold_conc_sugars
        
        self._simulate_upstream_units()
        fermentation_reactor.simulate()
        
        final_env_vol = te_r.env
        vol_spike_added = te_r.tot_vol_glu_feed_added
        initial_env_vol = final_env_vol - vol_spike_added
        sugars_in_initial_feed = initial_env_vol * self.target_conc_sugars
        sugars_in_spikes = vol_spike_added * self.conc_sugars_feed_spike
        self.splitter.split = sugars_in_initial_feed/(sugars_in_initial_feed+sugars_in_spikes) # split to initial feed
        
        # self._simulate_upstream_units()
        # fermentation_reactor.simulate()
        
    def _estimate_set_V_and_dil(self,
                               desired_conc_sugars,
                               evaporator, mixer, 
                               evaporator_V_lb, evaporator_V_ub,
                               mixer_dil_lb, mixer_dil_ub,
                               f_simulate_units_sequential,
                               terminal_stream):
        
        get_conc_sugars = self.get_conc_sugars
        f_simulate_units_sequential()
        # breakpoint()
        
        if mixer.outs[0].F_vol:
            def _evaporator_obj_f(V):
                evaporator.V = V
                f_simulate_units_sequential()
                return get_conc_sugars(terminal_stream) - desired_conc_sugars
                    
            def _mixer_obj_f(water_to_sugar_mol_ratio):
                mixer.water_to_sugar_mol_ratio = water_to_sugar_mol_ratio
                f_simulate_units_sequential()
                return get_conc_sugars(terminal_stream) - desired_conc_sugars
            
            _evaporator_obj_f(evaporator_V_lb)
            
            if _mixer_obj_f(mixer_dil_lb) < 0: # if there is too low a conc even with no dilution
                try:
                    IQ_interpolation(_evaporator_obj_f, evaporator_V_lb, evaporator_V_ub, ytol=1e-3)
                except:
                    breakpoint()
            elif _mixer_obj_f(mixer_dil_ub) > 0:
                try:
                    IQ_interpolation(_mixer_obj_f, mixer_dil_lb, mixer_dil_ub, ytol=1e-3)
                except:
                    breakpoint()
            else:
                _evaporator_obj_f(evaporator_V_lb)
                try:
                    IQ_interpolation(_mixer_obj_f, mixer_dil_lb, mixer_dil_ub, ytol=1e-3)
                except:
                    breakpoint()
    
    def _simulate_feed_units(self):
        self.splitter.simulate()
        for i in self.feed_units_sequential:
            i.simulate()
    
    def _simulate_spike_units(self):
        self.splitter.simulate()
        for i in self.spike_units_sequential:
            i.simulate()
    
    def _simulate_upstream_units(self):
        # self.splitter.simulate()
        self._simulate_feed_units()
        self._simulate_spike_units()
    
    def get_conc_sugars(self, stream):
        return stream.imass[self.sugar_IDs].sum()/stream.ivol['Water']
    
    def get_feed_conc_sugars(self):
        return self.get_conc_sugars(stream=self.feed_mixer.outs[0])
    
    def get_spike_conc_sugars(self):
        return self.get_conc_sugars(stream=self.spike_mixer.outs[0])
