# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import os

import nskinetics as nsk
import biosteam as bst
import tellurium as te
import simplesbml

__all__ = ('te_r', 'reset_kinetic_reaction_system',)

#%%
antimony_filename = 's_cerevisiae_ferm_fb_inhib_mod_ibo_antimony.txt'
r = te.loadAntimonyModel(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), antimony_filename))

#%% Create NSKinetics kinetic model
te_r = nsk.KineticModel(r)
te_r._units['time'] = 'h'
te_r._units['conc'] = 'g/L'
te_r.default_max_n_glu_spikes = 200

#%% Declare fed-batch + stage-1 events via the NSKinetics event API
glucose_feed_spike = nsk.FeedSpike(
    species='s_glu',
    when='s_glu <= threshold_conc_glu_spike',
    target='target_conc_glu_spike',
    feed_conc='conc_glu_feed_spike',
    volume_var='env',
    max_count='max_n_glu_spikes',
    count_var='n_glu_spikes',
    last_vol_var='last_vol_glu_feed_added',
    tot_vol_var='tot_vol_glu_feed_added',
    delay='glucose_feed_spikeDelay',
    priority=5,
    name='glucose_feed_spike',
)
for _event in glucose_feed_spike.expand():
    te_r.add_event(_event)

te_r.add_event(nsk.Event(when='time >= stage_1_max_time',
                         do={'is_aerobic': '0'},
                         name='stage_1_complete_max_time'))
te_r.add_event(nsk.Event(when='x >= stage_1_max_x',
                         do={'is_aerobic': '0'},
                         name='stage_1_complete_x_target'))

te_r.compile_events()

def reset_kinetic_reaction_system(km, reset_spike_cap=True, **kwargs):
    km._te.reset()
    r_te = km._te
    r_te.n_glu_spikes = 0
    r_te.last_vol_glu_feed_added = 0.
    r_te.tot_vol_glu_feed_added = 0.
    r_te.env = 1.
    r_te.is_aerobic = 1
    if reset_spike_cap: r_te.max_n_glu_spikes = km.default_max_n_glu_spikes

te_r.reset_func = reset_kinetic_reaction_system

#%% Set tolerances
integrator = r.integrator
integrator.absolute_tolerance = 1e-10
integrator.relative_tolerance = 1e-9
# integrator.variable_step_size = True

#%% Document references for IBO pathway kinetic parameters

# r.k_6 = 2.9 * 88.06 * 1e-6 * 60 # 2.9 uM/min # Zhao 2021 https://doi.org/10.3390/foods10051013
# r.K_6 = 31.8 * 88.06 * 1e-3 # 31.8 mM # Zhao 2021 https://doi.org/10.3390/foods10051013

#%%
simulate = False
if simulate:
    reset_kinetic_reaction_system(te_r)
    r.s_glu = 100 # initial glucose conc
    r.x = 1 # initial biomass conc
    print(r.s_glu, r.x, r.X_a, r.s_EtOH, r.s_acetate, r.s_acetald, r.s_AL, r.s_DHI, r.s_KIV, r.s_IBO)
    
    r.simulate(0, 300, 2000,
               # ['time', 's_glu', '[s_glu]'],
               # ['time', 'X_a', 'X_AcDH', 
               #  'a',
               #  ],
               )
    
    print(r.s_glu, r.x, r.X_a, r.s_EtOH, r.s_acetate, r.s_acetald, r.s_AL, r.s_DHI, r.s_KIV, r.s_IBO)
    r.plot()
