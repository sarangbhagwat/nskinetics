# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import nskinetics as nsk
import numpy as np
from warnings import filterwarnings

def test_simple_ESP_inverse_modeling_and_doe():
    plot = False
    show_progress = False
    show_output = False
    show_warnings = False
    
    if not show_warnings: filterwarnings("ignore")
    
    np.random.seed(3221)
    
    #%% Create a SpeciesSystem object with defined kinetic parameters
    sp_sys = nsk.SpeciesSystem('sp_sys', 
                           ['E', 'S', 'ES', 'P'], # enzyme, substrate, enzyme-substrate complex, product
                           concentrations=[1e-4, 1e-4, 0., 0.])
    
    # Describe reactions by writing chemical equations and kinetic parameter info
    reactions = [
                'E + S <-> ES; kf = 12.0, kb = 10.0', # kf = kon, kb = koff
                'ES -> E + P; kf = 32.0' # kf = kcat (enzyme turnover number)
                ]
    
    # Generate a ReactionSystem from strings
    rxn_sys = nsk.ReactionSystem(ID='rxn_sys', 
                                     reactions=reactions,
                                     species_system=sp_sys)
    
    KM = (rxn_sys.reactions[1].kf + rxn_sys.reactions[0].kb)/rxn_sys.reactions[0].kf
    
    #%% Define batches to run
    
    batch_concs = [
                   # np.array([1e-4, 1e-4, 0., 0.]),
                   np.array([0.08966552, 0.58624828, 0., 0.]),
                   # np.array([0.08966552, 0.62072759, 0., 0.]),
                   # np.array([0.0621069, 0.10353793, 0., 0.]),
                   ]
    batch_t_spans = [
                     # [0, 0.5*24*3600],
                     [0, 300],
                     # [0, 30],
                     # [0, 30],
                     ]
    
    n_batches = len(batch_concs)
    batch_indices = range(n_batches)
    
    #%% Simulate all batches and save results
    
    dfs_results = []
    
    for i, concentrations, t_span in zip(batch_indices, 
                                         batch_concs, 
                                         batch_t_spans):
        sp_sys.concentrations = concentrations
        # Simulate the ReactionSystem
        rxn_sys.solve(t_span=t_span,
                      sp_conc_for_events={'S':1e-6},
                      # filename=f'solution{i}',
                      save_events=False)                              
        dfs_results.append(rxn_sys._solution_dfs[0])
        # Plot results
        if plot:
            rxn_sys.plot_solution()
            rxn_sys.plot_solution(sps_to_include=['ES'])
    
    #%% Fit kinetic parameters to saved results from simulating with original set of kinetic parameters
    
    rxn_sys.reactions[1]._freeze_kf = True
    rxn_sys.reactions[0]._freeze_kb = True
    
    def set_koff(p, kcat=rxn_sys.reactions[1].kf, KM=KM):
        # With a known KM and known (frozen) kcat, 
        # set koff based on the kon value being used by the solver
        kon = p[0]
        rxn_sys.reactions[0]._kb = KM*kon - kcat
    
    
    rxn_sys.fit_reaction_kinetic_parameters_to_data(
                                                    # data=[f'solution{i}.xlsx'
                                                    #       for i in batch_indices],
                                                    data=dfs_results,
                                                    p0=1*np.ones(len(rxn_sys.reaction_kinetic_params)),
                                                    use_only=[
                                                              # 'E', 
                                                              'S', 
                                                              # 'P'
                                                              ],
                                                    call_before_each_solve=[set_koff],
                                                    options={'disp':show_progress},
                                                    # plot_during_fit=True,
                                                    show_progress=show_progress,
                                                    show_output=show_output,
                                                    n_minimize_runs=2,
                                                    n_de_runs=3,
                                                    differential_evolution_kwargs={'maxiter':100,
                                                                                   'bounds': [(0,5e4)],
                                                                                   'polish': False,
                                                                                   'disp': False,
                                                                                   'workers': 1,
                                                                                   'popsize': 10,
                                                                                   'tol': 1e-3,
                                                                                   }
                                                    )
    
    #%% Simulate batches with inverse-modeled rxn_sys
    for i, concentrations, t_span in zip(range(len(batch_concs)), 
                                         batch_concs, 
                                         batch_t_spans):
        sp_sys.concentrations = concentrations
        # Simulate the ReactionSystem
        rxn_sys.solve(t_span=t_span,
                      sp_conc_for_events={'S':1e-6},
                      # filename=f'solution{i}',
                      save_events=False,
                      )                  
        # Plot results            
        if plot:
            rxn_sys.plot_solution()
            rxn_sys.plot_solution(sps_to_include=['ES'])
        
    #%% Design of experiments
    param_keys = rxn_sys.reaction_kinetic_param_keys
    candidate_initials = {'E': np.linspace(1e-4, 0.1, 30), 'S': np.linspace(1e-4, 1., 30), 'ES': [0.0], 'P': [0.0]}
    spike_options = None
    t_eval = np.linspace(0, 30, 50)
    output_idx = [sp_sys.index('S'), 
                  sp_sys.index('P')]  # Measuring S and P
    
    best_expts = nsk.doe.optimal_design(
        rxn_sys,
        param_keys,
        candidate_initials,
        t_eval,
        spike_options=spike_options,
        output_idx=output_idx,
        epsilon=1e-4,
        top_n=3,
        show_output=show_output,
        )
    
    # Tests
    
    assert np.allclose(rxn_sys._fitsol[0], 
                       np.array([11.99998575]),
                       rtol=1e-3, atol=1e-8)
    
    assert np.allclose(rxn_sys._fitsol[1], 
                       0.9999999999992784,
                       rtol=1e-5, atol=1e-8)
    
    assert np.allclose(best_expts[0]['y0'],
                       np.array([0.08966552, 
                                 0.55176897, 
                                 0., 
                                 0.]),
                       rtol=1e-5, atol=1e-8)
    
    assert np.allclose(best_expts[0]['FIM'],
                       np.array([[5.12779117e-05]]),
                       rtol=1e-1, atol=1e-8)
    
    assert np.allclose(best_expts[0]['score'],
                       5.127791265124837e-05,
                       rtol=1e-1, atol=1e-8)
