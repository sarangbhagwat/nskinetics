# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np
import itertools

__all__ = ('design_experiments',
           'compute_sensitivities', 
           'd_optimality',
           'fisher_information_matrix',)

def compute_sensitivities(rxn_sys, t_eval, y0, param_inds, epsilon=1e-4, spikes=None, output_idx=None):
    """
    Compute finite-difference sensitivities of ReactionSystem outputs with respect to kinetic parameters.
    
    Parameters
    ----------
    rxn_sys : ReactionSystem
        The reaction system to be evaluated.
    t_eval : array-like
        Time points at which to evaluate the solution.
    y0 : array-like
        Initial concentrations of all species.
    param_inds : list of str
        List of parameter indices to compute sensitivities for.
    epsilon : float, optional
        Relative perturbation size for finite difference (default is 1e-4).
    spikes : dict or None, optional
        Dictionary specifying species and corresponding spike amounts (optional).
    output_idx : list of int or None, optional
        Indices of species whose concentrations are observed/measured. If None, all are used.
        
    Returns
    -------
    sensitivities : ndarray of shape (n_obs, n_params)
        Flattened sensitivity matrix, where n_obs = len(t_eval) * len(output_idx),
        and n_params = len(param_inds).
    
    """
    base_params = rxn_sys.reaction_kinetic_params.copy()
    t_span = (np.min(t_eval), np.max(t_eval))
    
    base_sol = rxn_sys.solve(t_span=t_span,
                             t_eval=t_eval, 
                             y0=y0, 
                             spikes=spikes)['y']
    
    if output_idx is None:
        output_idx = range(base_sol.shape[1])
    
    sensitivities = []
    
    for param in param_inds:
        perturbed_params = base_params.copy()
        delta = epsilon * abs(base_params[param]) if base_params[param] != 0 else epsilon
        perturbed_params[param] += delta
        
        rxn_sys.set_reaction_kinetic_params(perturbed_params)
        perturbed_sol = rxn_sys.solve(t_span=t_span,
                                 t_eval=t_eval, 
                                 y0=y0, 
                                 spikes=spikes)['y']
        
        # Compute finite difference sensitivity
        sens = (perturbed_sol[:, output_idx] - base_sol[:, output_idx]) / delta
        sensitivities.append(sens.reshape(-1))  # Flatten time × output dimension
        
        rxn_sys.set_reaction_kinetic_params(base_params)  # Reset to original parameters
        
    return np.stack(sensitivities, axis=1)  # shape: (n_obs, n_params)


def fisher_information_matrix(sensitivity_matrix, sigma=1.0):
    """
    Compute the Fisher Information Matrix (FIM) from a sensitivity matrix.
    
    Parameters
    ----------
    sensitivity_matrix : ndarray
        Sensitivity matrix of shape (n_obs, n_params).
    sigma : float, optional
        Assumed standard deviation of measurement noise (default is 1.0).
        
    Returns
    -------
    FIM : ndarray of shape (n_params, n_params)
        Fisher Information Matrix: F = Sᵀ S / σ²
    
    """
    return sensitivity_matrix.T @ sensitivity_matrix / sigma**2


def d_optimality(FIM):
    """
    Compute the D-optimality score of a Fisher Information Matrix.
    
    Parameters
    ----------
    FIM : ndarray
        Fisher Information Matrix.
        
    Returns
    -------
    score : float
        D-optimality criterion, defined as det(FIM).
        A small regularization term is added for numerical stability.
    
    """
    return np.linalg.det(FIM + 1e-12 * np.eye(FIM.shape[0]))


def design_experiments(rxn_sys, param_keys, candidate_initials, t_eval,
                       spike_options=None, output_idx=None, top_n=5,
                       epsilon=1e-4):
    """
    Perform design of experiments by exhaustively evaluating combinations of initial concentrations
    and optional spike conditions, returning those that maximize parameter identifiability.
    
    Parameters
    ----------
    rxn_sys : ReactionSystem
        The reaction system to be evaluated.
    param_keys : list of str
        List of parameter names to be identified.
    candidate_initials : dict
        Dictionary mapping species names to lists of candidate initial concentrations.
        Example: {'E': [0.1, 1.0], 'S': [0.5, 2.0]}
    t_eval : array-like
        Time points at which to simulate the system.
    spike_options : dict or None, optional
        Dictionary mapping species to lists of possible spike values (or None for no spikes).
        Example: {'S': [0.5, 1.0]}.
    output_idx : list of int or None, optional
        Indices of species that are experimentally measured (default is all).
    top_n : int, optional
        Number of top experimental designs to return based on D-optimality (default is 5).
        
    Returns
    -------
    top_designs : list of dict
        Each dict contains:
            'y0'     : ndarray, initial concentrations
            'spikes' : dict or None, spike specification
            'FIM'    : ndarray, Fisher Information Matrix
            'score'  : float, D-optimality score
    
    """
    designs = []
    
    # Create all combinations of initial concentrations
    species = list(candidate_initials.keys())
    initial_combos = list(itertools.product(*candidate_initials.values()))
    
    param_inds = [rxn_sys.reaction_kinetic_param_keys.index(param_key)
                  for param_key in param_keys]
    
    # Create all combinations of spike values
    if spike_options:
        spike_species = list(spike_options.keys())
        spike_combos = list(itertools.product(*spike_options.values()))
    else:
        spike_combos = [None]
        
    for y0_vals in initial_combos:
        y0 = np.array(y0_vals)
        
        for spike_setting in spike_combos:
            if spike_setting is not None:
                spikes = {name: [val] for name, val in zip(spike_options.keys(), spike_setting)}
            else:
                spikes = None
                
            try:
                # Compute sensitivity and FIM
                S = compute_sensitivities(rxn_sys, t_eval, np.array(y0), param_inds, spikes=spikes, output_idx=output_idx, epsilon=epsilon)
                FIM = fisher_information_matrix(S)
                score = d_optimality(FIM)
                
                # Store experiment design and score
                designs.append({
                    'y0': y0,
                    'spikes': spikes,
                    'FIM': FIM,
                    'score': score
                })
            except Exception as e:
                print(f"[Warning] Failed for y0={y0}, spikes={spikes}: {e}")
                
    # Sort by D-optimality and return top_n
    top_designs = sorted(designs, key=lambda d: d['score'], reverse=True)[:top_n]
    return top_designs
