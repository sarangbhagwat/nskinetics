# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""
Created on Sun Jan 11 18:10:51 2026

@author: sarangbhagwat
"""
from nskinetics.biorefineries import isobutanol
isobutanol.load()
from nskinetics.biorefineries.isobutanol import scenarios

model = isobutanol.models.models_EtOH_IBO_corn.model
plot_kinetic_results = isobutanol.models.plot_kinetic_results


def load_simulate_baseline(scenario='B', # 'A' or 'B'
                           plot=False,
                           stream_IDs=('ethanol', 'isobutanol'), # products whose MPSPs are solved
                           IRR_for_MPSP=0.15, # fixed IRR at which MPSPs are solved
                           n_sims=3, # simulations in this kernel (in-process stability check)
                           ):
    # Consolidated loader: workbook baselines + distributions + feeding
    # strategy + one baseline model_specification (scenarios.load_scenario).
    bundle = scenarios.load_scenario(scenario)
    model_specification = bundle['model_specification']
    solve_TEA = bundle['solve_TEA']
    feeding_kwargs = bundle['feeding_kwargs']
    expected_MPSPs = bundle['expected']

    all_results = []
    for i in range(n_sims):
        model_specification(**feeding_kwargs)
        results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
        scenarios.assert_MPSPs_pinned(expected_MPSPs, results, i+1)
        scenarios.assert_spike_feed_residual(i+1)
        if all_results:
            scenarios.assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)

    if plot:
        fig, ax = plot_kinetic_results(xlim=(0,80), ylim=(0,250))

    return all_results
