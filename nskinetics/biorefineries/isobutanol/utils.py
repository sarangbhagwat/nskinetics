# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.


import csv

__all__ = ('generate_save_kinetic_parameter_distributions',)


def generate_save_kinetic_parameter_distributions(nsk_kinetic_model, filename, 
                                                  fermentation_reactor_ID='V406',
                                                  shape='triangular', factors=(0.8, 1.0, 1.2)):
    r_te = nsk_kinetic_model._te
    all_params = r_te.getGlobalParameterIds()
    kinetic_params = [i for i in all_params if i[:2].lower()=='k_']
    
    data = []
    data.append(['Parameter name', 'Element', 'Kind', 'Units', 
                 'Baseline', 'Shape', 'Lower', 'Midpoint', 'Upper', 
                 'References', 'Load statement'])
    
    for k in kinetic_params:
        if 'k' in k:
            units = 'g_per_l_per_h' 
        elif 'K' in k:
            units = 'g_per_l'
        else:
            units = None
        baseline = r_te.__getattribute__(k)
        data.append([k, 'Fermentation Kinetics', 'coupled', units,
                     baseline, shape, baseline*factors[0], baseline*factors[1], baseline*factors[2], 
                     None, f'{fermentation_reactor_ID}.nsk_kinetic_model._te.{k} = x'])
    
    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
