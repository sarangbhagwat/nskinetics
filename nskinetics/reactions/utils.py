# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np

#%% Reading an equation string to get stoichiometries and kinetic parameters

def get_eqn_and_param_info_strs(equation_str, 
                                     delimiter,
                                     eqn_only_idfiers,
                                     kf_idfiers, kb_idfiers):
    
    eqn_only, param_info = None, None
    if delimiter in equation_str:
        str_parts = equation_str.split(delimiter)
        for str_part in str_parts:
            
            for eqn_only_idfier in eqn_only_idfiers:
                if eqn_only_idfier in str_part:
                    eqn_only = str_part
                    
            for param_idfier in kf_idfiers+kb_idfiers:
                if param_idfier in str_part:
                    param_info = str_part
    
    else:
        for eqn_only_idfier in eqn_only_idfiers:
            if eqn_only_idfier in equation_str:
                eqn_only = equation_str
        if eqn_only is None:
            raise ValueError(f'Could not parse equation string {equation_str}.')
    
    # print(equation_str)
    # print(eqn_only, '|', param_info)
    return eqn_only, param_info


def reformat_to_handle_numbers_in_sp_IDs(eqn_only,
                                         equation_str, 
                                         all_sp_IDs,
                                         conjugations):
    eqn_only_split = eqn_only.split(' ')
    # Handle instances of stoichiometry number immediately next to species name
    # e.g., 2A rather than 2 A; make sure 2A isn't a chemical
    replace={}
    last_str_was_float = False
    for str_ in eqn_only_split:
        if str_ not in all_sp_IDs + conjugations: # is it a sp_ID or conjugation
            resolved=False
            try:
                if last_str_was_float: 
                    raise ValueError(f'Equation string "{equation_str}" has at least numbers one after the other that with neither being part of a registered chemical name in the species system {all_sp_IDs}.')
                float(str_) # is it a number by itself (stoichiometry)
                last_str_was_float = True
                resolved = True
            except Exception as e:
                use_index_upto=1
                for i in range(1, len(str_)):
                    try:
                        float(str[:i])
                        use_index_upto+=1
                    except:
                        if i>1:
                            stoich_part = str_[:use_index_upto]
                            sp_ID_part = str_[use_index_upto:]
                            replace[str_] = stoich_part, sp_ID_part
                            resolved=True
                            if last_str_was_float:
                                raise e
                            break
                        else:
                            RuntimeError(f'Error reading equation string "{equation_str}"')
                if not use_index_upto==len(str_)-1:
                    last_str_was_float = False
            if not resolved:
                raise ValueError(f'Equation string "{equation_str}" contains element "{str_}" that was not identified as a species, stoichiometry number, conjugation')
    for k, v in replace.items():
        ind = eqn_only_split.index(k)
        eqn_only_split[ind] = v[0]
        eqn_only_split.insert(ind+1, v[1])
    
    # print(eqn_only_split)
    return eqn_only_split


def get_stoichiometry(eqn_only_split, all_sp_IDs, arrows):
    # Get stoichiometry array
    stoichiometry = []
    arrow_ind = None
    for arrow in arrows:
        if arrow in eqn_only_split:
            arrow_ind = eqn_only_split.index(arrow)
    
    for sp_ID in all_sp_IDs:
        if sp_ID in eqn_only_split:
            try:
                stoichiometry.append(float(eqn_only_split[eqn_only_split.index(sp_ID)-1]))
            except:
                stoichiometry.append(1.)
            if eqn_only_split.index(sp_ID)<arrow_ind:
                stoichiometry[-1] *= -1
        else:
            stoichiometry.append(0.)
            
    # print(stoichiometry)
    return stoichiometry


def get_kinetic_params(param_info, param_info_junk, kf_idfiers, kb_idfiers):
    if param_info is None: return None, None
    param_info = param_info.replace('=', ' ')
    for i in param_info_junk:
        param_info = param_info.replace(i, ' ')
    param_info = param_info.split(' ')
    # Get kinetic parameters, if any
    kf_kb = [None, None]
    if param_info is not None:
        param_info_clean = [i for i in param_info if not i in param_info_junk + ['']]

        for i in range(len(param_info_clean)-1):
            if param_info_clean[i] in kf_idfiers:
                kf_kb[0] = float(param_info_clean[i+1])
            elif param_info_clean[i] in kb_idfiers:
                kf_kb[1] = float(param_info_clean[i+1])
    return kf_kb[0], kf_kb[1]


def read_equation_str(equation_str, species_system):
    stoichiometry = []
    
    arrows = ['->', '<->']
    conjugations = ['+'] + arrows
    
    delimiter=';'
    
    eqn_only_idfiers = arrows
    kf_idfiers = ['kf',]
    kb_idfiers = ['kb',]
    
    param_info_junk = ['=', ',',]
    
    all_sp_IDs = [i.ID for i in species_system.all_sps]
    
    eqn_only, param_info = get_eqn_and_param_info_strs(
        equation_str=equation_str, 
        delimiter=delimiter,
        eqn_only_idfiers=eqn_only_idfiers,
        kf_idfiers=kf_idfiers, kb_idfiers=kb_idfiers)
    
    eqn_only_split = reformat_to_handle_numbers_in_sp_IDs(
        equation_str=equation_str, 
        eqn_only=eqn_only, 
        all_sp_IDs=all_sp_IDs,
        conjugations=conjugations)
    
    stoichiometry = get_stoichiometry(eqn_only_split=eqn_only_split, 
                                      all_sp_IDs=all_sp_IDs,
                                      arrows=arrows)
    
    kf, kb = get_kinetic_params(param_info=param_info, 
                       param_info_junk=param_info_junk,
                       kf_idfiers=kf_idfiers, kb_idfiers=kb_idfiers)
    
    return np.array(stoichiometry), kf, kb
