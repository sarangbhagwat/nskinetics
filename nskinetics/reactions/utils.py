# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 22:52:21 2025

@author: sarangbhagwat
"""
import numpy as np

#%% Reading an equation string to get stoichiometries and kinetic parameters

def get_eqn_and_param_info_split_str(equation_str, 
                                     param_info_start,
                                     param_info_junk):
    param_info_junk.append(param_info_start)
    split_str, param_info = None, None
    if param_info_start in equation_str:
        split_str = equation_str[:equation_str.index(param_info_start)].split(' ')
        param_info = equation_str[equation_str.index(param_info_start):].replace('=', ' ')
        for i in param_info_junk:
            param_info = param_info.replace(i, ' ')
        param_info = param_info.split(' ')
    else:
        split_str = equation_str.split(' ')
    
    return split_str, param_info

def reformat_to_handle_numbers_in_sp_IDs(equation_str, 
                                         split_str,
                                         all_sp_IDs,
                                         conjugations):
    # Handle instances of stoichiometry number immediately next to species name
    # e.g., 2A rather than 2 A; make sure 2A isn't a chemical
    replace={}
    last_str_was_float = False
    for str_ in split_str:
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
        ind = split_str.index(k)
        split_str[ind] = v[0]
        split_str.insert(ind+1, v[1])
        
    return split_str

def get_stoichiometry(split_str, all_sp_IDs, arrows):
    # Get stoichiometry array
    stoichiometry = []
    arrow_ind = None
    for arrow in arrows:
        if arrow in split_str:
            arrow_ind = split_str.index(arrow)
    
    for sp_ID in all_sp_IDs:
        if sp_ID in split_str:
            try:
                stoichiometry.append(float(split_str[split_str.index(sp_ID)-1]))
            except:
                stoichiometry.append(1.)
            if split_str.index(sp_ID)<arrow_ind:
                stoichiometry[-1] *= -1
        else:
            stoichiometry.append(0.)
    
    return stoichiometry

def get_kinetic_params(param_info, param_info_junk):
    # Get kinetic parameters, if any
    kf_kb = [None, None]
    if param_info is not None:
        param_info_clean = [i for i in param_info if not i in param_info_junk + ['']]
        kf_chars = ['kf',]
        kb_chars = ['kb',]
        for i in range(len(param_info_clean)-1):
            if param_info_clean[i] in kf_chars:
                kf_kb[0] = float(param_info_clean[i+1])
            elif param_info_clean[i] in kb_chars:
                kf_kb[1] = float(param_info_clean[i+1])
    return kf_kb[0], kf_kb[1]

def read_equation_str(equation_str, species_system):
    stoichiometry = []
    
    arrows = ['->', '<->']
    conjugations = ['+'] + arrows
    param_info_start=';'
    param_info_junk=['=', ',',]
    all_sp_IDs = [i.ID for i in species_system.all_sps]
    
    split_str, param_info = get_eqn_and_param_info_split_str(
        equation_str=equation_str, 
        param_info_start=param_info_start,
        param_info_junk=param_info_junk)
    
    split_str = reformat_to_handle_numbers_in_sp_IDs(
        equation_str=equation_str, 
        split_str=split_str, 
        all_sp_IDs=all_sp_IDs,
        conjugations=conjugations)
    
    stoichiometry = get_stoichiometry(split_str=split_str, 
                                      all_sp_IDs=all_sp_IDs,
                                      arrows=arrows)
    kf, kb = get_kinetic_params(param_info=param_info, 
                       param_info_junk=param_info_junk)

    return np.array(stoichiometry), kf, kb
