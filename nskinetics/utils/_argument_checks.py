# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np

__all__ = ('is_number', 
           'is_array_of_numbers',
           'is_list_of_strings')

def is_number(i):
    return isinstance(i,float)\
            or isinstance(i,int)\
            or isinstance(i, np.floating)\
            or isinstance(i, np.integer)
            
def is_array_of_numbers(i):
    return isinstance(i, np.ndarray) and\
        (np.issubdtype(i.dtype, np.floating) or
         np.issubdtype(i.dtype, np.integer))

def is_list_of_strings(i):
    for s in i:
        if not isinstance(s, str):
            return False
    return True
