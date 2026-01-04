# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import math
import numpy as np
__all__ = ('get_index_nearest_element_from_sorted_array',)

#%%
def get_index_nearest_element_from_sorted_array(array, value):
    """
    Finds the value in a sorted NumPy array closest to a given value 
    using binary search for efficiency.
    """
    # Find the insertion point that maintains the sorted order
    idx = np.searchsorted(array, value, side="left")
    
    # Handle edge cases
    if idx == len(array):
        return -1
    if idx == 0:
        return 0
    
    # Compare the values at the left and right insertion points
    if math.fabs(value - array[idx-1]) < math.fabs(value - array[idx]):
        return idx-1
    else:
        return idx