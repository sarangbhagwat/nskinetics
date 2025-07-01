# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np

__all__ = ('Species', 'SpeciesSystem',)

#%% Species and species system

class Species():
    """
    Abstract class for a single chemical species.
    
    Parameters
    ----------
    ID : str
    
    """
    def __init__(self, ID):
        self.ID = ID

class SpeciesSystem():
    """
    Abstract class for a system of chemical species.
    Note this class will, by convention,
    use "sp" to denote a single species and
    use "sps" to denote multiple species.
    
    Parameters
    ----------
    ID : str
        ID.
    all_sps : list
        List of Species objects or strings. 
        If the latter, new Species objects will be created.
    concentrations : list or np.ndarray, optional
        Concentrations of species, indexed the same way as
        all_sps. Defaults to zero-array of length equal to
        that of all_sps.
    
    """
    def __init__(self, ID, all_sps, concentrations=None):
        self.ID = ID
        processed_all_sps = []
        for i in all_sps:
            if isinstance(i, str):
                processed_all_sps.append(Species(ID=i))
            elif isinstance(i, Species):
                processed_all_sps.append(i)
            else:
                raise TypeError(f"\nProvided member of all_sps '{i}' must be of type str or Species.\n")
        self.all_sps = processed_all_sps
        
        if concentrations is None:
            concentrations = np.zeros(len(processed_all_sps))
        
        self.concentrations = np.array(concentrations)
        
    def indices(self, some_sps):
        all_sps = self.all_sps
        indices = [all_sps.index(i) for i in some_sps]
        return indices
    
    def index(self, sp):
        if isinstance(sp, Species):
            return self.all_sps.index(sp)
        elif isinstance(sp, str):
            return self.index_from_ID(sp)
    
    def index_from_ID(self, sp_ID):
        index = 0
        all_sps = self.all_sps
        for i in all_sps:
            if i.ID==sp_ID:
                return index
            index += 1