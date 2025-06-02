# -*- coding: utf-8 -*-
"""
Created on Thu May 29 17:37:32 2025

@author: sarangbhagwat
"""
import numpy as np

__all__ = ('Species, SpeciesSystem')

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
    """
    def __init__(self, ID, all_sps, concentrations):
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
        self.concentrations = np.array(concentrations)
        
    def get_indices(self, some_sps):
        all_sps = self.all_sps
        indices = [all_sps.index(i) for i in some_sps]
        return indices
    
        