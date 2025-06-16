# -*- coding: utf-8 -*-
"""
Created on Thu May 29 17:39:02 2025

@author: sarangbhagwat
"""

import numpy as np

__all__ = ('ReactionSystem', 'RxnSys')

#%% Reaction system

class ReactionSystem():
    """
    Abstract class for a system of reactions.
    
    Parameters
    ----------
    ID : str
        ID.
    reactions : list
        List of Rxn, ReversibleRxn, or RxnSystem objects.
    species_system : SpeciesSystem
        A SpeciesSystem object containing all species
        involved in this system of reactions.
    """
    def __init__(self, ID, reactions, species_system):
        self.ID = ID
        self.reactions = reactions
        self.species_system = species_system
    
    def get_dconcs_dt(self):
        reactions = self.reactions
        # species_concs_vector = self.species_system.all_sps
        # breakpoint()
        return np.sum([r.get_dconcs_dt() for r in reactions], axis=0)

RxnSys = ReactionSystem
