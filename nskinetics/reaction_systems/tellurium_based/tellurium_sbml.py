# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import tellurium as te

__all__ = ('TelluriumReactionSystem',)

#%% Species and species system

class TelluriumReactionSystem():
    """
    Abstract class for a chemical species.
    
    Parameters
    ----------
    te : tellurium.roadrunner.extended_roadrunner.ExtendedRoadRunner
        A Tellurium extended RoadRunner object.\
        
    """
    def __init__(self, te, units=None):
        self._te = te
        if units is not None:
            self._units = units
        else:
            self._units = {'time': 'min',
                           'conc': 'M'}
    @classmethod
    def from_sbml(filepath):
        return te.loadSBMLModel(filepath)
    
    def reset(self):
        self._te.reset()
