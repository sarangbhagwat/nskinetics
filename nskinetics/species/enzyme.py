# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

from . import Species
import numpy as np

__all__ = ('Enzyme', 'EnzymeComplex',
           'enzyme_complex_joiner')

enzyme_complex_joiner = '~'

#%% Enzyme and enzyme complex

class Enzyme(Species):
    """
    Class for an enzyme species.
    Can be used for simple or multi-module (assembly line)
    enzymes.
    Site availability properties can be used by ReactionSystem
    objects to come up with a list of reactions associated with
    this enzyme.
    
    Parameters
    ----------
    ID : str
    unavailable_sites: List or set, optional
        List or set of (integer/float, Species/string) tuples 
        (can be mixed freely) representing unavailable sites.
        
        Each tuple contains a Species (no Species IDs) that 
        cannot bind with the site corresponding to the contained 
        integer/float.
        The only allowable string is 'all', for which the site 
        is considered unavailable to all species.
        If the integer/float is np.inf, this is assumed to denote
        all possible sites.
        
        Defaults to an empty set (i.e., all sites are available).
    
    """
    def __init__(self, ID, 
                 unavailable_sites=None,):
        Species.__init__(self, ID=ID)
        self._unavailable_sites = set(unavailable_sites)\
            if unavailable_sites is not None\
            else set()
        
    @property
    def unavailable_sites(self):
        return self._unavailable_sites

    def is_available(self, site, binding_species=None):
        bsp = binding_species
        ua = self._unavailable_sites
        is_generally_unavailable = self.is_generally_unavailable
        is_explicitly_unavailable = self.is_explicitly_unavailable
        
        # if site unavailable
        if (np.inf, 'all') in ua\
            or is_generally_unavailable(site=site)\
            or is_generally_unavailable(binding_species=bsp)\
            or is_explicitly_unavailable(site=site, binding_species=bsp):
               return 0

        # if site not unavailable
        else:
            return 1
        
    def is_explicitly_unavailable(self, site=None, binding_species=None):
        if (site, binding_species) in self._unavailable_sites:
            return 1
        else:
            return 0
    
    def is_generally_unavailable(self, site=None, binding_species=None):
        ua = self._unavailable_sites
        bsp = binding_species
        if site is None and bsp is None:
            raise ValueError('Either site or binding_species must be provided.')
        elif site is not None and (site, 'all') in ua:
            return 1
        elif bsp is not None and (np.inf, bsp) in ua:
            return 1
        else:
            return 0
        
    def make_site_available(self, site, binding_species=None):
        ua = self._unavailable_sites
        bsp = binding_species if binding_species is not None else 'all'
        is_generally_unavailable = self.is_generally_unavailable
        is_explicitly_unavailable = self.is_explicitly_unavailable
        
        if is_explicitly_unavailable(site, bsp):
            ua.remove((site, bsp))
        if is_generally_unavailable(site=site):
            pass
            # !!! implementation in progress
            
    def make_site_unavailable(self, site, binding_species):
        ua = self._unavailable_sites
        bsp = binding_species
        if is_available(site=site, binding_species=bsp):
            ua.add((site, bsp))
            
class EnzymeComplex(Enzyme):
    """
    Class for an enzyme-(substrate/inhibitor) complex species.
    
    Parameters
    ----------
    ID : str
    unavailable_sites: List or set, optional
        List or set of (integer/float, Species/string) tuples 
        (can be mixed freely) representing unavailable sites.
        
        Each tuple contains a Species (no Species IDs) that 
        cannot bind with the site corresponding to the contained 
        integer/float.
        The only allowable string is 'all', for which the site 
        is considered unavailable to all species.
        If the integer/float is np.inf, this is assumed to denote
        all possible sites.
        
        Defaults to [(np.inf, 'all')] (i.e., no sites are available).
    
    """
    def __init__(self, ID, 
                 enzyme,
                 binding_species,
                 binding_site,
                 unavailable_sites=[(np.inf, 'all')],
                 sites_unavailable_to_specified_binding_species=None):
        
        assert isinstance(enzyme, Enzyme)
        assert isinstance(binding_species, Species)
        assert not isinstance(binding_species, Enzyme)
        
        Enzyme.__init__(self, ID=ID, 
                        unavailable_sites=unavailable_sites)
        
        self._enzyme = enzyme
        self._binding_site = binding_site
    
    def from_enzyme_and_binding_species(enzyme,
                                        binding_species,
                                        binding_site,
                                        unavailable_sites='all',):

        return Enzyme(ID=enzyme+enzyme_complex_joiner+binding_species,
                      unavailable_sites=unavailable_sites,
                      )
    