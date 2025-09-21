# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np

__all__ = ('Species', 
           'Compartment',
           'SpeciesSystem',
           'SimpleSpecies', 'ComplexSpecies')

#%% Species and species system

class Species():
    """
    Abstract class for a chemical species.
    
    Parameters
    ----------
    ID : str
        ID.
    MW: float or int
        Molar mass (molecular weight) of this species.
        
    """
    def __init__(self, ID, MW=1.):
        self._ID = ID
        self._MW = MW
    
    @property
    def ID(self): 
        return self._ID
    
    @property
    def MW(self): 
        return self._MW
    

class SimpleSpecies(Species):
    """
    Abstract class for a single chemical species that is not a complex
    of other chemical species participating in the expected species system.
    
    Parameters
    ----------
    ID : str
        ID.
    MW: float or int
        Molar mass (molecular weight) of this species.
        
    """
    def __init__(self, ID, MW=1.):
        Species.__init__(self, ID=ID, MW=MW)


class ComplexSpecies(Species):
    """
    Represents a complex formed by two or more constituent Species.
    Must provide exactly one combination of these parameters:
        - ID and simple_species_system; 
        - components
    Other parameters are optional.
    
    Parameters
    ----------

    ID : str, optional
        If provided, must be a string created by conjugating IDs of components 
        with '~'. For components with non-1 stoichiometries, ID must include these
        in the format f'number*species_ID' (as shown below), with or without parentheses
        (either is fine, although will be stored with parentheses). 
        Can be conjugated in any order, but will be stored in alphabetic order of 
        component Species IDs. Defaults to a string automatically
        created by this method from components.
        E.g., if components == {'E':1, 'S':1, 'I':2}, ID defaults to 'E~2*(I)~S'.
        
    components : dict, optional
        If provided, must be a dictionary mapping Species objects or Species IDs 
        to their stoichiometric coefficients in the complex.
        E.g., {speciesA: 1, speciesB: 2}.

    MW: float or int, optional
        Molar mass (molecular weight) of this species. Defaults to the stoichiometry-adjusted
        sum of its component species' molecular weights.
            
    simple_species_system: SimpleSpeciesSystem, optional
        Simple species system that includes all components of this complex species.
        
    Attributes
    ----------
    ID : str
        ID of the complex species.
    components : dict
        Species objects with their stoichiometric coefficients.
    MW : float
        Molecular weight of the complex, automatically computed from the sum of the components.
    """
    
    def __init__(self, ID=None, simple_species_system=None, components=None, MW=None):
        if ID is None and components is None:
            raise ValueError('\nEither ID or components must be provided, but both were None.')
        
        if ID is None:
            if isinstance(components, list):
                components = {i: 1 for i in components}
            
            for k, v in list(components.items()):
                if isinstance(k, str):
                    k_sp = simple_species_system.species(k)
                    components[k_sp] = v
                    del components[k]
                    
            self._components = components
            ID = self.ID_from_components(components)
            
        elif components is None:
            self._components = self.components_from_ID(ID)
        
        else:
            ValueError(f'\nOnly one of ID or components must be provided, but both were provided (ID={ID}, components={components}).')
        
        if MW is None:
            # calculate the MW as the sum of its parts
            MW = sum(species.MW * stoich for species, stoich in components.items())
        
        Species.__init__(self, ID=ID, MW=MW)
        
    def ID_from_components(self, components):
        components_sorted = sorted(components.items(), key=lambda i: i[0].ID)
        ID = ''
        for c in components_sorted:
            if c[1]==1:
                ID+='~'+c[0].ID
            else:
                ID+='~'+str(c[1])+'*('+c[0].ID+')'
        ID = ID[1:]
        return ID


#%%
class SpeciesSystem():
    """
    Abstract class for a system containing 
    defined chemical species.
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
    volume: float
        Volume of this system (used to calculate 'amounts'
                               by multiplying with 'concentrations')
    """
    def __init__(self, ID, all_sps, concentrations=None,
                 volume=1.0):
        self.ID = ID
        
        processed_all_sps = []
        for i in all_sps:
            if isinstance(i, str):
                processed_all_sps.append(Species(ID=i))
            elif isinstance(i, (Species, SpeciesSystem)):
                processed_all_sps.append(i)
            else:
                raise TypeError(f"\nProvided member of all_sps '{i}' must be of type str or Species.\n")
        self.all_sps = processed_all_sps
        
        if concentrations is None:
            concentrations = np.zeros(len(processed_all_sps))
        else:
            concentrations = np.array(concentrations)
        
        self._concentrations = concentrations
        
        self._volume = volume
        
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
    
    def contains(self, sp):
        if sp in self.all_sps + self.all_sp_IDs:
            return True
        else:
            return False
        
    @property
    def all_sp_IDs(self):
        all_sps = self.all_sps
        all_sp_IDs = []
        for sp in all_sps:
            if isinstance(sp, str):
                all_sp_IDs.append(sp)
            else:
                all_sp_IDs.append(sp.ID)
        return all_sp_IDs
    
    def species(self, ID_or_index):
        if isinstance(ID_or_index, int):
            return self.all_sps[ID_or_index]
        elif isinstance(ID_or_index, str):
            return self.all_sps[self.all_sp_IDs.index(ID_or_index)]
        else:
            raise ValueError(f'\nID_or_index must be int or str, but instead got {ID_or_index} of type {type(ID_or_index)}.')
    
    @property
    def concentrations(self):
        return self._concentrations
    
    @concentrations.setter
    def concentrations(self, concentrations):
        self._concentrations = concentrations
    
    def add_species(self, species, concentration=0.0):
        all_sps = self.all_sps
        if isinstance(species, str):
            species = Species(ID=species)
            all_sps.append(species)
        elif isinstance(species, Species):
            all_sps.append(species)
        else:
            raise TypeError(f"\nProvided member of all_sps '{species}' must be of type str or Species.\n")
        
        concs = list(self._concentrations)
        concs.insert(all_sps.index(species), concentration)
        self._concentrations = np.array(concs)
    
    def remove_species(self, species):
        all_sps, all_sp_IDs = self.all_sps, self.all_sp_IDs
        index = None
        if isinstance(species, str):
            index = all_sp_IDs.index(species)
        elif isinstance(species, Species):
            index = all_sps.index(species)
        else:
            raise TypeError(f"\nProvided member of all_sps '{species}' must be of type str or Species.\n")
        
        all_sps.pop(index)
        concs = list(self._concentrations)
        concs.pop(index)
        self._concentrations = np.array(concs)
    
    @property
    def volume(self):
        return self._volume
    
    @volume.setter
    def volume(self, volume):
        self._volume = volume
        
    @property
    def amounts(self):
        return self.volume * self.concentrations

Compartment = SpeciesSystem
