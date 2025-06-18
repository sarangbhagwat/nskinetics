# -*- coding: utf-8 -*-
"""
Created on Mon May 26 13:50:14 2025

@author: sarangbhagwat
"""
import numpy as np
from matplotlib import pyplot as plt
from ..species import SpeciesSystem
from ..reactions import Rxn, RevRxn
from .reaction_system import RxnSys

__all__ = ('EnzymeSubstrateProduct', 'ESP', 'CompetitiveInhibition',
           'NonCompetitiveInhibition', 
           'UnCompetitiveInhibition',
           'MechanismBasedInhibition')

#%% Core enzyme-catalyzed reaction system (can also be used for competitive inhibition)

class EnzymeSubstrateProduct(RxnSys):
    """
    Class for a system of reactions involving 
    a single enzyme, substrate, and product.
    
    Parameters
    ----------
    ID : str
        ID.
    enzyme : str
        Enzyme.
    substrate : str
        Substrate.
    es_complex : str
        Enzyme-Substrate complex.
    product : str
        Product.
    kon : float
        Forward reaction rate constant for the formation
        of enzyme-substrate complex from enzyme and
        substrate.
    koff : float
        Backward reaction rate constant for the formation
        of enzyme-substrate complex from enzyme and
        substrate.
    kcat : float
        Rate constant for the dissociation of
        enzyme-substrate complex to form enzyme
        and product.
    species_system : SpeciesSystem
        A SpeciesSystem object containing all species
        involved in this system of reactions.
    """
    def __init__(self, 
                 ID, 
                 kon, koff, kcat,
                 species_system,
                 enzyme=None, 
                 substrate=None, 
                 es_complex=None,
                 product=None, 
                 chem_equations=None,
                 stoichiometries=None):
        
        self.loading_rxn = loading_rxn = RevRxn(ID=ID+'_load', 
                                   reactants=[enzyme, substrate], 
                                   products=[es_complex],
                                   kf=kon, kb=koff,
                                   species_system=species_system)
        self.catalyzed_rxn = catalyzed_rxn = Rxn(ID=ID+'_cat',
                           reactants=[es_complex],
                           products=[enzyme, product],
                           kf=kcat,
                           species_system=species_system)
        
        RxnSys.__init__(self, ID=ID, 
                           reactions=[loading_rxn, catalyzed_rxn], 
                           species_system=species_system)
        
        self.enzyme = enzyme
        self.substrate = substrate
        self.es_complex = es_complex
        self.product = product
        
    @property
    def kon(self):
        return self.loading_rxn.kf
    @kon.setter
    def kon(self, kon):
        self.loading_rxn.kf = kon
    
    @property
    def koff(self):
        return self.loading_rxn.kb
    @koff.setter
    def koff(self, koff):
        self.loading_rxn.kb = koff
        
    @property
    def kcat(self):
        return self.catalyzed_rxn.kf
    @kcat.setter
    def kcat(self, kcat):
        self.catalyzed_rxn.kf = kcat
    
CompetitiveInhibition = ESP = EnzymeSubstrateProduct

#%% Non-competitive inhibition

class NonCompetitiveInhibition(RxnSys):
    """
    Class for a system of reactions involving 
    a single enzyme and an allosteric inhibitor.
    Consistent with non-competitive inhibition, these reversibly
    form an enzyme-inhibitor complex and can further reversibly 
    form an enzyme-inhibitor-substrate complex. The latter can
    further reversibly dissociate to form the enzyme-substrate 
    complex.
    
    Parameters
    ----------
    enzyme : str
        Enzyme.
    substrate : str
        Substrate.
    inhibitor : str
        Inhibitor.
    es_complex : str
        Enzyme-Substrate complex.
    ei_complex : str
        Enzyme-Inhibitor complex.
    esi_complex : str
        Enzyme-Substrate-Inhibitor complex.
    kon_ei : float
        Forward reaction rate constant for the formation
        of enzyme-inhibitor complex from enzyme and
        inhibitor.
    koff_ei : float
        Backward reaction rate constant for the formation
        of enzyme-inhibitor complex from enzyme and
        inhibitor.
    kon_ei_esi : float
        Forward reaction rate constant for the formation
        of enzyme-substrate-inhibitor complex from 
        enzyme-inhibitor complex and substrate.
    koff_ei_esi : float
        Backward reaction rate constant for the formation
        of enzyme-substrate-inhibitor complex from 
        enzyme-inhibitor complex and substrate.
    kon_es_esi : float
        Forward reaction rate constant for the formation
        of enzyme-substrate-inhibitor complex from 
        enzyme-substrate complex and inhibitor.
    koff_es_esi : float
        Backward reaction rate constant for the formation
        of enzyme-substrate-inhibitor complex from 
        enzyme-substrate complex and inhibitor.
    species_system : SpeciesSystem
        A SpeciesSystem object containing all species
        involved in this system of reactions.
    """
    def __init__(self, 
                 ID, 
                 enzyme, 
                 inhibitor,
                 substrate, 
                 ei_complex,
                 es_complex,
                 esi_complex,
                 kon_ei, koff_ei,
                 kon_ei_esi, koff_ei_esi,
                 kon_es_esi, koff_es_esi,
                 species_system):
        self.e_i_loading_rxn = e_i_loading_rxn =\
                                RevRxn(ID=ID+'_e_i_load', 
                                   reactants=[enzyme, inhibitor], 
                                   products=[ei_complex],
                                   kf=kon_ei, kb=koff_ei,
                                   species_system=species_system)
        self.ei_s_loading_rxn = ei_s_loading_rxn =\
                                RevRxn(ID=ID+'_ei_s_load', 
                                   reactants=[ei_complex, substrate], 
                                   products=[esi_complex],
                                   kf=kon_ei_esi, kb=koff_ei_esi,
                                   species_system=species_system)
        self.es_i_loading_rxn = es_i_loading_rxn =\
                                RevRxn(ID=ID+'_es_i_load', 
                                   reactants=[es_complex, inhibitor], 
                                   products=[esi_complex],
                                   kf=kon_es_esi, kb=koff_es_esi,
                                   species_system=species_system)
                                
        RxnSys.__init__(self, ID=ID, 
                           reactions=[e_i_loading_rxn, ei_s_loading_rxn, es_i_loading_rxn], 
                           species_system=species_system)
        
        self.enzyme = enzyme
        self.inhibitor = inhibitor
        self.substrate = substrate
        self.ei_complex = ei_complex
        self.es_complex = es_complex
        self.esi_complex = esi_complex
        
    @property
    def kon_ei(self):
        return self.e_i_loading_rxn.kf
    @kon_ei.setter
    def kon_ei(self, kon):
        self.e_i_loading_rxn.kf = kon
    
    @property
    def koff_ei(self):
        return self.e_i_loading_rxn.kb
    @koff_ei.setter
    def koff_ei(self, koff):
        self.e_i_loading_rxn.kb = koff
    
    @property
    def kon_ei_esi(self):
        return self.ei_s_loading_rxn.kf
    @kon_ei_esi.setter
    def kon_ei_esi(self, kon):
        self.ei_s_loading_rxn.kf = kon
    
    @property
    def koff_ei_esi(self):
        return self.ei_s_loading_rxn.kb
    @koff_ei_esi.setter
    def koff_ei_esi(self, koff):
        self.ei_s_loading_rxn.kb = koff
    
    @property
    def kon_es_esi(self):
        return self.es_i_loading_rxn.kf
    @kon_es_esi.setter
    def kon_es_esi(self, kon):
        self.es_i_loading_rxn.kf = kon
    
    @property
    def koff_es_esi(self):
        return self.es_i_loading_rxn.kb
    @koff_es_esi.setter
    def koff_es_esi(self, koff):
        self.es_i_loading_rxn.kb = koff


#%% Uncompetitive inhibition

class UnCompetitiveInhibition(RxnSys):
    """
    Class for a system of reactions involving 
    a single enzyme-substrate complex and an allosteric inhibitor.
    Consistent with uncompetitive inhibition, these reversibly
    form an enzyme-inhibitor-substrate complex.
    
    Parameters
    ----------
    es_complex : str
        Enzyme-Substrate complex.
    inhibitor : str
        Inhibitor.
    esi_complex : str
        Enzyme-Substrate-Inhibitor complex.
    kon_es_esi : float
        Forward reaction rate constant for the formation
        of enzyme-substrate-inhibitor complex from 
        enzyme-substrate complex and inhibitor.
    koff_es_esi : float
        Backward reaction rate constant for the formation
        of enzyme-substrate-inhibitor complex from 
        enzyme-substrate complex and inhibitor.
    """
    def __init__(self, 
                 ID, 
                 es_complex, 
                 inhibitor,
                 esi_complex,
                 kon_es_esi, koff_es_esi,
                 species_system):
        self.es_i_loading_rxn = es_i_loading_rxn =\
                                RevRxn(ID=ID+'_es_i_load', 
                                   reactants=[es_complex, inhibitor], 
                                   products=[esi_complex],
                                   kf=kon_es_esi, kb=koff_es_esi,
                                   species_system=species_system)
                                
        RxnSys.__init__(self, ID=ID, 
                           reactions=[es_i_loading_rxn], 
                           species_system=species_system)
        
        self.es_complex = es_complex
        self.inhibitor = inhibitor
        self.esi_complex = esi_complex
        
    @property
    def koff_es_esi(self):
        return self.es_i_loading_rxn.kb
    @koff_es_esi.setter
    def koff_es_esi(self, koff):
        self.es_i_loading_rxn.kb = koff


#%% Mechanism-based / "suicide" inhibition

class MechanismBasedInhibition(RxnSys):
    """
    Class for a system of reactions involving 
    a single enzyme, inhibitor, unstable enzyme-inhibitor
    complex, and stable enzyme-inhibitor complex (consistent
    with mechanism-based or "suicide" inhibition).
    
    Parameters
    ----------
    ID : str
        ID.
    enzyme : str
        Enzyme.
    substrate : str
        Substrate.
    ei_unstable_complex : str
        Unstable Enzyme-Inhibitor complex.
    ei_stable_complex : str
        Stable Enzyme-Inhibitor complex.
    product : str
        Product.
    kon : float
        Forward reaction rate constant for the formation
        of enzyme-substrate complex from enzyme and
        substrate.
    koff : float
        Backward reaction rate constant for the formation
        of enzyme-substrate complex from enzyme and
        substrate.
    kstabilize : float
        Rate constant for the conversion of
        unstable enzyme-inhibitor complex to 
        form stable enzyme-inhibitor complex.
    species_system : SpeciesSystem
        A SpeciesSystem object containing all species
        involved in this system of reactions.
    """
    def __init__(self, 
                 ID, 
                 enzyme, 
                 inhibitor, 
                 ei_unstable_complex,
                 ei_stable_complex,
                 kon, koff, kstabilize,
                 species_system):
        self.loading_rxn = loading_rxn = RevRxn(ID=ID+'_load', 
                                   reactants=[enzyme, inhibitor], 
                                   products=[ei_unstable_complex],
                                   kf=kon, kb=koff,
                                   species_system=species_system)
        self.stabilizing_rxn = stabilizing_rxn = Rxn(ID=ID+'_stabilize',
                           reactants=[ei_unstable_complex],
                           products=[ei_stable_complex],
                           kf=kstabilize,
                           species_system=species_system)
        
        RxnSys.__init__(self, ID=ID, 
                           reactions=[loading_rxn, stabilizing_rxn], 
                           species_system=species_system)
        
        self.enzyme = enzyme
        self.inhibitor = inhibitor
        self.ei_unstable_complex = ei_unstable_complex
        self.ei_stable_complex = ei_stable_complex
        
    @property
    def kon(self):
        return self.loading_rxn.kf
    @kon.setter
    def kon(self, kon):
        self.loading_rxn.kf = kon
    
    @property
    def koff(self):
        return self.loading_rxn.kb
    @koff.setter
    def koff(self, koff):
        self.loading_rxn.kb = koff
        
    @property
    def kstabilize(self):
        return self.catalyzed_rxn.kf
    @kstabilize.setter
    def kstabilize(self, kstabilize):
        self.catalyzed_rxn.kf = kstabilize
