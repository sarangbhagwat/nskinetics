# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import nskinetics as nsk
import numpy as np
from warnings import filterwarnings

def test_simple_ESP_inhib_predefined(plot=True, 
                                     show_progress=False,
                                     show_output=True,
                                     show_warnings=True):
    kcat = 32.
    kon = 1e2
    koff = 1e1
    
    # Default parameters for basic examples
    
    # Initial concentrations
    E_conc = 1e-3
    S_conc = 1e-1
    I_CI_conc = 1e-3
    I_NCI_conc = 0. #!!! debug NCI
    I_UCI_conc = 1e-3
    I_MBI_conc = 1e-3
    
    # Inhibition kinetics parameters
    kon_CI = 100
    koff_CI = 50
    kcat_CI = 30
    
    kon_ei=100
    koff_ei=90
    kon_es_esi_NCI=80
    koff_es_esi_NCI=70
    kon_ei_esi=75
    koff_ei_esi=60
    
    kon_es_esi_UCI = 60
    koff_es_esi_UCI = 50

    kon_MBI = 50
    koff_MBI = 20
    kstabilize_MBI = 40
    
    sp_sys = nsk.SpeciesSystem('multipurpose_sp_sys', ['E', 'S', 'ES', 'P', 
                                                   'I_CI', 'EI_CI', 'Q',
                                                   'I_NCI', 'EI_NCI', 'ESI_NCI',
                                                   'I_UCI', 'ESI_UCI',
                                                   'I_MBI', 'EI_MBI_unstable', 'EI_MBI_stable',
                                                   
                                                   ],
                           concentrations=[E_conc, S_conc, 0., 0., 
                                           I_CI_conc, 0., 0.,
                                           I_NCI_conc, 0., 0.,
                                           I_UCI_conc, 0.,
                                           I_MBI_conc, 0., 0.,
                                           ])
    
    MM_rxns = nsk.EnzymeSubstrateProduct(ID='MM', enzyme='E', substrate='S', es_complex='ES', product='P',
                              kon=kon, koff=koff, kcat=kcat, 
                              species_system=sp_sys)
    
    CI_rxns = nsk.CompetitiveInhibition(ID='CI', enzyme='E', substrate='I_CI', es_complex='EI_CI', product='Q',
                                    kon=kon_CI, koff=koff_CI, kcat=kcat_CI, 
                                    species_system=sp_sys)
    
    NCI_rxns = nsk.NonCompetitiveInhibition(ID='NCI', enzyme='E', inhibitor='I_NCI',
                                        substrate='S',
                                        ei_complex='EI_NCI',
                                        es_complex='ES',
                                        esi_complex='ESI_NCI',
                                        kon_ei=kon_ei, koff_ei=koff_ei,
                                        kon_es_esi=kon_es_esi_NCI, koff_es_esi=koff_es_esi_NCI,
                                        kon_ei_esi=kon_ei_esi, koff_ei_esi=koff_ei_esi,
                                        species_system=sp_sys)
    
    UCI_rxns = nsk.UnCompetitiveInhibition(ID='UCI', inhibitor='I_UCI', es_complex='ES',
                                       esi_complex='ESI_UCI',
                                       kon_es_esi=kon_es_esi_UCI, koff_es_esi=koff_es_esi_UCI, 
                                       species_system=sp_sys)
    
    MBI_rxns = nsk.MechanismBasedInhibition(ID='MBI', enzyme='E', inhibitor='I_MBI', 
                                            ei_unstable_complex='EI_MBI_unstable', 
                                            ei_stable_complex='EI_MBI_stable', 
                                            kon=kon_MBI, koff=koff_MBI, 
                                            kstabilize=kstabilize_MBI, 
                                            species_system=sp_sys)
    
    rxn_sys = nsk.RxnSys(ID='rxn_sys', 
                           reactions=[
                                      MM_rxns, 
                                      CI_rxns, 
                                      NCI_rxns, 
                                      UCI_rxns,
                                      MBI_rxns,
                                      ], 
                           species_system=sp_sys)
    

    # print(rxn_sys)
    # breakpoint()
    # Simulate the ReactionSystem
    
    rxn_sys.solve(t_span=[0, 10000],
                  sp_conc_for_events={'S':1e-6},)
    
    # Plot results
    if plot:
        rxn_sys.plot_solution()
    
    # Test absolute results
    assert np.allclose(rxn_sys._solution['t_events'], 
                       np.array([8382.63919662]), 
                       rtol=1e-3, atol=1e-3)
    
    assert np.allclose(rxn_sys._solution['y_events'], 
                       np.array([[4.82309338e-05],
                                 [1.00000000e-06],
                                 [4.38485153e-10],
                                 [1.00071398e-01],
                                 [8.30604981e-05],
                                 [7.85510180e-09],
                                 [9.17764773e-04],
                                 [1.00000000e-20],
                                 [1.00000000e-20],
                                 [1.00000000e-20],
                                 [9.99999745e-04],
                                 [5.26187882e-13],
                                 [4.87484236e-05],
                                 [5.59811388e-09],
                                 [9.52139030e-04]]), 
                       rtol=1e-3, atol=1e-8)
    
    # Make the same rxn_sys from strings
    
    reactions = [
                f'E + S <-> ES; kf={kon}, kb = {koff}',
                f'ES -> E + P; kf = {kcat}',
                
                f'E + I_CI <-> EI_CI; kf={kon_CI}, kb={koff_CI}',
                f'EI_CI -> E + Q; kf={kcat_CI}',
                
                f'E + I_NCI <-> EI_NCI; kf={kon_ei}, kb={koff_ei}',
                f'S + EI_NCI <-> ESI_NCI; kf={kon_ei_esi}, kb={koff_ei_esi}',
                f'ES + I_NCI <-> ESI_NCI; kf={kon_es_esi_NCI}, kb={koff_es_esi_NCI}',
                
                f'ES + I_UCI <-> ESI_UCI; kf={kon_es_esi_UCI}, kb={koff_es_esi_UCI}',
                
                f'E + I_MBI <-> EI_MBI_unstable; kf={kon_MBI}, kb={koff_MBI}',
                f'EI_MBI_unstable -> EI_MBI_stable; kf = {kstabilize_MBI}'
                ]
    
    # Generate a ReactionSystem from strings
    rxn_sys_from_str = nsk.ReactionSystem(ID='rxn_sys_from_str', 
                                     reactions=reactions,
                                     species_system=sp_sys)
    
    # print(rxn_sys_from_str)
    
    # Simulate
    sp_sys.concentrations = np.array([E_conc, S_conc, 0., 0., 
                    I_CI_conc, 0., 0.,
                    I_NCI_conc, 0., 0.,
                    I_UCI_conc, 0.,
                    I_MBI_conc, 0., 0.,
                    ])
    
    rxn_sys_from_str.solve(t_span=[0, 10000],
                  sp_conc_for_events={'S':1e-6},)
    
    # Plot results
    if plot:
        rxn_sys_from_str.plot_solution()
        
    # Test comparative results

    assert np.allclose(rxn_sys._solution['t_events'], 
                       rxn_sys_from_str._solution['t_events'], 
                       rtol=1e-3, atol=1e-3)
    
    assert np.allclose(rxn_sys._solution['y_events'], 
                       rxn_sys_from_str._solution['y_events'], 
                       rtol=1e-3, atol=1e-8)
    
    # Test comparative reaction strings
    for (r1, r2) in zip(rxn_sys.reactions_flattened, rxn_sys_from_str.reactions_flattened):
        r1_str = r1.__str__()
        r2_str = r2.__str__()
        assert r1_str[r1_str.index(':')+1:] == r2_str[r2_str.index(':')+1:]
        
    return rxn_sys