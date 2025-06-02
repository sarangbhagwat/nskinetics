from nskinetics.species import SpeciesSystem
from nskinetics.reactions import RxnSys
from nskinetics.nonsteady import EnzymeSubstrateProduct,\
    CompetitiveInhibition, NonCompetitiveInhibition, UnCompetitiveInhibition,\
    MechanismBasedInhibition
from nskinetics.steady import MichaelisMenten

import numpy as np
from matplotlib import pyplot as plt 

#%% Inhibition examples initialization

# Competitive inhibition
# S0 < KM, S0 ~= E0, high KD
# -----------------
inhib_label = 'competitive'
# Initial concentrations
E_conc = 1e-4
S_conc = 1e-4
I_CI_conc = 1e-2
I_NCI_conc = 0
I_UCI_conc = 0
I_MBI_conc = 0

# Core enzyme parameters
kcat = 32.
kon = 1e4
koff = 1e1
KM = (koff+kcat)/kon

# Inhibition kinetics parameters
kon_CI = 1e4
koff_CI = 1e2
kcat_CI = 1e4

kon_ei=1e4
koff_ei=1e2
kon_es_esi=1e4
koff_es_esi=1e2
kon_ei_esi=1e4
koff_ei_esi=1e2

kon_es_esi = 1e4
koff_es_esi = 1e2

kon_MBI = 1e4
koff_MBI = 1e2
kstabilize_MBI = 32.

# Simulation parameters
t0 = 0.
tmax = 7200
dt = 5e-2
max_abs_remaining_substrate = 1e-6
max_rel_substrate_depletion = None
include_substrate_in_complex_for_max = True


# # Non-Competitive inhibition
# # S0 < KM, S0 ~= E0, high KD
# # -----------------
# inhib_label = 'non-competitive'
# # Initial concentrations
# E_conc = 1e-4
# S_conc = 1e-4
# I_CI_conc = 0
# I_NCI_conc = 1e-3
# I_UCI_conc = 0
# I_MBI_conc = 0

# # Core enzyme parameters
# kcat = 32.
# kon = 1e4
# koff = 1e1
# KM = (koff+kcat)/kon

# # Inhibition kinetics parameters
# kon_CI = 1e4
# koff_CI = 1e2
# kcat_CI = 1e4

# kon_ei=1e4
# koff_ei=1e2
# kon_es_esi=1e4
# koff_es_esi=1e2
# kon_ei_esi=1e4
# koff_ei_esi=1e2

# kon_es_esi = 1e4
# koff_es_esi = 1e2

# kon_MBI = 1e4
# koff_MBI = 1e2
# kstabilize_MBI = 32.

# # Simulation parameters
# t0 = 0.
# tmax = 7200
# dt = 5e-2
# max_abs_remaining_substrate = 1e-6
# max_rel_substrate_depletion = None
# include_substrate_in_complex_for_max = True


# # Un-Competitive inhibition
# # S0 < KM, S0 ~= E0, high KD
# # -----------------
# inhib_label = 'uncompetitive'
# # Initial concentrations
# E_conc = 1e-4
# S_conc = 1e-4
# I_CI_conc = 0
# I_NCI_conc = 0
# I_UCI_conc = 1e-3
# I_MBI_conc = 0

# # Core enzyme parameters
# kcat = 32.
# kon = 1e4
# koff = 1e1
# KM = (koff+kcat)/kon

# # Inhibition kinetics parameters
# kon_CI = 1e4
# koff_CI = 1e2
# kcat_CI = 1e4

# kon_ei=1e4
# koff_ei=1e2
# kon_es_esi=1e4
# koff_es_esi=1e2
# kon_ei_esi=1e4
# koff_ei_esi=1e2

# kon_es_esi = 1e4
# koff_es_esi = 1e2

# kon_MBI = 1e4
# koff_MBI = 1e2
# kstabilize_MBI = 32.

# # Simulation parameters
# t0 = 0.
# tmax = 7200
# dt = 5e-2
# max_abs_remaining_substrate = 1e-6
# max_rel_substrate_depletion = None
# include_substrate_in_complex_for_max = True


# # "Mechanism-Based" inhibition
# # S0 < KM, S0 ~= E0, high KD
# # -----------------
# inhib_label = '"mechanism-based"'
# # Initial concentrations
# E_conc = 1e-4
# S_conc = 1e-4
# I_CI_conc = 0
# I_NCI_conc = 0
# I_UCI_conc = 0
# I_MBI_conc = 1e-2

# # Core enzyme parameters
# kcat = 32.
# kon = 1e4
# koff = 1e1
# KM = (koff+kcat)/kon

# # Inhibition kinetics parameters
# kon_CI = 1e4
# koff_CI = 1e2
# kcat_CI = 1e4

# kon_ei=1e4
# koff_ei=1e2
# kon_es_esi=1e4
# koff_es_esi=1e2
# kon_ei_esi=1e4
# koff_ei_esi=1e2

# kon_es_esi = 1e4
# koff_es_esi = 1e2

# kon_MBI = 1e4
# koff_MBI = 1e2
# kstabilize_MBI = 32.

# # Simulation parameters
# t0 = 0.
# tmax = 7200
# dt = 5e-2
# max_abs_remaining_substrate = 1e-6
# max_rel_substrate_depletion = None
# include_substrate_in_complex_for_max = True

#%%
print(S_conc, KM)

#%% Example - Multipurpose system WITH inhibition

sp_sys = SpeciesSystem('multipurpose_sp_sys', ['E', 'S', 'ES', 'P', 
                                               'I_CI', 'EI_CI', 'Q',
                                               'I_NCI', 'EI_NCI', 'ESI_NCI',
                                               'I_UCI', 'ESI_UCI',
                                               'I_MBI', 'EI_MCI_unstable', 'EI_MCI_stable',
                                               
                                               ],
                       concentrations=[E_conc, S_conc, 0., 0., 
                                       I_CI_conc, 0., 0.,
                                       I_NCI_conc, 0., 0.,
                                       I_UCI_conc, 0.,
                                       I_MBI_conc, 0., 0.,
                                       ])

MM_rxns = EnzymeSubstrateProduct(ID='MM', enzyme='E', substrate='S', es_complex='ES', product='P',
                          kon=kon, koff=koff, kcat=kcat, 
                          species_system=sp_sys)

CI_rxns = CompetitiveInhibition(ID='CI', enzyme='E', substrate='I_CI', es_complex='EI_CI', product='Q',
                                kon=kon_CI, koff=koff_CI, kcat=kcat_CI, 
                                species_system=sp_sys)

NCI_rxns = NonCompetitiveInhibition(ID='NCI', enzyme='E', inhibitor='I_NCI',
                                    substrate='S',
                                    ei_complex='EI_NCI',
                                    es_complex='ES',
                                    esi_complex='ESI_NCI',
                                    kon_ei=kon_ei, koff_ei=koff_ei,
                                    kon_es_esi=kon_es_esi, koff_es_esi=koff_es_esi,
                                    kon_ei_esi=kon_ei_esi, koff_ei_esi=koff_ei_esi,
                                    species_system=sp_sys)

UCI_rxns = UnCompetitiveInhibition(ID='UCI', inhibitor='I_UCI', es_complex='ES',
                                   esi_complex='ESI_UCI',
                                   kon_es_esi=kon_es_esi, koff_es_esi=koff_es_esi, 
                                   species_system=sp_sys)

MBI_rxns = MechanismBasedInhibition(ID='MBI', enzyme='E', inhibitor='I_MBI',
                                    ei_unstable_complex='EI_MBI_unstable',
                                    ei_stable_complex='EI_MBI_stable',
                                    kon=kon_MBI, koff=koff_MBI,
                                    kstabilize=kstabilize_MBI,
                                    species_system=sp_sys
                                    )
multipurpose_rxn_sys = RxnSys(ID='multipurpose_rxn_sys', 
                       reactions=[
                                  MM_rxns, 
                                  CI_rxns, 
                                  NCI_rxns, 
                                  UCI_rxns,
                                  ], 
                       species_system=sp_sys)

def simulate(t0=t0, tmax=tmax, dt=dt,
             max_rel_substrate_depletion=max_rel_substrate_depletion,
             max_abs_remaining_substrate=max_abs_remaining_substrate,
             include_substrate_in_complex_for_max=include_substrate_in_complex_for_max):
    ts = [t0]
    Es = [sp_sys.concentrations[0]]
    Ss = [sp_sys.concentrations[1]]
    ESs = [sp_sys.concentrations[2]]
    Ps = [sp_sys.concentrations[3]]
    Qs = [sp_sys.concentrations[6]]
    ESI_NCIs = [sp_sys.concentrations[9]]
    
    max_rem_S = 0.
    if max_rel_substrate_depletion is not None:
        max_rem_S = (1.-max_rel_substrate_depletion)*Ss[0]
    if max_abs_remaining_substrate is not None:
        max_rem_S=max(max_rem_S, max_abs_remaining_substrate)
    
    while ts[-1]<tmax and Ss[-1]+include_substrate_in_complex_for_max*ESs[-1]>max_rem_S:
        # if ts[-1]%5==0: breakpoint()
        sp_sys.concentrations += multipurpose_rxn_sys.get_dconcs_dt()*dt
        # sp_sys.concentrations[np.where(sp_sys.concentrations<0)] = 0.
        
        Es.append(sp_sys.concentrations[0])
        Ss.append(sp_sys.concentrations[1])
        ESs.append(sp_sys.concentrations[2])
        Ps.append(sp_sys.concentrations[3])
        Qs.append(sp_sys.concentrations[6])
        ESI_NCIs.append(sp_sys.concentrations[9])
        ts.append(ts[-1] + dt)
    return np.array(ts), np.array(Es), np.array(Ss),\
           np.array(ESs), np.array(Ps), np.array(Qs),\
           np.array(ESI_NCIs)

ts, Es, Ss, ESs, Ps, Qs, ESI_NCIs = simulate()

# # enzyme material balance check
# np.allclose(np.ones(len(Es))*Es[0], Es+ESs)
# # substrate and product material balance check
# np.allclose(np.ones(len(Ss))*Ss[0], Ss+Ps)

plt.plot(ts, Ss, label=f'w/ {inhib_label} inhibition')
plt.xlabel('time [s]')
plt.ylabel('Substrate concentration [mol/L]')
plt.legend()


#%% Example - Multipurpose system WITHOUT inhibition
# Set inhibitor concentrations to zero
I_CI_conc = 0
I_NCI_conc = 0
I_UCI_conc = 0
I_MBI_conc = 0

sp_sys = SpeciesSystem('multipurpose_sp_sys', ['E', 'S', 'ES', 'P', 
                                               'I_CI', 'EI_CI', 'Q',
                                               'I_NCI', 'EI_NCI', 'ESI_NCI',
                                               'I_UCI', 'ESI_UCI',
                                               'I_MBI', 'EI_MCI_unstable', 'EI_MCI_stable',
                                               
                                               ],
                       concentrations=[E_conc, S_conc, 0., 0., 
                                       I_CI_conc, 0., 0.,
                                       I_NCI_conc, 0., 0.,
                                       I_UCI_conc, 0.,
                                       I_MBI_conc, 0., 0.,
                                       ])

MM_rxns = EnzymeSubstrateProduct(ID='MM', enzyme='E', substrate='S', es_complex='ES', product='P',
                          kon=kon, koff=koff, kcat=kcat, 
                          species_system=sp_sys)

CI_rxns = CompetitiveInhibition(ID='CI', enzyme='E', substrate='I_CI', es_complex='EI_CI', product='Q',
                                kon=kon_CI, koff=koff_CI, kcat=kcat_CI, 
                                species_system=sp_sys)

NCI_rxns = NonCompetitiveInhibition(ID='NCI', enzyme='E', inhibitor='I_NCI',
                                    substrate='S',
                                    ei_complex='EI_NCI',
                                    es_complex='ES',
                                    esi_complex='ESI_NCI',
                                    kon_ei=kon_ei, koff_ei=koff_ei,
                                    kon_es_esi=kon_es_esi, koff_es_esi=koff_es_esi,
                                    kon_ei_esi=kon_ei_esi, koff_ei_esi=koff_ei_esi,
                                    species_system=sp_sys)

UCI_rxns = UnCompetitiveInhibition(ID='UCI', inhibitor='I_UCI', es_complex='ES',
                                   esi_complex='ESI_UCI',
                                   kon_es_esi=kon_es_esi, koff_es_esi=koff_es_esi, 
                                   species_system=sp_sys)

multipurpose_rxn_sys = RxnSys(ID='multipurpose_rxn_sys', 
                       reactions=[
                                  MM_rxns, 
                                  CI_rxns, 
                                  NCI_rxns, 
                                  UCI_rxns,
                                  ], 
                       species_system=sp_sys)

def simulate(t0=t0, tmax=tmax, dt=dt,
             max_rel_substrate_depletion=max_rel_substrate_depletion,
             max_abs_remaining_substrate=max_abs_remaining_substrate,
             include_substrate_in_complex_for_max=include_substrate_in_complex_for_max):
    ts = [t0]
    Es = [sp_sys.concentrations[0]]
    Ss = [sp_sys.concentrations[1]]
    ESs = [sp_sys.concentrations[2]]
    Ps = [sp_sys.concentrations[3]]
    Qs = [sp_sys.concentrations[6]]
    ESI_NCIs = [sp_sys.concentrations[9]]
    
    max_rem_S = 0.
    if max_rel_substrate_depletion is not None:
        max_rem_S = (1.-max_rel_substrate_depletion)*Ss[0]
    if max_abs_remaining_substrate is not None:
        max_rem_S=max(max_rem_S, max_abs_remaining_substrate)
    
    while ts[-1]<tmax and Ss[-1]+include_substrate_in_complex_for_max*ESs[-1]>max_rem_S:
        # if ts[-1]%5==0: breakpoint()
        sp_sys.concentrations += multipurpose_rxn_sys.get_dconcs_dt()*dt
        # sp_sys.concentrations[np.where(sp_sys.concentrations<0)] = 0.
        
        Es.append(sp_sys.concentrations[0])
        Ss.append(sp_sys.concentrations[1])
        ESs.append(sp_sys.concentrations[2])
        Ps.append(sp_sys.concentrations[3])
        Qs.append(sp_sys.concentrations[6])
        ESI_NCIs.append(sp_sys.concentrations[9])
        ts.append(ts[-1] + dt)
    return np.array(ts), np.array(Es), np.array(Ss),\
           np.array(ESs), np.array(Ps), np.array(Qs),\
           np.array(ESI_NCIs)

ts, Es, Ss, ESs, Ps, Qs, ESI_NCIs = simulate()

# # enzyme material balance check
# np.allclose(np.ones(len(Es))*Es[0], Es+ESs)
# # substrate and product material balance check
# np.allclose(np.ones(len(Ss))*Ss[0], Ss+Ps)

plt.plot(ts, Ss, label='w/o inhibition')
plt.xlabel('time [s]')
plt.ylabel('Substrate concentration [mol/L]')
plt.legend()


