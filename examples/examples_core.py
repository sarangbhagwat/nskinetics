from nskinetics import SpeciesSystem
from nskinetics import RxnSys
from nskinetics import EnzymeSubstrateProduct,\
    CompetitiveInhibition, NonCompetitiveInhibition, UnCompetitiveInhibition
from nskinetics import MichaelisMenten

import numpy as np
from matplotlib import pyplot as plt 

#%% Core examples initialization

# # S0 < KM, S0 ~= E0, high KD
# # -----------------
# # Initial concentrations
# E_conc = 1e-4
# S_conc = 1e-4

# # Core enzyme parameters
# kcat = 32.
# kon = 1e4
# koff = 1e1
# KM = (koff+kcat)/kon

# # Simulation parameters
# t0 = 0.
# tmax = 7200
# dt = 5e-2
# max_abs_remaining_substrate = 1e-6
# max_rel_substrate_depletion = None
# include_substrate_in_complex_for_max = True

# # S0 >> KM, S0 > E0, high KD
# # -----------------
# # Initial concentrations
# E_conc = 1e-2
# S_conc = 1e-1

# # Core enzyme parameters
# kcat = 32.
# kon = 1e4
# koff = 1e1
# KM = (koff+kcat)/kon

# # Simulation parameters
# t0 = 0.
# tmax = 500
# dt = 5e-2
# max_abs_remaining_substrate = 1e-6
# max_rel_substrate_depletion = None
# include_substrate_in_complex_for_max = True

# S0 >> KM, S0 >> E0, high KD
# -----------------
# Initial concentrations
E_conc = 1e-4
S_conc = 1e-4

# Core enzyme parameters
kcat = 32.
kon = 12
koff = 1e1
KM = (koff+kcat)/kon

# Simulation parameters
t0 = 0.
tmax = 30*3600
dt = 0.1
max_abs_remaining_substrate = 1e-6
max_rel_substrate_depletion = None
include_substrate_in_complex_for_max = True

#%% Default parameters for basic examples

# Initial concentrations
I_CI_conc = 0
I_NCI_conc = 0
I_UCI_conc = 0
I_MBI_conc = 0

# Inhibition kinetics parameters
kon_CI = 0
koff_CI = 0
kcat_CI = 0


kon_ei=0
koff_ei=0
kon_es_esi=0
koff_es_esi=0
kon_ei_esi=0
koff_ei_esi=0

kon_es_esi = 0
koff_es_esi = 0

#%%
print(S_conc, KM)

#%% Example - Multipurpose system

sp_sys = SpeciesSystem('multipurpose_sp_sys', ['E', 'S', 'ES', 'P', 
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

plt.plot(ts, Ss, label='nonsteady', linestyle='dashed',
         zorder=2, color='blue', 
         # alpha=0.4,
         )
plt.ylim(0, 1.1*S_conc)
plt.xlim(0, ts[-1])
plt.xlabel('time [s]')
plt.ylabel('Substrate concentration [mol/L]')
plt.legend()


#%% Example - STEADY state Michaelis Menten used for non-steady system

SS_MM_sp_sys = SpeciesSystem('SS_M', ['E', 'S','P',],
                       concentrations=[E_conc, S_conc, 0.,
                                       ])

SS_MM_rxns = MichaelisMenten(ID='SS_MM', 
                                        enzyme='E', substrate='S',product='P',
                                        kcat=kcat, KM=KM,
                                        species_system=SS_MM_sp_sys)

def simulate_SS_MM(t0=t0, tmax=tmax, dt=dt, 
             max_rel_substrate_depletion=max_rel_substrate_depletion,
             max_abs_remaining_substrate=max_abs_remaining_substrate):
    ts = [t0]
    Es = [SS_MM_sp_sys.concentrations[0]]
    Ss = [SS_MM_sp_sys.concentrations[1]]
    Ps = [SS_MM_sp_sys.concentrations[2]]
    
    max_rem_S = 0.
    if max_rel_substrate_depletion is not None:
        max_rem_S = (1.-max_rel_substrate_depletion)*Ss[0]
    if max_abs_remaining_substrate is not None:
        max_rem_S=max(max_rem_S, max_abs_remaining_substrate)
    
    while ts[-1]<tmax and Ss[-1]>max_rem_S:
        # if ts[-1]%5==0: breakpoint()
        SS_MM_sp_sys.concentrations += SS_MM_rxns.get_dconcs_dt()*dt
        # SS_MM_sp_sys.concentrations[np.where(SS_MM_sp_sys.concentrations<0)] = 0.
        
        Es.append(SS_MM_sp_sys.concentrations[0])
        Ss.append(SS_MM_sp_sys.concentrations[1])
        Ps.append(SS_MM_sp_sys.concentrations[2])
        ts.append(ts[-1] + dt)
    return ts, Es, Ss, Ps

ts_steady, Es_steady, Ss_steady, Ps_steady = simulate_SS_MM()

plt.plot(ts_steady, Ss_steady, label='steady', linestyle='dashed',
         zorder=2, color='orange')
plt.ylim(0, 1.1*S_conc)
plt.xlim(0, ts[-1])
plt.xlabel('time [s]')
plt.ylabel('Substrate concentration [mol/L]')

#%%
from scipy.integrate import solve_ivp

sp_sys.concentrations = np.array([E_conc, S_conc, 0., 0., 
                I_CI_conc, 0., 0.,
                I_NCI_conc, 0., 0.,
                I_UCI_conc, 0.,
                I_MBI_conc, 0., 0.,
                ])
    
def ode_system_RHS(t, concs):
    concs[np.where(concs<0)] = 0. #  not needed with a low enough atol
    sp_sys.concentrations = concs
    # dconcs_dt = MM_rxns.get_dconcs_dt()
    # print(t)
    # print(concs[:4])
    # print(dconcs_dt[:4])
    # print('\n')
    # if np.any(concs<0): breakpoint()
    return multipurpose_rxn_sys.get_dconcs_dt()

def S_obj_f(t, concs, S=max_abs_remaining_substrate):
    return concs[1]-S

sol = solve_ivp(ode_system_RHS, 
                t_span=[t0, tmax], 
                y0=sp_sys.concentrations,
                # t_eval=np.arange(0., 5.1, 0.1)*1000,
                # t_eval=[1000, 5000],
                atol=1e-12, # <= 1e-6*max(sp_sys.concentrations)
                rtol=1e-6, # 1e-6
                # the solver keeps the local error estimates less than atol + rtol * abs(y)
                events=S_obj_f,
                # method='LSODA',
                dense_output=False)

plt.plot(sol.t, sol.y[1, :], label='nonsteady_ode_RK', linestyle='solid',
         zorder=2, color='gray', 
         # linewidth=0.5,
         )
plt.ylim(0, 1.1*S_conc)
plt.xlim(0, ts[-1])
plt.xlabel('time [s]')
plt.ylabel('Substrate concentration [mol/L]')

#%%
from scipy.integrate import solve_ivp

sp_sys.concentrations = np.array([E_conc, S_conc, 0., 0., 
                I_CI_conc, 0., 0.,
                I_NCI_conc, 0., 0.,
                I_UCI_conc, 0.,
                I_MBI_conc, 0., 0.,
                ])
    
# def ode_system_RHS(t, concs):
#     concs[np.where(concs<0)] = 0. #  not needed with a low enough atol
#     sp_sys.concentrations = concs
#     # dconcs_dt = MM_rxns.get_dconcs_dt()
#     # print(t)
#     # print(concs[:4])
#     # print(dconcs_dt[:4])
#     # print('\n')
#     # if np.any(concs<0): breakpoint()
#     return multipurpose_rxn_sys.get_dconcs_dt()

# def S_obj_f(t, concs, S=max_abs_remaining_substrate):
#     return concs[1]-S

# sol = solve_ivp(ode_system_RHS, 
#                 t_span=[t0, tmax], 
#                 y0=sp_sys.concentrations,
#                 # t_eval=np.arange(0., 5.1, 0.1)*1000,
#                 # t_eval=[1000, 5000],
#                 atol=1e-12, # <= 1e-6*max(sp_sys.concentrations)
#                 rtol=1e-6, # 1e-6
#                 # the solver keeps the local error estimates less than atol + rtol * abs(y)
#                 events=S_obj_f,
#                 method='LSODA',
#                 dense_output=False)

sol = multipurpose_rxn_sys.solve(t_span=[t0, tmax], 
                                 sp_conc_for_events={'S':1e-6})

plt.plot(sol.t, sol.y[1, :], label='nonsteady_ode_BDF-stiffness', linestyle='solid',
         zorder=2, color='green', 
         # linewidth=0.5,
         )
plt.ylim(0, 1.1*S_conc)
plt.xlim(0, ts[-1])
plt.xlabel('time [s]')
plt.ylabel('Substrate concentration [mol/L]')

plt.legend()
plt.show()

#%%
def f_ode_RK_nonsteady():
    sp_sys.concentrations = np.array([E_conc, S_conc, 0., 0., 
                    I_CI_conc, 0., 0.,
                    I_NCI_conc, 0., 0.,
                    I_UCI_conc, 0.,
                    I_MBI_conc, 0., 0.,
                    ])
    
    sol = solve_ivp(ode_system_RHS, 
                    t_span=[t0, tmax], 
                    y0=sp_sys.concentrations,
                    # t_eval=np.arange(0., 5.1, 0.1)*1000,
                    # t_eval=[1000, 5000],
                    atol=1e-12, # <= 1e-6*max(sp_sys.concentrations)
                    rtol=1e-6, # 1e-6
                    # the solver keeps the local error estimates less than atol + rtol * abs(y)
                    events=S_obj_f,
                    # method='LSODA',
                    # method='BDF',
                    dense_output=False)

def f_ode_BDF_S_nonsteady():
    sp_sys.concentrations = np.array([E_conc, S_conc, 0., 0., 
                    I_CI_conc, 0., 0.,
                    I_NCI_conc, 0., 0.,
                    I_UCI_conc, 0.,
                    I_MBI_conc, 0., 0.,
                    ])
    
    sol = solve_ivp(ode_system_RHS, 
                    t_span=[t0, tmax], 
                    y0=sp_sys.concentrations,
                    # t_eval=np.arange(0., 5.1, 0.1)*1000,
                    # t_eval=[1000, 5000],
                    atol=1e-12, # <= 1e-6*max(sp_sys.concentrations)
                    rtol=1e-6, # 1e-6
                    # the solver keeps the local error estimates less than atol + rtol * abs(y)
                    events=S_obj_f,
                    method='LSODA',
                    # method='BDF',
                    dense_output=False)
    
def f_sim_nonsteady():
    sp_sys.concentrations = np.array([E_conc, S_conc, 0., 0., 
                    I_CI_conc, 0., 0.,
                    I_NCI_conc, 0., 0.,
                    I_UCI_conc, 0.,
                    I_MBI_conc, 0., 0.,
                    ])
    simulate()

def f_sim_steady():
    sp_sys.concentrations = np.array([E_conc, S_conc, 0., 0., 
                    I_CI_conc, 0., 0.,
                    I_NCI_conc, 0., 0.,
                    I_UCI_conc, 0.,
                    I_MBI_conc, 0., 0.,
                    ])
    simulate_SS_MM()
