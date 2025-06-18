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


#%%
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
                                 sp_conc_for_events={'S':max_abs_remaining_substrate})

#%%
# plt.plot(sol.t, sol.y[1, :], label='nonsteady_ODE_LSODA', linestyle='solid',
#          zorder=2, color='green', 
#          # linewidth=0.5,
#          )
# plt.ylim(0, 1.1*S_conc)
# # plt.xlim(0, ts[-1])
# plt.xlabel('time [s]')
# plt.ylabel('Substrate concentration [mol/L]')

# plt.legend()
# plt.show()
multipurpose_rxn_sys.plot_solution(sps_to_include=['E', 'S', 'ES', 'P'])
multipurpose_rxn_sys.plot_solution(sps_to_include=['E'])
multipurpose_rxn_sys.plot_solution(sps_to_include=['ES']) 
                                     
#%%
def f():
    sp_sys.concentrations = np.array([E_conc, S_conc, 0., 0., 
                I_CI_conc, 0., 0.,
                I_NCI_conc, 0., 0.,
                I_UCI_conc, 0.,
                I_MBI_conc, 0., 0.,
                ])
    sol = multipurpose_rxn_sys.solve(t_span=[t0, tmax], 
                                 sp_conc_for_events={'S':1e-6})