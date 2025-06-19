from nskinetics import SpeciesSystem
from nskinetics import RxnSys
from nskinetics import EnzymeSubstrateProduct, ReactionSystem,\
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



#%%
print(S_conc, KM)

#%% Example - Multipurpose system

sp_sys = SpeciesSystem('multipurpose_sp_sys', ['E', 'S', 'ES', 'P', 
                                               ],
                       concentrations=[E_conc, S_conc, 0., 0., 
                                       ])

reactions = reactions = [
    # EnzymeSubstrateProduct
    'E + S <-> ES; kf = 12.0, kb = 10.0',
    'ES -> E + P; kf = 32.0'
    ]

MM_rxns = ReactionSystem(ID='MM', 
                         reactions=reactions,
                          species_system=sp_sys)



#%%
sp_sys.concentrations = np.array([E_conc, S_conc, 0., 0., 
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

sol = MM_rxns.solve(t_span=[t0, tmax], 
                                 sp_conc_for_events={'S':max_abs_remaining_substrate})

#%%
MM_rxns.plot_solution(sps_to_include=['E', 'S', 'ES', 'P',])
MM_rxns.plot_solution(sps_to_include=['E'])
MM_rxns.plot_solution(sps_to_include=['ES']) 
