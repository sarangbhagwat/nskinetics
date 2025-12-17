# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import nskinetics as nsk
import biosteam as bst
import tellurium as te
import simplesbml

__all__ = ('R302',)

#%% Import SBML
# filename = 'BIOMD0000000245_url.xml' # Lei et al., 2001
# filepath = f'{filename}'
# model = simplesbml.loadSBMLFile(filepath)

sbml_url = 'https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1003250000/3/BIOMD0000000245_url.xml' # Lei et al., 2001
r = te.loadSBMLModel(sbml_url)
model = simplesbml.loadSBMLStr(r.getSBML())

#%% Replace biomass growth reaction(s) to include ethanol and acetate inhibition

replace_BM_growth_gluc = True
replace_BM_growth_acetate = False

### --- 1. Biomass growth on glucose --- ###
if replace_BM_growth_gluc:
    # remove existing reaction
    original_index_BM_growth_gluc = model.getListOfReactionIds().index('r7')
    assert 'BM growth (gluc)' in str(model.model.getReaction(original_index_BM_growth_gluc))
    model.model.removeReaction(original_index_BM_growth_gluc)
    
    # add new reaction and associated parameters
    # ethanol inhibition coefficient
    model.addParameter(param_id='k_7ie', val=0.040, # Schneiderman et al., 2015 https://doi.org/10.1016/j.biortech.2014.11.087
                       units='g_per_l_per_h', # g_per_g_per_h per Lei et al., but SBML model has g_per_l_per_h
                       ) 
    
    # acetate inhibition coefficient
    model.addParameter(param_id='k_7ia', val=0.120, # Schneiderman et al., 2015 https://doi.org/10.1016/j.biortech.2014.11.087
                       units='g_per_l_per_h', # g_per_g_per_h per Lei et al., but SBML model has g_per_l_per_h
                       ) 
    
    model.addReaction(rxn_id='r7', 
                      reactants=['s_glu', 
                                 'x', 's_EtOH', 's_acetate', # added modifiers to both reactants and products since modifiers can't be passed as argument
                                 ],
                      products=['0.732 x', '0.127 CO2', '0.063 Red', 
                                'x', 's_EtOH', 's_acetate',  # added modifiers to both reactants and products since modifiers can't be passed as argument
                                ],
                      expression='k_7 * s_glu * (2.71828 ^ (- k_7ie * s_EtOH)) * (2.71828 ^ (- k_7ia * s_acetate)) / (s_glu + K_7) * x * X_a * env',
                      # modifiers=['x', 's_EtOH'], # can't define here
                      )
### ------------------------------------ ###

### --- 2. Biomass growth on acetate --- ###
if replace_BM_growth_acetate:
    # remove existing reaction
    original_index_BM_growth_acetate = model.getListOfReactionIds().index('r8')
    assert 'BM growth (acetate)' in str(model.model.getReaction(original_index_BM_growth_acetate))
    model.model.removeReaction(original_index_BM_growth_acetate)
    
    # add new reaction and associated parameters
    # ethanol inhibition coefficient
    model.addParameter(param_id='k_8ie', val=0.020, # Schneiderman et al., 2015 https://doi.org/10.1016/j.biortech.2014.11.087
                       units='g_per_l_per_h', # g_per_g_per_h per Lei et al., but SBML model has g_per_l_per_h
                       ) 
    
    # acetate inhibition coefficient
    model.addParameter(param_id='k_8ia', val=0.120, # Schneiderman et al., 2015 https://doi.org/10.1016/j.biortech.2014.11.087
                       units='g_per_l_per_h', # g_per_g_per_h per Lei et al., but SBML model has g_per_l_per_h
                       ) 
    
    model.addReaction(rxn_id='r8', 
                      reactants=['s_acetate',
                                 'x', 's_glu', 's_EtOH', # added modifiers to both reactants and products since modifiers can't be passed as argument
                                 ],
                      products=['0.619 x', '0.325 CO2', '0.214 Red',
                                'x', 's_glu', 's_EtOH', # added modifiers to both reactants and products since modifiers can't be passed as argument
                                ],
                      expression='k_8 * s_acetate * (2.71828 ^ (- k_8ie * s_EtOH)) * (2.71828 ^ (- k_8ia * s_acetate)) / ((s_acetate + K_5e) * (1 + K_5i * s_glu)) * x * X_a * env',
                      # modifiers=['x', 's_glu', 's_EtOH'], # can't define here
                      )
### ------------------------------------ ###

#%% Replace metabolic reaction(s) to include ethanol and acetate inhibition

replace_glycolysis = True
replace_Acetald_dehydrogenase = True
replace_ADH = False # !!! not implemented

### --- 1. Glycolysis --- ###
if replace_glycolysis:
    # remove existing reaction
    original_index_glycolysis = model.getListOfReactionIds().index('r1')
    assert 'glycolysis' in str(model.model.getReaction(original_index_glycolysis))
    model.model.removeReaction(original_index_glycolysis)
    
    # add new reaction and associated parameters
    # ethanol inhibition coefficient
    model.addParameter(param_id='k_1ie', val=0.040, # Schneiderman et al., 2015 https://doi.org/10.1016/j.biortech.2014.11.087
                       units='g_per_l_per_h', # g_per_g_per_h per Lei et al., but SBML model has g_per_l_per_h
                       ) 
    
    # acetate inhibition coefficient
    model.addParameter(param_id='k_1ia', val=0.120, # Schneiderman et al., 2015 https://doi.org/10.1016/j.biortech.2014.11.087
                       units='g_per_l_per_h', # g_per_g_per_h per Lei et al., but SBML model has g_per_l_per_h
                       ) 
    
    model.addReaction(rxn_id='r1', 
                      reactants=['s_glu', 
                                 's_acetald', 'x', 's_EtOH', 's_acetate', # added modifiers to both reactants and products since modifiers can't be passed as argument
                                 ],
                      products=['0.978 s_pyr', '0.178 Red', 
                                's_acetald', 'x', 's_EtOH', 's_acetate',  # added modifiers to both reactants and products since modifiers can't be passed as argument
                                ],
                      expression='(k_1l * s_glu / (s_glu + K_1l) + k_1h * s_glu / (s_glu + K_1h) + k_1e * s_acetald * s_glu / (s_glu * (K_1i * s_acetald + 1) + K_1e)) * x * X_a * env * (2.71828 ^ (- k_1ie * s_EtOH)) * (2.71828 ^ (- k_1ia * s_acetate))',
                      # modifiers=['x', 's_EtOH'], # can't define here
                      )
### ------------------------------------ ###


### --- 2. Acetate production by Acetaldehyde dehydrogenase --- ###
if replace_Acetald_dehydrogenase:
    # remove existing reaction
    original_index_Acetald_dehydrogenase = model.getListOfReactionIds().index('r4')
    assert 'Acetald. dehydrogenase' in str(model.model.getReaction(original_index_Acetald_dehydrogenase))
    model.model.removeReaction(original_index_Acetald_dehydrogenase)
    
    # add new reaction and associated parameters
    # ethanol inhibition coefficient
    model.addParameter(param_id='k_4ie', val=0.040, # Schneiderman et al., 2015 https://doi.org/10.1016/j.biortech.2014.11.087
                       units='g_per_l_per_h', # g_per_g_per_h per Lei et al., but SBML model has g_per_l_per_h
                       ) 
    
    # acetate inhibition coefficient
    model.addParameter(param_id='k_4ia', val=0.120, # Schneiderman et al., 2015 https://doi.org/10.1016/j.biortech.2014.11.087
                       units='g_per_l_per_h', # g_per_g_per_h per Lei et al., but SBML model has g_per_l_per_h
                       ) 
    
    model.addReaction(rxn_id='r4', 
                      reactants=['s_acetald', 
                                 'x', 's_EtOH', 's_acetate', # added modifiers to both reactants and products since modifiers can't be passed as argument
                                 ],
                      products=['1.363 s_acetate', '0.363 Red', 
                                'x', 's_EtOH', 's_acetate',# added modifiers to both reactants and products since modifiers can't be passed as argument
                                ],
                      expression='k_4 * s_acetald * (2.71828 ^ (- k_4ie * s_EtOH)) * (2.71828 ^ (- k_4ia * s_acetate)) / (s_acetald + K_4) * x * X_a * X_AcDH * env'
                      # modifiers=['x', 's_EtOH'], # can't define here
                      )
    
### ------------------------------------ ###

### --- 3. Ethanol production by Alcohol Dehydrogenase --- ###
if replace_ADH:
    pass # !!! not implemented
### ------------------------------------ ###

#%% Replace cell deactivation reaction(s) 
replace_active_BM_degradation = True
### --- 1. Degradation of active biomass --- ###
if replace_active_BM_degradation:
    # remove existing reaction
    original_index_active_BM_degradation = model.getListOfReactionIds().index('r10')
    assert 'active BM degradation' in str(model.model.getReaction(original_index_active_BM_degradation))
    model.model.removeReaction(original_index_active_BM_degradation)
    
    # remove existing parameters
    model.model.removeParameter('k_10')
    model.model.removeParameter('K_10')
    model.model.removeParameter('k_10e')
    model.model.removeParameter('K_10e')
    
    # add new reaction and associated parameters
    # ethanol inhibition coefficient
    model.addParameter(param_id='k_10', val=0.06, # Atitallah et al., 2020 https://doi.org/10.1016/j.renene.2020.03.010
                       units='g_per_l_per_h',
                       ) 
    
    
    model.addReaction(rxn_id='r10', 
                      reactants=['a',
                                 # 'x', # added modifiers to both reactants and products since modifiers can't be passed as argument
                                 ],
                      products=[
                              # 'x', # added modifiers to both reactants and products since modifiers can't be passed as argument
                              ],
                      expression='k_10 * x * X_a * env'
                      )
    
### ------------------------------------ ###

#%% Add glucose feed spike event

model.addParameter(param_id='n_glu_spikes', val=0, 
                   units='dimensionless',
                   ) 

model.addEvent(event_id='glucose_feed_spike', 
               trigger='s_glu <= -1 && n_glu_spikes<2', 
               assignments={'s_glu': '100', "n_glu_spikes": "n_glu_spikes + 1"}, )

def update_threshold_glucose_conc_for_feed_spike(model, new_val):
    event_str_identifier_1 = '<event id="glucose_feed_spike" useValuesFromTriggerTime="true">\n        <trigger>\n          <math xmlns="http://www.w3.org/1998/Math/MathML">\n            <apply>\n              <leq/>\n              <ci> s_glu </ci>\n              <cn type="integer"> '
    sbml_str = model.toSBML()
    old_val_str = sbml_str[len(event_str_identifier_1) + sbml_str.index(event_str_identifier_1)]
    old_event_partial_str = f'<event id="glucose_feed_spike" useValuesFromTriggerTime="true">\n        <trigger>\n          <math xmlns="http://www.w3.org/1998/Math/MathML">\n            <apply>\n              <leq/>\n              <ci> s_glu </ci>\n              <cn type="integer"> {old_val_str} </cn>\n            </apply>'
    new_event_partial_str = f'<event id="glucose_feed_spike" useValuesFromTriggerTime="true">\n        <trigger>\n          <math xmlns="http://www.w3.org/1998/Math/MathML">\n            <apply>\n              <leq/>\n              <ci> s_glu </ci>\n              <cn type="integer"> {new_val} </cn>\n            </apply>'
    new_sbml_str = sbml_str.replace(old_event_partial_str, new_event_partial_str)
    return simplesbml.loadSBMLStr(new_sbml_str)

# model = update_threshold_glucose_conc_for_feed_spike(model, 10) # -1 => no spiking

#%% Load updated model SBML in Tellurium
r = te.loadSBMLModel(model.toSBML())

#%% Update parameter values through Tellurium model
r.D = 0 # dilution rate # note r.S_f (inlet glucose conc) does not matter when D=0 

# r.k_6 = 2.82 # ADH ethanol production

# r.k_2 = 0 # TCA pyruvate
# r.k_5 = 0 # TCA acetate
# r.k_5e = 0 # also TCA acetate

# r.k_6r = 0 # use of ethanol as substrate # keep this, as possible even in anaerobic conditions

# r.k_10 = 0 # glucose-based active biomass degradation
# r.K_10 = 102.5 # glucose-based active biomass degradation # updated based on Atitallah et al., 2020 https://doi.org/10.1016/j.renene.2020.03.010

# r.k_10e = 0 # ethanol-based active biomass degradation
# r.K_10e = 42.2 # ethanol-based active biomass degradation # updated based on Atitallah et al., 2020 https://doi.org/10.1016/j.renene.2020.03.010

#%%  simulations
r.reset()
r.s_glu = 100 # initial glucose conc
r.x = 2 # initial biomass conc
print(r.s_glu, r.x, r.X_a, r.s_EtOH, r.s_acetate, r.s_acetald)

r.simulate(0, 300, 300, 
           # ['time', 'X_a', 'X_AcDH', 
            # 'a',
            # ],
           )

print(r.s_glu, r.x, r.X_a, r.s_EtOH, r.s_acetate, r.s_acetald)
r.plot()

#%% Create NSKinetics Reaction System
te_r = nsk.TelluriumReactionSystem(r)
te_r._units['time'] = 'h'
te_r._units['conc'] = 'g/L'

#%%

bst.settings.set_thermo(['Water', 'Carbon', 'Glucose', 'Ethanol', 'Sucrose'])

feed = bst.Stream('feed', Water=100, Carbon=1, Glucose=1)

map_chemicals_nsk_to_bst = {'s_glu': 'Glucose',
                            'x': 'Carbon',
                            's_EtOH': 'Ethanol',}

def reset_kinetic_reaction_system(r):
    r.reset()
    r._te.n_glu_spikes = 0
    
R302 = nsk.units.NSKFermentation('R302', 
                                 ins=feed, 
                                 kinetic_reaction_system=te_r,
                                 map_chemicals_nsk_to_bst=map_chemicals_nsk_to_bst,
                                 n_simulation_steps=None,
                                 f_reset_kinetic_reaction_system=reset_kinetic_reaction_system,
                                 tau=3*24)

R302.simulate()

