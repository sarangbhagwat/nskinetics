# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import sympy
from sympy import symbols, Function, Eq, dsolve, diff


# Define constants
kon = 1e2
koff = 1e1
kcat = 32.

# Define symbols
t = symbols('t', real=True)
E = Function('E')(t)
S = Function('S')(t)
ES = Function('ES')(t)
# P = Function('P')(t)

# Define the system of ODEs
eq1 = Eq(diff(E, t), -kon*E*S + (koff+kcat)*ES)
eq2 = Eq(diff(S, t), -kon*E*S + koff*ES)
eq3 = Eq(diff(ES, t), kon*E*S - (koff+kcat)*ES)
# eq4 = Eq(diff(P, t), kcat*ES)

# Initial conditions
E_conc = 1e-4
S_conc = 1e-4
ics = {E.subs(t, 0): E_conc, 
       S.subs(t, 0): S_conc,
       ES.subs(t, 0): 0.,
       # P.subs(t, 0): 0.,
       }

# Solve the system
solution = dsolve([eq1, eq2, eq3, 
                   # eq4,
                   ], 
                  [E, S, ES, 
                   # P,
                   ], ics=ics)

# Print the solution
print(solution)