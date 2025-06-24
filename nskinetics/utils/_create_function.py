# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

__all__ = ('create_function', 
           'codify')

def codify(statement):
    statement = replace_apostrophes(statement)
    statement = replace_newline(statement)
    return statement

def replace_newline(statement):
    statement = statement.replace('\n', ';')
    return statement

def replace_apostrophes(statement):
    statement = statement.replace('’', "'").replace('‘', "'").replace('“', '"').replace('”', '"')
    return statement

def create_function(code, namespace):
    def wrapper_fn(statement):
        def f(t, concs):
            namespace['t'] = t
            namespace['concs'] = concs
            exec(codify(statement), namespace)
            return namespace['y']
        return f
    function = wrapper_fn(code)
    return function
