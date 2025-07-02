# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

from sklearn.metrics import r2_score
from scipy.optimize import minimize

__all__ = ('fit_multiple_dependent_variables',)

def fit_multiple_dependent_variables(f, 
                                     xdata, ydata,
                                     p0=None,
                                     fit_method='mean R^2',
                                     **kwargs):
    
    implemented_fit_methods = ('mean r^2',)
    
    if fit_method.lower()=='mean r^2':
        
        def load_get_mean_r2_score(p):
            ypred = f(xdata, p)
            return - r2_score(ypred, ydata)
        
        res = minimize(fun=load_get_mean_r2_score,
                 x0=p0,
                 **kwargs)
        
        return res.x, load_get_mean_r2_score(res.x), res.success
    
    else:
        raise ValueError(f'Method {fit_method} not implemented; must be one of {implemented_fit_methods}\n')