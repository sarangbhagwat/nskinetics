# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

from sklearn.metrics import r2_score
from scipy.optimize import minimize
import numpy as np

__all__ = ('fit_multiple_dependent_variables',)

def fit_multiple_dependent_variables(f, 
                                     xdata, ydata,
                                     p0,
                                     fit_method='mean R^2',
                                     r2_score_multioutput='uniform_average',
                                     n_minimize_runs=2,
                                     random_param_bound=1000.,
                                     show_progress=False,
                                     **kwargs):
    """
    Fit a model function to multiple dependent variables using a shared set of parameters.
    
    Parameters
    ----------
    f : callable
        A model function that takes two arguments: an array of independent variables `xdata`
        and a parameter array `p`, and returns an array of shape (N, M) or (M, N), where N is
        the number of dependent variables and M is the number of data points.
    xdata : array-like
        A 1D array of independent variable values of shape (M,).
    ydata : array-like
        A 2D array of observed dependent variable values of shape (N, M), where each row corresponds
        to a dependent variable.
    p0 : array-like, optional
        Initial guess for the parameters to be optimized.
        Only used in the first minimization run (thereafter random).
    fit_method : str, optional
        Method used for fitting. Currently, only 'mean R^2' is implemented (case-insensitive), which
        minimizes the negative mean R² score across all dependent variables.
    **kwargs : dict
        Additional keyword arguments passed to `scipy.optimize.minimize`.
        
    Returns
    -------
    p_opt : ndarray
        Optimized parameter values that best fit the model.
    score : float
        The mean R² score achieved by the optimized parameters.
    success : bool
        Whether the optimization was successful.
        
    Raises
    ------
    ValueError
        If an unsupported `fit_method` is provided.
        
    Notes
    -----
    - The function uses `scipy.optimize.minimize` for optimization.
    - R² scores are computed between the predicted and observed dependent variable arrays.
    - The loss minimized is the negative mean R² score.
    
    """
    
    implemented_fit_methods = ('mean r^2',)
    ydata_transpose = ydata.transpose()
    
    if fit_method.lower()=='mean r^2':
        
        def load_get_mean_r2_score(p):
            ypred = f(xdata, p).transpose()
            return 1 - r2_score(ypred, ydata_transpose,
                              multioutput=r2_score_multioutput)
        
        best_result = None
        for i in range(n_minimize_runs):
            if i>0:
                p0 = np.random.rand(1)*random_param_bound*np.random.rand(*p0.shape)
            
            if show_progress:
                print(f'\n\nOptimization run {i+1}:')
                print('------------------\n')
            result = minimize(fun=load_get_mean_r2_score,
                     x0=p0,
                     **kwargs)
            if show_progress:
                print('res.x =', result.x)
                print('R^2 =', 1. - result.fun)
                print('Success =', result.success)
            if best_result is None or result.fun < best_result.fun:
                best_result = result
        
        load_get_mean_r2_score(best_result.x)
        
        return best_result.x, 1. - best_result.fun, best_result.success
    
    else:
        raise ValueError(f'Method {fit_method} not implemented; must be one of {implemented_fit_methods}\n')