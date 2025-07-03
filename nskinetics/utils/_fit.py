# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

from sklearn.metrics import r2_score
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import numpy as np

__all__ = ('fit_multiple_dependent_variables',)

def fit_multiple_dependent_variables(f, 
                                     xdata, ydata,
                                     p0,
                                     r2_score_multioutput='uniform_average',
                                     n_minimize_runs=2,
                                     n_de_runs=5,
                                     minimize_kwargs=None,
                                     differential_evolution_kwargs=None,
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
    n_minimize_runs : int, optional
        Number of local optimization (minimization) runs to perform, each with a different
        starting point (first uses `p0`, others use DE output).
    n_de_runs : int, optional
        Number of differential evolution runs per minimization run to generate good starting
        guesses; the best DE result is used to initialize the local minimization.
    minimize_kwargs : dict
        Additional keyword arguments passed to `scipy.optimize.minimize`.
    differential_evolution_kwargs : dict
        Additional keyword arguments passed to `scipy.optimize.differential_evolution_kwargs`.
        
    Returns
    -------
    p_opt : ndarray
        Optimized parameter values that best fit the model.
    score : float
        The mean R² score achieved by the optimized parameters.
    success : bool
        Whether the finally used optimization result was flagged as successful.
        
    Notes
    -----
    - This function uses `scipy.optimize.minimize` and `scipy.optimize.differential_evolution` for optimization.
    - R² scores are computed between the predicted and observed dependent variable arrays.
    - The loss minimized is <one minus the mean R² score>.
    
    See Also
        --------
        `scipy.optimize.minimize <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`_
        `scipy.optimize.differential_evolution <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html>`_
        
    """
    
    implemented_fit_methods = ('mean r^2',)
    ydata_transpose = ydata.transpose()
    
    def load_get_obj_f(p):
        ypred = f(xdata, p).transpose()
        return 1 - r2_score(ypred, ydata_transpose,
                          multioutput=r2_score_multioutput)
    
    best_result = None
    for i in range(n_minimize_runs):
        if show_progress:
            print(f'\n\nMinimization run {i+1}:')
            print('------------------\n')
        
        best_de_result = None
        
        if i>0:
            # p0 = np.random.rand(1)*random_param_bound*np.random.rand(*p0.shape)
            for j in range(n_de_runs):
                bounds_de = differential_evolution_kwargs['bounds']
                bounds_de_to_use = []
                for (b1, b2) in bounds_de:
                    bounds_de_to_use.append((b1, b2/(10**j)))
                differential_evolution_kwargs_to_use = differential_evolution_kwargs.copy()
                differential_evolution_kwargs_to_use['bounds'] = bounds_de_to_use
                if show_progress:
                    print(f'\n\tDifferential evolution run {i+1}.{j+1}:')
                    print('\t------------------------------')
                    print('\n\t\tRunning differential evolution (DE) to get the initial guess (x0) for this minimization run ...\n')
                    print('\t\t', differential_evolution_kwargs_to_use)
                result_de = differential_evolution(load_get_obj_f,
                                                   **differential_evolution_kwargs_to_use)
                
                if show_progress:
                    print('\n\t\tDE run complete.')
                    print('\t\tres.x =', result_de.x)
                    print('\t\tR^2 =', 1. - result_de.fun)
                    print('\t\tSuccess =', result_de.success)
                
                if best_de_result is None or result_de.fun < best_de_result.fun:
                    best_de_result = result_de
        
        p0 = p0 if best_de_result is None else best_de_result.x
        best_result = best_de_result
        if show_progress:
            print('\n\tRunning minimization to get the final set of parameters ...')
            print(f'\tx0 = {p0}')
            if best_de_result is not None:
                print(f'\tInitial R^2 = {1. - best_de_result.fun}')
            print('\t\t', differential_evolution_kwargs)
        result = minimize(fun=load_get_obj_f,
                 x0=p0,
                 **minimize_kwargs)
        if show_progress:
            print('\n\tMinimization run complete.')
            print('\tres.x =', result.x)
            print('\tR^2 =', 1. - result.fun)
            print('\tSuccess =', result.success)
        if best_result is None or result.fun < best_result.fun:
            best_result = result
    
    load_get_obj_f(best_result.x)
    
    return best_result.x, 1. - best_result.fun, best_result.success
    