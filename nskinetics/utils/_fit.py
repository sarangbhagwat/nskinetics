# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np
from scipy.optimize import curve_fit

__all__ = ('fit_vector_output_model',)

def fit_vector_output_model(func, x_data, y_data, 
                            p0=None, 
                            bounds=(-np.inf, np.inf),
                            **kwargs):
    """
    Fit a model that returns a vector-valued output to multivariate y_data.

    Parameters:
    -----------
    func : callable
        A function f(x, *params) that returns a vector of length N for each scalar x.
    x_data : array_like, shape (M,)
        Independent variable values.
    y_data : array_like, shape (M, N)
        Observed dependent variable values.
    p0 : array_like, optional
        Initial guess for parameters.
    bounds : 2-tuple of array_like, optional
        Lower and upper bounds for parameters.

    Returns:
    --------
    popt : array
        Optimal parameters.
    pcov : 2D array
        Covariance matrix of the parameters.
    """
    x_data = np.asarray(x_data)
    y_data = np.asarray(y_data)

    M, N = y_data.shape
    if x_data.shape[0] != M:
        raise ValueError("x_data must have the same number of entries as rows in y_data")

    def wrapped_func(x_repeat, *params):
        # x_repeat is flattened x_data repeated N times for curve_fit
        x_vals = np.asarray(x_repeat)
        x_vals = x_vals.reshape(-1)  # (M * N,)
        outputs = [func(x_data[i], *params) for i in range(M)]  # shape (M, N)
        return np.array(outputs).reshape(-1)  # Flatten to (M*N,)

    # Flatten y_data to a vector
    y_flat = y_data.reshape(-1)

    # Create a dummy x vector to match shape (M*N,)
    x_dummy = np.repeat(x_data, y_data.shape[1])

    popt, pcov = curve_fit(wrapped_func, x_dummy, y_flat, p0=p0, bounds=bounds)
    return popt, pcov