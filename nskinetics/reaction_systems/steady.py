# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import numpy as np
from numba import njit
from scipy.optimize import curve_fit
from pandas import DataFrame, read_excel
from typing import Union

__all__ = ('michaelis_menten', 'fit_michaelis_menten')

#%% Michaelis Menten (STEADY-state only)

@njit(cache=True)
def michaelis_menten(S, Vmax, KM):
    return Vmax * S / (KM + S)

def fit_michaelis_menten(data: Union[DataFrame, dict, str],
                         E0: float,
                         s_col: str = 'S',
                         v_col: str = 'V',
                         sheet_name: Union[str, int] = 0,
                         normalize=False,
                         S_data=None,
                         V_data=None):
    """
    Fit Michaelis-Menten parameters (kcat, KM) from data in a DataFrame, 
    dict, or Excel file.
    
    Parameters:
    -----------
    data : DataFrame, dict, or str
        Input data:
        - pandas DataFrame with columns s_col and v_col
        - dict with keys s_col and v_col
        - Excel file path with columns s_col and v_col
    
    E0 : float
        Total enzyme concentration [E]_0
        
    s_col : str
        Column name for substrate concentration [S]
    
    v_col : str
        Column name for reaction rate V
        
    sheet_name : str or int
        Sheet name/index to read if loading from Excel
    
    Returns:
    --------
    kcat : float
        Turnover number (units of V / E0)
        
    KM : float
        Michaelis constant (same units as S)
    
    """
    if S_data is None:
        if isinstance(data, DataFrame):
            df = data
        elif isinstance(data, dict):
            df = DataFrame(data)
        elif isinstance(data, str):
            df = read_excel(data, sheet_name=sheet_name)
        else:
            raise TypeError("Data must be a DataFrame, dict, or Excel file path.")
    
        if s_col not in df or v_col not in df:
            raise ValueError(f"Columns '{s_col}' and '{v_col}' must be in the data.")
    
        S_data = df[s_col].to_numpy(dtype=float)
        V_data = df[v_col].to_numpy(dtype=float)
    
    S_data_norm, V_data_norm = None, None
    if normalize:
        S_scale = np.max(S_data)
        S_data_norm = S_data/S_scale
        V_scale = np.max(V_data)
        V_data_norm = V_data/V_scale
        
    
    S_data_for_fit = S_data if not normalize else S_data_norm
    V_data_for_fit = V_data if not normalize else V_data_norm
    
    # Initial guesses: Vmax=max(V), KM=median(S)
    p0 = [np.max(V_data_for_fit), np.median(S_data_for_fit)]
    
    popt, _ = curve_fit(michaelis_menten, 
                        S_data_for_fit, 
                        V_data_for_fit, 
                        p0=p0)
    
    Vmax_fit, KM_fit = popt
    
    if normalize:
        Vmax_fit *= V_scale
        KM_fit *= S_scale
    
    kcat_fit = Vmax_fit / E0
    
    return kcat_fit, KM_fit
