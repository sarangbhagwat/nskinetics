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

__all__ = ('michaelis_menten', 'fit_michaelis_menten',
           'tQSSA', 'fit_tQSSA',
           'michaelis_menten_rsb', 'fit_michaelis_menten_rsb')

#%% Michaelis Menten

@njit(cache=True)
def michaelis_menten(S, Vmax, KM):
    return Vmax * S / (KM + S)

briggs_haldane = michaelis_menten

# @njit(cache=True)
# def michaelis_menten_two_substrates_sequential_random(S, Vmax, ):
    
# @njit(cache=True)
# def michaelis_menten_two_substrates_sequential_ping_pong(S, Vmax, ):
    
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
    
    bounds = [(0,0),(np.inf, np.inf)]
    
    popt, _ = curve_fit(michaelis_menten, 
                        S_data_for_fit, 
                        V_data_for_fit, 
                        p0=p0,
                        bounds=bounds)
    
    Vmax_fit, KM_fit = popt
    
    if normalize:
        Vmax_fit *= V_scale
        KM_fit *= S_scale
    
    kcat_fit = Vmax_fit / E0
    
    return kcat_fit, KM_fit


#%% Michaelis-Menten with random sequential binding of two substrates

@njit(cache=True)
def michaelis_menten_rsb(A, B, E, kcat, KiA, KA, KB):
    return kcat*E*A*B / (KiA*KB + KB*A + KA*B + A*B)

briggs_haldane = michaelis_menten

# @njit(cache=True)
# def michaelis_menten_two_substrates_sequential_random(S, Vmax, ):
    
# @njit(cache=True)
# def michaelis_menten_two_substrates_sequential_ping_pong(S, Vmax, ):
    
def fit_michaelis_menten_rsb(data: Union[DataFrame, dict, str],
                             E0: float,
                             B0: float,
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
    p0 = [B0, 
          E0, 
          np.max(V_data_for_fit)/E0, 
          np.median(S_data_for_fit), 
          np.median(S_data_for_fit), 
          np.median(S_data_for_fit)]
    
    bounds = [(B0*(1-1e-6), 
               E0*(1-1e-6),
               0,
               0,
               0,
               0),
              
              (B0*(1+1e-6), 
               E0*(1+1e-6),
               np.inf, 
               np.inf, 
               np.inf, 
               np.inf)]
    
    popt, _ = curve_fit(michaelis_menten_rsb, 
                        S_data_for_fit, 
                        V_data_for_fit, 
                        p0=p0,
                        bounds=bounds)
    
    B, E, kcat_fit, KiA_fit, KA_fit, KB_fit = popt
    
    if normalize:
        kcat_fit *= V_scale
        KA_fit *= S_scale
        KB_fit *= S_scale
    
    # kcat_fit = Vmax_fit / E0
    
    return kcat_fit, KiA_fit, KA_fit, KB_fit


#%% Total quasi-steady-state approximation model (tQSSA)

@njit(cache=True)
def tQSSA(S, E, kcat, KM):
    return (kcat/2.0) * (E+S+KM - ((E+S+KM)**2 - 4*E*S)**0.5)

total_quasi_steady_state_approximation_model = tQSSA

def fit_tQSSA(data: Union[DataFrame, dict, str],
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
    p0 = [E0, 
          np.max(V_data_for_fit)/E0, 
          np.median(S_data_for_fit)]
    
    bounds = [(E0*(1-1e-6),0,0),(E0*(1+1e-6),np.inf, np.inf)]
    
    popt, _ = curve_fit(tQSSA, 
                        S_data_for_fit, 
                        V_data_for_fit, 
                        p0=p0,
                        bounds=bounds)
    
    E, kcat_fit, KM_fit = popt
    
    if normalize:
        kcat_fit *= V_scale
        KM_fit *= S_scale
    
    # kcat_fit = Vmax_fit / E0
    
    return kcat_fit, KM_fit

