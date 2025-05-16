import numpy as np
from scipy.interpolate import interp1d

# Data extracted from the table for saturated water (liquid only)
data = {
    "Temp_C": [
        0.01, 10, 20, 30, 40, 50, 60, 70, 80, 90,
        100, 125, 150, 175, 200, 225, 250, 275, 300, 325,
        350, 360, 374.15
    ],
    "vf_m3_per_kg": [
        0.001000, 0.001000, 0.001000, 0.001000, 0.00101, 0.00101, 0.00102, 0.00105, 0.00103, 0.00104,
        0.00104, 0.00107, 0.00109, 0.00112, 0.00116, 0.00120, 0.00125, 0.00132, 0.00140, 0.00153,
        0.00174, 0.00190, 0.00317
    ],
    "cp_kJ_per_kgK": [
        4.217, 4.193, 4.182, 4.179, 4.179, 4.181, 4.185, 4.190, 4.197, 4.205,
        4.216, 4.254, 4.310, 4.389, 4.497, 4.648, 4.867, 5.202, 5.762, 6.861,
        10.10, 14.6, 14.6
    ],
    "lambda_W_per_mK": [
        0.569, 0.587, 0.603, 0.618, 0.632, 0.643, 0.653, 0.662, 0.670, 0.676,
        0.681, 0.687, 0.687, 0.679, 0.665, 0.644, 0.616, 0.582, 0.541, 0.493,
        0.437, 0.400, 0.24
    ],
    "mu_kg_per_ms": [
        0.001755, 0.001301, 0.001002, 0.000797, 0.000651, 0.000544, 0.000462, 0.000400, 0.000350, 0.000311,
        0.000278, 0.000219, 0.000180, 0.000153, 0.000133, 0.0001182, 0.0001065, 0.0000972, 0.0000897, 0.0000790,
        0.0000648, 0.0000582, 0.000045
    ],
    "Pr": [
        8.8, 9.29, 6.95, 5.39, 4.31, 3.53, 2.96, 2.53, 2.19, 1.93,
        1.723, 1.358, 1.133, 0.990, 0.902, 0.853, 0.841, 0.869, 0.955, 1.100,
        1.50, 2.11, float('inf')
    ]
}

# Create interpolation functions with linear extrapolation
_interp_funcs = {
    "vf": interp1d(data["Temp_C"], data["vf_m3_per_kg"], kind='linear', fill_value='extrapolate'),
    "cp": interp1d(data["Temp_C"], data["cp_kJ_per_kgK"], kind='linear', fill_value='extrapolate'),
    "lambda": interp1d(data["Temp_C"], data["lambda_W_per_mK"], kind='linear', fill_value='extrapolate'),
    "mu": interp1d(data["Temp_C"], data["mu_kg_per_ms"], kind='linear', fill_value='extrapolate'),
    "Pr": interp1d(data["Temp_C"], data["Pr"], kind='linear', fill_value='extrapolate'),
}

def rho(T_C):
    """Returns density [kg/m^3] of saturated liquid water at T [Celsius]"""
    vf = _interp_funcs["vf"](T_C)
    return 1 / vf

def cp(T_C):
    """Returns isobaric specific heat capacity [kJ/kg-K] of saturated liquid water at T [Celsius]"""
    return _interp_funcs["cp"](T_C)

def thermal_conductivity(T_C):
    """Returns thermal conductivity [W/m-K] of saturated liquid water at T [Celsius]"""
    return _interp_funcs["lambda"](T_C)

def dynamic_viscosity(T_C):
    """Returns dynamic viscosity [kg/(m-s)] of saturated liquid water at T [Celsius]"""
    return _interp_funcs["mu"](T_C)

def prandtl_number(T_C):
    """Returns Prandtl number of saturated liquid water at T [Celsius]"""
    return _interp_funcs["Pr"](T_C)
