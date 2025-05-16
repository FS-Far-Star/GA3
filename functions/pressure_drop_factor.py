import numpy as np
import matplotlib.pyplot as plt

sigma = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
Kc = np.array([0.5,0.46,0.42,0.38,0.34,0.30,0.26,0.22,0.18,0.14,0.10])
Ke = np.array([1.0, 0.8, 0.62,0.47,0.32,0.20,0.10,0.03,-0.03,-0.08,-0.1])


Kc_coefficients = np.polyfit(sigma, Kc, 2)
Ke_coefficients = np.polyfit(sigma, Ke, 2)

def Kc_value(area_ratio):
    a,b,c = Kc_coefficients
    return c + b*area_ratio + a*area_ratio**2

def Ke_value(area_ratio):
    d,e,f = Ke_coefficients
    return f + e*area_ratio + d*area_ratio**2