import numpy as np

# Raw data for cold side (Pressure in bar, Flowrate in L/s)
pressure_cold = np.array([0.1584, 0.1958, 0.2493, 0.3127, 0.3723, 0.4436, 0.4950, 0.5318, 0.5739, 0.7077])
flowrate_cold = np.array([0.6580, 0.6290, 0.5830, 0.5380, 0.4670, 0.3920, 0.3210, 0.2790, 0.2210, 0.0000])

# Raw data for hot side (Pressure in bar, Flowrate in L/s)
pressure_hot = np.array([0.0932, 0.1688, 0.2209, 0.2871, 0.3554, 0.4041, 0.4853, 0.5260, 0.5665, 0.6239])
flowrate_hot = np.array([0.4360, 0.3870, 0.3520, 0.3110, 0.2600, 0.2290, 0.1670, 0.1180, 0.0690, 0.0010])

# Polynomial coefficients from curve fitting
# Cold side: Q = a*P^2 + b*P + c
cold_coeffs = np.polyfit(pressure_cold, flowrate_cold, 2)

# Hot side: Q = a*P^2 + b*P + c
hot_coeffs = np.polyfit(pressure_hot, flowrate_hot, 2)

def flowrate_cold_side(delta_p_bar):
    """Returns flowrate [litres/s] for a given pressure rise [bar] on the cold side using 2nd degree polynomial fit."""
    a, b, c = cold_coeffs
    return a * delta_p_bar**2 + b * delta_p_bar + c

def flowrate_hot_side(delta_p_bar):
    """Returns flowrate [litres/s] for a given pressure rise [bar] on the hot side using 2nd degree polynomial fit."""
    a, b, c = hot_coeffs
    return a * delta_p_bar**2 + b * delta_p_bar + c

