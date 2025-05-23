import numpy as np

# NEW raw data for cold side (Flowrate in L/s, Pressure in bar)
flowrate_cold = np.array([0.6333, 0.6083, 0.5750, 0.5083, 0.4250, 0.3583, 0.3083, 0.2417, 0.1917, 0.1583])
pressure_cold = np.array([0.1024, 0.1444, 0.1870, 0.2717, 0.3568, 0.4203, 0.4626, 0.5152, 0.5597, 0.5776])

# NEW raw data for hot side (Flowrate in L/s, Pressure in bar)
flowrate_hot = np.array([0.4826, 0.4340, 0.3924, 0.3507, 0.3021, 0.2535, 0.1979, 0.1493, 0.1111, 0.0694])
pressure_hot = np.array([0.0944, 0.1662, 0.2297, 0.2820, 0.3294, 0.3856, 0.4447, 0.5006, 0.5311, 0.5615])

# Polynomial coefficients from curve fitting (2nd degree)
cold_coeffs = np.polyfit(pressure_cold, flowrate_cold, 2)
hot_coeffs = np.polyfit(pressure_hot, flowrate_hot, 2)

def flowrate_cold_side(delta_p_bar):
    """Returns flowrate [litres/s] for a given pressure rise [bar] on the cold side using 2nd degree polynomial fit."""
    a, b, c = cold_coeffs
    return a * delta_p_bar**2 + b * delta_p_bar + c

def flowrate_hot_side(delta_p_bar):
    """Returns flowrate [litres/s] for a given pressure rise [bar] on the hot side using 2nd degree polynomial fit."""
    a, b, c = hot_coeffs
    return a * delta_p_bar**2 + b * delta_p_bar + c
