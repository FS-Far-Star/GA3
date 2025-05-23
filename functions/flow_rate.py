import numpy as np
from matplotlib import pyplot as plt

########################## 2024 ##########################
# Raw data for cold side (Flowrate in L/s, Pressure in bar)
flowrate_cold = np.array([0.6333, 0.6083, 0.5750, 0.5083, 0.4250, 0.3583, 0.3083, 0.2417, 0.1917, 0.1583])
pressure_cold = np.array([0.1024, 0.1444, 0.1870, 0.2717, 0.3568, 0.4203, 0.4626, 0.5152, 0.5597, 0.5776])

# Raw data for hot side (Flowrate in L/s, Pressure in bar)
flowrate_hot = np.array([0.4826, 0.4340, 0.3924, 0.3507, 0.3021, 0.2535, 0.1979, 0.1493, 0.1111, 0.0694])
pressure_hot = np.array([0.0944, 0.1662, 0.2297, 0.2820, 0.3294, 0.3856, 0.4447, 0.5006, 0.5311, 0.5615])

# Polynomial coefficients from curve fitting (2nd degree)
cold_coeffs_2024 = np.polyfit(pressure_cold, flowrate_cold, 2)
hot_coeffs_2024 = np.polyfit(pressure_hot, flowrate_hot, 2)

# plt.plot(pressure_cold,flowrate_cold,'-*')
# plt.plot(pressure_hot,flowrate_hot,'-.')

########################## 2025 ##########################
# Raw data for cold side (Pressure in bar, Flowrate in L/s)
pressure_cold = np.array([0.1584, 0.1958, 0.2493, 0.3127, 0.3723, 0.4436, 0.4950, 0.5318, 0.5739, 0.7077])
flowrate_cold = np.array([0.6580, 0.6290, 0.5830, 0.5380, 0.4670, 0.3920, 0.3210, 0.2790, 0.2210, 0.0000])

# Raw data for hot side (Pressure in bar, Flowrate in L/s)
pressure_hot = np.array([0.0932, 0.1688, 0.2209, 0.2871, 0.3554, 0.4041, 0.4853, 0.5260, 0.5665, 0.6239])
flowrate_hot = np.array([0.4360, 0.3870, 0.3520, 0.3110, 0.2600, 0.2290, 0.1670, 0.1180, 0.0690, 0.0010])

# Polynomial coefficients from curve fitting
cold_coeffs_2025 = np.polyfit(pressure_cold, flowrate_cold, 2)
hot_coeffs_2025 = np.polyfit(pressure_hot, flowrate_hot, 2)

# plt.plot(pressure_cold,flowrate_cold,'*')
# plt.plot(pressure_hot,flowrate_hot,'.')
# plt.ylabel('flow rate L/s')
# plt.xlabel('dP')
# plt.show()

def flowrate_cold_side(delta_p_bar,flag='2025'):
    if flag == '2024':
        """Returns flowrate [litres/s] for a given pressure rise [bar] on the cold side using 2nd degree polynomial fit."""
        a, b, c = cold_coeffs_2024
        return a * delta_p_bar**2 + b * delta_p_bar + c
    if flag == '2025':
        """Returns flowrate [litres/s] for a given pressure rise [bar] on the cold side using 2nd degree polynomial fit."""
        a, b, c = cold_coeffs_2025
        return a * delta_p_bar**2 + b * delta_p_bar + c

def flowrate_hot_side(delta_p_bar,flag='2025'):
    """Returns flowrate [litres/s] for a given pressure rise [bar] on the hot side using 2nd degree polynomial fit."""
    if flag == '2024':
        """Returns flowrate [litres/s] for a given pressure rise [bar] on the cold side using 2nd degree polynomial fit."""
        a, b, c = hot_coeffs_2024
        return a * delta_p_bar**2 + b * delta_p_bar + c
    if flag == '2025':
        """Returns flowrate [litres/s] for a given pressure rise [bar] on the cold side using 2nd degree polynomial fit."""
        a, b, c = hot_coeffs_2025
        return a * delta_p_bar**2 + b * delta_p_bar + c
