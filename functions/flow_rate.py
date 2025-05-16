import numpy as np
from scipy.interpolate import interp1d

# COLD SIDE data (Pressure in bar, Flowrate in L/s)
pressure_cold = [0.1584, 0.1958, 0.2493, 0.3127, 0.3723, 0.4436, 0.4950, 0.5318, 0.5739, 0.7077]
flowrate_cold = [0.6580, 0.6290, 0.5830, 0.5380, 0.4670, 0.3920, 0.3210, 0.2790, 0.2210, 0.0000]

# HOT SIDE data (Pressure in bar, Flowrate in L/s)
pressure_hot = [0.0932, 0.1688, 0.2209, 0.2871, 0.3554, 0.4041, 0.4853, 0.5260, 0.5665, 0.6239]
flowrate_hot = [0.4360, 0.3870, 0.3520, 0.3110, 0.2600, 0.2290, 0.1670, 0.1180, 0.0690, 0.0010]

# Create interpolation functions with linear extrapolation
_cold_interp = interp1d(pressure_cold, flowrate_cold, kind='linear', fill_value='extrapolate')
_hot_interp = interp1d(pressure_hot, flowrate_hot, kind='linear', fill_value='extrapolate')

def flowrate_cold_side(delta_p_bar):
    """Returns flowrate [litres/s] for a given pressure rise [bar] on the cold side."""
    return float(_cold_interp(delta_p_bar))

def flowrate_hot_side(delta_p_bar):
    """Returns flowrate [litres/s] for a given pressure rise [bar] on the hot side."""
    return float(_hot_interp(delta_p_bar))