import numpy as np
import pandas as pd
from functions import NTU_function as NTU

# input
T1_i = 20
T2_i = 60

# design constants
L = 0.23        # m     length 350mm limit - 50*2 ends for nozzle - 20mm for end plates = 230mm
d_sh = 0.064    # m
d_noz = 0.02    # m
d_i = 0.006     # m
d_o = 0.008     # m

results = []

n=0
for tube_passes in range(1, 4):  # Passes: 1 to 3
    for N_b in range(0, 15):  # Baffles
        for N in range(1, 20//tube_passes+1):  # Tubes
            # Packing geometry logic
            if N * tube_passes == 1:
                Y = (d_sh - d_o)/2
            elif N * tube_passes>1 and N * tube_passes<=7:
                Y = (d_sh - 2*d_o)/4
            elif N * tube_passes>7 and N * tube_passes<=19: 
                Y = (d_sh - 3*d_o)/6 

            try:
                m_1, m_2, T1_out, T2_out, effectiveness, Q_t = NTU.NTU_analysis(
                    T1_i, T2_i, L, d_sh, d_noz, d_i, d_o, N, Y, N_b, tube_passes
                )
                results.append({
                    'N': N,
                    'N_b': N_b,
                    'tube_passes': tube_passes,
                    'm_dot_1 (kg/s)': round(m_1, 3),
                    'm_dot_2 (kg/s)': round(m_2, 3),
                    'T1_out (°C)': round(T1_out, 3),
                    'T2_out (°C)': round(T2_out, 3),
                    'ε': round(effectiveness, 3),
                    'Q_t': round(Q_t, 3)
                })
            except Exception as e:
                # Skip configurations that cause errors
                continue
            # print(tube_passes,N_b,N)

# print(results)
# Create DataFrame
df_results = pd.DataFrame(results)
# Save to CSV if needed
df_results.to_csv("ntu_simulation_results.csv", index=False)
print("Simulation complete. Results saved to 'ntu_simulation_results.csv'.")
