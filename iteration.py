import numpy as np
import pandas as pd
from functions import NTU_function as NTU

# input
T1_i = 20
T2_i = 60

# design constants
L = 0.236        # m     length 350mm limit - 50*2 ends for nozzle - 20mm for end plates = 230mm
d_sh = 0.064    # m
d_noz = 0.02    # m
d_i = 0.006     # m
d_o = 0.008     # m
#shell_passes = 1
cut_percent = 0.2

results = []

for tube_passes in range(1,4):  # Passes: 1 to 3
    for N_b in range(0, 20):  # Baffles
        for N in range(1, 20//tube_passes):  # Tubes
            for shell_passes in range(1,3):
                m_1, Re_s, h_s, Delta_P1, T1_out, m_2, Re_tube, h_i, Delta_P2, T2_out, effectiveness, Q_t,stability = NTU.NTU_analysis(
                    T1_i, T2_i, L, d_sh, d_noz, d_i, d_o, N, N_b, tube_passes, shell_passes,cut_percent,flag = '2025'
                )
                results.append({
                    'N': N,
                    'N_b': N_b,
                    'tube_passes': tube_passes,
                    'shell_passes': shell_passes,
                    'm_dot_1 (kg/s)': round(m_1, 3),
                    'Re_sh': round(Re_s, 1),
                    'h_s': round(h_s, 2),
                    'ΔP1 (Pa)': round(Delta_P1, 1),
                    'T1_out (°C)': round(T1_out, 3),
                    'm_dot_2 (kg/s)': round(m_2, 3),
                    'Re_tube': round(Re_tube, 1),
                    'h_i': round(h_i, 2),
                    'ΔP2 (Pa)': round(Delta_P2, 1),
                    'T2_out (°C)': round(T2_out, 3),
                    'ε': round(effectiveness, 3),
                    'Q_t (kW)': round(Q_t/1000, 3),
                    'stability': stability,
                })

                # print(tube_passes,N_b,N)

# print(results)
# Create DataFrame
df_results = pd.DataFrame(results)
# Save to CSV if needed
df_results.to_csv("ntu_simulation_results.csv", index=False,encoding='utf-8-sig')
print("Simulation complete. Results saved to 'ntu_simulation_results.csv'.")
