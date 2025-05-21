import numpy as np
import pandas as pd
from functions import NTU_function as NTU

# input
T1_i = 20
T2_i = 60

# design constants
L = 0.278        # m     length 350mm limit - 50*2 ends for nozzle - 20mm for end plates = 230mm
d_sh = 0.064    # m
d_noz = 0.02    # m
d_i = 0.006     # m
d_o = 0.008     # m

results = []

for tube_passes in range(1,4):  # Passes: 1 to 3
    for N_b in range(0, 10):  # Baffles
        # for N in range(1, 20//tube_passes+1):  # Tubes
        for N in range (6,7):

            # try:
                m_1, Re_s, h_s, Delta_P1, T1_out, m_2, Re_tube, Nu_i, Delta_P2, T2_out, effectiveness, Q_t = NTU.NTU_analysis(
                    T1_i, T2_i, L, d_sh, d_noz, d_i, d_o, N, N_b, tube_passes
                )
                results.append({
                    'N': N,
                    'N_b': N_b,
                    'tube_passes': tube_passes,
                    'm_dot_1 (kg/s)': round(m_1, 3),
                    'Re_sh': round(Re_s, 1),
                    'h_s': round(h_s, 2),
                    'ΔP1 (Pa)': round(Delta_P1, 1),
                    'T1_out (°C)': round(T1_out, 3),
                    'm_dot_2 (kg/s)': round(m_2, 3),
                    'Re_tube': round(Re_tube, 1),
                    'Nu_i': round(Nu_i, 2),
                    'ΔP2 (Pa)': round(Delta_P2, 1),
                    'T2_out (°C)': round(T2_out, 3),
                    'ε': round(effectiveness, 3),
                    'Q_t (W)': round(Q_t, 3)
                })

            # except Exception as e:
            #     # Skip configurations that cause errors
            #     print('error')
            #     continue
            # print(tube_passes,N_b,N)

# print(results)
# Create DataFrame
df_results = pd.DataFrame(results)
# Save to CSV if needed
df_results.to_csv("ntu_simulation_results.csv", index=False)
print("Simulation complete. Results saved to 'ntu_simulation_results.csv'.")
