import iteration

import pandas as pd
import numpy as np
# design constants
L = 0.236        # m     length 350mm limit - 50*2 ends for nozzle - 20mm for end plates = 230mm
d_sh = 0.064    # m
d_noz = 0.02    # m
d_i = 0.006     # m
d_o = 0.008     # m
Bc_percent = 0.2

def baffle_cut_area(N,tube_passes,d_o,d_sh, Bc_percent):
    h = (Bc_percent / 100) * d_sh
    R = d_sh / 2

    # Clamp for safety
    if h > d_sh or h < 0:
        raise ValueError("Invalid baffle cut height.")

    theta = 2 * np.arccos((R - h) / R)  # radians
    A_cut = (R**2 / 2) * (theta - np.sin(theta))
    A_cut -= N*tube_passes*0.25*np.pi*d_o**2*0.8
    return A_cut

def mass_calc(d_sh, Bc_percent, N, N_b, tube_passes):
    tubes = 0.2 * L * N * tube_passes
    shell = 0.650 * L
    baffles =2.39* N_b * baffle_cut_area(N,tube_passes,d_o,d_sh, Bc_percent)
    nozzles = 0.025*4
    mass = tubes + shell + baffles + nozzles
    return mass

# Load the data
df = pd.read_csv("ntu_simulation_results.csv")

initial_count = len(df)

# --- Step 1: Filter out unstable runs ---
df_filtered = df[df['stability'] == True]
stability_removed = initial_count - len(df_filtered)

# --- Step 2: Filter m_dot bounds ---
m1_min, m1_max = 0 , 0.6843 * 0.99  # kg/s
m2_min, m2_max = 0 , 0.4536 * 0.99  # kg/s

before_flow_filter = len(df_filtered)
df_filtered = df_filtered[
    (df_filtered['m_dot_1 (kg/s)'].between(m1_min, m1_max)) &
    (df_filtered['m_dot_2 (kg/s)'].between(m2_min, m2_max))
]
flow_removed = before_flow_filter - len(df_filtered)

# --- Step 3: Compute and filter mass ---
df_filtered['mass (kg)'] = df_filtered.apply(
    lambda row: mass_calc(d_sh, Bc_percent, row['N'], row['N_b'], row['tube_passes']),
    axis=1
)

# Apply upper bound filter on mass
mass_max = 1.2  # kg — update as needed

before_mass_filter = len(df_filtered)
df_filtered = df_filtered[df_filtered['mass (kg)'] <= mass_max]
mass_removed = before_mass_filter - len(df_filtered)

# --- Step 4: Reorder columns ---
df_filtered.drop(columns=['stability'], inplace=True)

desired_order = ['mass (kg)', 'N', 'N_b', 'tube_passes', 'T1_out (°C)', 'T2_out (°C)', 'ε']
remaining_cols = [col for col in df_filtered.columns if col not in desired_order]
final_order = desired_order + remaining_cols

df_filtered = df_filtered[final_order]

# Save the cleaned data
df_filtered.to_csv("post_screening.csv", index=False, encoding='utf-8-sig')

# --- Summary report ---
print("Filtered and reordered results saved to 'post_screening.csv'.")
print(f"Initial entries:         {initial_count}")
print(f"Removed for instability: {stability_removed}")
print(f"Removed for flow bounds: {flow_removed}")
print(f"Removed for mass bounds: {mass_removed}")
print(f"Final entries:           {len(df_filtered)}")
