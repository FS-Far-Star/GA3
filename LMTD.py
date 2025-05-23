from functions import transport_properties,flow_rate,pressure_drop_factor, b_coefficients, a_coefficients
import numpy as np

# input
T1_i = 20
T2_i = 60

# transport properties
T_mean = (T1_i + T2_i)/2
cp = transport_properties.cp(T_mean) *1000
rho = transport_properties.rho(T_mean)
mu = transport_properties.dynamic_viscosity(T_mean)
Pr = transport_properties.prandtl_number(T_mean)
k_w = transport_properties.thermal_conductivity(T_mean)
k_tube = 386    # copper tube

# design
L = 0.35        #m length
d_sh = 0.064    #m shell diamter
d_noz = 0.02    #m nozzle diameter
d_i = 0.006     #m tube ID
d_o = 0.008     #m tube OD
d_h = 0.025     #m hose diameter
N = 13          # tubes
N_b = 9         # baffles
B = L /(N_b+1)  #m baffle spacing
arrangement = 'square'  # 'triangular'
tube_passes = 1

# Packing geometry logic
if N * tube_passes == 1:
    Y = (d_sh - d_o)/2
    d_otl = d_o
elif N * tube_passes>1 and N * tube_passes<=7:
    Y = (d_sh - 2*d_o)/4
    d_otl = 2*Y+d_o
elif N * tube_passes>7 and N * tube_passes<=19: 
    Y = (d_sh - 3*d_o)/6 
    d_otl = 4*Y+d_o
else: 
    Y = (d_sh - 4*d_o)/8 
    d_otl = 5*Y+d_o
# Y = 0.012 # fixed values

# Areas
A_noz = 0.25 * np.pi * d_noz**2     # m^2
A_tube = 0.25 * np.pi * d_i**2      # m^2
A_pipe = 0.25 * np.pi * d_sh**2     # m^2
A_sh = d_sh/Y*(Y-d_o)*B             # m^2   
A_i = np.pi*d_i*L                   # m^2
A_o = np.pi*d_o*L                   # m^2
A_ht = N*np.pi*d_i*L                # m^2
A_hose = 0.25 * np.pi * d_h**2      # m^2

# mass flow initial guess
m_1 = 0.5   #kg/s
m_2 = 0.45  #kg/s

################# hydraulic analysis #################
error = 2000
counter = 0
while error > 0.001 and counter < 20: 
    m_tube = m_2/N                  # kg/s, mass flow per tubee
    v_tube = m_tube/(rho*A_tube)    # m/s
    Re_tube = rho*v_tube*d_i/mu     # tube Reynold's number
    # print(rho,v_tube,d_i,mu)

    # d_sh_adjusted = d_sh*A_sh/A_pipe 
    # v_sh = m_1/(rho*A_sh)
    # Re_sh = rho * v_sh*d_sh_adjusted/mu
    # print(rho,v_sh,d_sh_adjusted,mu)

    # friction loss 2
    # print('Re_tube',np.round(Re_tube,0))
    f = (1.82*np.log10(Re_tube)-1.64)**(-2)   # friction factor
    friction_loss2 = 0.5 * rho * v_tube**2 * (f*L*tube_passes/d_i)
    # print('friction loss 2',np.round(friction_loss2,1))

    # entry exit loss 2
    sigma = N * A_tube/A_pipe
    Ke1 = pressure_drop_factor.Ke_value(tube_passes*sigma) # tube exit
    Kc1 = pressure_drop_factor.Kc_value(tube_passes*sigma) # tube entrance
    Ke2 = pressure_drop_factor.Ke_value(tube_passes*sigma/2) # tube exit
    Kc2 = pressure_drop_factor.Kc_value(tube_passes*sigma/2) # tube entrance
    end_loss2 = 0.5 * rho * v_tube**2 * ((Ke1 + Kc1) + (tube_passes-1)*(Ke2+Kc2))
    # print('end loss 2',np.round(end_loss2,1))

    # nozzle loss 2
    v_noz2 = m_2/(rho*A_noz)    # m/s nozzle speed
    nozzle_loss2 = 2*0.5*rho*v_noz2**2
    # print('nozzle loss 2',np.round(nozzle_loss2,1))

    # shell loss 1
    # CUED method
    # if arrangement == 'square':
    #     a = 0.34     
    # elif arrangement == 'triangular':
    #     a = 0.2
    # print('Re_sh',np.round(Re_sh,0))
    # shell_loss = 4*a*Re_sh**-0.15*N*rho*v_sh**2
    # print('shell loss 1',np.round(shell_loss,1))

    S_m = B * ((d_sh - d_otl)+ (d_otl - d_o)*(Y-d_o)/Y) # valid for triangular only
    G_s = m_1/S_m
    Re_s = d_o * G_s / mu
    b = b_coefficients.b3()/(1 + 0.14*Re_s**b_coefficients.b4())
    f = b_coefficients.b1(Re_s) * (1.33/(Y/d_o))**b * Re_s**b_coefficients.b2(Re_s)
    P_p = Y * 3**0.5/2
    N_tcc = (d_sh/P_p)*(1 - 2*0.2)
    shell_loss = 2*f*(G_s**2/rho)*N_tcc


    # nozzle loss 1
    v_noz1 = m_1/(rho*A_noz)    # m/s nozzle speed
    nozzle_loss1 = 2*0.5*rho*v_noz1**2
    # print('nozzle loss 1',np.round(nozzle_loss1,1))

    # hose loss 1 
    v_hose1 = m_1/(rho*A_hose)
    hose_loss1 = 22.26*0.5*rho*v_hose1**2

    # hose loss 2 
    v_hose2 = m_2/(rho*A_hose)
    hose_loss2 = 23.86*0.5*rho*v_hose2**2

    # total dP
    Delta_P2 = friction_loss2 + end_loss2 + nozzle_loss2 + hose_loss2    # hot side
    Delta_P1 = shell_loss + nozzle_loss1 + hose_loss1                    # cold side
    #print('dP1:',np.round(Delta_P1,1),'dP2:',np.round(Delta_P2,1))

    # check flow rate
    m_1_calculated = flow_rate.flowrate_cold_side(Delta_P1/10**5) * rho/1000    # dP must be converted to bar; Q must be converted to m dot
    m_2_calculated = flow_rate.flowrate_hot_side(Delta_P2/10**5) * rho/1000
    # print(m_1_calculated,m_2_calculated)
    error = max(abs(m_1 - m_1_calculated),abs(m_2 - m_2_calculated))
    
    m_1_calculated = max(m_1_calculated ,0)
    m_2_calculated = max(m_2_calculated ,0)

    m_1 = (m_1_calculated-m_1)*0.25+m_1
    m_2 = (m_2_calculated-m_2)*0.25+m_2
    counter +=1

print('m_dot_1 = ',np.round(m_1,3),'kg/s')
print('m_dot_2 = ',np.round(m_2,3),'kg/s')
# print(counter)

m_tube = m_2/N                  # kg/s, mass flow per tubee
v_tube = m_tube/(rho*A_tube)    # m/s
Re_tube = rho*v_tube*d_i/mu     # tube Reynold's number

################# thermal analysis - LMTD #################
Nu_i = 0.023 * Re_tube **0.8 * Pr **0.3
h_i = Nu_i*k_w/d_i

# CUED method
# if arrangement == 'square':
#     c = 0.15     
# elif arrangement == 'triangular':
#     c = 0.2
# Nu_o = c * Re_sh **0.6 * Pr **0.3
# h_o = Nu_o*k_w/d_o
# H = 1/(1/h_i+1/h_o*A_i/A_o+ A_i*np.log(d_o/d_i)/(2*np.pi*k_tube*L*tube_passes))

cp1 = transport_properties.cp(T1_i)*1000
cp2 = transport_properties.cp(T2_i)*1000

a_1, a_2 , a_3 , a_4 = a_coefficients.a_table(Re_s)
j = a_1 * (1.33/(Y/d_o))**(a_3/(1+0.14*(Re_s)**a_4)) * (Re_s**a_2) 
h_s = j * cp1 * G_s * (Pr ** (-2/3))
H = 1/(1/h_i+1/h_s*A_i/A_o+ A_i*np.log(d_o/d_i)/(2*np.pi*k_tube*L*tube_passes))

# Define nonlinear equations to solve for T1_out and T2_out
from scipy.optimize import fsolve

def safe_LMTD(dT1, dT2):
    if dT1 <= 0 or dT2 <= 0:
        return 1e-6  # prevents log of 0 or negative temps
    elif abs(dT1 - dT2) < 1e-6:
        return (dT1 + dT2) / 2
    else:
        return (dT1 - dT2) / np.log(dT1 / dT2)

def equations(vars):
    T1_out, T2_out = vars

    Q1 = m_1 * cp1 * (T1_out - T1_i)
    Q2 = m_2 * cp2 * (T2_i - T2_out)

    deltaT1 = T2_i - T1_out
    deltaT2 = T2_out - T1_i
    # print(deltaT1,deltaT2)

    LMTD = safe_LMTD(deltaT1, deltaT2)
    # print(LMTD)
    F = 1
    Q_HA = H * A_ht * LMTD * F

    return [Q1 - Q2, Q1 - Q_HA]

guess = [35, 45]
T1_out, T2_out = fsolve(equations, guess)
print('T_1_out = ', np.round(T1_out, 3), '°C')
print('T_2_out = ', np.round(T2_out, 3), '°C')
a = min(m_1*cp1,m_2*cp2)
print(a)
effectiveness = m_1 * cp1 * (T1_out - T1_i)/(min(m_1*cp1,m_2*cp2)*max(T2_i - T1_out,T2_out - T1_i))
print('ε =',np.round(effectiveness,3))