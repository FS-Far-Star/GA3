from functions import transport_properties,flow_rate
import numpy as np

#### LMTD approach ####
# input
T1_i = 20
T2_i = 60

# transport properties
T_mean = (T1_i + T2_i)/2
cp = transport_properties.cp(T_mean)
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
N = 13          # tubes
N_b = 9         # baffles
Y = 0.014       # m tube center-center distance
B = L /(N_b+1)  #m baffle spacing

# Areas
A_noz = 0.25 * np.pi * d_noz**2    # m^2
A_tube = 0.25 * np.pi * d_i**2     # m^2
A_pipe = 0.25 * np.pi * d_sh**2    # m^2
A_sh = d_sh/Y*(Y-d_o)*B         # m^2   

# mass flow initial guess
m_1 = 0.5   #kg/s
m_2 = 0.45  #kg/s

#################hydraulic analysis#################
error = 2000
counter = 0
while error > 1000: 
    m_tube = m_2/N                  # kg/s, mass flow per tubee
    v_tube = m_tube/(rho*A_tube)    # m/s
    Re_tube = rho*v_tube*d_i/mu     # tube Reynold's number
    # print(rho,v_tube,d_i,mu)

    d_sh_adjusted = d_sh*A_sh/A_pipe 
    v_sh = m_1/(rho*A_sh)
    Re_sh = rho * v_sh*d_sh_adjusted/mu
    # print(rho,v_sh,d_sh_adjusted,mu)

    # friction loss 2
    # print('Re_tube',np.round(Re_tube,0))
    f = (1.82*np.log10(Re_tube)-1.64)**(-2)   # friction factor
    friction_loss2 = 0.5 * rho * v_tube**2 * (f*L/d_i)
    # print('friction loss 2',np.round(friction_loss2,1))

    # entry exit loss 2
    sigma = N * A_tube/A_pipe
    Ke = 0.45
    Kc = 0.8
    end_loss2 = 0.5 * rho * v_tube**2 * (Ke + Kc)
    # print('end loss 2',np.round(end_loss2,1))

    # nozzle loss 2
    v_noz2 = m_2/(rho*A_noz)    # m/s nozzle speed
    nozzle_loss2 = 2*0.5*rho*v_noz2**2
    # print('nozzle loss 2',np.round(nozzle_loss2,1))

    # shell loss 1
    # print('Re_sh',np.round(Re_sh,0))
    a = 0.34     # 0.34 for square, 0.2 for triangular
    shell_loss = 4*a*Re_sh**-0.15*N*rho*v_sh**2
    # print('shell loss 1',np.round(shell_loss,1))

    # nozzle loss 1
    v_noz1 = m_1/(rho*A_noz)    # m/s nozzle speed
    nozzle_loss1 = 2*0.5*rho*v_noz1**2
    # print('nozzle loss 1',np.round(nozzle_loss1,1))

    # total dP
    Delta_P2 = friction_loss2 + end_loss2 + nozzle_loss2    # hot side
    Delta_P1 = shell_loss + nozzle_loss1                    # cold side
    # print('dP1:',np.round(Delta_P1,1),'dP2:',np.round(Delta_P2,1))

    # check flow rate
    m_1_calculated = flow_rate.flowrate_cold_side(Delta_P1/10**5) * rho/1000    # dP must be converted to bar; Q must be converted to m dot
    m_2_calculated = flow_rate.flowrate_hot_side(Delta_P2/10**5) * rho/1000
    # print(m_1_calculated,m_2_calculated)
    error = max(abs(m_1 - m_1_calculated),abs(m_2 - m_2_calculated))

    m_1 = m_1_calculated
    m_2 = m_2_calculated
    counter +=1

print('m_dot_1 = ',np.round(m_1,3),'kg/s')
print('m_dot_2 = ',np.round(m_2,3),'kg/s')
# print(counter)