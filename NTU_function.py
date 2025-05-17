from functions import transport_properties,flow_rate,pressure_drop_factor
import numpy as np

def NTU_analysis(T1_i,T2_i,L,d_sh,d_noz,d_i,d_o,N,Y,N_b,arrangement = 'triangular'):
    # transport properties
    T_mean = (T1_i + T2_i)/2
    cp = transport_properties.cp(T_mean) *1000
    rho = transport_properties.rho(T_mean)
    mu = transport_properties.dynamic_viscosity(T_mean)
    Pr = transport_properties.prandtl_number(T_mean)
    k_w = transport_properties.thermal_conductivity(T_mean)
    k_tube = 386    # copper tube
    B = L /(N_b+1)  #m baffle spacing
    d_h = 0.025     #m hose diameter

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
    while error > 0.001: 
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
        Ke = pressure_drop_factor.Ke_value(sigma)
        Kc = pressure_drop_factor.Kc_value(sigma)
        end_loss2 = 0.5 * rho * v_tube**2 * (Ke + Kc)
        # print('end loss 2',np.round(end_loss2,1))

        # nozzle loss 2
        v_noz2 = m_2/(rho*A_noz)    # m/s nozzle speed
        nozzle_loss2 = 2*0.5*rho*v_noz2**2
        # print('nozzle loss 2',np.round(nozzle_loss2,1))

        # shell loss 1
        # print('Re_sh',np.round(Re_sh,0))
        if arrangement == 'square':
            a = 0.34     
        elif arrangement == 'triangular':
            a = 0.2
        shell_loss = 4*a*Re_sh**-0.15*N*rho*v_sh**2
        # print('shell loss 1',np.round(shell_loss,1))

        # nozzle loss 1
        v_noz1 = m_1/(rho*A_noz)    # m/s nozzle speed
        nozzle_loss1 = 2*0.5*rho*v_noz1**2
        # print('nozzle loss 1',np.round(nozzle_loss1,1))

            # hose loss 1 
        v_hose1 = m_1/(rho*A_hose)
        hose_loss1 = 22.26*0.5*rho*v_hose1

        # hose loss 2 
        v_hose2 = m_2/(rho*A_hose)
        hose_loss2 = 23.86*0.5*rho*v_hose2

        # total dP
        Delta_P2 = friction_loss2 + end_loss2 + nozzle_loss2 + hose_loss2    # hot side
        Delta_P1 = shell_loss + nozzle_loss1 + hose_loss1                    # cold side
        #print('dP1:',np.round(Delta_P1,1),'dP2:',np.round(Delta_P2,1))

        # check flow rate
        m_1_calculated = flow_rate.flowrate_cold_side(Delta_P1/10**5) * rho/1000    # dP must be converted to bar; Q must be converted to m dot
        m_2_calculated = flow_rate.flowrate_hot_side(Delta_P2/10**5) * rho/1000
        # print(m_1_calculated,m_2_calculated)
        error = max(abs(m_1 - m_1_calculated),abs(m_2 - m_2_calculated))

        m_1 = m_1_calculated
        m_2 = m_2_calculated
        counter +=1
    ################# thermal analysis - NTU #################
    Nu_i = 0.023 * Re_tube **0.8 * Pr **0.3
    if arrangement == 'square':
        c = 0.15     
    elif arrangement == 'triangular':
        c = 0.2
    Nu_o = c * Re_sh **0.6 * Pr **0.3
    h_i = Nu_i*k_w/d_i
    h_o = Nu_o*k_w/d_o
    H = 1/(1/h_i+1/h_o*A_i/A_o+ A_i*np.log(d_o/d_i)/(2*np.pi*k_tube*L))

    cp1 = transport_properties.cp(T1_i)*1000
    cp2 = transport_properties.cp(T2_i)*1000

    def effectiveness_ntu_counterflow(m_1, cp1, T1_i, m_2, cp2, T2_i, H, A_ht):
        # Thermal capacities
        C1 = m_1 * cp1
        C2 = m_2 * cp2
        C_min = min(C1, C2)
        C_max = max(C1, C2)
        C_r = C_min / C_max

        # NTU
        NTU = H * A_ht / C_min

        # Effectiveness (ε) for counterflow
        if C_r != 1:
            epsilon = (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))
        else:
            epsilon = NTU / (1 + NTU)

        # Heat transfer
        Q = epsilon * C_min * (T2_i - T1_i)  # assumes T2 is hot, T1 is cold

        # Outlet temperatures
        T1_out = T1_i + Q / C1
        T2_out = T2_i - Q / C2

        return T1_out, T2_out

    T1_out, T2_out = effectiveness_ntu_counterflow(
        m_1, cp1, T1_i, m_2, cp2, T2_i, H, A_ht
    )
    effectiveness = m_1 * cp1 * (T1_out - T1_i)/(min(m_1*cp1,m_2*cp2)*max(T2_i - T1_out,T2_out - T1_i))
    return [m_1,m_2,T1_out,T2_out,effectiveness]


# # input
# T1_i = 20
# T2_i = 60

# # design
# L = 0.35        #m length
# d_sh = 0.064    #m shell diamter
# d_noz = 0.02    #m nozzle diameter
# d_i = 0.006     #m tube ID
# d_o = 0.008     #m tube OD
# N = 13          # tubes
# N_b = 9         # baffles
# Y = 0.014       # m tube center-center distance

# m_1,m_2,T1_out,T2_out,effectiveness = NTU_analysis(T1_i,T2_i,L,d_sh,d_noz,d_i,d_o,N,Y,N_b)
# print('m_dot_1 = ',np.round(m_1,3),'kg/s')
# print('m_dot_2 = ',np.round(m_2,3),'kg/s')
# print('T_1_out = ', np.round(T1_out, 3), '°C')
# print('T_2_out = ', np.round(T2_out, 3), '°C')
# print('ε =',np.round(effectiveness,3))