from functions import transport_properties,flow_rate,pressure_drop_factor, b_coefficients, a_coefficients
import numpy as np

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

def NTU_analysis(T1_i,T2_i,L,d_sh,d_noz,d_i,d_o,N,N_b,tube_passes,arrangement = 'triangular'):
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
    # Y = 0.012 # fixed values

    # transport properties
    T_mean = (T1_i + T2_i)/2
    # cp = transport_properties.cp(T_mean) *1000
    rho = transport_properties.rho(T_mean)
    mu = transport_properties.dynamic_viscosity(T_mean)
    Pr = transport_properties.prandtl_number(T_mean)
    k_w = transport_properties.thermal_conductivity(T_mean)
    k_tube = 386    # copper tube
    B = L /(N_b+1)  #m baffle spacing
    # print(B)
    d_h = 0.025     #m hose diameter

    # Areas
    A_noz = 0.25 * np.pi * d_noz**2     # m^2
    A_tube = 0.25 * np.pi * d_i**2      # m^2
    A_pipe = 0.25 * np.pi * d_sh**2     # m^2
    A_sh = d_sh/Y*(Y-d_o)*B             # m^2   
    A_i = np.pi*d_i*L*tube_passes       # m^2
    A_o = np.pi*d_o*L*tube_passes       # m^2
    A_ht = N*np.pi*d_i*L*tube_passes    # m^2
    A_hose = 0.25 * np.pi * d_h**2      # m^2

    # mass flow initial guess
    m_1 = 0.5   #kg/s
    m_2 = 0.45  #kg/s

    ################# hydraulic analysis #################
    error = 2000
    counter = 0
    while error > 0.001: 
        m_tube = m_2/N                  # kg/s, mass flow per tube
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
        S_m = B * ((d_sh - d_otl)+ (d_otl - d_o)*(Y-d_o)/Y) # valid for triangular only
        G_s = m_1/S_m
        Re_s = d_o * G_s / mu

        # print('Re_sh',np.round(Re_sh,0))
        # if arrangement == 'square':
        #     a = 0.34     
        # elif arrangement == 'triangular':
        #     a = 0.2
        b = b_coefficients.b3()/(1 + 0.14*Re_s**b_coefficients.b4())
        f = b_coefficients.b1(Re_s) * (1.33/(Y/d_o))**b * Re_s**b_coefficients.b2(Re_s)
        P_p = Y * 3**0.5/2
        N_tcc = (d_sh/P_p)*(1 - 2*0.2)
        shell_loss = 2*f*(G_s**2/rho)*N_tcc*(N_b-1)
        # shell_loss = 4*a*Re_sh**-0.15*N*rho*v_sh**2
        # print('shell loss 1',np.round(shell_loss,1))

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
        # print('dP1:',np.round(Delta_P1,1),'dP2:',np.round(Delta_P2,1))

        # check flow rate
        m_1_calculated = flow_rate.flowrate_cold_side(Delta_P1/10**5) * rho/1000    # dP must be converted to bar; Q must be converted to m dot
        m_2_calculated = flow_rate.flowrate_hot_side(Delta_P2/10**5) * rho/1000
        # print(m_1_calculated,m_2_calculated)
        error = max(abs(m_1 - m_1_calculated),abs(m_2 - m_2_calculated))

        m_1 = (m_1_calculated-m_1)*0.25+m_1
        m_2 = (m_2_calculated-m_2)*0.25+m_2
        if m_1 <0 or np.isnan(m_1) or m_2 <0 or np.isnan(m_2):
            print('invalid m_dot')
        counter +=1
    ################# thermal analysis - NTU #################

    if m_1 <0 or np.isnan(m_1) or m_2 <0 or np.isnan(m_2):
        m_1 = []
        m_2 = []
        effectiveness = []
        T1_out = []
        T2_out = []
    else:
        m_tube = m_2/N                  # kg/s, mass flow per tubee
        v_tube = m_tube/(rho*A_tube)    # m/s
        Re_tube = rho*v_tube*d_i/mu     # tube Reynold's number
        Nu_i = 0.023 * Re_tube **0.8 * Pr **0.3
        h_i = Nu_i*k_w/d_i

        # CUED method
        if arrangement == 'square':
            c = 0.15     
        elif arrangement == 'triangular':
            c = 0.2
        Nu_o = c * Re_sh **0.6 * Pr **0.3
        h_o = Nu_o*k_w/d_o
        # H = 1/(1/h_i+1/h_o*A_i/A_o+ A_i*np.log(d_o/d_i)/(2*np.pi*k_tube*L*tube_passes))

        cp1 = transport_properties.cp(T1_i)*1000
        cp2 = transport_properties.cp(T2_i)*1000

        a_1, a_2 , a_3 , a_4 = a_coefficients.a_table(Re_s)
        j = a_1 * (1.33/(Y/d_o))**(a_3/(1+0.14*(Re_s)**a_4)) * (Re_s**a_2) 
        h_s = j * cp1 * G_s * (Pr ** (-2/3))
        H = 1/(1/h_i+1/h_s*A_i/A_o+ A_i*np.log(d_o/d_i)/(2*np.pi*k_tube*L*tube_passes))

        print(h_o,h_s)
        # print('H = ',H)
        # print('A_ht = ',A_ht)

       
        T1_out, T2_out = effectiveness_ntu_counterflow(
            m_1, cp1, T1_i, m_2, cp2, T2_i, H, A_ht
        )
        effectiveness = m_1 * cp1 * (T1_out - T1_i)/(min(m_1*cp1,m_2*cp2)*max(T2_i - T1_out,T2_out - T1_i))
        Q_t = m_1 * (T1_out-T1_i) * cp1
    return [m_1,Re_s,h_s, Delta_P1,T1_out,m_2, Re_tube, Nu_i, Delta_P2,T2_out,effectiveness,Q_t]


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