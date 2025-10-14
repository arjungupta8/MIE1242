import matplotlib.pyplot as plt
import numpy as np


# Fixed Physical Parameters
length = 0.1 # Length of heat sink, [m]
width = 0.05 # Width of heat sink, [m]
height = 0.02 # Height of heat sink, [m]
thickness = 0.0002 # Fin thickness, [m]

# Fixed Fluid Parameters
rho = 1 # [kg/m3]
mu = 0.0001 # [Pa*s = kg/ms]
k_air = 0.02 # thermal conductivity, [W/mk]
prandtl = 0.7 # Prandtl Number = viscous diffusion rate / thermal diffusion rate
sp_heat = 1000 # specific heat, [J/kgK], amount of energy required to raise temperature of 1kg of fluid by 1K
const_A = 0.40377 # Numeric Value for Coefficient A for constant wall temperature

# Changes in some scenarios
flow_def = 5 # Default Flow Rate [cfm]
pitch_def = 0.002 # Default Fin Pitch, [m]

# Other
cfm_to_m3s = 0.00047194745
Q_heater = 50 # [W]
T_amb = 25 # [degC], baseline ambient temperature for scenario 3


def calculate_Nu(Dh, Re, f):
    Gz = (Dh*Re*prandtl)/length # Graetz Number
    Nu_L_Lev = const_A * ((f*Re*Gz)**(1/3))
    Nu_inf = 7.541 # fully developed transport constant for flow between parallel plates with constant wall temperature
    n = (Nu_inf + 45.5)/14.5
    O_lev = (Nu_inf - 7.16)/5
    Nu = ((((Nu_L_Lev)**n) + (Nu_inf - O_lev)**n)**(1/n)) + O_lev

    return Nu


def pressure_drop(f, rho, vel, Dh):
    return f * length / Dh * rho * (vel ** 2) * 0.5



def compute_results(pitch, flow):
    n_chan = int(np.floor(width/pitch)) # Calculates the number of channels given the fin width and fin pitch
    n_fins = n_chan + 1
    gap = (width - (thickness * n_fins)) / n_chan # gap between each channel. Heat sink width - combined fins / num channels
    A_chan = gap * height # Cross-section of a single channel
    Q_chan = flow / n_chan
    vel = Q_chan / A_chan # velocity of flow
    Dh = (4*gap*height)/((2*gap)+(2*height))
    # print ('Pitch and Dh', pitch, Dh)
    # print ('Gap: ', gap)
    # print ('Number of chan: ', n_chan)
    Re = rho * vel * Dh / mu
    # print (f'Pitch: {pitch}, Gap: {gap}, Num Chan: {n_chan}, Num Fins: {n_fins}, Re: {Re}')
    f = 96/Re # Darcy Friction Factor for parallel plates
    Nu = calculate_Nu(Dh, Re, f)
    h = Nu * k_air / Dh # Calculates the convective heat transfer coefficient
    fin_area = (n_fins-1) * 2 * height * length # Surface area per fin = height * length * 2 (2 sides). Outside surface of outside fins do not have the same HT
    base_area = length * (gap*n_chan) # Surface area of heat sink base
    A_total = fin_area + (2*base_area) # Inside of Top and Bottom Surfaces also considered
    hA = h * A_total

    delta_p = f * length / Dh * rho * (vel ** 2) * 0.5 # Darcy - Weisbach for major losses, no minor losses
    fan_power = delta_p * flow # Is this needed?

    return Nu, hA, delta_p


def fan_pressure(flow):
    return -3.6012 * (flow ** 2) - 6.26 * flow + 597.62



def scenario_1():
    flow = flow_def * cfm_to_m3s
    fin_pitch = np.linspace(0.001, 0.005, 9)
    Nu_all = [] # array to hold the resulting Nu's.
    hA_all = [] # array to hold the resulting convective heat transfer coefficients * Surface Areas (asks for hA)
    R_conv_all = [] # Array to hold the resulting thermal resistances. = 1/ (h*SurfaceArea)
    delta_p_all = []
    for pitch in fin_pitch:
        Nu, hA, delta_p = compute_results(pitch, flow)
        Nu_all.append(Nu)
        hA_all.append(hA)
        R_conv = 1/hA
        R_conv_all.append(R_conv)
        delta_p_all.append(delta_p)

    plt.figure()
    plt.plot(fin_pitch, Nu_all, marker='o')
    plt.xlabel('Fin pitch (m)')
    plt.ylabel('Average Nusselt number per fin')
    plt.title('Scenario 1: Nu vs Fin Pitch (Q = 5 CFM)')
    plt.savefig('a2_s1_nu_vs_pitch.png', dpi=300)
    plt.close()

    plt.figure()
    plt.plot(fin_pitch, hA_all, marker='o')
    plt.xlabel('Fin pitch (m)')
    plt.ylabel('Overall Heat Transfer Coefficient (h*A)')
    plt.title('Scenario 1: Heat Transfer Coefficient vs Fin Pitch (Q = 5 CFM)')
    plt.savefig('a2_s1_hA_vs_pitch.png', dpi=300)
    plt.close()

    plt.figure()
    plt.plot(fin_pitch, R_conv_all, marker='o')
    plt.xlabel('Fin pitch (m)')
    plt.ylabel('Convective Thermal Resistance')
    plt.title('Scenario 1: Convective Thermal Resistance vs Fin Pitch (Q = 5 CFM)')
    plt.savefig('a2_s1_res_vs_pitch.png', dpi=300)
    plt.close()

    plt.figure()
    plt.plot(fin_pitch, delta_p_all, marker='o')
    plt.xlabel('Fin pitch (m)')
    plt.ylabel('Pressure Drop (Pascals)')
    plt.title('Scenario 1: Overall Pressure Drop vs Fin Pitch (Q = 5 CFM)')
    plt.savefig('a2_s1_deltaP_vs_pitch.png', dpi=300)
    plt.close()

def scenario_2():
    flow_range = np.linspace(1, 10, 10)
    Nu_all = []
    hA_all = []
    R_conv_all = []
    delta_p_all = []
    for flows in flow_range:
        flow = flows * cfm_to_m3s
        Nu, hA, delta_p = compute_results(pitch_def, flow)
        Nu_all.append(Nu)
        hA_all.append(hA)
        R_conv = 1/hA
        R_conv_all.append(R_conv)
        delta_p_all.append(delta_p)

    plt.figure()
    plt.plot(flow_range, hA_all, marker='o')
    plt.xlabel('Flow Rate (CFM)')
    plt.ylabel('Heat Transfer Coefficient (hA)')
    plt.title('Scenario 2: Heat Transfer Coefficient vs Flow Rate (Fin Pitch = 2mm)')
    plt.savefig('a2_s2_hA_vs_flow_rate.png', dpi=300)
    plt.close()

    plt.figure()
    plt.plot(flow_range, R_conv_all, marker='o')
    plt.xlabel('Flow Rate (CFM)')
    plt.ylabel('Convective Thermal Resistance (1/hA)')
    plt.title('Scenario 2: Convective Thermal Resistance vs Flow Rate (Fin Pitch = 2mm)')
    plt.savefig('a2_s2_res_vs_flow_rate.png', dpi=300)
    plt.close()

    plt.figure()
    plt.plot(flow_range, Nu_all, marker='o')
    plt.xlabel('Flow Rate (CFM)')
    plt.ylabel('Average Nusselt number per fin')
    plt.title('Scenario 2: Average Nusselt Number vs Flow Rate (Fin Pitch = 2mm)')
    plt.savefig('a2_s2_Nu_vs_flow_rate.png', dpi=300)
    plt.close()

    plt.figure()
    plt.plot(flow_range, delta_p_all, marker='o')
    plt.xlabel('Flow Rate (CFM)')
    plt.ylabel('Pressure Drop (Pascals)')
    plt.title('Scenario 2: Pressure Drop vs Flow Rate (Fin Pitch = 2mm)')
    plt.savefig('a2_s2_deltaP_vs_flow_rate.png', dpi=300)
    plt.close()


def scenario_3():
    fin_pitches = np.linspace(0.001, 0.005, 5)
    flow_range = np.linspace(1, 20, 20)

    delta_T_all = []
    delta_P_all = []
    out_T_all = []

    for pitch in fin_pitches:
        deltaP_allPitches = []
        Nu_all = []
        hA_all = []
        R_conv_all = []
        for flows in flow_range:
            flow = flows * cfm_to_m3s
            Nu, hA, delta_p = compute_results(pitch, flow)
            Nu_all.append(Nu)
            hA_all.append(hA)
            R_conv = 1 / hA
            R_conv_all.append(R_conv)
            deltaP_allPitches.append(delta_p)
        deltaP_array = np.array(deltaP_allPitches) # convert to NumPy array to make it easier to work with
        deltaP_fan = fan_pressure(flow_range) # can send the whole array as it is NumPy and will return array
        diff = np.abs(deltaP_fan - deltaP_array) # NumPy arrays be like MATLAB arrays. So this works
        idx = np.nanargmin(diff) # Looks for index where difference is the smallest - operating point
        flow_operating = (flow_range[idx]) * cfm_to_m3s
        deltaP_operating = deltaP_allPitches[idx]
        deltaT_air = Q_heater / (rho * flow_operating * sp_heat)
        delta_T_all.append(deltaT_air)
        delta_P_all.append(deltaP_operating)
        out_T = T_amb + deltaT_air
        out_T_all.append(out_T)


        plt.figure()
        plt.plot(flow_range, deltaP_array, marker = 'o')
        plt.xlabel('Flow Rate (CFM)')
        plt.ylabel('Pressure Drop (Pascals)')
        plt.title(f'Flow Rate vs Overall Pressure Drop through Heat Sink, Fin Pitch = {pitch}m')
        plt.savefig(f'a2_s3_flow_vs_deltaP_{pitch}.png', dpi=300)
        plt.close()

    plt.figure()
    plt.plot(fin_pitches, delta_T_all, marker = 'o')
    plt.xlabel('Fin Pitch (m)')
    plt.ylabel('Air Temperature Rise (degK)')
    plt.title(f'Air Temperature Rise vs Fin Pitch, for Operating Flow Rate')
    plt.savefig('a2_s3_pitch_vs_temprise.png', dpi=300)
    plt.close()

    plt.figure()
    plt.plot(fin_pitches, delta_P_all, marker = 'o')
    plt.xlabel('Fin Pitch (m)')
    plt.ylabel('Pressure Drop (Pascals)')
    plt.title(f'Operating Pressure Drop vs Fin Pitch, for Operating Flow Rate')
    plt.savefig('a2_s3_pitch_vs_pressureDrop.png', dpi=300)
    plt.close()

    plt.figure()
    plt.plot(fin_pitches, out_T_all, marker = 'o')
    plt.xlabel('Fin Pitch (m)')
    plt.ylabel('Outlet Air Temperature (degC)')
    plt.title(f'Outlet Air Temperature vs Fin Pitch, for Operating Flow Rate')
    plt.savefig('a2_s3_pitch_vs_outAirTemp.png', dpi=300)
    plt.close()

if __name__ == '__main__':
    scenario_1()
    scenario_2()
    scenario_3()









