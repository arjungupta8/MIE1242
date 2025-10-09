import numpy as np


# Fixed Physical Parameters
length = 0.1 # Length of heat sink, [m]
width = 0.05 # Width of heat sink, [m]
height = 0.02 # Height of heat sink, [m]
thickness = 0.0002 # Fin thickness, [m]

# Fixed Fluid Parameters
density = 1 # [kg/m3]
viscosity = 1 # [Pa*s = kg/ms]
k = 0.02 # thermal conductivity, [W/mk]
prandtl = 0.7 # Prandtl Number = viscous diffusion rate / thermal diffusion rate
sp_heat = 1000 # specific heat, [J/kgK], amount of energy required to raise temperature of 1kg of fluid by 1K
# Nusselt # is the convective heat transfer / conductive heat transfer

# Changes in some scenarios
flow_def = 5 # [cfm]
pitch = 0.002 # Fin Pitch, [m]

# Other
cfm_to_m3s = 0.00047194745
Q_heater = 50 # [W]
T_amb = 25 # [degC], baseline ambient temperature for scenario 3

def hydraulic_diameter(gap):
    a = gap
    b = height
    return ((4*a*b)/((2*a)+(2*b)))





def scenario_1():
    flow = flow_def * cfm_to_m3s
    fin_pitch = np.linspace(0.001, 0.005, 9)


