import numpy as np
import matplotlib.pyplot as plt
from numpy.ma.mrecords import addfield

# Variables that do not change
g = 9.81
eps = 2e-6 # SS Episilon value in meters, 0.002 in mm, so divide by 1000 for meters.
rho = 1000  # Default Density: [kg/m^3]
pipes = {
        2: (0.2, 0.02),
        3: (0.3, 0.02),
        4: (1.0, 0.02),
        5: (0.2, 0.02),
        6: (2.0, 0.01),  # heat exchanger
        7: (1.0, 0.02),
        8: (0.2, 0.02),
    } # Pipe Sections: (L, D) in meters

def reynolds_number(V, D, mu, rho=rho):
    return rho * V * D / mu



def haaland_friction(Re, D): # Returns the friction factor
    if Re < 2000:
        return 64.0 / Re  # laminar
    return 1.0 / (-1.8 * np.log10((eps / (3.7 * D)) ** 1.1 + 6.9 / Re)) ** 2

def pressure_drop_friction(f, L, D, Q, rho=rho):
    A = np.pi*((D/2)**2) # Cross-sectional area of the pipe
    V = Q/A # Velocity = Flow Rate / Area

    return f * (L/D) * 0.5 * rho * V**2 # Head Loss Equation

def pressure_drop_minor(K, d, Q, rho=rho):
    A = np.pi*((d/2)**2) # Cross-sectional area of the pipe
    V = Q/A # Velocity = Flow Rate / Area
    return K * 0.5 * rho * V**2

def pressure_drop_branch1(pipes=pipes):



def main():
    scenario_number = 1 # Options: 1, 2, 3

    # Default variables, but change in scenarios
    Qt = 5 # Default Flow Rate, LPM - Scenario 1
    K_valve = 5 # Default Minor Loss Coeff for Valve - Scenario 2
    mu = 0.001 # Default Viscosity: [Pa s], Pascal Seconds = kg * m^-1 * s^-1 - Scenario 3







if __name__ == '__main__':
    main()