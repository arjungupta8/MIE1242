import numpy as np
import matplotlib.pyplot as plt

# Global variables that do not change
g = 9.81
eps = 2e-6 # SS Episilon value in meters, 0.002 in mm, so divide by 1000 for meters.
rho = 1000  # Default Density: [kg/m^3]
K_elbow = 0.40 # 90 Degree, Long Radius Elbow
K_SE = 0.5625 # Sudden Expansion
K_SC = 0.315 # Sudden Contraction
K_tee_branch = 1 # Tee Connection, Branch (90 Degree flow)
K_tee_line = 0.24 # Tee Connection, Line (Straight)

# Default global variables that change with scenarios
Qt_Default = 5e-3/60 # Default Flow Rate, m3/s (Converted from LPM) - Scenario 1
K_valve_default = 5 # Default Minor Loss Coeff for Valve - Scenario 2
mu_default = 0.001 # Default Viscosity: [Pa s], Pascal Seconds = kg * m^-1 * s^-1 - Scenario 3

pipes = {
        2: (0.2, 0.02),
        3: (0.3, 0.02),
        4: (1.0, 0.02),
        5: (0.2, 0.02),
        6: (2.0, 0.01),  # heat exchanger
        7: (1.0, 0.02),
        8: (0.2, 0.02),
    } # Pipe Sections: (L, D) in meters

def reynolds_number(V, D, mu):
    return rho * V * D / mu

def haaland_friction(Re, D): # Returns the friction factor
    if Re < 2000:
        return 64.0 / Re  # laminar
    return 1.0 / (-1.8 * np.log10((eps / (3.7 * D)) ** 1.1 + 6.9 / Re)) ** 2

def pressure_drop_friction(Q, L, D, mu):
    A = np.pi*((D/2)**2) # Cross-sectional area of the pipe
    V = Q/A # Velocity = Flow Rate / Area
    Re = reynolds_number(V, D, mu)
    f = haaland_friction(Re, D)

    return f * (L/D) * 0.5 * rho * V**2 # Head Loss Equation

def pressure_drop_minor(K, D, Q):
    A = np.pi*((D/2)**2) # Cross-sectional area of the pipe
    V = Q/A # Velocity = Flow Rate / Area
    return K * 0.5 * rho * V**2

def pressure_drop_branch1(Q1, K_valve, mu):
    dp = 0
    for sec in [2, 3, 4]:
        L, D = pipes[sec]
        dp += pressure_drop_friction(Q1, L, D, mu)
    dp += pressure_drop_minor(K_valve, D, Q1) # Valve Minor Loss
    dp += (2*(pressure_drop_minor(K_tee_branch, D, Q1)))
    dp += (2*(pressure_drop_minor(K_elbow, D, Q1)))
    return dp


def pressure_drop_branch2(Q2, K_valve, mu):
    dp=0
    for sec in [5, 6, 7, 8]:
        L, D = pipes[sec]
        dp += pressure_drop_friction(Q2, L, D, mu)
    dp += (2*(pressure_drop_minor(K_tee_line, D, Q2)))
    dp += (pressure_drop_minor(K_SE, D, Q2))
    dp += (pressure_drop_minor(K_SC, D, Q2))
    return dp


def solver(Qt, K_valve, mu):
    # Iterative solver
    def residual(Q1):
        Q2 = Qt - Q1
        if Q2 <= 0 or Q1 <= 0:
            return 1e6 # This is a failsafe - If the solution does not make physical sense/if backflow is present
        dp1 = pressure_drop_branch1(Q1, K_valve, mu)
        dp2 = pressure_drop_branch2(Q2, K_valve, mu)
        return dp1 - dp2

    Q1_low, Q1_high = 1e-8, Qt-1e-8
    for _ in range(100):
        Q1_mid = 0.5*(Q1_low + Q1_high)
        if residual(Q1_low)*residual(Q1_mid) < 0:
            Q1_high = Q1_mid
        else:
            Q1_low = Q1_mid
    Q1 = 0.5*(Q1_low + Q1_high)
    Q2 = Qt - Q1
    return Q1, Q2


def scenario_1():
    Qt_values = np.linspace(1e-3/60, 10e-3/60, 100)
    K_valve = K_valve_default
    mu = mu_default
    Q1_list, Q2_list, ratio = [], [], []
    for Qt in Qt_values:
        Q1, Q2  = solver(Qt, K_valve, mu)
        Q1_list.append(Q1*60000) # Convert Q1 and Q2 back to LPM by * 60000
        Q2_list.append(Q2*60000)
        ratio.append(Q1/Q2)

    plt.figure()
    plt.plot(Qt_values * 60000, ratio, 'o-')
    plt.xlabel("Qt [LPM]")
    plt.ylabel("Q1/Q2")
    plt.title("Scenario 1: Flow split ratio vs Qt")
    plt.grid(True)
    plt.savefig("scenario1_ratio.png", dpi=300)
    plt.show()

    # Plot 2: Q1, Q2, Qt
    plt.figure()
    plt.plot(Qt_values * 60000, Q1_list, 'o-', label="Q1")
    plt.plot(Qt_values * 60000, Q2_list, 'o-', label="Q2")
    plt.plot(Qt_values * 60000, Qt_values * 60000, 'k--', label="Qt")
    plt.xlabel("Qt [LPM]")
    plt.ylabel("Flowrate [LPM]")
    plt.title("Scenario 1: Q1, Q2, Qt vs Qt")
    plt.legend()
    plt.grid(True)
    plt.savefig("scenario1_flows.png", dpi=300)
    plt.show()

def scenario_2():
    Qt = Qt_Default
    K_valve_values = np.linspace(0, 10, 100)
    mu = mu_default
    Q1_list, Q2_list, ratio = [], [], []
    for K_valve in K_valve_values:
        Q1, Q2  = solver(Qt, K_valve, mu)
        Q1_list.append(Q1*60000) # Convert Q1 and Q2 back to LPM by * 60000
        Q2_list.append(Q2*60000)
        ratio.append(Q1/Q2)

    plt.figure()
    plt.plot(K_valve_values, ratio, 'o-')
    plt.xlabel("Valve K")
    plt.ylabel("Q1/Q2")
    plt.title("Scenario 2: Flow split ratio vs Valve Loss")
    plt.grid(True)
    plt.savefig("scenario2_ratio.png", dpi=300)
    plt.show()

    # Plot 2: Q1, Q2, Qt
    plt.figure()
    plt.plot(K_valve_values, Q1_list, 'o-', label="Q1")
    plt.plot(K_valve_values, Q2_list, 'o-', label="Q2")
    plt.plot(K_valve_values, [Qt * 60000] * len(K_valve_values), 'k--', label="Qt")
    plt.xlabel("Valve K")
    plt.ylabel("Flowrate [LPM]")
    plt.title("Scenario 2: Q1, Q2, Qt vs Valve Loss")
    plt.legend()
    plt.grid(True)
    plt.savefig("scenario2_flows.png", dpi=300)
    plt.show()

def scenario_3():
    Qt = Qt_Default
    K_valve = K_valve_default
    mu_values = np.linspace(0.001, 0.01, 100)
    Q1_list, Q2_list, ratio = [], [], []
    for mu in mu_values:
        Q1, Q2  = solver(Qt, K_valve, mu)
        Q1_list.append(Q1*60000) # Convert Q1 and Q2 back to LPM by * 60000
        Q2_list.append(Q2*60000)
        ratio.append(Q1/Q2)

    # Plot 1: ratio
    plt.figure()
    plt.plot(mu_values, ratio, 'o-')
    plt.xlabel("Viscosity [Pa.s]")
    plt.ylabel("Q1/Q2")
    plt.title("Scenario 3: Flow split ratio vs Viscosity")
    plt.grid(True)
    plt.savefig("scenario3_ratio.png", dpi=300)
    plt.show()

    # Plot 2: Q1, Q2, Qt
    plt.figure()
    plt.plot(mu_values, Q1_list, 'o-', label="Q1")
    plt.plot(mu_values, Q2_list, 'o-', label="Q2")
    plt.plot(mu_values, [Qt * 60000] * len(mu_values), 'k--', label="Qt")
    plt.xlabel("Viscosity [Pa.s]")
    plt.ylabel("Flowrate [LPM]")
    plt.title("Scenario 3: Q1, Q2, Qt vs Viscosity")
    plt.legend()
    plt.grid(True)
    plt.savefig("scenario3_flows.png", dpi=300)
    plt.show()



def main():


    scenario_1()
    scenario_2()
    scenario_3()








if __name__ == '__main__':
    main()