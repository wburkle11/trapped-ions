# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:24:16 2025

@author: iontrap
"""

#4 LEVEL OPTICAL BLOCH EQUATIONS

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants (set to reasonable values for your system)
Gamma_1 = 1e6  # decay rate for state 1
Gamma_2 = 1e6  # decay rate for state 2
Gamma_3 = 1e6  # decay rate for state 3
Gamma_4 = 1e6  # decay rate for state 4

# Rabi frequencies for transitions between states
Omega_12 = 1e6  # Rabi frequency for transition 1 -> 2
Omega_13 = 1e6  # Rabi frequency for transition 1 -> 3
Omega_24 = 1e6  # Rabi frequency for transition 2 -> 4

# Detunings
Delta_12 = 0  # detuning for transition 1 -> 2
Delta_13 = 0  # detuning for transition 1 -> 3
Delta_23 = 0  # detuning for transition 2 -> 3

# System of differential equations
def rho_system(t, y):
    # Extract real and imaginary parts of the coherences
    rho_11, rho_12_real, rho_12_imag, rho_13_real, rho_13_imag, \
    rho_22, rho_23_real, rho_23_imag, rho_33, rho_24_real, rho_24_imag, rho_44 = y

    # Population equations (diagonal elements)
    d_rho_11 = -Gamma_1 * rho_11 + Omega_12 * (rho_12_real + rho_12_imag) + Omega_13 * (rho_13_real + rho_13_imag)
    d_rho_22 = -Gamma_2 * rho_22 + Omega_12 * (rho_12_real - rho_12_imag)
    d_rho_33 = -Gamma_3 * rho_33 + Omega_13 * (rho_13_real - rho_13_imag)
    d_rho_44 = -Gamma_4 * rho_44 + Omega_24 * (rho_23_real - rho_23_imag)

    # Coherence equations (off-diagonal elements)
    gamma_12 = 1e6  # Decay rate for coherence between |1> and |2>
    gamma_13 = 1e6  # Decay rate for coherence between |1> and |3>
    gamma_23 = 1e6  # Decay rate for coherence between |2> and |3|
    gamma_24 = 1e6  # Decay rate for coherence between |2> and |4|

    # rho_12: Coherence between states |1> and |2>
    d_rho_12_real = -1j * Delta_12 * rho_12_imag - gamma_12 * rho_12_real + Omega_12 * (rho_11 - rho_22)
    d_rho_12_imag = 1j * Delta_12 * rho_12_real - gamma_12 * rho_12_imag + Omega_12 * (rho_11 - rho_22)

    # rho_13: Coherence between states |1> and |3>
    d_rho_13_real = -1j * Delta_13 * rho_13_imag - gamma_13 * rho_13_real + Omega_13 * (rho_11 - rho_33)
    d_rho_13_imag = 1j * Delta_13 * rho_13_real - gamma_13 * rho_13_imag + Omega_13 * (rho_11 - rho_33)

    # rho_23: Coherence between states |2> and |3>
    d_rho_23_real = -1j * Delta_23 * rho_23_imag - gamma_23 * rho_23_real + Omega_24 * (rho_22 - rho_44)
    d_rho_23_imag = 1j * Delta_23 * rho_23_real - gamma_23 * rho_23_imag + Omega_24 * (rho_22 - rho_44)

    # rho_24: Coherence between states |2> and |4>
    d_rho_24_real = -gamma_24 * rho_24_real + Omega_12 * (rho_22 - rho_44)
    d_rho_24_imag = -gamma_24 * rho_24_imag + Omega_12 * (rho_22 - rho_44)

    return [d_rho_11, d_rho_12_real, d_rho_12_imag, d_rho_13_real, d_rho_13_imag,
            d_rho_22, d_rho_23_real, d_rho_23_imag, d_rho_33, d_rho_24_real, d_rho_24_imag, d_rho_44]

# Initial conditions (assuming the system starts fully in state |1⟩)
initial_conditions = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Time span and evaluation points
t_span = (0, 1e-6)  # time range for the solution (adjust as needed)
t_eval = np.linspace(0, 1e-6, 1000)

# Solve the differential equations
solution = solve_ivp(rho_system, t_span, initial_conditions, t_eval=t_eval)

# Extract the population of the excited state |2⟩ over time
rho_22 = solution.y[5]  # Populations of state 2

# Plot the population of the excited state
plt.plot(solution.t, rho_22)
plt.xlabel('Time (s)')
plt.ylabel('Population of |2⟩')
plt.title('Excited State Population in a 4-Level System')
plt.show()