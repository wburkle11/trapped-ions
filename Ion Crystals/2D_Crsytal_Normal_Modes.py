# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 12:00:55 2025

@author: iontrap
"""

#Code for calculating Equilibrium positions and normal modes for 2D Ion Crystals: 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sp
import scipy.constants as scipy
import qutip as qt
import pandas as pd
import math 
from scipy.constants import epsilon_0, pi, elementary_charge
from scipy.linalg import eigh


#IMPORT 2D Equil Positions from ILYOUNGS Code
file_path = r'Z:\Users\Wes\Equilibrium_Position_Data\2DPos2.csv'
# Load the CSV file into a DataFrame
df = pd.read_csv(file_path, header=None)

# Convert DataFrame to a NumPy array (matrix)
matrix = df.to_numpy()

#print(matrix)

def axial_modes(omega_z, final_positions, m_Yb):
    
    # Define constants
    qe = elementary_charge  # Charge of an electron (same as in Mathematica)
    ell = (qe**2 / (8 * pi * epsilon_0 * m_Yb * omega_z**2))**(1/3)
    
    # Scale positions
    scaled_pos = np.array(final_positions) / ell
    scaled_pos = np.where(np.abs(scaled_pos) < 1e-8, 0, scaled_pos)
    #print(scaled_pos)
    n_ions = len(final_positions)
    
    # Construct Kzz matrix
    Kzz = np.zeros((n_ions, n_ions))
    for i in range(n_ions):
        for j in range(n_ions):
            if i != j:
                dx = scaled_pos[i, 0] - scaled_pos[j, 0]
                dy = scaled_pos[i, 2] - scaled_pos[j, 2]
                Kzz[i, j] = 1 / (dx**2 + dy**2)**(1.5)
            else:
                Kzz[i, i] = 1 - sum(
                    1 / ( (scaled_pos[k, 0] - scaled_pos[j, 0])**2 +
                           (scaled_pos[k, 2] - scaled_pos[j, 2])**2 )**(1.5)
                    for k in range(n_ions) if k != j
                )
    
    # Compute eigenvalues and eigenvectors
    eigvals, eigvecs = eigh(Kzz)
    
    # Convert eigenvalues to frequencies (MHz)
    z_mode_freqs = eigvals * omega_z / (2 * pi) * 1e-6
    
    return z_mode_freqs, eigvecs



# Example Usage:
M_Yb = 171*1.67262192e-27
frad = 1.5e6 * 2*math.pi  # Trap Frequency
positions = matrix  # Example particle positions
N = len(positions)

evals, evecs = axial_modes(frad, positions, M_Yb)
evals_rounded = np.round(evals, 2)


# Plot the data
# plt.scatter(x_values, y_values, s=5)
plt.figure(figsize=(14, 4))  # Set figure size
# Plot the data
plt.plot(evals_rounded, np.zeros_like(evals_rounded), 'bo', markersize=2)  # 1D plot using zeros for y-values
x_values_rounded = np.round(evals_rounded, 2)
plt.xticks(evals_rounded, rotation=90)
plt.tick_params(axis='x', labelsize=10)  # Adjust the fontsize for x-axis numbers

#Normal Mode Bandwidth
bandwidth = evals_rounded[-1] - evals_rounded[0]
plt.annotate(f'Normal Mode Bandwidth: {bandwidth:.2f} MHz', 
             xy=(0.5 * (evals_rounded[-1] + evals_rounded[0]), 0.04), 
             fontsize=10, ha='center', va='bottom', color='blue')

plt.annotate(f'Number of Ions = {N} ', 
             xy=(0.5 * (evals_rounded[-1] + evals_rounded[0]), 0.03), 
             fontsize=10, ha='center', va='bottom', color='blue')

# Hide y-axis
plt.yticks([])  # Remove y-axis ticks
plt.gca().spines['left'].set_visible(False)  # Optionally hide the left spine
plt.axvline(x=evals_rounded[-1], color='blue', label=f'COM = {evals[-1]:.2f} MHz', linewidth=0.8, linestyle='--')
# Add labels and title (optional)
plt.xlabel('Frequency (MHz)')
plt.title('2D-Crystal Normal Modes')
plt.legend()

# Show the plot
plt.show()    