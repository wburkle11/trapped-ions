# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 10:46:38 2025

@author: iontrap
"""

#Code for calculating Equilibrium positions and normal modes for 1D Ion chains: 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sp
import scipy.constants as scipy
import qutip as qt
import pandas as pd

#definitons
q_e = scipy.elementary_charge
eps_o = scipy.epsilon_0
M_Yb = 171*1.67262192e-27
delta_k = np.sqrt(2)*2*np.pi/(355e-9)
hbar = scipy.hbar

fz = .58 #MHz   #Trap Freq in Axial direction
fx = 4 #MHz  #Trap Freq in Radial direction
omega_z = 2*np.pi*(fz)*(10**6)  #Angular Frequencies
omega_x = 2*np.pi*(fx)*(10**6)


lc = (q_e**2/(4*np.pi*eps_o*M_Yb*omega_z**2))**(1/3)
recoil = (hbar*(delta_k)**2 ) / (2*M_Yb)


N = 2 #number of ions


def equil_positions(N):

  """
  Calculates equilibrium positions for N ions under harmonic trapping potentials \
  along all spatial axis.

  Parameters:
  N (int): The number of ions in the system.

  Returns:
  chopped_solution (array): Array of equilibrium positions.

  """
  def coupled_equations(u):
      equations = []
      for m in range(N):
          sum1 = sum(1 / (u[m] - u[n])**2 for n in range(m))
          sum2 = sum(1 / (u[m] - u[n])**2 for n in range(m+1, N))
          equations.append(u[m] - sum1 + sum2)
      return equations

  # Generate initial conditions
  initial_guess = [m / 10 for m in range(1, N+1)]

  # Solve the coupled equations
  solution = fsolve(coupled_equations, initial_guess)

  # Chop (close to zero precision) and return the solution
  chopped_solution = np.round(solution, decimals=4)

  return chopped_solution


def scaled_positions(u,lc):
  """
  Calculates the true equilibrium positions of ions in meters.

  Parameters:
  u (array): Array of equilibrium positions.
  lc (float): The lattice constant in meters.

  Returns:
  u*lc (array): Array of equilibrium positions in meters.

  """
  pos = equil_positions(N)

  return pos*lc

def eigensystem(frad,fax,N):
  """
  Calculates the mode spectrum for a given set of ions at specific trap frequencies.

  Parameters:
  frad (float): Radial Trap frequency in Hz.
  fax (float): Axial trap frequency in Hz.
  num_ions (int): The number of ions in the system.
  lc (float): The lattice constant in meters.

  Returns:
  eigenvalues (array): Normal mode frequencies in MHz. Need to be
  multiplied by 10^6 to get in Hz.

  eigenvectors (array): Array of eigenvectors.
  """
  u = equil_positions(N)
  #u = threes
  print(u)


  Anm = np.empty((N,N))

  for n in range(N):
    for m in range(N):
      if n == m:
        sum1 = 0.0
        for p in range(N):
          if(p!=m):
            sum1 += 1/abs(u[m] - u[p])**3
        Anm[n][n]= (frad/fax)**2-sum1
      elif n!=m:
        sum2 = 1/ abs(u[m] - u[n])**3
        Anm[n][m]= sum2
  
  eigenvalues, eigenvectors = np.linalg.eig(Anm)
  eigenvalues = np.sqrt(eigenvalues)*fz
  eigenvectors = eigenvectors.T

  return eigenvalues, eigenvectors


#Compute Evals / Evects 

evals, evects = eigensystem(fx,fz,N)

evals = evals*10**6 #convert to Hz
#print(evals/10e6)
sorted_evals = sorted(evals/10**6)
#choose mode to detune from
COM_freq = np.max(evals) #COM frequency
COM_freq_index = np.argmax(evals) #COM frequency index

# Example data
x_values = sorted_evals

# Create x-values, which will be the indices (0, 1, 2, ...)
y_values = [1]*len(x_values)

# Plot the data
# plt.scatter(x_values, y_values, s=5)
plt.figure(figsize=(14, 4))  # Set figure size
# Plot the data
plt.plot(x_values, np.zeros_like(x_values), 'bo', markersize=2)  # 1D plot using zeros for y-values
x_values_rounded = np.round(x_values, 2)
plt.xticks(x_values_rounded, rotation=90)
plt.tick_params(axis='x', labelsize=10)  # Adjust the fontsize for x-axis numbers

#Normal Mode Bandwidth
bandwidth = sorted_evals[-1] - sorted_evals[0]
plt.annotate(f'Normal Mode Bandwidth: {bandwidth:.2f} MHz', 
             xy=(0.5 * (sorted_evals[-1] + sorted_evals[0]), 0.04), 
             fontsize=10, ha='center', va='bottom', color='blue')

plt.annotate(f'Number of Ions = {N} ', 
             xy=(0.5 * (sorted_evals[-1] + sorted_evals[0]), 0.03), 
             fontsize=10, ha='center', va='bottom', color='blue')

# Hide y-axis
plt.yticks([])  # Remove y-axis ticks
plt.gca().spines['left'].set_visible(False)  # Optionally hide the left spine
plt.axvline(x=x_values[-1], color='blue', label=f'COM = {x_values[-1]:.2f} MHz', linewidth=0.8, linestyle='--')
# Add labels and title (optional)
plt.xlabel('Frequency (MHz)')
plt.title('1D-Chain Normal Modes')
plt.legend()

# Show the plot
plt.show()
#print(sorted_evals[-1])
#print("Normal Mode Freq from CAR:", sorted_evals)

#print("Normal Mode Bandwidth",(sorted_evals[-1]-sorted_evals[0]), "MHz")