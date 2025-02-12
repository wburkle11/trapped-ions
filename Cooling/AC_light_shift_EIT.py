# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 21:43:22 2025

@author: iontrap
"""

#CALCULATING AC LIGHT SHIFT FROM KIM SUPPLEMENTAL

import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math
from sympy import symbols, Eq, solve

#DEFINITIONS

delpi = np.linspace(-30*1e6, 80*1e6, 800)* 2*math.pi  #PROBE DETUNING
delpi2 = 2*math.pi*50*1e6


zeeman = 2*math.pi*5*1e6  #Zeeman Splitting 
del_d = 2*math.pi*55*1e6  #DRIVE DETUNING
#print(del_d / (2*math.pi))
 
gamma = 2*math.pi * 19.6*1e6          #total decay rate out of the excited state
delsig_p = del_d + zeeman           #total detuning of sigma_plus component of driving beam
delsig_m = del_d - zeeman           #total detuning of sigma_minus component of driving beam
omega_sig_m = 2*math.pi * 18*1e6      #Rabi of sigma_minus component of driving beam
omega_sig_p = 2*math.pi * 18*1e6      #Rabi of sigma_plus component of driving beam
omega_pi = 2*math.pi * 6*1e6          #Rabi of probe beam
nu = 2*math.pi*2*1e6                  #Trap Freq. we wish to cool
nu2 = np.linspace(1*1e6, 4*1e6, 1000)* 2*math.pi  #PROBE DETUNING

omega = omega_sig_m                 #Parameter for computing pee in Monroe paper (..questionable..)

lambda_ = symbols('lambda')

# Define the equation
equation = Eq(1/4 * lambda_ * (4*lambda_*(lambda_ - delsig_p) * (lambda_ - delsig_m) - 
              (lambda_ - delsig_p) * omega_sig_m**2 - 
              (lambda_ - delsig_m) * omega_sig_p**2), 0)

# Solve the equation for lambda
solutions = solve(equation, lambda_) 

# Extract the second-to-last solution
second_to_last_solution = solutions[-2]

# Print the absolute value of the second-to-last solution
adjusted_sol = abs(second_to_last_solution - delsig_m)  # Subtract del_d and take the absolute value
print(adjusted_sol / (2*math.pi*1e6), 'MHz')  # Print the result in MHz (divide by 2*pi)