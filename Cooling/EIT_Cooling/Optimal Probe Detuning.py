# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:10:03 2024

@author: iontrap
"""
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
#Script for finding the optimal probe detuning: 
    
#Based on the relationship between eigenvalues and trap freq found in the Monroe Paper: 
    
nu = 2*(np.pi)* 1.5e6


omega_sig_m = 2*(np.pi) * 17e6
omega_sig_p = 2*(np.pi) * 17e6
omega_pi = 2*(np.pi) * 4e6

delta = np.linspace(10e6, 70e6, 500)* 2*(np.pi)
#delta = 2*(np.pi) * 30e6


nu1 = (np.sqrt(delta**2 + omega_pi**2 + omega_sig_m**2))/2 - delta/2

#print(nu1)


plt.plot(delta/ 2*(np.pi), nu1/ 2*(np.pi))

plt.axhline(nu/ 2*(np.pi), color='black', linestyle='--')  


# Add labels and title
plt.xlabel('delta')
plt.ylabel('Trap Freq')
plt.title('Optimal Probe Detuning')


# Add a grid for better readability
plt.grid(True)

plt.show() 