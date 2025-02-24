# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 11:43:00 2024

@author: iontrap
"""

#optimal probe detuning from ROOS paper: 
    
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

delta = np.linspace(10e6, 100e6, 500)* 2*(np.pi)

rabi_sigma = 2*(np.pi) * 17e6


sig = 1/2 * (np.sqrt(rabi_sigma**2 + delta**2) - np.abs(delta))

plt.plot(delta/ 2*(np.pi), sig/ 2*(np.pi))
plt.axhline(4500000, color = "black")

#plt.axhline(nu/ 2*(np.pi), color='black', linestyle='--')  


# Add labels and title
plt.xlabel('delta')
plt.ylabel('Light Shift')
plt.title('Optimal Probe Detuning')


# Add a grid for better readability
plt.grid(True)

plt.show() 

zed = 70e6
print(zed)