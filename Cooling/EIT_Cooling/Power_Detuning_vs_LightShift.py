# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 13:55:07 2025

@author: iontrap
"""

# CODE for my quick study of the effects on the light shift and cooliing bandwidth as a function of our main exp parameters (detuning, power) 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sp
import scipy.constants as scipy
import qutip as qt
import pandas as pd

x1_det_values = np.arange(15, 76, 5)

#print(x1_det_values)

y1_det_values = [2.79, 2.44, 2.16, 1.93, 1.74, 1.58, 1.45, 1.33, 1.24, 1.15, 1.08, 1.01, 0.95]

y2_det_values = [1.65, 1.76, 1.77, 1.74, 1.68, 1.61, 1.53, 1.46, 1.39, 1.32, 1.26, 1.20, 1.14]

fig = plt.figure(figsize=(12, 5))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))  # figsize adjusts the overall size

# Plot on the first subplot (ax1)
ax1.plot(x1_det_values, y1_det_values)
ax1.set_title('Detuning vs. Light Shift')
ax1.set_xlabel(r'Overlapped Detuning (MHz)')
ax1.set_ylabel('Light Shift (MHz)')

# Plot on the second subplot (ax2)
ax2.plot(x1_det_values, y2_det_values)
ax2.set_title('Detuning vs Cooling Bandwidth')
ax2.set_xlabel('Overlapped Detuning (MHz)')
ax2.set_ylabel('Cooling Bandwidth (MHz)')

# Adjust spacing between the plots
plt.tight_layout()


x2_pow_values = np.arange(15, 55, 3)

#print(x2_pow_values)

y1_pow_values = [1.07, 1.45, 1.83, 2.19, 2.53, 2.84, 3.11, 3.34, 3.54, 3.71, 3.86, 3.98, 4.09, 4.18]

y2_pow_values = [1.18, 1.53, 1.84, 2.08, 2.27, 2.39, 2.47, 2.51, 2.52, 2.52, 2.51, 2.49, 2.47, 2.45]

fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))  # figsize adjusts the overall size

# Plot on the first subplot (ax1)
ax3.plot(x2_pow_values, y1_pow_values)
ax3.set_title('Drive Power vs. Light Shift')
ax3.set_xlabel(r'Drive Rabi Freq (MHz)')
ax3.set_ylabel('Light Shift (MHz)')

# Plot on the second subplot (ax2)
ax4.plot(x2_pow_values, y2_pow_values)
ax4.set_title('Drive Power vs Cooling Bandwidth')
ax4.set_xlabel('Drive Rabi Freq (MHz)')
ax4.set_ylabel('Cooling Bandwidth (MHz)')

# Adjust spacing between the plots
plt.tight_layout()

## Extra Computations 

y1_det_span = np.max(y1_det_values) - np.min(y1_det_values)
y2_pow_span = np.max(y2_pow_values) - np.min(y2_pow_values)
y2_det_span = np.max(y2_det_values) - np.min(y2_det_values)
print(y2_pow_span)