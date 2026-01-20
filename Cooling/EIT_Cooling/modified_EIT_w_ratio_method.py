# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 14:11:21 2025

@author: iontrap
"""

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


delpi = np.linspace(-30e6, 80e6, 500)* 2*(np.pi)
#delpi = 2*(np.pi)*60e6

#print(delpi)
zeeman = 2*(np.pi)*5e6
del_d = 2*(np.pi)*55e6

#gamma is total decay rate out of the excited state

#we have 3 detunings set here, 2 for coupling beams, 1 for probe. 
gamma = 2*(np.pi) * 19.6e6
delsig_p = del_d + zeeman
delsig_m = del_d - zeeman
omega_sig_m = 2*(np.pi) * 17e6
omega_sig_p = 2*(np.pi) * 17e6
omega_pi = 2*(np.pi) * 4e6
nu = 2*(np.pi)*2.5e6
#print(del_d/2*(np.pi))

x1 = delpi
x2 = delpi + nu
x3 = delpi - nu

def compute_peekim(x):

    Z = 4*gamma**2 * (x - delsig_m)**2*(x - delsig_p)**2 + (4*x*(x - delsig_p)*(x - delsig_m) - (x - delsig_p)* omega_sig_m**2 - (x - delsig_m)* omega_sig_p**2)**2

    W = 16*(x - delsig_m)**2 * (x - delsig_p)**2/Z
    
    return W

result_car = compute_peekim(x1)
#print(result_car)
result_rsb = compute_peekim(x2)
#print(result_rsb)
result_bsb = compute_peekim(x3)
#print(result_bsb)

#PLOTTING FANO ABS PROFILE

plt.plot(delpi / 2*(np.pi), result_car)

plt.axvline(delsig_p/2*(np.pi), color='black', linestyle='--')  # Red dashed line
plt.axvline(delsig_m/2*(np.pi), color='black', linestyle='--')  # Blue dash-dot line

# Add labels and title
plt.xlabel('delpi')
plt.ylabel('ABS.')
plt.title('EIT Fano Profile')

# Add a grid for better readability
plt.grid(True)

plt.show() 


#PLOTTING NBAR

nbar = (result_car + result_bsb)/(result_rsb - result_bsb)

# Normalize y-values to be between 0.1 and 10
#a, b = 0.1, 10
#nbar_rescaled = a + ((nbar - np.min(nbar)) / (np.max(nbar) - np.min(nbar))) * (b - a)

plt.plot(delpi / 2*(np.pi), nbar)
plt.axhline(0, color='black', linestyle='--')  # Red dashed line
#plt.xlim(40e6 , 64e6)
# Set both x and y axes to a logarithmic scale
#plt.xscale('log')
#plt.yscale('log')

plt.show()


#LIGHT SHIFT (Roos EQ 1)

sig = 1/2*(np.sqrt(omega_sig_p**2 + del_d**2) - np.abs(del_d))

#print(sig / 2*(np.pi))

#plt.plot(sig / 2*(np.pi), nbar)
#plt.show()

