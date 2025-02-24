# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 15:34:04 2025

@author: iontrap
"""
#EIT ABSORBTION SPECTRUM

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


delpi = np.linspace(-30e6, 70e6, 500)* 2*(np.pi)
#delpi = 2*(np.pi)*63e6

#print(delpi)
zeeman = 2*(np.pi)*5e6
del_d = 2*(np.pi)*60e6

#gamma is total decay rate out of the excited state

#we have 3 detunings set here, 2 for coupling beams, 1 for probe. 
gamma = 2*(np.pi) * 19.6e6
delsig_p = del_d + zeeman
delsig_m = del_d - zeeman
omega_sig_m = 2*(np.pi) * 15e6
omega_sig_p = 2*(np.pi) * 15e6
omega_pi = 2*(np.pi) * 6e6
nu = 2*(np.pi)*4e6
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
