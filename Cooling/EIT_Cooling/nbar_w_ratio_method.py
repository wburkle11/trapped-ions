# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 16:51:28 2025

@author: iontrap
"""

## FINAL optimal detuning code: 
    
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

delpi = np.linspace(50e6, 120e6, 40)* 2*(np.pi)
#delpi = 2*(np.pi)*60e6

#print(delpi)
zeeman = 2*(np.pi)*5e6
del_d = 2*(np.pi)*55e6

diff = (delpi - del_d) / 2*(np.pi)

#gamma is total decay rate out of the excited state

#we have 3 detunings set here, 2 for coupling beams, 1 for probe. 
gamma = 2*(np.pi) * 19.6e6
delsig_p = del_d + zeeman
delsig_m = del_d - zeeman
omega_sig_m = 2*(np.pi) * 17e6
omega_sig_p = 2*(np.pi) * 17e6
omega_pi = 2*(np.pi) * 4e6

#secular frequency
nu = 5e6
rsb = (delpi + nu)
bsb = (delpi - nu)
car = delpi

def scattering_amp(x):
    #EQ 14 (Kim, Supplemental) 

    Z = 4*gamma**2 * (x - delsig_m)**2*(x - delsig_p)**2 + (4*x*(x - delsig_p)*(x - delsig_m) - (x - delsig_p)* omega_sig_m**2 - (x - delsig_m)* omega_sig_p**2)**2

    W = 16*(x - delsig_m)**2 * (x - delsig_p)**2/Z
    
    return W

result_bsb = scattering_amp(bsb)
result_rsb = scattering_amp(rsb)
result_car = scattering_amp(car)


#avg phonon number

nbar = (result_car + result_bsb)/(result_rsb - result_bsb)

#print(nbar)

#nbar_normalized = (nbar - np.min(nbar)) / (np.max(nbar) - np.min(nbar))

plt.plot(diff, nbar)

#plt.yscale('log')

# Add a grid for better readability
plt.grid(True)

plt.show() 

