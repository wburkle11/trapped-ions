# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 16:34:07 2025

@author: iontrap
"""

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.ticker import ScalarFormatter


delpi = np.linspace(20e6, 65e6, 1000)* 2*(np.pi)
#delpi = 2*(np.pi)*60e6
#print(delpi / 2*np.pi)

#print(delpi)
zeeman = 2*(np.pi)*5e6
del_d = 2*(np.pi)*65e6
#del_d = np.linspace(38.0e6, 42.0e6, 500)* 2*(np.pi)

#gamma is total decay rate out of the excited state

#we have 3 detunings set here, 2 for coupling beams, 1 for probe. 
gamma = 2*(np.pi) * 19.6e6
delsig_p = del_d + zeeman
delsig_m = del_d - zeeman
omega_sig_m = 2*(np.pi) * 8e6
omega_sig_p = 2*(np.pi) * 8e6
omega_pi = 2*(np.pi) * 6e6
nu = 2*(np.pi)*3e6
#print(del_d/2*(np.pi))

x1 = delpi
x2 = delpi + nu
x3 = delpi - nu

Z1 = 4*gamma**2 * (x1 - delsig_m)**2*(x1 - delsig_p)**2 + (4*x1*(x1 - delsig_p)*(x1 - delsig_m) - (x1 - delsig_p)* omega_sig_m**2 - (x1 - delsig_m)* omega_sig_p**2)**2

W1 = 16*(x1 - delsig_m)**2 * (x1 - delsig_p)**2/Z1

Z2 = 4*gamma**2 * (x2 - delsig_m)**2*(x2 - delsig_p)**2 + (4*x2*(x2 - delsig_p)*(x2 - delsig_m) - (x2 - delsig_p)* omega_sig_m**2 - (x2 - delsig_m)* omega_sig_p**2)**2

W2 = 16*(x2 - delsig_m)**2 * (x2 - delsig_p)**2/Z2

Z3 = 4*gamma**2 * (x3 - delsig_m)**2*(x3 - delsig_p)**2 + (4*x3*(x3 - delsig_p)*(x3 - delsig_m) - (x3 - delsig_p)* omega_sig_m**2 - (x3 - delsig_m)* omega_sig_p**2)**2

W3 = 16*(x3 - delsig_m)**2 * (x3 - delsig_p)**2/Z3

#print(W1)
#print(W2)
#print(W3)

nbar = W1 + W3 / W2 - W3

#print(nbar)

# Normalize y-values to be between 0.1 and 10
a, b = 0.1, 10
nbar_rescaled = a + ((nbar - np.min(nbar)) / (np.max(nbar) - np.min(nbar))) * (b - a)

min_index = np.argmin(nbar)
x_min = delpi[min_index]

#print(min_index)
print(nbar[min_index] , "Min Nbar")
diff = (delpi - del_d / 2*np.pi)
#print(diff)

plt.plot((delpi / 2*np.pi), nbar)
plt.axvline(x= x_min / 2*np.pi, color = 'r')
plt.xlabel('Delpi / 2pi (MHz)')
plt.ylabel('nbar')

optimal_del = x_min / 6.28
print(optimal_del, "Optimal Detuning")

#print(2*math.pi)

# Set both x and y axes to a logarithmic scale
#plt.xscale('log')
#plt.yscale('log')

plt.show()



