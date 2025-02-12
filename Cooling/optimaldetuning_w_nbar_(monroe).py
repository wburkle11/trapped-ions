# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 18:27:53 2025

@author: iontrap
"""

#Lets try recreating exactly what monroe did.... 

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

gamma = 2*(math.pi) * 19.6e6

del0 = np.linspace(4.45, 4.65, 300)*gamma
#delpi = 4.54 * gamma
#print(delpi)

#print(delpi)
zeeman = 2*(math.pi)*7.7e6
del1 = 4.5*gamma
#del_d = np.linspace(4.45, 4.65, 300)*gamma
#print(del_d)
#calc = 2*(np.pi)*80e6
#print(calc + zeeman)

#gamma is total decay rate out of the excited state

#we have 3 detunings set here, 2 for coupling beams, 1 for probe. 
#delsig_p = del_d + zeeman
#delsig_m = del_d - zeeman
omegam1 = 0.8 * gamma
omega1 = 2.0*gamma
omega0 = 0.76*gamma
#omega = omega_sig_p 
omega = omega1

#Excited State Population
nu = 2*(math.pi)*4.45e6
#nu = np.linspace(0, 0.25, 300)*gamma

#print(nu)
x1 = del0 
x2 = del0 + nu
x3 = del0 - nu

def calculate_diff(x):
    
    diff = (x - del1)
    
    return diff

diff_r = calculate_diff(x1)
diff_r2 = calculate_diff(x2)
diff_r3 = calculate_diff(x3)
# print(diff_r)


def compute_pee(x, y):

    Z = 8*y**2*omega1**2*omega0**2*gamma + 2*y**2*gamma**3*omega**2 - 4*x*y*omega1**4*gamma + 1/2*omega**6*gamma + 8*y**2*gamma*(del1**2*omega0**2 + x**2*omega1**2) + 4*del1*y*omega0**2*gamma

    pee = (4*y**2*omega0**2*omega1**2*gamma)/Z
    
    return pee
    

result_car = compute_pee(x1, diff_r)
#print(result_car)
result_rsb = compute_pee(x2, diff_r2)
#print(result_rsb)
result_bsb = compute_pee(x3, diff_r3)
#print(result_bsb)

#result_c1 = compute_pee(nu, diff_r)
#result_c2 = compute_pee(-nu, diff_r)


#Define nbar

#nbar = result_rsb/result_bsb / (1 - result_rsb/result_bsb)
#print(nbar)

nbar = (result_car + result_bsb) / result_rsb - result_bsb
#print(nbar2)


#PLOT OF ABSORBTION PROFILE

#plt.plot(delpi / gamma, result_car)
#plt.axvline(vline/gamma, color='black', linestyle='--')  # Blue dash-dot line
#plt.axvline(del_d/gamma, color='black', linestyle='--')  # Blue dash-dot line
#plt.axvline(del_d/gamma + nu/gamma, color='red', linestyle='--')  # Blue dash-dot line
#plt.axvline(del_d/gamma - nu/gamma, color='blue', linestyle='--')  # Blue dash-dot line
# Add a grid for better readability
#plt.grid(True)

#plt.show() 

#Cooling Bandwidth (Monroe Supp. Eq 16)
#W_c = (1 + np.sqrt(2))*omega_sig_p**2/(1.5 + np.sqrt(2))*del_d
#print(W_c)
#PLOT OF NBAR

# Normalize y-values to be between 0.1 and 10
a, b = 0.1, 10
nbar_rescaled = a + ((nbar - np.min(nbar)) / (np.max(nbar) - np.min(nbar))) * (b - a)

plt.plot(del0 / gamma , nbar_rescaled)
min_index = np.argmin(nbar)
x_min = del0[min_index]
print(x_min / 6.28)
plt.xlabel('Delpi / Gamma')
plt.ylabel('nbar')
#plt.axvline(x= x_min , color = 'r')
# Set both x and y axes to a logarithmic scale
#plt.xscale('log')
plt.yscale('log')

#plt.plot(delpi / gamma , nbar_rescaled)

#plt.axvline(W_c / gamma, color='black', linestyle='--')  # Blue dash-dot line
plt.show()


#plot of trap freq vs nbar
#plt.plot(nu / gamma , nbar_rescaled)
#plt.yscale('log')
#plt.show()

