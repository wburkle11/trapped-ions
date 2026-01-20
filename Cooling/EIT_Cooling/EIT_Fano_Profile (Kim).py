# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 11:17:55 2024

@author: iontrap
"""

### This script is for calculating a Fano-Like absorbtion profile for EIT cooling:
    
### MAIN SOURCE: Kim, Supplemental Material (Double EIT Cooling)

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


delpi = np.linspace(-70e6, 30e6, 500)* 2*(np.pi)

#print(delpi)
zeeman = 2*(np.pi)*5e6
del_d = -1*2*(np.pi)*40e6

#gamma is total decay rate out of the excited state

#we have 3 detunings set here, 2 for coupling beams, 1 for probe. 
gamma = 2*(np.pi) * 19.6e6
delsig_p = del_d + zeeman
delsig_m = del_d - zeeman
omega_sig_m = 2*(np.pi) * 13e6
omega_sig_p = 2*(np.pi) * 13e6
omega_pi = 2*(np.pi) * 4e6
print(delsig_p / 2*(np.pi))


#EQ 14 (Kim, Supplemental) 

Z = 4*gamma**2 * (delpi - delsig_m)**2*(delpi - delsig_p)**2 + (4*delpi*(delpi - delsig_p)*(delpi - delsig_m) - (delpi - delsig_p)* omega_sig_m**2 - (delpi - delsig_m)* omega_sig_p**2)**2

W = 16*(delpi - delsig_m)**2 * (delpi - delsig_p)**2/Z

#print(W)

plt.plot(delpi / 2*(np.pi), W)

plt.axvline(delsig_p/2*(np.pi), color='black', linestyle='--')  # Red dashed line
plt.axvline(delsig_m/2*(np.pi), color='black', linestyle='--')  # Blue dash-dot line


# Add labels and title
plt.xlabel('delpi')
plt.ylabel('ABS.')
plt.title('EIT Fano Profile')


# Add a grid for better readability
plt.grid(True)

plt.show() 