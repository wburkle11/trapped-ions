# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 14:06:31 2024

@author: iontrap
"""

### This script is for calculating a Fano-Like absorbtion profile for EIT cooling:
    
### MAIN SOURCE: Monroe, Supplemental Material (Effecient Single EIT Cooling)

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


delpi = np.linspace(-50e6, 90e6, 1000)* 2*(np.pi)

#print(delpi)


gamma = 2*(np.pi) * 19.6e6
delsig_p = gamma*3.69
delsig_m = gamma*4.47
omega_sig_m = gamma*0.2
omega_sig_p = gamma*0.7
omega_pi = gamma*0.35


#EQ 14 (Kim, Supplemental) 

Z = 4*gamma**2 * (delpi - delsig_m)**2*(delpi - delsig_p)**2 + (4*delpi*(delpi - delsig_p)*(delpi - delsig_m) - (delpi - delsig_p)* omega_sig_m**2 - (delpi - delsig_m)* omega_sig_p**2)**2

W = 16*(delpi - delsig_m)**2 * (delpi - delsig_p)**2/Z

print(W)

plt.plot(delpi / gamma, W)

# Add labels and title
plt.xlabel('delpi')
plt.ylabel('ARB.')
plt.title('EIT Fano Profile')

# Add a grid for better readability
plt.grid(True)

plt.show() 