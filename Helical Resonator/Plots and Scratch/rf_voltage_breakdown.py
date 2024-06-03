# This notebook is for caluclating the resonant frequency in which we should expect some voltage breakdown.

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#Constants

omega_r = 1e6 #1 MHz
Q = 1.602*10**-19
m = 171*1.67*10**-27
U_0 = 1
kappa = 1
resonant = 2.4e7
V_0 = 100

#From Olmshenck Thesis
#d_0 = 0.00021  #310 micron
#z_0 = 0.0013    #100 micron

#For FiberTrap
d_0 = 0.00031
z_0 = 1e-4

#From Phil
a = Q**2/(2*m**2*d_0**4)
b = (kappa*U_0)/z_0**2
beta = math.sqrt(a-b)
#print(beta)

#From Olmshenck
beta_0 = Q/(2**(1/2)*m*d_0**2)
#print(beta_0)

omega_t = (V_0/omega_r)*beta
#print(V_0, '(Volts) -Break Down Vrf')
#print(omega_t*1e6, '(MHz)')

god = resonant*omega_r/beta
#print(god)


