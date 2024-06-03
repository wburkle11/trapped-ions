import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


def compute_res_freq(C_t, b, D, d, d_0):

    # Helical Resonator Dimensions
    # Sheild Height:
    #B = 0.189  # meters
    # Sheild Diameter
    #D = 0.09997  # meters
    # Main Coil Diameter
    #d = 0.020  # meters
    # Diameter of Main Coil Wire
    #d_0 = 0.003
    
    rat3 = d / D
    # Winding Pitch
    t = 2*d_0 + 0.001  # meters
    # Coil Height
    B = b + D / 2  # meters
    # Number of turns in COIL
    N_c = b/t
    l_c = (np.pi)*d*N_c
    # Resistivity of Copper(Au)
    p = 1.7e-8  # ohm-meter
    # Capacitance of connecting wires
    C_w = 0.0001e-12  # Farads


    # Capactiance between Main Coil and Sheild [Siverns Eq 26]
    rat1 = D / d
    K_Cs = 39.37 * (0.75 / (np.log(D / d))) * 1e-12  # farads/meter
    C_s = b * K_Cs  # farads

    # Coil Self Capacitance    [Siverns Eq 25]
    H = 11.26 * (b / d) + 8 + (27 / np.sqrt(b / d))  # farads/meter
    C_C = H * d * 1e-12  # farads

    # Coil Inductance [Siverns Eq 27]
    K_Lc = 39.37 * (0.025 * (d ** 2) * (1 - ((rat3) ** 2)) / (t ** 2)) * 1e-6  # henrys/meter
    # Inductance
    L_C = b * K_Lc

    # Predicted Resonant Frequency
    omega_0 = 1 / np.sqrt((C_s + C_t + C_w + C_C) * L_C)
    omega_0t = omega_0/(2*(np.pi))

    return omega_0t


#C_t = 1e-12 * np.linspace(0.0001, 100, 100)
#res_freq = compute_res_freq(C_t)

#fig, ax = plt.subplots(figsize = (10, 8))
#ax.plot(C_t, res_freq, '-r')
#plt.ylim(0, 80000000)
#ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.yaxis.set_minor_locator(AutoMinorLocator())
#plt.xlabel('Trap Capacitance [F]')
#plt.ylabel('Resonant Frequency (Hz)')

#value = compute_res_freq(0.3e-10)
#print('For 30pF -> Resonant Frequency =', value*10**-6, 'MHz')

# Adding gridlines
#plt.grid()
#plt.show()