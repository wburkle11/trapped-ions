import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


def compute_res_freq(C_t):
    # Helical Resonator Dimensions
    # Sheild Height:
    B = 0.190  # meters
    # Sheild Diameter
    D = 0.160  # meters
    # Main Coil Diameter
    d = 0.085  # meters
    rat3 = d / D
    # Winding Pitch
    t = 0.014  # meters
    # Diameter of Main Coil Wire
    d_0 = 0.004
    # Uncoiled Main Coil Length
    l_c = 2.1  # meters
    # Number of turns in COIL
    N_c = 7.9
    # Coil Height
    b = B - D / 2  # meters
    # print(b)
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


C_t = 1e-12 * np.linspace(0.0001, 100, 100)
res_freq = compute_res_freq(C_t)

fig, ax = plt.subplots(figsize = (10, 8))
ax.plot(C_t, res_freq, '-r')
plt.ylim(0, 80000000)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.xlabel('Trap Capacitance [F]')
plt.ylabel('Resonant Frequency (Hz)')

# Coordinates of the point (x, y)
#x = 0.3e-10
#y = 1.5e7
# Plot the point
#plt.plot(x, y, 'bo')  # 'bo' stands for blue circle marker
#c = 0.5e7
#plt.errorbar(x, y, yerr=c, fmt="o", color="b")

# Adding gridlines
plt.grid()
plt.show()