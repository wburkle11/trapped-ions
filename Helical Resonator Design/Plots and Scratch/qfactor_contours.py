#Code used for making contour plots (Siverns Figure 10)

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import plotly.graph_objects as gols

from matplotlib.colors import LinearSegmentedColormap


# def the_ratio(d):
    
#     #sheild diameter
#     D = 0.108
#     ratio = d/D
#     return ratio

# d = np.linspace(0.001, 0.20, 20)
# big = the_ratio(d)
# big1 = big/2

def custom_colormap():
    magma_cmap = plt.get_cmap('gist_rainbow')
    colors = magma_cmap(np.linspace(0, 1, 256))

    # Change color of the first index to white
    colors[0] = [0, 0, 0, 1]  # RGB values with alpha = 1

    return LinearSegmentedColormap.from_list('custom_colormap', colors)

def compute_q_factor_2(C_t, d, D, R_t, omega_0):

    # Helical Resonator Dimensions
    # Sheild Height:
    B = 0.12  # meters
    # Sheild Diameter
    #D = 0.108  # meters
    # Main Coil Diameter
    #d = 0.01  # meters
    rat3 = d / D
    #print(rat3)
    # Winding Pitch
    t = 0.009  # meters
    # Diameter of Main Coil Wire
    d_0 = 0.005
    # Uncoiled Main Coil Length
    l_c = 0.88  # meters
    # Number of turns in COIL
    N_c = 6.75
    # Coil Height
    # b = B - D / 2  # meters
    b = 0.066
    # print(b)
    # Resistivity of Copper(Au)
    p = 1.7e-8  # ohm-meter
    # Capacitance of connecting wires
    C_w = 0.0001e-12  # Farads

    # EXTRA CAPACITANCES FROM CPCTV. DIVIDER
    C_d1 = 2e-13
    C_d2 = 2e-11

    C_pick = 1 / C_d1 + 1 / C_d2
    C_pick1 = 1 / C_pick
    # print(C_pick1)

    # Capactiance between Main Coil and Sheild [Siverns Eq 26]
    #big2 = 1/big1
    rat2 = D/d
    K_Cs = 39.37 * (0.75 / (np.log(rat2))) * 1e-12  # farads/meter
    C_s = b * K_Cs  # farads

    # Coil Self Capacitance    [Siverns Eq 25]
    H = 11.26 * (b / d) + 8 + (27 / np.sqrt(b / d))  # farads/meter
    C_C = H * d * 1e-12  # farads

    # Coil Inductance [Siverns Eq 27]
    K_Lc = 39.37 * (0.025 * (d ** 2) * (1 - ((rat3) ** 2)) / (t ** 2)) * 1e-6  # henrys/meter
    # Inductance
    L_C = b * K_Lc
    #print(L_C)
    #jay = (C_s + C_t + C_w + C_C) * L_C
    #print (jay)

    # Predicted Resonant Frequency
    #omega_0 = 1 / np.sqrt((C_s + C_t + C_w + C_C) * L_C)

    # Predicted Value of b:
    C_gamma = C_t + C_w
    K_cd = 35 * d * 1e-12
    K_cb = 11.26 * 1e-12
    b1 = (C_gamma + K_cd) / (K_Cs + K_cb)
    b2 = np.sqrt((K_Cs + K_cb) / ((C_gamma + K_cd) ** 2 * K_Lc * omega_0 ** 2) + 1 / 4) - 1 / 2
    b3 = b1 * b2

    # Checking Validity of [Siverns Eq 24]
    X_Cc = 1 / (omega_0 * C_C)
    X_Ct = 1 / (omega_0 * C_t)
    X_Cw = 1 / (omega_0 * C_w)
    X_Cs = 1 / (omega_0 * C_s)
    # Reactance of Coil Inductance
    X_Lc = omega_0 * L_C
    X_R = 1 / (omega_0 * (C_s + C_w))
    X_T = 1 / (omega_0 * C_t)
    X1 = X_Lc + X_Cc

    # Skin Depth of Copper (function of omega_0 in Hz)
    sigma = np.sqrt((2 * p) / ((omega_0) * 4 * (np.pi) * 1e-7))  # meters

    # Estimate of Coil-Sheild Junction Resistance [Siverns Eq 36]
    # Length of Solder Junction
    l_j = 0.003  # meters
    # Diameter of Solder Junction
    d_j = 0.003  # meters
    R_j = (p * l_j) / ((np.pi) * d_j * sigma)  # ohms

    # Sheild Resistance     [Siverns Eq 35]
    # Number of turn current undergoes in sheild:
    N_s = (b * l_c) / (4 * (np.pi) * (D - d) ** 2)
    # Distance current will travel from bottom to top of sheild
    l_s = N_s * np.sqrt(((D ** 2) * ((np.pi) ** 2)) + (b / N_s) ** 2)
    R_s = (N_s * p * l_s) / (b * sigma)  # ohms

    # Coil Resistance  [Sivers Eq 34]
    R_c = (p * l_c) / (d_0 * (np.pi) * sigma)

    # TOTAL RESISTANCE [Siverns Eq 23]
    R_ESRt1 = ((R_c * (X_Cc ** 2)) / (R_c ** 2 + (X1) ** 2))
    R_ESRt2 = ((R_t * (X_R ** 2)) / (R_t ** 2 + (X_R + X_T) ** 2))
    R_ESR = R_ESRt1 + R_ESRt2 + R_s + R_j

    # Impedance Route
    #Z_coil = (1 / (1j * X_Lc + R_c) + 1 / (1j * X_Cc)) ** -1
    #Z_E = (1 / (1j * X_Ct + R_t) + 1 / (1j * X_Cw) + 1 / (1j * X_Cs)) ** -1
    #Z_total = Z_coil + Z_E + R_s + R_j
    #Z_total_real = np.real(Z_total)
    #Z_coil_real = np.real(Z_coil)
    #Z_E_real = np.real(Z_E)

    # Q FACTOR!!!  [Siverns Eq 22]
    Q = X_Lc / R_ESR
    return Q
   

#C_t = 1e-12 * np.linspace(0.0001, 100, 100)
#d=0.01
d = np.linspace(0.001, 0.2, 20)
D = np.linspace(1, d, 20)
#print(D)
q_factors_r1 = compute_q_factor_2(5e-12, d, D, 1.0, 26e6)
total_q = np.abs(q_factors_r1)
x_ax = d/D
print(d)
print(x_ax)
print(total_q)
#print(big1)
#print(d)
# ndim = total_q.ndim
# print(ndim)
#print(len(total_q))
#print(len(big1))
#print(len(d))
#print(q_factors_r1)
#q_factors_r3 = compute_q_factor_2(20e-12, d, 3.0, 1.256e8)
#q_factors_r5 = compute_q_factor_2(50e-12, d, 5.0, 1.256e8)
# Create a meshgrid for x and y
rev = d[::-1]

# Define the constant z value
#Z = np.ones_like(X) * 250  # Setting z to 10

#print(Z)

# Plot the heatmap
# plt.figure(figsize=(8, 6))
# plt.imshow(total_q, cmap='viridis', interpolation='nearest')
# plt.colorbar(label='Value')
# plt.title('Heatmap Plot')
# plt.xlabel('X')
# d = plt.xticks()
# big1 = plt.yticks()
# plt.ylabel('Y')
# plt.show()

# Plotting
# Labels


# Create custom colormap
cmap = custom_colormap()

# Heat map
fig, ax = plt.subplots(figsize = (8, 8))
im = ax.imshow(total_q, cmap=cmap)

# Add the color bar
cbar = ax.figure.colorbar(im, ticks=np.arange(0, np.max(total_q)+1, 20), ax = ax)
cbar.ax.set_ylabel("Radial Secular Frequency (MHz)", rotation = -90, va = "bottom")

# Add the labels
ax.set_xticks(np.arange(len(D)), labels=D)
ax.set_yticks(np.arange(len(rev)), labels=rev)

# naming the x axis
plt.xlabel('RF Voltage Amplitude (V)')
# naming the y axis
plt.ylabel('Resonant Frequency (MHz)')

# Rotate the labels of the X-axis
plt.setp(ax.get_xticklabels(), rotation=40,
         ha="right", rotation_mode="anchor")
plt.show()

# fig, ax = plt.subplots(1, 1)

# ax.contour(X, Y, total_q)
 
# ax.set_title('Contour Plot') 
# ax.set_xlabel('feature_x') 
# ax.set_ylabel('feature_y') 
  
# plt.show() 