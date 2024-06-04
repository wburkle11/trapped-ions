
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


def compute_q_factor(C_t, R_t, b, D, d, d_0):

    ''' Input:  Primary Resonator Parameters; b, D, d, d_0
                Trap Capacitace: C_t
                Trap Resistance: R_t
        
        Output: Qfactor; Q 
    '''
    
    
    #These parameters below are calculated as depending on arguments to the function:
    
    rat3 = d/D
    # Winding Pitch
    t = 2*d_0 + 0.001  
    # Coil Height
    B = b + D / 2  
    # Number of turns in COIL
    N_c = b/t
    l_c = (np.pi)*d*N_c
    # Resistivity of Copper(Au)
    p = 1.7e-8  
    # Capacitance of connecting wires
    C_w = 0.0001e-12  

    

    # Capactiance between Main Coil and Sheild [Siverns Eq 26]
    rat1 = D / d
    K_Cs = 39.37 * (0.75 / (np.log(D / d))) * 1e-12  
    C_s = b * K_Cs  # farads


    # Coil Self Capacitance    [Siverns Eq 25]
    H = 11.26 * (b / d) + 8 + (27 / np.sqrt(b / d))  
    C_C = H * d * 1e-12  # farads


    # Coil Inductance [Siverns Eq 27]
    K_Lc = 39.37 * (0.025 * (d ** 2) * (1 - ((rat3) ** 2)) / (t ** 2)) * 1e-6  
    # Inductance
    L_C = b * K_Lc


    # Predicted Resonant Frequency
    omega_0 = 1 / np.sqrt((C_s + C_t + C_w + C_C) * L_C)



    # Circuit Reactances [Siverns Eq 24]
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
    sigma = np.sqrt((2 * p) / ((omega_0) * 4 * (np.pi) * 1e-7)) 

    # Estimate of Coil-Sheild Junction Resistance [Siverns Eq 36]
    # Length of Solder Junction
    l_j = 0.003  # meters
    # Diameter of Solder Junction
    d_j = 0.003  # meters
    R_j = (p * l_j) / ((np.pi) * d_j * sigma)  


    # Sheild Resistance     [Siverns Eq 35]
    # Number of turn current undergoes in sheild:
    N_s = (b * l_c) / (4 * (np.pi) * (D - d) ** 2)
    # Distance current will travel from bottom to top of sheild
    l_s = N_s * np.sqrt(((D ** 2) * ((np.pi) ** 2)) + (b / N_s) ** 2)
    R_s = (N_s * p * l_s) / (b * sigma)  


    # Coil Resistance  [Sivers Eq 34]
    R_c = (p * l_c) / (d_0 * (np.pi) * sigma)


    # TOTAL RESISTANCE [Siverns Eq 23]
    R_ESRt1 = ((R_c * (X_Cc ** 2)) / (R_c ** 2 + (X1) ** 2))
    R_ESRt2 = ((R_t * (X_R ** 2)) / (R_t ** 2 + (X_R + X_T) ** 2))
    R_ESR = R_ESRt1 + R_ESRt2 + R_s + R_j


    # Impedance Route (should agree with TOTAL RESISTANCE calculation above)
    Z_coil = (1 / (1j * X_Lc + R_c) + 1 / (1j * X_Cc)) ** -1
    Z_E = (1 / (1j * X_Ct + R_t) + 1 / (1j * X_Cw) + 1 / (1j * X_Cs)) ** -1
    Z_total = Z_coil + Z_E + R_s + R_j
    Z_total_real = np.real(Z_total)
    Z_coil_real = np.real(Z_coil)
    Z_E_real = np.real(Z_E)

    # Q FACTOR!!!  [Siverns Eq 22]
    Q = X_Lc / R_ESR
    return Q


