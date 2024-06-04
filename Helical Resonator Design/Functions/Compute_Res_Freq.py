
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


def compute_res_freq(C_t, b, D, d, d_0):
    
    ''' Input:  Primary Resonator Parameters; b, D, d, d_0
                Trap Capacitace: C_t
        
        Output: Resonant Frequency; omega_0 
    '''
    
    #These parameters below are calculated as depending on arguments to the function: 
        
    rat3 = d / D
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
    C_s = b * K_Cs  

   
    # Coil Self Capacitance    [Siverns Eq 25]
    H = 11.26 * (b / d) + 8 + (27 / np.sqrt(b / d))  
    C_C = H * d * 1e-12  

    
    # Coil Inductance [Siverns Eq 27]
    K_Lc = 39.37 * (0.025 * (d ** 2) * (1 - ((rat3) ** 2)) / (t ** 2)) * 1e-6  
    # Inductance
    L_C = b * K_Lc

    
    # Predicted Resonant Frequency
    omega_0 = 1 / np.sqrt((C_s + C_t + C_w + C_C) * L_C)
    omega_0t = omega_0/(2*(np.pi))

    return omega_0t


