## CODE WRITTEN FOR MAXIMIZING OUR QFACTOR 

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def max_qfactor3(C_t, R_t, t, D, b, d, d_0):
    
    # Helical Resonator Dimensions
    # Sheild Height:
    #B = 0.189  # meters
    # Sheild Diameter
    #D = 0.09997  # meters
    # Main Coil Diameter
    #d = 0.070  # meters
    
    # Diameter of Main Coil Wire
    #d_0 = 0.003
    
    
    rat3 = d / D
    # Winding Pitch
    #t = 2*d_0  # meters
    # Coil Height
    B = b + D / 2  # meters
    # Number of turns in COIL
    N_c = b/t
    l_c = (np.pi)*d*N_c
    # Resistivity of Copper(Au)
    p = 1.7e-8  # ohm-meter
    # Capacitance of connecting wires
    C_w = 0.0001e-12  # Farads
    
    
    
   
   
    
   #Length of rods in Vaccum
    l = 0.20   #meters
    #Radius of rods
    r = 5e-5   #meters
    #Resistivity of Copper(Au)
    p = 1.7e-8   #ohm-meter
    #Trap Capacitance
    #C_t = 35e-12   #Farads
    #Capacitance of connecting wires
    C_w = 0.0001e-12   #Farads

    #Capactiance between Main Coil and Sheild [Siverns Eq 26]
    rat1 = D/d
    K_Cs = 39.37*(0.75/(np.log(D/d)))*1e-12        #farads/meter
    C_s = b*K_Cs                                     #farads
    #print(C_s)


    #Coil Self Capacitance    [Siverns Eq 25]
    H = 11.26*(b/d)+8+(27/np.sqrt(b/d))   #farads/meter
    C_C = H*d*1e-12   #farads


    #Coil Inductance [Siverns Eq 27]
    K_Lc = 39.37*(0.025*(d**2)*(1-((rat3)**2))/(t**2))*1e-6  #henrys/meter
    #Inductance
    L_C = b*K_Lc
    #print(L_C)


    #Predicted Resonant Frequency
    omega_0 = 1/np.sqrt((C_s + C_t + C_w + C_C)*L_C)
    testing = omega_0/(2*(np.pi))
    #print(testing, 'Predicted Resonant Frequency Hz [Siverns Eq 21]')

    #Predicted Value of b:
    C_gamma = C_t + C_w
    K_cd = 35*d*1e-12
    K_cb = 11.26*1e-12
    b1 = (C_gamma + K_cd)/(K_Cs + K_cb)
    b2 = np.sqrt((K_Cs + K_cb)/((C_gamma + K_cd)**2*K_Lc*omega_0**2) +1/4) - 1/2
    b3 = b1*b2


    #Checking Validity of [Siverns Eq 24]
    X_Cc = 1/(omega_0*C_C)
    X_Ct = 1/(omega_0*C_t)
    X_Cw = 1/(omega_0*C_w)
    X_Cs = 1/(omega_0*C_s)
    #Reactance of Coil Inductance
    X_Lc = omega_0*L_C
    X_R = 1/(omega_0*(C_s+C_w))
    X_T = 1/(omega_0*C_t)
    X1 = X_Lc + X_Cc


    #Skin Depth of Copper (function of omega_0 in Hz)
    sigma = 2*np.sqrt((p)/((omega_0)*4*(np.pi)*1e-7))   #meters


    #Trap Resistance [Siverns Eq 33]
    A = (np.pi)*r**2 
    R = (p*l)/A                          #ohms
    #R_t = 3                           #ohms


    #Estimate of Coil-Sheild Junction Resistance [Siverns Eq 36]
    #Length of Solder Junction
    l_j = 0.003                        #meters
    #Diameter of Solder Junction
    d_j = 0.003                         #meters
    R_j = (p*l_j)/((np.pi)*d_j*sigma)  #ohms


    #Sheild Resistance     [Siverns Eq 35]
    #Number of turn current undergoes in sheild: 
    N_s = (b*l_c)/(4*(np.pi)*(D-d)**2)
    #Distance current will travel from bottom to top of sheild
    l_s = N_s*np.sqrt(((D**2)*((np.pi)**2))+(b/N_s)**2)
    R_s = (N_s*p*l_s)/(b*sigma)   #ohms


    #Coil Resistance  [Sivers Eq 34]
    R_c = (p*l_c)/(d_0*(np.pi)*sigma)


    #TOTAL RESISTANCE [Siverns Eq 23]
    R_ESRt1 = ((R_c*(X_Cc**2))/(R_c**2 + (X1)**2))
    R_ESRt2 = ((R_t*(X_R**2))/(R_t**2 + (X_R + X_T)**2))
    R_ESR = R_ESRt1 + R_ESRt2 + R_s + R_j 
    #print(R_ESRt1, R_ESRt2)


    # #Impedance Route (this agrees with equations above)
    # Z_coil = (1/(1j*X_Lc + R_c) + 1/(1j*X_Cc))**-1
    # Z_E = (1/(1j*X_Ct + R_t) + 1/(1j*X_Cw) + 1/(1j*X_Cs))**-1
    # Z_total = Z_coil + Z_E + R_s + R_j
    # Z_tot_real = np.real(Z_total)
    # Z_tot_coil = np.real(Z_coil)
    # Z_tot_E = np.real(Z_E)


    #Constants for Approx. R_ESR
    a = C_t/(C_s+C_w)
    alpha = a/(a+1)
    #APPROXIMATED Total Resistance [Siverns Eq 24]
    R_ESRa = R_j + R_c +R_s +R_t*(alpha**2)


    #Q FACTOR!!!  [Siverns Eq 22]
    Q = X_Lc/R_ESR
    #print(Q, 'Q Factor [Siverns Eq 22]')
    
    return Q 

#C_t = 30e-12
#R_t = 3 
#omega_0 = 14164762.44
#t = np.linspace(0.006, 0.030, 200)

#computed_qs = max_qfactor3(C_t, R_t, t)
#print(computed_qs)

#plt.plot(t, computed_qs)

# fail = d > 0.070
# filtered_fail = d[fail]
# #print(filtered_fail)
# first_element = filtered_fail[0]
# # plt.fill_between(d, computed_qs, where=( d >= first_element), color='lightblue')

#plt.axvline(x = 0.088, color='red', linestyle='--')




#plt.fill_between(d, computed_qs, where=( d >= first_element), color='lightblue')


# Adding labels and title
#plt.xlabel('tau (meters)')
#plt.ylabel('Q')


# Displaying the plot
#plt.show()