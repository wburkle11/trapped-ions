# -*- coding: utf-8 -*-
#Notebook for Helical Resonator Calculations: 
#Equations in the file are exact equations taken from Sivers 2012 "Application of RF voltages to ion traps via helical resonators" 

#Importing Libraries
import math
import matplotlib.pyplot as plt
import numpy as np


#Helical Resonator Dimensions
#Sheild Height:
B = 0.120   #meters
#Sheild Diameter
D = 0.108   #meters
#Main Coil Diameter
d = 0.042   #meters
rat3 = d/D
#Winding Pitch
t = 0.009   #meters
#Diameter of Main Coil Wire
d_0 = 0.005
#Uncoiled Main Coil Length
l_c = 0.88   #meters
#Number of turns in COIL
N_c = 6.75
#Coil Height
b = B - D/2    #meters
#print(b)

#Length of rods in Vaccum
l = 0.20   #meters
#Radius of rods
r = 5e-5   #meters
#Resistivity of Copper(Au)
p = 1.7e-8   #ohm-meter
#Trap Capacitance
C_t = 40e-12   #Farads
#print(C_t, 'Trap Capacitance')
#Capacitance of connecting wires
C_w = 0.0001e-12   #Farads

#Capactiance between Main Coil and Sheild [Siverns Eq 26]
rat1 = D/d
K_Cs = 39.37*(0.75/(math.log(D/d)))*1e-12        #farads/meter
C_s = b*K_Cs                                     #farads
#print(C_s)   #farads

#Coil Self Capacitance    [Siverns Eq 25]
H = 11.26*(b/d)+8+(27/math.sqrt(b/d))   #farads/meter
C_C = H*d*1e-12   #farads
#print(C_C)   #farads

#Coil Inductance [Siverns Eq 27]
K_Lc = 39.37*(0.025*(d**2)*(1-((rat3)**2))/(t**2))*1e-6  #henrys/meter
#Inductance
L_C = b*K_Lc
#print(L_C)

#Predicted Resonant Frequency
omega_0 = 1/math.sqrt((C_s + C_t + C_w + C_C)*L_C)
#print(omega_0/(2*(math.pi)), 'Predicted Resonant Frequency Hz [Siverns Eq 21]')

#Predicted Value of b:
C_gamma = C_t + C_w
K_cd = 35*d*1e-12
K_cb = 11.26*1e-12
b1 = (C_gamma + K_cd)/(K_Cs + K_cb)
#print(b1)
b2 = math.sqrt((K_Cs + K_cb)/((C_gamma + K_cd)**2*K_Lc*omega_0**2) +1/4) - 1/2
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
sigma = math.sqrt((2*p)/((omega_0)*4*(math.pi)*1e-7))   #meters
#print(sigma)


#Theory (Siverns Eq 27) Limited Too:
if b/d >= 1:
    print('Geometric Constraint is Valid')
else:
    print('Geometric Constraint is Violeted')


#Trap Resistance [Siverns Eq 33]
A = (math.pi)*r**2 
R = (p*l)/A                          #ohms
R_t = 0.1                               #ohms
#print(R_t, 'Ohms, Trap Resistance')


#Estimate of Coil-Sheild Junction Resistance [Siverns Eq 36]
#Length of Solder Junction
l_j = 0.003                          #meters
#Diameter of Solder Junction
d_j = 0.003                         #meters
R_j = (p*l_j)/((math.pi)*d_j*sigma)  #ohms
#print(R_j)   #ohms


#Sheild Resistance     [Siverns Eq 35]
#Number of turn current undergoes in sheild: 
N_s = (b*l_c)/(4*(math.pi)*(D-d)**2)
#Distance current will travel from bottom to top of sheild
l_s = N_s*math.sqrt(((D**2)*((math.pi)**2))+(b/N_s)**2)
R_s = (N_s*p*l_s)/(b*sigma)   #ohms
#print(R_s)   #ohms


#Coil Resistance  [Sivers Eq 34]
R_c = (p*l_c)/(d_0*(math.pi)*sigma)
#print(R_c)   #ohms


#TOTAL RESISTANCE [Siverns Eq 23]
#tester = (X_Cc*(X_Lc**2 + R_c**2)*R_c*X_Cc)/((X_Lc*X_Cc + X_Lc**2 + R_c**2)**2 + (R_c*X_Cc)**2)
#tester2 = (R_t*X_Cs*X_Cw)/(R_t**2*(X_Cs + X_Cw)**2 + (X_Cs*(X_Ct + X_Cw)) + X_Ct*X_Cw)**2
R_ESRt1 = ((R_c*(X_Cc**2))/(R_c**2 + (X1)**2))
R_ESRt2 = ((R_t*(X_R**2))/(R_t**2 + (X_R + X_T)**2))
#R_TEST = (X_R*(X_Ct**2 + R_t**2)*R_t*X_R)/((X_Ct*X_R + X_Ct**2 + R_t**2)**2 + (R_t*X_R)**2)
R_ESR = R_ESRt1 + R_ESRt2 + R_s + R_j
#print(R_ESRt1, R_ESRt2, R_s, R_j)
#print(R_ESRt1, tester)
#print(R_ESRt2, R_TEST)

#Impedance Route
Z_coil = (1/(1j*X_Lc + R_c) + 1/(1j*X_Cc))**-1
Z_E = (1/(1j*X_Ct + R_t) + 1/(1j*X_Cw) + 1/(1j*X_Cs))**-1
Z_total = Z_coil + Z_E + R_s + R_j
god = np.real(Z_total)
#print(god, R_ESR, '---Real part of Z_total / Total Circuit Resistance')
god2 = np.real(Z_coil)
#print(god2, R_ESRt1, '---Real part of Z_coil / 1st term of Total Circuit Resistance')
god3 = np.real(Z_E)
#print(god3, R_ESRt2, '---Real part of Z_E / 2nd term of Total Circuit Resistance')


#Constants for Approx. R_ESR
a = C_t/(C_s+C_w + C_C)
alpha = a/(a+1)
#print(a, alpha)


#APPROXIMATED Total Resistance [Siverns Eq 24]
R_ESRa = R_j + R_c +R_s +R_t*(alpha**2)


#Q FACTOR!!!  [Siverns Eq 22]
Q = X_Lc/R_ESR
print(Q, 'Q Factor [Siverns Eq 22]')

