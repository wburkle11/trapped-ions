

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from Helical_Resonator_Functions import *


#FUNCTIONS THAT CREATE PLOTS AND STORE VALUES FROM PREVIOUS FUNCTIONS

def init_res_freq(b, D, d, d_0):
    
    ''' Returns a plot of Resonant Freq vs. C_t, and Resonant Freq value '''
    
    C_t = 1e-12 * np.linspace(0.0001, 100, 100)
    omegas = compute_res_freq(C_t, b, D, d, d_0)
    value = compute_res_freq(C_t1, b, D, d, d_0)
    
    
    fig1, ax1 = plt.subplots(figsize = (10, 8))
    ax1.plot(C_t, omegas, '-r')
    plt.ylim(0, 80000000)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    plt.xlabel('Trap Capacitance [F]')
    plt.ylabel('Resonant Frequency (Hz)')
    plt.grid()
    plt.show()
    
    return fig1, ax1, value

def init_qfactor(b, D, d, d_0):
    
    '''returns computed resonator qfactor values (depends on R_t), and plot of qfactor vs. C_t '''
    
    C_t = 1e-12 * np.linspace(0.0001, 100, 100) 
    qs = compute_q_factor(C_t, R_t, b, D, d, d_0)
    qss = compute_q_factor(C_t1, R_t, b, D, d, d_0)
    
    
    fig2, ax2 = plt.subplots(figsize = (10, 8))
    ax2.plot(C_t, qs, '-r')
    plt.legend(["R_t = 3 Ohm", "R_t = 5 Ohm", "R_t = 7 Ohm", "R_t = 10 Ohm"], prop = { "size": 15 }, loc="upper right")
    plt.ylim(0, 1000)
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    plt.xlabel('Trap Capacitance [F]')
    plt.ylabel('Q')
    plt.grid()
    plt.show()
    
    return fig2, ax2, qss


def optimize_D(b, d, d_0): 
    
    ''' returns plot of optimal D s.t. qfactor is maximized, as well as associated value '''
    
    D = np.linspace(0.075, 0.09997, 200)
    Dqs = max_qfactor2(C_t1, R_t, D, b, d, d_0)
    
    fig3 = plt.plot(D, Dqs)
    index_max = np.argmax(Dqs)
    best_D = D[index_max]
    
    plt.axvline(x = best_D, color='green', linestyle='--')
    plt.xlabel('D (meters)')
    plt.ylabel('Q')
    plt.show()
    
    return fig3, best_D

def optimize_d(D, b, d_0):
    
    ''' returns plot of optimal d s.t. qfactor is maximized, as well as associated value '''
    
    d = np.linspace(0.0001, 0.07, 200)
    dqs = max_qfactor2(C_t1, R_t, D, b, d, d_0)

    fig4 = plt.plot(d, dqs)
    index_qmax = np.argmax(dqs)
    best_d = d[index_qmax]
    
    plt.axvline(x = best_d, color='green', linestyle='--')
    plt.xlabel('d (meters)')
    plt.ylabel('Q')
    plt.show()
    
    return fig4, best_d

def optimize_b(D, d, d_0):
    
    ''' returns plot of optimal b s.t. qfactor is maximized, as well as associated value '''
    
    b = np.linspace(0.04, 0.2, 200)
    bqs = max_qfactor2(C_t1, R_t, D, b, d, d_0)
    
    fig5 = plt.plot(b, bqs)
    index_max = np.argmax(bqs)
    best_b = b[index_max]
    
    plt.axvline(x = best_b, color='green', linestyle='--')
    plt.xlabel('b (meters)')
    plt.ylabel('Q')
    plt.show()
    
    return fig5, best_b


def optimize_d_0(D, b, d):
    
    ''' returns plot of optimal d_0 s.t. qfactor is maximized, as well as associated value '''
    
    d_0 = np.linspace(0.002, 0.017, 200)
    d_0qs = max_qfactor2(C_t1, R_t, D, b, d, d_0)
    
    fig6 = plt.plot(d_0, d_0qs)
    
    index_max = np.argmax(d_0qs)
    best_d_0 = d_0[index_max]
    
    
    plt.axvline(x = best_d_0, color='green', linestyle='--')
    plt.xlabel('d_0 (meters)')
    plt.ylabel('Q')
    plt.show()
    
    return fig6, best_d_0

def optimize_tau(D, b, d, d_0):
    
    ''' returns plot of optimal tau s.t. qfactor is maximized, as well as associated value '''
    
    t = np.linspace(0.004, 0.030, 200)
    tqs = max_qfactor3(C_t1, R_t, t, D, b, d, d_0)
    
    fig7 = plt.plot(t, tqs)
    
    index_max = np.argmax(tqs)
    best_t = t[index_max]
    
    
    plt.axvline(x = best_t, color='green', linestyle='--')
    plt.xlabel('tau (meters)')
    plt.ylabel('Q')
    plt.show()
    
    return fig7, best_t


# Helical Resonator Dimensions
# Sheild Height:
b = 0.178  
# Sheild Diameter
D = 0.0997  
# Main Coil Diameter
d = 0.055  
# Diameter of Main Coil Wire
d_0 = 0.003


#Trap Specific Parameters
R_t = 3 
C_t1 = 0.3e-10


fig1, ax1, value = init_res_freq(b, D, d, d_0)
fig2, ax2, qss = init_qfactor(b, D, d, d_0)
fig3, best_D = optimize_D(b, d, d_0)
fig4, best_d = optimize_d(D, b, d_0)
fig4, best_b = optimize_b(D, d, d_0)
fig6, best_d_0 = optimize_d_0(D, b, d)
fig6, best_t = optimize_tau(D, b, d, d_0)

print('─' * 60)
print('For 30pF -> Q Factor =', qss)
print('For 30pF -> Resonant Frequency =', value*10**-6, 'MHz')
print('--' * 30)
print('Optimal D =', best_D)
print('Optimal d =', best_d)
print('Optimal b =', best_b)
print('Optimal d_0 =', best_d_0)
print('Optimal tau =', best_t)
