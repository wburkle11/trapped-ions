# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:36:52 2025

@author: iontrap
"""

#EIT COOLING FUNCTIONS

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.ticker import ScalarFormatter
from sympy import symbols, Eq, solve
from scipy.ndimage import gaussian_filter
from matplotlib.lines import Line2D


#DEFINITIONS

delpi = np.linspace(-30*1e6, 80*1e6, 50000)* 2*math.pi  #PROBE DETUNING

zeeman = 2*math.pi*5*1e6  #Zeeman Splitting 
del_d = 2*math.pi*60*1e6  #DRIVE DETUNING
delpi3 = np.linspace(del_d / (2*math.pi) + 2*1e6, del_d / (2*math.pi) + 20*1e6, 50000)* 2*math.pi  #PROBE DETUNING
#print(delpi3 / (2*math.pi))
delpi2 = del_d - zeeman
 
gamma = 2*math.pi * 19.6*1e6          #total decay rate out of the excited state
delsig_p = del_d + zeeman           #total detuning of sigma_plus component of driving beam
delsig_m = del_d - zeeman           #total detuning of sigma_minus component of driving beam
omega_sig_m = 2*math.pi * 19.5*1e6      #Rabi of sigma_minus component of driving beam
omega_sig_p = omega_sig_m      #Rabi of sigma_plus component of driving beam
omega_pi = 2*math.pi * 5*1e6          #Rabi of probe beam
nu = 2*math.pi*1.3*1e6                  #Trap Freq. we wish to cool
nu2 = np.linspace(0.5*1e6, 5*1e6, 909)* 2*math.pi  #PROBE DETUNING

omega = omega_sig_m                 #Parameter for computing pee in Monroe paper (..questionable..)


#CALCULATIONS

def calculate_diff(x):
    
    ''' This function computes a simple difference between the detunings of the probe and driving beams respectively, 
    What goes into this function as the detuning of the probe beam, can be the probe beam detuning itself, or the 
    probe beam detuning +/- the trap freq we wish to cool'''
    
    diff = (x - delsig_m)
    
    return diff

x1 = delpi 
x2 = delpi + nu
x3 = delpi - nu
x4 = delpi2
x5 = delpi2 + nu2
x6 = delpi2 - nu2
x7 = delpi3
x8 = delpi3 + nu
x9 = delpi3 - nu

diff_car = calculate_diff(x1)
diff_red = calculate_diff(x2)
diff_blue = calculate_diff(x3)

diff_omega2_mon_car = calculate_diff(x4)
diff_omega2_mon_red = calculate_diff(x5)
diff_omega2_mon_blue = calculate_diff(x6)

diff_car2 = calculate_diff(x7)
diff_red2 = calculate_diff(x8)
diff_blue2 = calculate_diff(x9)

def compute_peekim(x):
    
    ''' 
    
    This function computes the resulting "square of scatterring amplitude", with expressions provided
    in the Kim paper. It is a function of "x" and other defined variables. "x" represents the
    argument for probe detuning, as this argument will take on 3 different values 
    
    -delpi
    -delpi + nu
    -delpi - nu
    
    '''

    Z = 4*gamma**2 * (x - delsig_m)**2*(x - delsig_p)**2 + (4*x*(x - delsig_p)*(x - delsig_m) - (x - delsig_p)* omega_sig_m**2 - (x - delsig_m)* omega_sig_p**2)**2

    W = 16*(x - delsig_m)**2 * (x - delsig_p)**2/Z
    
    return W/gamma

kim_car_result = compute_peekim(x1)
kim_red_result = compute_peekim(x2)
kim_blue_result = compute_peekim(x3)

kim_omega2_car = compute_peekim(x4)
kim_omega2_red = compute_peekim(x5)
kim_omega2_blue = compute_peekim(x6)

kim_car_result2 = compute_peekim(x7)
kim_red_result2 = compute_peekim(x8)
kim_blue_result2 = compute_peekim(x9)


def compute_peemonroe(x, y):
    
    ''' 
    
    Computes resulting excited state population based on expressions in the Monroe paper. The "x" argument 
    is the input for probe beam detuning, which is itself, and itself +/- trap frequency. The "y" argument represents the difference 
    between probe and drive beam detunings, corrresponding to the choice of argument "x"
    
    '''

    Z = 8*y**2*omega_sig_m**2*omega_pi**2*gamma + 2*y**2*gamma**3*omega**2 - 4*x*y*omega_sig_m**4*gamma + 1/2*omega**6*gamma + 8*y**2*gamma*(delsig_m**2*omega_pi**2 + x**2*omega_sig_m**2) + 4*delsig_m*y*omega_pi**2*gamma

    pee = (4*y**2*omega_pi**2*omega_sig_m**2*gamma)/Z
    
    return pee

mon_car_result = compute_peemonroe(x1, diff_car)
mon_red_result = compute_peemonroe(x2, diff_red)
mon_blue_result = compute_peemonroe(x3, diff_blue)

mon_omega2_car = compute_peemonroe(x4, diff_omega2_mon_car)
mon_omega2_red = compute_peemonroe(x5, diff_omega2_mon_red)
mon_omega2_blue = compute_peemonroe(x6, diff_omega2_mon_blue)



def compute_nbar_ratio(x, y, z):
    
    '''
    
    This function will compute nbar using the ratio method with carrier excited state population, red excited state population, 
    and blue excited state population as inputs. "x" = Carrier, "y" = RSB, "z" = BSB
    
    '''
    
    nbar = (x + z) / (y - z)
    
    return nbar

nbar_kim = compute_nbar_ratio(kim_car_result, kim_red_result, kim_blue_result)
nbar_mon = compute_nbar_ratio(mon_car_result, mon_red_result, mon_blue_result)

nbar_kim_2 = compute_nbar_ratio(kim_omega2_car, kim_omega2_red, kim_omega2_blue)
nbar_mon_2 = compute_nbar_ratio(mon_omega2_car, mon_omega2_red, mon_omega2_blue)

nbar_kim_3 = compute_nbar_ratio(kim_car_result2, kim_red_result2, kim_blue_result2)

# Apply Gaussian smoothing
nbar_kim_smooth = gaussian_filter(nbar_kim, sigma=8)  # Adjust sigma for smoothing strength
nbar_mon_smooth = gaussian_filter(nbar_mon, sigma=8)
nbar_kim_2_smooth = gaussian_filter(nbar_kim_2, sigma=8)
nbar_mon_2_smooth = gaussian_filter(nbar_mon_2, sigma=8)


#AC LIGHT SHIFT (KIM)
lambda_ = symbols('lambda')

    
# Define the equation
equation = Eq(1/4 * lambda_ * (4*lambda_*(lambda_ - delsig_p) * (lambda_ - delsig_m) - 
                  (lambda_ - delsig_p) * omega_sig_m**2 - 
                  (lambda_ - delsig_m) * omega_sig_p**2), 0)

 # Solve the equation for lambda
solutions = solve(equation, lambda_) 

    # Extract the second-to-last solution
second_to_last_solution = solutions[-2]

    # Print the absolute value of the second-to-last solution
adjusted_sol = abs(second_to_last_solution - delsig_m)  # Subtract del_d and take the absolute value
print('AC Light Shift =',adjusted_sol / (2*math.pi*1e6), 'MHz')  # Print the result in MHz (divide by 2*pi)


#Define a cut off threshhold and masks

delpi_lower_threshold2 = delsig_m - 0.5*1e6*(2*math.pi) 
delpi_upper_threshold2 = delsig_m + 1.5*1e6*(2*math.pi) 
mask3 = (delpi >= delpi_lower_threshold2) & (delpi <= delpi_upper_threshold2)
mask = (delpi >= delpi_lower_threshold2) & (delpi <= delpi_upper_threshold2)
filtered_delpi = delpi[mask]
filty2_kim = nbar_kim_smooth[mask]
delpi_lower_threshold = delsig_m - 2*math.pi*5*1e6
delpi_upper_threshold = delsig_m + 2*math.pi*10*1e6
mask2 = (delpi >= delpi_lower_threshold) & (delpi <= delpi_upper_threshold)
zoomed_delpi = delpi[mask2]
zoomed_kim_result = kim_car_result[mask2]
zoomed_mon_result = mon_car_result[mask2]
delpi_lower_threshold3 = del_d + 2*1e6*(2*math.pi) 
delpi_upper_threshold3 = del_d + 7*1e6*(2*math.pi)
extra_mask = (delpi >= delpi_lower_threshold3) & (delpi <= delpi_upper_threshold3)
zeeman_delpi = delpi[extra_mask]


#zeeman_delpi2 = zeeman_delpi[:len(nbar_kim_2)]  # Trim x to match the length of y
#print(len(zeeman_delpi))
#print(zeeman_delpi / (2*math.pi))
#nbar_kim_2_interpolated = np.interp(zeeman_delpi, delpi, nbar_kim_2)
#print(len(nbar_kim_2))
#print(nbar_kim_2)


def find_min(x):
    
    '''
    
    This function is used for finding the minimum value of nbar given a scan of probe detunings
    
    '''
    
    min_index = np.argmin(x)
    x_min = delpi[min_index]
    
    return x_min


min_index_kim = find_min(filty2_kim)
min_index_kim_MHZ = min_index_kim / 1e6
minval_nbarkim = np.min(filty2_kim)


target_nbar_kim = minval_nbarkim



#  Find the indices where y == target_y
optimal_detunings = []  # To store the optimal detuning values

# Define a function to find the closest match to a target value
def find_closest(arr, target):
    # Calculate the absolute differences between the target and all elements in the array
    differences = np.abs(arr - target)
    # Find the index of the minimum difference
    closest_index = np.argmin(differences)
    # Return the closest value and its index
    return closest_index, arr[closest_index]

# Get the closest matches for Kim
min_index_kim, closest_nbar_kim = find_closest(filty2_kim, target_nbar_kim)
min_index_kim_MHZ = min_index_kim / 1e6  # Convert to MHz if needed
corresponding_x_value_kim = delpi[min_index_kim] / (2 * np.pi)


def find_opt_detuning(x, y, z, delpi, tolerance=1e-20):
    '''
    This function finds the optimal detuning corresponding to the minimum nbar value.
    
    x: matching indices (from np.where)
    y: filtered nbar (nbar_rescaled_kim or nbar_rescaled_mon)
    z: target nbar (minval_nbarkim or minval_nbarmon)
    delpi: probe detuning (x-values)
    tolerance: tolerance for floating-point comparison (default is 1e-6)
    '''
    
    # Find indices where y is close to target nbar (within the given tolerance)
    matching_indices = np.where(np.isclose(y, z, atol=tolerance))[0]
    
    if matching_indices.size > 0:
        for index in matching_indices:
            optimal_detuning = delpi[index] / (2*math.pi)
            optimal_detunings.append(optimal_detuning)  # Store the value
            #print(f"Exact match: Optimal Detuning = {optimal_detuning} Hz, Minimum Nbar = {y[index]}")
        
    return optimal_detunings


# Step 1: Find the x-value corresponding to the minimum value of filty2_kim
min_nbar_index = np.argmin(filty2_kim)  # Minimum value of the nbar data
min_nbar_value = filty2_kim[min_nbar_index]  # Index of the minimum value
min_nbar_x_value = filtered_delpi[min_nbar_index] / (2 * math.pi)  # Corresponding x-value
overlap_index = np.argmin(np.abs(filtered_delpi - delsig_m))
closest_x = filty2_kim[overlap_index]

# PLOTS 

#FULL-FANO-PROFILE

# Create a 1x2 grid of subplots (1 row, 2 columns)
#fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))  # figsize adjusts the overall size
fig = plt.figure(figsize=(12, 5))

# Add a single axis
ax1 = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
label1 = r'$\Delta_{d}$' + f' = {del_d / (2*math.pi*1e6):.2f} MHz'
label2 = r'$\delta_{B}$' + f' = {zeeman / (2*math.pi*1e6):.2f} MHz'
label3 = r'$\Omega_{+} \equiv \Omega_{-}$' + f' = {omega_sig_p / (2*math.pi*1e6):.2f} MHz'
label4 = r'$\Omega_{\pi}$' + f' = {omega_pi / (2*math.pi*1e6):.2f} MHz'
combo_label1 = f'{label1} \n {label2} \n {label3} \n {label4}'
# Plot on the first subplot (ax1)
ax1.plot(delpi / (2*math.pi), kim_car_result, label=combo_label1)
ax1.set_title('FANO-Profile (Double-EIT)')
ax1.set_xlabel(r'$\Delta_{\pi} / 2\pi$ (Hz)', fontsize = 12)
ax1.set_ylabel('Fluorescence')
liner1 = ax1.axvline(x=delsig_m / (2*math.pi), color='black', linestyle='--', label=r'$\Delta_{\pi} \equiv \Delta = \Delta_{d} - \delta_{B}$' + f' = {delsig_m / (2*math.pi*1e6):.2f} MHz')
liner2 = ax1.axvline(x=delsig_p / (2*math.pi), color='black', linestyle='--', label=r'$\Delta_{\pi} \equiv \Delta = \Delta_{d} + \delta_{B}$' + f' = {delsig_p / (2*math.pi*1e6):.2f} MHz')
custom_line = Line2D([0], [0], color='none')  # Line with no color (invisible line)
ax1.legend([custom_line, liner1, liner2], [combo_label1, r'$\Delta_{\pi} \equiv \Delta_{-} = \Delta_{d} - \delta_{B}$' + f' = {delsig_m / (2*math.pi*1e6):.2f} MHz', r'$\Delta_{\pi} \equiv \Delta_{+} = \Delta_{d} + \delta_{B}$' + f' = {delsig_p / (2*math.pi*1e6):.2f} MHz' ], loc='best', fontsize=12)

# Plot on the second subplot (ax2)
# ax2.plot(delpi / (2*math.pi), mon_car_result)
# ax2.set_title('FANO-Profile (Monroe, Single-EIT)')
# ax2.set_xlabel('delpi / 2pi')
# ax2.set_ylabel('Excited State Population')

# Adjust spacing between the plots
plt.tight_layout()



#ZOOMED IN FANO PROFILE 

# Create a 1x2 grid of subplots (1 row, 2 columns)
#fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))  # figsize adjusts the overall size
fig2 = plt.figure(figsize=(12, 5))

# Add a single axis
ax3 = fig2.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot

# Plot on the first subplot (ax1)
ax3.plot(zoomed_delpi / (2*math.pi), zoomed_kim_result)
ax3.set_title(r'FANO-Profile $(\Delta_{\pi} \equiv \Delta_{-} = \Delta_{d} - \delta_{B})$')
ax3.set_xlabel(r'$\Delta_{\pi} / 2\pi$ (Hz)', fontsize = 12)
ax3.set_ylabel('Fluorescence')
ax3.axvline(x=delsig_m / (2*math.pi), color='black', linestyle='--', label=r'$\Delta_{-}$' + f' = {delsig_m / (2*math.pi*1e6):.2f} MHz')
ax3.axvline(x= (delsig_m + adjusted_sol) / (2*math.pi), color='green', label=r'$\delta_{+}$' + f' = {adjusted_sol / (2*math.pi*1e6):.2f} MHz')
ax3.axvline(x=delsig_m / (2*math.pi) + nu / (2*math.pi), color='red', label=r'$\Delta_{-} + \nu$' + f' = {delsig_m / (2*math.pi*1e6) + nu / (2*math.pi*1e6):.2f} MHz')
ax3.axvline(x=delsig_m / (2*math.pi) - nu / (2*math.pi), color='blue', label=r'$\Delta_{-} - \nu$' + f' = {delsig_m / (2*math.pi*1e6) - nu / (2*math.pi*1e6):.2f} MHz')
ax3.axvline(x=delsig_m / (2*math.pi) + nu / (2*math.pi), color='red', label=r'$\nu$' + f' = {nu / (2*math.pi*1e6):.2f} MHz', alpha=0)
ax3.legend(loc='best', fontsize=12)

# Plot on the second subplot (ax2)
# ax4.plot(zoomed_delpi / (2*math.pi), zoomed_mon_result)
# ax4.set_title('FANO-Profile (Monroe, Single-EIT)')
# ax4.set_xlabel('delpi / 2pi')
# ax4.set_ylabel('Excited State Population')
# ax4.axvline(x=delsig_m / (2*math.pi), color='black', label=f'delsig_m = {delsig_m / (2*math.pi*1e6):.2f} MHz')
# ax4.axvline(x=delsig_m / (2*math.pi) + nu / (2*math.pi), color='red', label=f'delsig_m + nu = {delsig_m / (2*math.pi*1e6) + nu / (2*math.pi*1e6):.2f} MHz')
# ax4.axvline(x=delsig_m / (2*math.pi) - nu / (2*math.pi), color='blue', label=f'delsig_m - nu = {delsig_m / (2*math.pi*1e6) - nu / (2*math.pi*1e6):.2f} MHz')
# ax4.legend(loc='best', fontsize=10)
# Adjust spacing between the plots
plt.tight_layout()


#NBAR

#fig2, (ax5, ax6) = plt.subplots(1, 2, figsize = (12, 5))
fig3 = plt.figure(figsize=(12, 5))

# Add a single axis
ax5 = fig3.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot

# Plot on the first subplot (ax1)
combined_label1 = (r'Minimum $\bar{n}$' + f' = {closest_nbar_kim:.3f} \n' 
                  r'Optimal $\Delta_{\pi}$' + f' = {min_nbar_x_value / 1e6:.2f} MHz \n' 
                  r'$\bar{n}\;(\Delta_{-})$' + f' = {closest_x:.3f}')

ax5.plot(filtered_delpi / (2*math.pi), filty2_kim, label=combined_label1)
ax5.set_title('Minimum Nbar')
ax5.set_yscale('log')
ax5.set_xlabel(r'$\Delta_{\pi} / 2\pi$ (Hz)', fontsize = 12)
ax5.set_ylabel(r'$\bar{n}$', fontsize = 12)
line2 = ax5.axvline(x=delsig_m / (2*math.pi), color='black', linestyle='--', label=r'$\Delta_{-}$' + f' = {delsig_m / (2*math.pi*1e6):.2f} MHz')

# Combine all legends into one
ax5.legend([custom_line, line2], [combined_label1, r'$\Delta_{-}$' + f' = {delsig_m / (2*math.pi*1e6):.2f} MHz'], loc='best', fontsize=12)

# Plot on the second subplot (ax2)
# combined_label2 = f'min nbar = {minval_nbarmon:.2f} , optimal delpi = {optimal_detuning_mon[1]/1e6:.2f} MHz'
# ax6.plot(filtered_delpi / (2*math.pi), filtered_nbar_mon, label=combined_label2)
# ax6.set_title('nbar (Monroe)')
# ax6.set_yscale('log')
# ax6.set_xlabel('delpi / 2pi')
# ax6.set_ylabel('nbar')
# line4 = ax6.axvline(x=delsig_m / (2*math.pi), color='black', label=f'delsig_m = {delsig_m / (2*math.pi*1e6):.2f} MHz')

# # Combine all legends into one
# ax6.legend(loc='best', fontsize=10)

# Adjust spacing between the plots
plt.tight_layout()


#COOLING BANDWIDTH

# Create a 1x2 grid of subplots (1 row, 2 columns)
#fig, (ax7, ax8) = plt.subplots(1, 2, figsize=(12, 5))  # figsize adjusts the overall size
fig4 = plt.figure(figsize=(12, 5))

# Add a single axis
ax7 = fig4.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
# Plot on the first subplot (ax1)
ax7.plot(nu2 / (2*math.pi*1e6), nbar_kim_2)
#print(nbar_kim_2)
ax7.set_title('Nbar vs Trap Frequency')
ax7.set_xlabel(r'$\nu / 2\pi$ (MHz)', fontsize = 12)
ax7.set_ylabel(r'$\bar{n}$', fontsize = 12)
ax7.set_yscale('log')
ax7.axvline(x=adjusted_sol / (2*math.pi*1e6), color='green', label=r'$\delta_{+}$' + f' = {adjusted_sol / (2*math.pi*1e6):.2f} MHz')
ax7.legend(loc='best', fontsize=12)
# Plot on the second subplot (ax2)
# ax8.plot(nu2 / (2*math.pi), nbar_rescaled2_mon)
# ax8.set_title('Nbar vs Omega (Monroe)')
# ax8.set_xlabel('nu / 2pi')
# ax8.set_ylabel('nbar')

# Adjust spacing between the plots
plt.tight_layout()


#DETUNING DIFFERENCE (Supplemental Plot)

# fig6 = plt.figure(figsize=(12, 5))
# ax6 = fig6.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot

# ult_diff = (delpi3 - del_d) / (2*math.pi)
# #print(len(ult_diff))

# ax6.plot(ult_diff, nbar_kim_3)

# ax6.set_title('Relative Detuning vs. Nbar')
# ax6.set_xlabel('delpi - deld')
# ax6.set_ylabel('nbar')
# ax6.set_yscale('log')
# line8 = ax6.axvline(x= zeeman / (2*math.pi), color='red', linestyle='--', label=f'Zeeman Shift = {zeeman / (2*math.pi*1e6):.2f} MHz')
# ax6.legend(loc='best', fontsize=10)

# # Adjust spacing between the plots
# plt.tight_layout()

# # Show the plots
# plt.show()