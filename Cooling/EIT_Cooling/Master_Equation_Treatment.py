# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 14:13:23 2025

@author: iontrap
"""
import numpy as np
import qutip as qt
import math
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve
from matplotlib.lines import Line2D

# ----------------------------------------
# Section 1: Definitions
# ----------------------------------------

# Define the parameters
Gamma = 2 * math.pi * 19.6 * 1e6  # total decay rate out of the excited state
Omega_s_minus = 2 * math.pi * 22 * 1e6  # Rabi frequency for sigma- transition
Omega_r = 2 * math.pi * 2 * 1e6  # Rabi frequency for r transition
Omega_s_plus = 2 * math.pi * 22 * 1e6  # Rabi frequency for sigma+ transition
Delta_d = 2 * math.pi * 55 * 1e6  # DRIVE DETUNING
delta_B = 2 * math.pi * 5 * 1e6  # Zeeman Splitting
delsig_m = Delta_d - delta_B
delsig_p = Delta_d + delta_B
nu = 2*math.pi*2*1e6             #Trap Freq. we wish to cool
nu2 = np.linspace(0.5*1e6, 6*1e6, 909)* 2*math.pi  
delpi2 = Delta_d - delta_B

# Define the range of Delta_p values
Delta_p_values = np.linspace(-30 * 1e6, 80 * 1e6, 1000) * 2 * math.pi  # Range of Delta_p


# ----------------------------------------
# Section 2: Calculations
# ----------------------------------------

def calculate_diff(x):
    
    ''' This function computes a simple difference between the detunings of the probe and driving beams respectively, 
    What goes into this function as the detuning of the probe beam, can be the probe beam detuning itself, or the 
    probe beam detuning +/- the trap freq we wish to cool'''
    
    diff = (x - delsig_m)
    
    return diff

x1 = Delta_p_values 
x2 = Delta_p_values + nu
x3 = Delta_p_values - nu
x4 = delpi2 
x5 = delpi2 + nu2
x6 = delpi2 - nu2


diff_car = calculate_diff(x1)
diff_red = calculate_diff(x2)
diff_blue = calculate_diff(x3)

diff_car2 = calculate_diff(x4)
diff_red2 = calculate_diff(x5)
diff_blue2 = calculate_diff(x6)

#------------------------------------------------------------------------

def ss_master_equation_sol(delpi, del_d, omega_sig_m, omega_sig_p, omega_pi, zeeman):
    ''' This function calculates the steady state solution to the master equation'''

    # Define the basis states
    N = 4
    psi_0 = qt.basis(N, 2)  # |0>
    psi_plus = qt.basis(N, 1)  # |+>
    psi_minus = qt.basis(N, 3)  # |->  
    psi_e = qt.basis(N, 0)  # |e>

    # Define the collapse operators
    c1 = np.sqrt(Gamma / 3) * psi_plus * psi_e.dag()  # Decay from |e> to |+>
    c2 = np.sqrt(Gamma / 3) * psi_0 * psi_e.dag()  # Decay from |e> to |0>
    c3 = np.sqrt(Gamma / 3) * psi_minus * psi_e.dag()  # Decay from |e> to |->  
    c_ops = [c1, c2, c3]

    # Ensure delpi is iterable
    if np.isscalar(delpi):  
        delpi = [delpi]  # Convert single value into a list

    # Initialize an array to store the population in |e|
    population_e = np.zeros(len(delpi), dtype=complex)

    # Loop over each value of Delta_p
    for i, delta_p in enumerate(delpi):
        # Define the Hamiltonian for the current Delta_p
        H_s = qt.Qobj([[0, omega_sig_m / 2, -omega_pi / 2, omega_sig_p / 2],
                       [omega_sig_m / 2, del_d + zeeman, 0, 0],
                       [-omega_pi / 2, 0, delta_p, 0],
                       [omega_sig_p / 2, 0, 0, del_d - zeeman]])

        # Solve for the steady state
        rho_ss = qt.steadystate(H_s, c_ops, method='direct', tol=1e-10)

        # Calculate the population in the excited state |e>
        population_e[i] = (psi_e.dag() * rho_ss * psi_e).tr()

    # Return a single value if input was a scalar, otherwise return an array
    return population_e[0] if len(population_e) == 1 else population_e

population_e1 = ss_master_equation_sol(Delta_p_values, Delta_d, Omega_s_minus, Omega_s_plus, Omega_r, delta_B)

car_result = ss_master_equation_sol(x1, Delta_d, Omega_s_minus, Omega_s_plus, Omega_r, delta_B)
red_result = ss_master_equation_sol(x2, Delta_d, Omega_s_minus, Omega_s_plus, Omega_r, delta_B)
blue_result = ss_master_equation_sol(x3, Delta_d, Omega_s_minus, Omega_s_plus, Omega_r, delta_B)

car_result2 = ss_master_equation_sol(x4, Delta_d, Omega_s_minus, Omega_s_plus, Omega_r, delta_B)
red_result2 = ss_master_equation_sol(x5, Delta_d, Omega_s_minus, Omega_s_plus, Omega_r, delta_B)
blue_result2 = ss_master_equation_sol(x6, Delta_d, Omega_s_minus, Omega_s_plus, Omega_r, delta_B)

#------------------------------------------------------------------------

def compute_nbar_ratio(x, y, z):
    
    '''
    
    This function will compute nbar using the ratio method with carrier excited state population, red excited state population, 
    and blue excited state population as inputs. "x" = Carrier, "y" = RSB, "z" = BSB
    
    '''
    
    nbar = (x + z) / (y - z)
    
    return nbar


nbar_result = compute_nbar_ratio(car_result, red_result, blue_result)
nbar_result2 = compute_nbar_ratio(car_result2, red_result2, blue_result2)

#------------------------------------------------------------------------

#AC LIGHT SHIFT (KIM)
lambda_ = symbols('lambda')

    
# Define the equation
equation = Eq(1/4 * lambda_ * (4*lambda_*(lambda_ - delsig_p) * (lambda_ - delsig_m) - 
                  (lambda_ - delsig_p) * Omega_s_minus**2 - 
                  (lambda_ - delsig_m) * Omega_s_plus**2), 0)

 # Solve the equation for lambda
solutions = solve(equation, lambda_) 

    # Extract the second-to-last solution
second_to_last_solution = solutions[-2]

    # Print the absolute value of the second-to-last solution
adjusted_sol = abs(second_to_last_solution - delsig_m)  # Subtract del_d and take the absolute value
print('AC Light Shift =',adjusted_sol / (2*math.pi*1e6), 'MHz')  # Print the result in MHz (divide by 2*pi)

#------------------------------------------------------------------------

#Creating Mask for 2nd Plot 
delpi_lower_threshold = delsig_m - 2*math.pi*5*1e6
delpi_upper_threshold = delsig_m + 2*math.pi*10*1e6
mask2 = (Delta_p_values >= delpi_lower_threshold) & (Delta_p_values <= delpi_upper_threshold)
zoomed_delpi = Delta_p_values[mask2]
zoomed_result = population_e1[mask2]

#Creating Mask for 3rd Plot
delpi_lower_threshold2 = delsig_m - 0.5*1e6*(2*math.pi) 
delpi_upper_threshold2 = delsig_m + 1.5*1e6*(2*math.pi) 
mask = (Delta_p_values >= delpi_lower_threshold2) & (Delta_p_values <= delpi_upper_threshold2)
filtered_delpi = Delta_p_values[mask]
filty_result = nbar_result[mask]

#------------------------------------------------------------------------

def find_min(x):
    
    '''
    
    This function is used for finding the minimum value of nbar given a scan of probe detunings
    
    '''
    
    min_index = np.argmin(x)
    x_min = Delta_p_values[min_index]
    
    return x_min


min_index_kim = find_min(filty_result)
min_index_kim_MHZ = min_index_kim / 1e6
minval_nbarkim = np.min(filty_result)


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
min_index_kim, closest_nbar_kim = find_closest(filty_result, target_nbar_kim)
min_index_kim_MHZ = min_index_kim / 1e6  # Convert to MHz if needed
corresponding_x_value_kim = Delta_p_values[min_index_kim] / (2 * np.pi)


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
min_nbar_index = np.argmin(filty_result)  # Minimum value of the nbar data
min_nbar_value = filty_result[min_nbar_index]  # Index of the minimum value
min_nbar_x_value = filtered_delpi[min_nbar_index] / (2 * math.pi)  # Corresponding x-value
overlap_index = np.argmin(np.abs(filtered_delpi - delsig_m))
closest_x = filty_result[overlap_index]


# ----------------------------------------
# Section 3: Plots
# ----------------------------------------


# Plot the results
fig1 = plt.figure(figsize=(12, 5))

# Add a single axis
ax1 = fig1.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot

label1 = r'$\Delta_{d}$' + f' = {Delta_d / (2*math.pi*1e6):.2f} MHz'
label2 = r'$\delta_{z}$' + f' = {delta_B / (2*math.pi*1e6):.2f} MHz'
label3 = r'$\Omega_{\sigma_{+}} \equiv \Omega_{\sigma_{-}}$' + f' = {Omega_s_plus / (2*math.pi*1e6):.2f} MHz'
label4 = r'$\Omega_{\pi}$' + f' = {Omega_r / (2*math.pi*1e6):.2f} MHz'
combo_label1 = f'{label1} \n {label2} \n {label3} \n {label4}'
ax1.plot(Delta_p_values.real / (2 * math.pi * 1e6), population_e1.real, 'b-', label=combo_label1)
ax1.set_xlabel('Probe Detuning $\\Delta_p$ (MHz)')
ax1.set_ylabel('Population in $|e\\rangle$')
liner1 = ax1.axvline(x=delsig_m / (2*math.pi*1e6), color='black', linestyle='--', label=r'$\Delta_{\pi} \equiv \Delta_{-} = \Delta_{d} - \delta_{B}$' + f' = {delsig_m / (2*math.pi*1e6):.2f} MHz')
liner2 = ax1.axvline(x=delsig_p / (2*math.pi*1e6), color='black', linestyle='--', label=r'$\Delta_{\pi} \equiv \Delta_{+} = \Delta_{d} + \delta_{B}$' + f' = {delsig_p / (2*math.pi*1e6):.2f} MHz')
custom_line = Line2D([0], [0], color='none')  # Line with no color (invisible line)
ax1.legend([custom_line, liner1, liner2], [combo_label1, r'$\Delta_{\pi} \equiv \Delta_{-} = \Delta_{d} - \delta_{B}$' + f' = {delsig_m / (2*math.pi*1e6):.2f} MHz', r'$\Delta_{\pi} \equiv \Delta_{+} = \Delta_{d} + \delta_{B}$' + f' = {delsig_p / (2*math.pi*1e6):.2f} MHz' ], loc='best', fontsize=12)
ax1.set_title('Population in $|e\\rangle$ vs Probe Detuning $\\Delta_p$')

ax1.legend()

#-----------------------------------------------------------------------------------------

# Plot the results
fig2 = plt.figure(figsize=(12, 5))

# Add a single axis
ax1 = fig2.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot

ax1.plot(zoomed_delpi.real / (2 * math.pi*1e6), zoomed_result.real, 'b-')
ax1.set_xlabel('Probe Detuning $\\Delta_p$ (MHz)')
ax1.set_ylabel('Population in $|e\\rangle$')
ax1.set_title('Population in $|e\\rangle$ vs Probe Detuning $\\Delta_p$')
ax1.axvline(x=delsig_m / (2*math.pi*1e6), color='black', linestyle='--', label=r'$\Delta_{-}$' + f' = {delsig_m / (2*math.pi*1e6):.2f} MHz')
ax1.axvline(x= (delsig_m + adjusted_sol) / (2*math.pi*1e6), color='green', label=r'$\delta_{+}$' + f' = {adjusted_sol / (2*math.pi*1e6):.2f} MHz')
ax1.axvline(x=delsig_m / (2*math.pi*1e6) + nu / (2*math.pi*1e6), color='red', label=r'$\Delta_{-} + \nu$' + f' = {delsig_m / (2*math.pi*1e6) + nu / (2*math.pi*1e6):.2f} MHz')
ax1.axvline(x=delsig_m / (2*math.pi*1e6) - nu / (2*math.pi*1e6), color='blue', label=r'$\Delta_{-} - \nu$' + f' = {delsig_m / (2*math.pi*1e6) - nu / (2*math.pi*1e6):.2f} MHz')
ax1.axvline(x=delsig_m / (2*math.pi*1e6) + nu / (2*math.pi*1e6), color='red', label=r'$\nu$' + f' = {nu / (2*math.pi*1e6):.2f} MHz', alpha=0)
ax1.legend(loc='best', fontsize=12)


#-----------------------------------------------------------------------------------------

fig3 = plt.figure(figsize=(12, 5))

# Add a single axis
ax1 = fig3.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot

# Plot on the first subplot (ax1)
combined_label1 = (r'Minimum $\bar{n}$' + f' = {closest_nbar_kim.real:.3f} \n' 
                  r'Optimal $\Delta_{\pi}$' + f' = {min_nbar_x_value.real / 1e6:.2f} MHz \n' 
                  r'$\bar{n}\;(\Delta_{-})$' + f' = {closest_x.real:.3f}')

ax1.plot(filtered_delpi.real / (2*math.pi*1e6), filty_result.real, 'b-', label=combined_label1)
ax1.set_title('Minimum Nbar')
ax1.set_yscale('log')
ax1.set_xlabel(r'$\Delta_{\pi} / 2\pi$ (MHz)', fontsize = 12)
ax1.set_ylabel(r'$\bar{n}$', fontsize = 12)
line2 = ax1.axvline(x=delsig_m / (2*math.pi*1e6), color='black', linestyle='--', label=r'$\Delta_{-}$' + f' = {delsig_m / (2*math.pi*1e6):.2f} MHz')

# Combine all legends into one
ax1.legend([custom_line, line2], [combined_label1, r'$\Delta_{-}$' + f' = {delsig_m / (2*math.pi*1e6):.2f} MHz'], loc='best', fontsize=12)

#-----------------------------------------------------------------------------------------

#COOLING BANDWIDTH

# Create a 1x2 grid of subplots (1 row, 2 columns)
#fig, (ax7, ax8) = plt.subplots(1, 2, figsize=(12, 5))  # figsize adjusts the overall size
fig5 = plt.figure(figsize=(12, 5))

# Add a single axis
ax7 = fig5.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
# Plot on the first subplot (ax1)
ax7.plot(nu2.real / (2*math.pi*1e6), nbar_result2.real,'-b')
#print(nbar_kim_2)
ax7.set_title('Nbar vs Trap Frequency')
ax7.set_xlabel(r'$\nu / 2\pi$ (MHz)', fontsize = 12)
ax7.set_ylabel(r'$\bar{n}$', fontsize = 12)
ax7.set_yscale('log')
liner = ax7.axvline(x=adjusted_sol / (2*math.pi*1e6), color='green', label=r'$\delta_{+}$' + f' = {adjusted_sol / (2*math.pi*1e6):.2f} MHz')
ax7.axhline(y=0.1, color='black', linestyle='--')

intersection_x_values = []
for i in range(1, len(nu2)):
    # Find where the y-values cross 0.1 by checking when we go from below to above (or above to below)
    if (nbar_result2[i-1] < 0.1 and nbar_result2[i] > 0.1) or (nbar_result2[i-1] > 0.1 and nbar_result2[i] < 0.1):
        # Linear interpolation to find the exact x where y=0.1
        x1, x2 = nu2[i-1], nu2[i]
        y1, y2 = nbar_result2[i-1], nbar_result2[i]
        
        # Linear interpolation for the intersection point
        x_intersection = x1 + (0.1 - y1) * (x2 - x1) / (y2 - y1)
        intersection_x_values.append(x_intersection)
        
    
# Calculate and print the differences between consecutive intersections
if len(intersection_x_values) > 1:
    for i in range(1, len(intersection_x_values)):
        diff = intersection_x_values[i] - intersection_x_values[i - 1]
        

custom_label = (r' Bandwidth $(\bar{n} \leq 0.1)$' f' = {diff.real/(2*math.pi*1e6):.2f} MHz \n'
                r'$\delta_{+}$' + f' = {adjusted_sol / (2*math.pi*1e6):.2f} MHz')
ax7.legend([liner] ,[custom_label], loc='best', fontsize=12)
