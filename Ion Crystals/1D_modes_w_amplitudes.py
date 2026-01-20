# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 17:30:24 2025

@author: iontrap
"""

"""
Transverse (radial) normal modes for a 1D trapped-ion chain.

This version plots each mode as a delta (stem) at its frequency,
with height equal to a chosen "normal-mode amplitude":

- Default: participation of a single ion i, height = |V_{i,k}|^2
- Optional: weighted participation, height = sum_i w_i * |V_{i,k}|^2
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.constants as scipy

# ----------------------------
# Physical constants & inputs
# ----------------------------
q_e   = scipy.elementary_charge
eps0  = scipy.epsilon_0
m_Yb  = 171 * 1.67262192e-27
hbar  = scipy.hbar

# Trap frequencies (in MHz)
fz_MHz = 0.3   # axial
fx_MHz = 3.0   # radial (transverse, e.g., x)

# Convert to angular frequencies (rad/s)
omega_z = 2*np.pi * fz_MHz * 1e6
omega_x = 2*np.pi * fx_MHz * 1e6

# Characteristic length (Coulomb/axial balance)
lc = (q_e**2 / (4*np.pi*eps0*m_Yb*omega_z**2))**(1/3)

# Chain size
N = 8

# ----------------------------
# Equilibrium positions (dimensionless)
# ----------------------------
def equilibrium_positions_dimless(N):
    """
    Solve for dimensionless equilibrium positions u_m (along z/axial),
    defined by: u_m - Σ_{n<m} 1/(u_m - u_n)^2 + Σ_{n>m} 1/(u_m - u_n)^2 = 0
    """
    def equations(u):
        eqs = []
        for m in range(N):
            s1 = sum(1.0/(u[m] - u[n])**2 for n in range(m))
            s2 = sum(1.0/(u[m] - u[n])**2 for n in range(m+1, N))
            eqs.append(u[m] - s1 + s2)
        return eqs

    # symmetric initial guess around 0
    a = 1.5
    init = [(k - (N-1)/2)*a for k in range(N)]
    sol = fsolve(equations, init)
    return np.array(sol, dtype=float)

def equilibrium_positions_meters(u_dimless, lc):
    return u_dimless * lc

# ----------------------------
# Radial (transverse) mode matrix and eigensystem
# ----------------------------
def radial_mode_matrix(u):
    """
    Build the transverse (radial) dynamical matrix B_nm in the standard
    dimensionless form where eigenvalues β_k give ω_k by ω_k = ω_z * sqrt(β_k).
    """
    N = len(u)
    B = np.zeros((N, N), dtype=float)
    for m in range(N):
        s = 0.0
        for p in range(N):
            if p != m:
                s += 1.0 / abs(u[m] - u[p])**3
        B[m, m] = (omega_x/omega_z)**2 - s
        for n in range(N):
            if n != m:
                B[m, n] = 1.0 / abs(u[m] - u[n])**3
    return B

def radial_modes(u):
    """
    Returns:
        freqs_MHz_sorted: eigenfrequencies (MHz), sorted ascending
        eigvecs_sorted: eigenvectors (columns correspond to sorted freqs)
        freqs_MHz_all: unsorted (for indexing if needed)
    """
    B = radial_mode_matrix(u)
    beta, V = np.linalg.eig(B)   # β_k and eigenvectors (columns)
    idx = np.argsort(beta.real)
    beta_sorted = beta.real[idx]
    V_sorted = V[:, idx]
    omega_k = omega_z * np.sqrt(np.clip(beta_sorted, a_min=0.0, a_max=None))
    freqs_Hz  = omega_k / (2*np.pi)
    freqs_MHz = freqs_Hz / 1e6
    return freqs_MHz, V_sorted, ( (omega_z*np.sqrt(np.clip(beta.real,0,None)))/(2*np.pi)/1e6 )

# ----------------------------
# Run: positions, modes, plot (delta functions)
# ----------------------------
u = equilibrium_positions_dimless(N)
z_positions_m = equilibrium_positions_meters(u, lc)

freqs_MHz, eigvecs, _ = radial_modes(u)

# ----------------------------
# Choose how to set delta heights
# ----------------------------
# Option A (default): pick a single ion and show its participation across modes
ion_index = N // 2            # center ion by default; change as desired
heights = np.abs(eigvecs[ion_index, :])**2  # |V_{i,k}|^2 per mode k

# Option B (alternative): weighted participation (e.g. LD/beam projection)
# Provide a weight vector of length N (nonnegative), then uncomment:
# weight_vector = np.ones(N) / N  # example: average over all ions
# weight_vector = np.zeros(N); weight_vector[ion_index] = 1.0  # same as Option A
# heights = (weight_vector[:, None] * (np.abs(eigvecs)**2)).sum(axis=0)

# (Optional) normalize heights to max=1 for nicer plotting
if heights.max() > 0:
    heights = heights / heights.max()

# ----------------------------
# Plot: delta (stem) spectrum
# ----------------------------
plt.figure(figsize=(14, 4))

markerline, stemlines, baseline = plt.stem(freqs_MHz, heights, use_line_collection=True)
plt.setp(markerline, 'marker', 'o', 'markersize', 4)
plt.setp(stemlines,  'linewidth', 1.5)
plt.setp(baseline,   'linewidth', 0.8)

plt.xlabel('Transverse mode frequency (MHz)')
plt.ylabel(r'Mode amplitude (normalized)')
plt.title(f'Radial Modes (N={N}) — delta heights = participation of ion {ion_index}')

# Bandwidth and COM guide
bandwidth = freqs_MHz.max() - freqs_MHz.min()
plt.annotate(f'Bandwidth: {bandwidth:.3f} MHz',
             xy=(0.5*(freqs_MHz.max()+freqs_MHz.min()), 0.9),
             ha='center', va='bottom')
plt.axvline(fx_MHz, linestyle='--', linewidth=0.8, label=f'COM ≈ {fx_MHz:.2f} MHz')

plt.ylim(0, 1.05)
plt.legend()
plt.tight_layout()
plt.show()

print("z positions (µm):", z_positions_m*1e6)
