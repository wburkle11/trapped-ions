# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 16:05:57 2025

@author: iontrap
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0, pi, elementary_charge
from scipy.linalg import eigh

def axial_modes(omega_z, positions_xyz, m_ion, plane_axes=(0,2)):
    """
    Axial (out-of-plane) modes for a 2D crystal.

    positions_xyz : (N,>=3) array of equilibrium positions [m]
    plane_axes    : which two columns form the crystal plane (e.g., (0,2) for x–z plane)
    Returns:
        f_MHz : (N,) axial mode frequencies in MHz (COM will be ~ omega_z/2π)
        modes : (N,N) eigenvectors (columns)
        lam   : (N,) eigenvalues of Kzz (dimensionless) = (ω_n/ω_z)^2
    """
    q = elementary_charge
    # Correct scaling length: ℓ^3 = q^2/(4π ε0 m ω_z^2)
    ell = (q*q / (4*pi*epsilon_0 * m_ion * omega_z**2))**(1/3)

    # Use the two coordinates that lie in the plane and scale by ℓ
    xy = np.asarray(positions_xyz, dtype=float)[:, plane_axes] / ell
    N = xy.shape[0]

    # Pairwise in-plane distances
    dx = xy[:, 0][:, None] - xy[:, 0][None, :]
    dy = xy[:, 1][:, None] - xy[:, 1][None, :]
    rho2 = dx*dx + dy*dy
    np.fill_diagonal(rho2, np.inf)           # avoid self-term
    C = 1.0 / (rho2 ** 1.5)                   # C_ij = 1/ρ_ij^3
    np.fill_diagonal(C, 0.0)

    # Dimensionless axial matrix Kzz:
    #   K_ii = 1 - Σ_{j≠i} C_ij
    #   K_ij = + C_ij  (i≠j)
    Kzz = C.copy()
    np.fill_diagonal(Kzz, 1.0 - C.sum(axis=1))

    lam, modes = eigh(Kzz)                    # lam = (ω_n/ω_z)^2
    lam = np.clip(lam, 0.0, None)             # guard tiny negatives
    omega = omega_z * np.sqrt(lam)            # ω_n = ω_z * sqrt(lam)
    f_MHz = omega / (2*np.pi*1e6)
    return f_MHz, modes, lam

# -------- Example usage with your file --------
import pandas as pd, math
file_path = r'Z:\Users\Wes\Equilibrium_Position_Data\2DPos3.csv'
positions = pd.read_csv(file_path, header=None).to_numpy()

m_Yb = 171 * 1.67262192e-27
omega_z = 2*math.pi*3.2e6                   # 1.50 MHz COM (rad/s)

# If your crystal plane is x–z, keep plane_axes=(0,2). If x–y, change to (0,1).
f_MHz, modes, lam = axial_modes(omega_z, positions, m_Yb, plane_axes=(0,2))

# ------ Plot in your requested style ------
N = len(f_MHz)
f_sorted = np.sort(f_MHz)
f_com = omega_z / (2*np.pi*1e6)
bandwidth = f_sorted[-1] - f_sorted[0]

plt.figure(figsize=(10, 2.2))
plt.scatter(f_sorted, np.zeros_like(f_sorted), s=8, c='blue')
plt.axvline(f_com, ls='--', c='blue', label=f'COM = {f_com:.2f} MHz')
plt.yticks([])
plt.xlabel('Frequency (MHz)')
plt.title('2D-Crystal Normal Modes')

cx = 0.5*(f_sorted[-1] + f_sorted[0])
plt.text(0.5, 0.72, f'Normal Mode Bandwidth: {bandwidth:.2f} MHz',
         transform=plt.gca().transAxes, ha='center', color='royalblue')
plt.text(0.5, 0.25, f'Number of Ions = {N}',
         transform=plt.gca().transAxes, ha='center', color='royalblue')

xpad = 0.03*(f_sorted[-1] - f_sorted[0] + 1e-12)
plt.xlim(f_sorted[0]-xpad, f_sorted[-1]+xpad)
plt.legend(loc='upper left', framealpha=0.85)
plt.grid(axis='x', alpha=0.2)
plt.tight_layout()
plt.show()

