# -*- coding: utf-8 -*-
"""
Transverse (radial) normal modes for a 1D trapped-ion chain.

- Solves for equilibrium positions along z (dimensionless u_m)
- Builds the radial dynamical matrix B_nm:
    B_nm = [ (ωx/ωz)^2 - Σ_{p≠m} 1/|u_m - u_p|^3 ] δ_nm + (1-δ_nm) / |u_m - u_n|^3
- Eigenvalues β_k give ω_k via:  ω_k = ω_z * sqrt(β_k)
  (So the COM transverse mode has ω_COM = ω_x, as expected.)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.constants as scipy
from matplotlib.lines import Line2D

# ----------------------------
# Physical constants & inputs
# ----------------------------
q_e   = scipy.elementary_charge
eps0  = scipy.epsilon_0
m_Yb  = 171 * 1.67262192e-27
hbar  = scipy.hbar

# Trap frequencies (in MHz)
fz_MHz = 0.6  # axial
fx_MHz = 2.0   # radial (transverse, e.g., x)

# Convert to angular frequencies (rad/s)
omega_z = 2*np.pi * fz_MHz * 1e6
omega_x = 2*np.pi * fx_MHz * 1e6

# Characteristic length (Coulomb/axial balance)
lc = (q_e**2 / (4*np.pi*eps0*m_Yb*omega_z**2))**(1/3)

# Chain size
N = 5

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

    # Symmetric initial guess around 0
    # spacing heuristic ~ (3/2) for N up to ~10; not sensitive
    a = 1.5
    init = [(k - (N-1)/2)*a for k in range(N)]

    sol = fsolve(equations, init)
    return np.array(sol, dtype=float)

def equilibrium_positions_meters(u_dimless, lc):
    """Convert dimensionless u_m to meters."""
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
        # diagonal
        s = 0.0
        for p in range(N):
            if p != m:
                s += 1.0 / abs(u[m] - u[p])**3
        B[m, m] = (omega_x/omega_z)**2 - s

        # off-diagonals
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
    beta, V = np.linalg.eigh(B)       # β_k and eigenvectors (columns)
    # Sort by frequency
    idx = np.argsort(beta.real)
    beta_sorted = beta.real[idx]
    V_sorted = V[:, idx]

    # ω_k = ω_z * sqrt(β_k)
    # Note: For stability, β_k >= 0. If trap is too weak, some β could go negative → imaginary ω.
    omega_k = omega_z * np.sqrt(np.clip(beta_sorted, a_min=0.0, a_max=None))

    # Convert to Hz and MHz
    freqs_Hz  = omega_k / (2*np.pi)
    freqs_MHz = freqs_Hz / 1e6
    return freqs_MHz, V_sorted, ( (omega_z*np.sqrt(np.clip(beta.real,0,None)))/(2*np.pi)/1e6 )

# ----------------------------
# Run: positions, modes, plot
# ----------------------------
u = equilibrium_positions_dimless(N)
z_positions_m = equilibrium_positions_meters(u, lc)

freqs_MHz, eigvecs, _ = radial_modes(u)

# Exact frequencies
xvals = freqs_MHz

plt.figure(figsize=(14, 4))

# --- Vertical lines at each mode frequency ---
for f in xvals:
    plt.axvline(f, ymin=0.1, ymax=0.9, color='black', linewidth=2)

# --- Format tick labels with 2 decimal places ---
tick_labels = [f"{f:.2f}" for f in xvals]
plt.xticks(xvals, tick_labels, rotation=90)

plt.yticks([])
plt.gca().spines['left'].set_visible(False)
plt.xlabel('Mode frequency (MHz)')
plt.title('Motional Modes (Radial)')

bandwidth = xvals.max() - xvals.min()
plt.annotate(f'Bandwidth: {bandwidth:.2f} MHz',
             xy=(0.5*(xvals.max()+xvals.min()), 0.04),
             ha='center', va='bottom')

plt.tight_layout()
plt.show()

print("z positions (µm):", z_positions_m*1e6)

