# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 11:40:36 2026

@author: iontrap
"""

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# 0) Physical constants / inputs
# ----------------------------
e = 1.602176634e-19
eps0 = 8.8541878128e-12
kC = 1.0 / (4*np.pi*eps0)
hbar = 1.054571817e-34
amu = 1.66053906660e-27

mYb171 = 171 * amu

N = 5

# Trap frequencies (Hz)
wz_Hz = 0.4e6      # axial
wx_Hz = 1.5e6      # radial (x). assumes y similar
wz = 2*np.pi*wz_Hz
wx = 2*np.pi*wx_Hz

# Raman geometry -> effective delta-k magnitude
lam = 355e-9
k = 2*np.pi/lam
theta_deg = 90.0                 # angle between arms
dk = 2*k*np.sin(np.deg2rad(theta_deg)/2)

# ----------------------------
# 1) Equilibrium positions along z (1D chain)
#    Solve dU/dz_i = 0 for harmonic + Coulomb
# ----------------------------
def equilibrium_positions_1d(N, wz, m, z0_guess_scale=5e-6, max_iter=200, tol=1e-12):
    """
    Returns equilibrium positions z_i (meters), centered near 0.
    Uses Newton iteration on forces.
    """
    # initial guess: equally spaced
    z = (np.arange(N) - (N-1)/2) * z0_guess_scale

    def forces(z):
        # Force = -dU/dz
        F = -m*(wz**2)*z
        for i in range(N):
            for j in range(N):
                if j == i:
                    continue
                dz = z[i] - z[j]
                F[i] += kC * e*e * dz / (np.abs(dz)**3)
        return F

    def jacobian(z):
        # J_ij = dF_i/dz_j
        J = np.zeros((N, N))
        # harmonic part
        J += np.diag(-m*wz**2 * np.ones(N))
        # Coulomb part
        for i in range(N):
            for j in range(N):
                if j == i:
                    continue
                dz = z[i] - z[j]
                r = np.abs(dz)
                # d/dz_i of (dz/r^3) = (1/r^3) - 3 dz^2 / r^5 = (1 - 3)/r^3 for 1D? careful:
                # Since r=|dz|, dz^2 = r^2, so:
                # d/dz_i (dz/r^3) = 1/r^3 - 3*dz*(d r / dz_i)/r^4
                # In 1D with r=|dz| and dz_i derivative gives sign(dz), simplifies to -2/r^3
                # The exact 1D second derivative of 1/|dz| is 2/|dz|^3.
                # For force form, the Jacobian elements become:
                coef = kC*e*e
                # contribution to J_ii
                J[i, i] += -2*coef / (r**3)
                # contribution to J_ij
                J[i, j] += +2*coef / (r**3)
        return J

    for _ in range(max_iter):
        F = forces(z)
        if np.linalg.norm(F, ord=np.inf) < tol:
            break
        J = jacobian(z)
        dz = np.linalg.solve(J, -F)
        z = z + dz

    # center to remove floating drift
    z = z - np.mean(z)
    return z

z_eq = equilibrium_positions_1d(N, wz, mYb171)

# ----------------------------
# 2) Radial normal modes (transverse x)
#    Small oscillations in x around x=0:
#    U ~ (1/2)m wx^2 sum x_i^2  + (1/2) sum_{i<j} kC e^2 (x_i - x_j)^2 / |z_ij|^3
#    => dynamical matrix:
#       D_ii = wx^2 - sum_{j!=i} (2*kC e^2)/(m*|z_ij|^3)
#       D_ij = + (2*kC e^2)/(m*|z_ij|^3)  for i!=j
#    Eigenvalues give mode frequencies^2 in rad/s.
# ----------------------------
def radial_mode_matrix(z, wx, m):
    N = len(z)
    D = np.zeros((N, N), dtype=float)
    for i in range(N):
        s = 0.0
        for j in range(N):
            if j == i:
                continue
            rij = np.abs(z[i]-z[j])
            c = 2*kC*e*e/(m*(rij**3))
            s += c
            D[i, j] += c
        D[i, i] += wx**2 - s
    return D

D = radial_mode_matrix(z_eq, wx, mYb171)
w2, B = np.linalg.eigh(D)   # B columns are eigenvectors
w_modes = np.sqrt(np.maximum(w2, 0.0))  # rad/s

# Sort by frequency (lowest to highest)
idx = np.argsort(w_modes)
w_modes = w_modes[idx]
B = B[:, idx]

# Normalize eigenvectors (should already be orthonormal from eigh)
# B_{i m} is participation of ion i in mode m.

# ----------------------------
# 3) Lamb–Dicke parameters for each ion/mode
#    eta_{i m} = dk * b_{i m} * sqrt(hbar/(2 m w_m))
# ----------------------------
eta = np.zeros((N, N))
for m in range(N):
    zpf = np.sqrt(hbar/(2*mYb171*w_modes[m]))
    eta[:, m] = dk * B[:, m] * zpf

# ----------------------------
# 4) Smooth MS-style gate model (global drive)
#    Assume a bichromatic beatnote with instantaneous frequency:
#        mu(t) = w_target + delta_profile(t)
#    For each mode m, detuning in the rotating picture:
#        delta_m(t) = mu(t) - w_m
#    Define phase:
#        phi_m(t) = ∫_0^t delta_m(t') dt'
#
#    With global Rabi envelope Omega(t), define:
#        A_m = ∫_0^{tg} Omega(t) * exp(i phi_m(t)) dt
#    Residual displacement:
#        alpha_{i m}(tg) = -i * eta_{i m} * A_m
#
#    Spin-spin phase kernel (global Omega):
#        K_m = ∫ dt1 Omega(t1) e^{i phi_m(t1)} * [∫_0^{t1} Omega(t2) e^{-i phi_m(t2)} dt2]
#        chi_{ij} = sum_m eta_{i m} eta_{j m} * Im(K_m)
#
#    NOTE: This captures the standard MS geometric-phase structure in a way
#    that generalizes cleanly to time-dependent detuning.
# ----------------------------


def aese_delta_ramp(t, delta_max, delta_min, tau_d, j=3):
    """
    Implements Eqs. (18)-(19) from the paper.

    Parameters
    ----------
    t : array
        Time array, 0 <= t <= tau_d
    delta_max, delta_min : float
        Detuning endpoints (rad/s)
    tau_d : float
        Duration of the detuning ramp
    j : int
        Power-law exponent (paper uses j=3)

    Returns
    -------
    delta(t) : array
        Time-dependent detuning
    """
    b = delta_max**(-j)
    c = (2 / tau_d) * (delta_min**(-j) - delta_max**(-j))

    g = t/2 - (tau_d/(4*np.pi))*np.sin(2*np.pi*t/tau_d)

    return (b + c*g)**(-1/j)

def simulate_smooth_gate(w_modes, eta, target_mode=0,
                         tg=250e-6, dt=50e-9,
                         Omega0=(2*np.pi*80e3)/eta,
                         delta_start=2*np.pi*400e3,
                         delta_min=2*np.pi*30e3,
                         tau_d=60e-6,
                         t_edge=20e-6,
                         j=3):
    """
    Returns:
      t (array)
      Omega(t) (array)
      mu(t) (array)
      chi_ij (NxN) accumulated spin-spin phases (radians)
      alpha_im (NxN) residual displacements (complex)
      A_m (N,) mode integrals (complex)
      K_m (N,) mode phase kernels (complex)
    """
    N = eta.shape[0]
    t = np.arange(0, tg, dt)

    # Global amplitude envelope (smooth on/off)
    Omega = Omega0

    # AESE detuning profile: ramp down for tau_d, hold, ramp up for tau_d
    delta_t = np.full_like(t, delta_min, dtype=float)

    mask1 = t < tau_d
    mask2 = t > (tg - tau_d)

    # ramp down
    delta_t[mask1] = aese_delta_ramp(
        t[mask1],
        delta_max=delta_start,
        delta_min=delta_min,
        tau_d=tau_d,
        j=j
    )

    # ramp up (time-reversed)
    delta_t[mask2] = aese_delta_ramp(
        tg - t[mask2],
        delta_max=delta_start,
        delta_min=delta_min,
        tau_d=tau_d,
        j=j
    )

    mu_t = w_modes[target_mode] + delta_t  # instantaneous beatnote frequency

    # Mode integrals
    A_m = np.zeros(N, dtype=np.complex128)
    K_m = np.zeros(N, dtype=np.complex128)

    for m in range(N):
        delta_m = mu_t - w_modes[m]
        phi_m = np.cumsum(delta_m) * dt

        exp_iphi = np.exp(1j * phi_m)
        exp_imphi = np.conj(exp_iphi)

        # A_m = ∫ Omega e^{i phi} dt
        A_m[m] = np.sum(Omega * exp_iphi) * dt

        # K_m = ∫ dt1 Omega(t1) e^{i phi(t1)} [∫_0^{t1} Omega(t2)e^{-i phi(t2)} dt2]
        b = np.cumsum(Omega * exp_imphi) * dt
        K_m[m] = np.sum(Omega * exp_iphi * b) * dt

    # Residual displacement: alpha_{i m} = -i eta_{i m} A_m
    alpha_im = -1j * eta * A_m[None, :]

    # Spin-spin phases: chi_ij = sum_m eta_im eta_jm Im(K_m)
    ImK = np.imag(K_m)
    chi = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            chi[i, j] = np.sum(eta[i, :] * eta[j, :] * ImK)

    return t, Omega, mu_t, chi, alpha_im, A_m, K_m


def plot_J_couplings(J, title=r"$J_{ij}$", normalize=True, show_diag_white=True):
    """
    Plot J_ij as a heatmap..
    """
    Jp = np.array(J, dtype=float, copy=True)

    # optionally whiten diagonal
    if show_diag_white:
        np.fill_diagonal(Jp, 0.0)

    # normalize off-diagonal entries to [-1, 1]
    if normalize:
        mask = ~np.eye(Jp.shape[0], dtype=bool)
        max_abs = np.max(np.abs(Jp[mask])) if np.any(mask) else np.max(np.abs(Jp))
        if max_abs > 0:
            Jp = Jp / max_abs

    n = Jp.shape[0]

    fig, ax = plt.subplots(figsize=(3.6, 3.2), dpi=160)

    im = ax.imshow(Jp, vmin=-1, vmax=1, cmap="RdBu_r", interpolation="nearest")

    # Remove tick labels (or keep if you want)
    ax.set_xticks([])
    ax.set_yticks([])

    # Gridlines like a matrix
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", linestyle="-", linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Bold border
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

    ax.set_title(title, pad=8)

    # Colorbar with -1,0,1 ticks
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_ticks([-1, 0, 1])
    cbar.ax.set_ylabel(r"$J_{ij}$", rotation=90)

    plt.tight_layout()
    plt.show()


# 5) Run 
# ----------------------------
t, Omega_t, mu_t, chi, alpha_im, A_m, K_m = simulate_smooth_gate(
    w_modes=w_modes, eta=eta,
    target_mode=0,
    tg=250e-6,
    dt=100e-9,
    Omega0=2*np.pi*35e3,
    delta_start=2*np.pi*400e3,
    delta_min=2*np.pi*50e3,
    tau_d=80e-6,     # <-- AESE ramp duration
    t_edge=20e-6,    # <-- amplitude turn on/off time
    j=3
)

alpha_rms_per_mode = np.sqrt(np.mean(np.abs(alpha_im)**2, axis=0))
alpha_rms_per_ion  = np.sqrt(np.mean(np.abs(alpha_im)**2, axis=1))

print("=== Equilibrium z positions (um) ===")
print(np.round(z_eq*1e6, 3))

print("\n=== Radial mode frequencies (MHz) ===")
print(np.round(w_modes/(2*np.pi*1e6), 6))

print("\n=== RMS residual displacement |alpha_im| (per mode) ===")
print(alpha_rms_per_mode)

print("\n=== Spin-spin phase matrix chi_ij (radians) ===")
np.set_printoptions(precision=4, suppress=True)
print(chi)

plot_J_couplings(chi, title=r"$J_{ij}$ (normalized)")
