# -*- coding: utf-8 -*-
"""
EIT-cooling steady-state scans — simplified, readable version
Keeps your physics & plots, but removes duplicate code and clarifies flow.
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import qutip as qt
from sympy import symbols, Eq, solve

# ----------------------------------------
# 0) Constants & knobs (edit here)
# ----------------------------------------
TWOPI = 2 * math.pi

Gamma         = TWOPI * 19.6e6        # total decay rate out of |e>
Omega_s_minus = TWOPI * 22e6          # Rabi for σ−
Omega_s_plus  = TWOPI * 22e6          # Rabi for σ+
Omega_pi      = TWOPI * 1.5e6         # Rabi for π
Delta_d       = TWOPI * 45e6          # drive detuning
delta_z       = TWOPI * 5e6           # Zeeman splitting

delsig_m = Delta_d - delta_z          # Δ−
delsig_p = Delta_d + delta_z          # Δ+
delpi2   = delsig_p                   # convenience alias (your code used this)

nu       = TWOPI * 3.2e6              # single trap freq for sideband points
nu2      = np.linspace(0.5e6, 6e6, 909) * TWOPI  # scan of trap freqs (rad/s)

Delta_p_values = np.linspace(-30e6, 80e6, 1000) * TWOPI  # probe detuning scan (rad/s)
Delta_p_values2 = np.linspace(delsig_m - 4e6, delsig_m + 4e6, 1000)

# ----------------------------------------
# 1) Steady-state solver → pop in |e>
# ----------------------------------------
def pop_e_steady(Delta_pi):
    """
    Returns excited-state population vs probe detuning Δπ.
    Accepts scalar or array (rad/s). Output is np.ndarray (complex).
    Basis ordering: |e>, |+>, |0>, |->
    """
    Delta_pi = np.atleast_1d(Delta_pi)

    N = 4
    psi_e    = qt.basis(N, 0)
    psi_plus = qt.basis(N, 1)
    psi_0    = qt.basis(N, 2)
    psi_minus= qt.basis(N, 3)

    # Decay channels: e -> {+, 0, -} with equal branching
    c_ops = [
        np.sqrt(Gamma/3) * psi_plus  * psi_e.dag(),
        np.sqrt(Gamma/3) * psi_0     * psi_e.dag(),
        np.sqrt(Gamma/3) * psi_minus * psi_e.dag(),
    ]

    pop = np.empty(Delta_pi.size, dtype=complex)

    for i, dpi in enumerate(Delta_pi):
        H = qt.Qobj([
            [0,                   Omega_s_minus/2,  -Omega_pi/2,       Omega_s_plus/2],
            [Omega_s_minus/2,     Delta_d + delta_z, 0,                0],
            [-Omega_pi/2,         0,                 dpi,              0],
            [Omega_s_plus/2,      0,                 0,       Delta_d - delta_z],
        ])
        rho_ss = qt.steadystate(H, c_ops, method='direct', tol=1e-10)
        pop[i] = (psi_e.dag() * rho_ss * psi_e).tr()

    return pop


# ----------------------------------------
# 2) nbar estimator from carrier/sidebands
# ----------------------------------------
def nbar_ratio(car, rsb, bsb):
    """ nbar = (carrier + blue) / (red - blue) """
    car, rsb, bsb = np.asarray(car), np.asarray(rsb), np.asarray(bsb)
    return (car + bsb) / (rsb - bsb)


# ----------------------------------------
# 3) AC light shift
# ----------------------------------------
def cooling_bright_line_from_eigs(Delta_plus, Delta_minus, Omega_sigma_plus, Omega_sigma_minus):
    H = np.array([[0.0,                     Omega_sigma_minus/2.0,  Omega_sigma_plus/2.0],
                  [Omega_sigma_minus/2.0,   Delta_plus,             0.0],
                  [Omega_sigma_plus/2.0,    0.0,                    Delta_minus]], dtype=float)
    # Hermitian → use eigvalsh for sorted real eigenvalues
    lam = np.linalg.eigvalsh(H)
    return lam[-1]  # largest eigenvalue = right bright branch

# compute right bright position (in rad/s)
x_cool = cooling_bright_line_from_eigs(delsig_p, delsig_m, Omega_s_plus, Omega_s_minus)


# ----------------------------------------
# 4) Compute all curves
# ----------------------------------------
# Carrier / sidebands vs Δp (using single ν)
car   = pop_e_steady(Delta_p_values)
rsb   = pop_e_steady(Delta_p_values + nu)
bsb   = pop_e_steady(Delta_p_values - nu)
nbar1 = nbar_ratio(car, rsb, bsb)

# Carrier (scalar) & sidebands vs ν (for bandwidth)
car2  = pop_e_steady(delpi2)           # scalar → 1-element array
rsb2  = pop_e_steady(delpi2 + nu2)     # array
bsb2  = pop_e_steady(delpi2 - nu2)     # array
nbar2 = nbar_ratio(car2, rsb2, bsb2)




# ----------------------------------------
# 5) Convenience helpers (minima, masks)
# ----------------------------------------
def argmin_xy(x, y):
    i = int(np.argmin(y))
    return i, float(x[i]), float(y[i])

def mask_between(x, a, b):
    return (x >= a) & (x <= b)


# ----------------------------------------
# 6) Plots (three main figures + bandwidth)
# ----------------------------------------
# Fig 1: Pop in |e> vs Δπ
plt.figure(figsize=(10, 4))
plt.plot(Delta_p_values/(TWOPI*1e6), car.real, 'b-')
plt.axvline(delsig_m/(TWOPI*1e6), color='k', ls='--', label=r'$\Delta_{-}$')
plt.axvline(delsig_p/(TWOPI*1e6), color='k', ls=':',  label=r'$\Delta_{+}$')
plt.xlabel(r'$\Delta_{p}/2\pi$ (MHz)')
plt.ylabel(r'Population in $|e\rangle$')
plt.title(r'Absorption Profile')
#plt.legend()
plt.tight_layout()

# --- Fig 2: Zoom near dark resonance; mark null + eigenvalue-based right bright line ---
# choose a window that covers Δ+ and the right bright line
lo = min(delsig_p, x_cool) - TWOPI*5e6
hi = max(delsig_p, x_cool) + TWOPI*10e6
zoom_mask = mask_between(Delta_p_values, lo, hi)
x_zoom = Delta_p_values[zoom_mask]
y_zoom = car[zoom_mask].real

# dark null (minimum)
i_null = int(np.argmin(y_zoom))
x_null = float(x_zoom[i_null])

# separation from null to the right bright line (MHz)
sep_MHz = abs((x_cool - x_null) / (TWOPI*1e6))

plt.figure(figsize=(10, 4))
plt.plot(x_zoom/(TWOPI*1e6), y_zoom, 'b-')
plt.axvline(x_null/(TWOPI*1e6), color='k', ls='--')
#plt.axvline(x_cool/(TWOPI*1e6), color='g', ls='-',
            #label=fr'AC Light Shift ≈ {sep_MHz:.2f} MHz)')
plt.axvline((delsig_p + nu)/(TWOPI*1e6), color='r', ls='-',
            label=fr'COM RSB ≈ {nu/(TWOPI*1e6):.2f} MHz)')
plt.axvline((delsig_p - nu)/(TWOPI*1e6), color='b', ls='-',
            label=fr'COM BSB ≈ {nu/(TWOPI*1e6):.2f} MHz)')
plt.xlabel(r'$\Delta_{p}/2\pi$ (MHz)')
plt.ylabel(r'Population in $|e\rangle$')
plt.title('Absorption Profile')
plt.legend()
plt.tight_layout()

# Fig 3: nbar vs Δπ (log scale) + optimum marker
zoom_mask2 = mask_between(Delta_p_values, delsig_p - TWOPI*0.5e6, delsig_p + TWOPI*2.5e6)
i_min, x_min, y_min = argmin_xy(Delta_p_values[zoom_mask2], nbar1.real[zoom_mask2])

# Evaluate nbar at Δp = Δd + δz (≈ 50 MHz)
i_delp = np.argmin(np.abs(Delta_p_values - delsig_p))
nbar_at_delp = nbar1.real[i_delp]

plt.figure(figsize=(10, 4))
plt.plot(Delta_p_values[zoom_mask2]/(TWOPI*1e6), nbar1.real[zoom_mask2], 'b-')
plt.yscale('log')
plt.axvline(delsig_p/(TWOPI*1e6), color='k', ls='--',
            label=fr'$\Delta_p = \Delta_d + \delta_z$' 
                  f'\n$\overline{{n}}$ = {nbar_at_delp:.3g}')
#plt.axvline(x_min/(TWOPI*1e6), color='g', ls=':',  label=f'opt ≈ {sep_MHz:.2f} MHz)')
plt.xlabel(r'$\Delta_{p}/2\pi$ (MHz)')
plt.ylabel(r'$\bar{n}$')
plt.title(r'$\bar{n}$ vs $\Delta_{p}$')
plt.legend()
plt.tight_layout()

# Estimate bandwidth where nbar <= 0.1
mask_bw = nbar2.real <= 0.1
nu_ok = nu2[mask_bw]
bandwidth_mhz = (nu_ok.max() - nu_ok.min()) / (TWOPI*1e6)
shift = x_cool - delsig_p

fig, ax = plt.subplots(figsize=(10,4))
ax.plot(nu2/(TWOPI*1e6), nbar2.real, 'b-')
ax.set_yscale('log')
ax.axhline(0.1, color='k', ls='--', label='threshold = 0.1')

if np.any(mask_bw):
    x1, x2 = nu_ok.min()/(TWOPI*1e6), nu_ok.max()/(TWOPI*1e6)
    ax.axvspan(x1, x2, alpha=0.15, color='C1',
               label=rf'bandwidth: {bandwidth_mhz:.2f} MHz')

ax.set_xlabel(r'$\nu/2\pi$ (MHz)')
ax.set_ylabel(r'$\bar{n}$')
ax.set_title('Cooling bandwidth')
ax.legend()
fig.tight_layout()


# Print a couple of useful numbers
print(f"Optimum Δπ/2π ≈ {x_min/(TWOPI*1e6):.3f} MHz, min nbar ≈ {y_min.real:.3g}")
print(f"AC light shift δ+/2π ≈ {shift/(TWOPI*1e6):.3f} MHz")


