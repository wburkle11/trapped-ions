# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 15:41:45 2025

@author: iontrap
"""

# yb_readout_bright_and_dark_with_leaks_opt_threshold_twoNA.py
# Adds a second bright distribution for NA = 0.35 (η = 0.032).
# Threshold is still optimized using the original NA = 0.28 (η = 0.02) case.

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.special import gammainc, factorial

# =========================
# EDIT THESE PARAMETERS
# =========================
eta = 0.02               # collection efficiency (NA = 0.28 -> ~2%)
eta_hi = 0.032           # collection efficiency for NA = 0.35
lambda_bright_total = 11.3   # desired *total* bright mean counts (signal + background) for NA=0.28
lambda_bg = 0.59             # constant Poisson background mean (sets dark ~1.4)

alpha2 = 2.0e-4          # bright -> dark leak per emitted photon
alpha1 = 0.4e-5          # dark -> bright leak per emitted photon 

p_bright_prior = 0.5     # prior for bright (equal priors by default)

# Derived: the "signal" part of the bright expectation (no-leak λ0 = η R_sc t)
lambda0_signal = max(lambda_bright_total - lambda_bg, 0.0)          # for NA=0.28
lambda0_signal_hi = (eta_hi/eta) * lambda0_signal                   # scale for NA=0.35 (3.2%)

# =========================
# PMFs WITHOUT BACKGROUND (signal channel only)
# =========================
def p_bright_with_leak_signal(n, lambda0, eta, alpha2):
    """
    p_bright(n) = exp[-(1+α2/η)λ0] λ0^n / n!  +  (α2/η)/(1+α2/η)^(n+1) * P(n+1, (1+α2/η)λ0)
    where P is the regularized lower incomplete gamma (scipy.special.gammainc).
    """
    n = np.asarray(n, dtype=np.int64)
    r = alpha2 / eta
    x = (1.0 + r) * lambda0
    term1 = np.exp(-x) * (lambda0 ** n) / factorial(n, exact=False)
    term2 = (r / ((1.0 + r) ** (n + 1))) * gammainc(n + 1, x)
    return term1 + term2

def p_dark_with_leak_signal(n, lambda0, eta, alpha1):
    """
    p_dark(n) = e^{-α1 λ0 / η} * [ δ_{n0} + (α1/η)/(1-α1/η)^{n+1} * P(n+1, (1-α1/η)λ0) ]
    (valid for α1/η < 1; typical in experiments)
    """
    n = np.asarray(n, dtype=np.int64)
    r = alpha1 / eta
    if r >= 1.0:
        raise ValueError("This model assumes α1/η < 1.")
    x = (1.0 - r) * lambda0
    pref = np.exp(-alpha1 * lambda0 / eta)
    delta = (n == 0).astype(float)
    term2 = (r / ((1.0 - r) ** (n + 1))) * gammainc(n + 1, x)
    return pref * (delta + term2)

# =========================
# Convolution with Poisson background
# =========================
def convolve_with_background(pmf_signal, lambda_bg):
    N_bg_max = int(np.ceil(lambda_bg + 8 * np.sqrt(max(lambda_bg, 1.0))))
    k = np.arange(0, max(40, N_bg_max) + 1)
    pmf_bg = poisson.pmf(k, lambda_bg)
    pmf = np.convolve(pmf_signal, pmf_bg)
    cutoff = np.searchsorted(np.cumsum(pmf), 0.999999)
    return pmf[:max(cutoff+1, 60)]

# =========================
# Build PMFs on a sufficiently wide range (based on NA=0.28 signal)
# =========================
N0_max = int(np.ceil(lambda0_signal + 8*np.sqrt(max(lambda0_signal, 1.0))))
N = np.arange(0, max(60, N0_max) + 1)

# Signal-only PMFs
pmf_b_sig      = p_bright_with_leak_signal(N, lambda0_signal,    eta,    alpha2)  # NA=0.28
pmf_b_sig_hiNA = p_bright_with_leak_signal(N, lambda0_signal_hi, eta_hi, alpha2)  # NA=0.35
pmf_d_sig      = p_dark_with_leak_signal(  N, lambda0_signal,    eta,    alpha1)

# Convolve with background to get the *total* PMFs
pmf_bright      = convolve_with_background(pmf_b_sig,      lambda_bg)   # NA=0.28
pmf_bright_hiNA = convolve_with_background(pmf_b_sig_hiNA, lambda_bg)   # NA=0.35
pmf_dark        = convolve_with_background(pmf_d_sig,      lambda_bg)

# Resize N to fit the convolved PMFs (use longest length)
L = max(len(pmf_bright), len(pmf_bright_hiNA), len(pmf_dark))
def pad(p): 
    return np.pad(p, (0, L-len(p)), constant_values=0.0)
pmf_bright, pmf_bright_hiNA, pmf_dark = pad(pmf_bright), pad(pmf_bright_hiNA), pad(pmf_dark)
N = np.arange(0, L)

# =========================
# Find optimal threshold k* (maximize average fidelity using NA=0.28 bright)
# =========================
cdf_b = np.cumsum(pmf_bright)
cdf_d = np.cumsum(pmf_dark)
tail_b = 1.0 - cdf_b

def F_b_for_k_from(pmf_tail, k):      # P(N >= k)
    return 1.0 if k <= 0 else pmf_tail[k-1]

def F_d_for_k_from(cdf, k):           # P(N < k)
    return 0.0 if k <= 0 else cdf[k-1]

def tail_from_pmf(p): return 1.0 - np.cumsum(p)

Ks = np.arange(0, N[-1] + 1, dtype=int)
F_b = np.array([F_b_for_k_from(tail_b, k)     for k in Ks])   # using NA=0.28 bright
F_d = np.array([F_d_for_k_from(cdf_d,  k)     for k in Ks])

pri_b = float(p_bright_prior); pri_d = 1.0 - pri_b
F_avg = pri_b * F_b + pri_d * F_d

k_opt = int(Ks[np.argmax(F_avg)])
Fb_opt_028 = F_b_for_k_from(tail_b, k_opt)
Fd_opt     = F_d_for_k_from(cdf_d,  k_opt)
Favg_opt   = pri_b * Fb_opt_028 + pri_d * Fd_opt

# Also compute bright fidelity for the NA=0.35 curve at the SAME k*
tail_b_hi = tail_from_pmf(pmf_bright_hiNA)
Fb_opt_035 = F_b_for_k_from(tail_b_hi, k_opt)
Favg_opt_035 = pri_b * Fb_opt_035 + pri_d * Fd_opt

print(f"\n=== Optimal threshold using NA=0.28 ===")
print(f"k*_opt = {k_opt}")
print(f"F_b(NA=0.28, k*) = {Fb_opt_028:.6f}   F_d(k*) = {Fd_opt:.6f}   F_avg(k*) = {Favg_opt:.6f}")
print(f"F_b(NA=0.35, k*) = {Fb_opt_035:.6f}   SAME AS ABOVE        F_avg(k*) = {Favg_opt_035:.6f}")

# =========================
# Optimal threshold for NA=0.35 bright vs dark (minimize overlap)
# =========================
cdf_b_hi = np.cumsum(pmf_bright_hiNA)
tail_b_hi = 1.0 - cdf_b_hi

# Fidelity arrays vs threshold k, using NA=0.35 bright + the same dark PMF
F_b_035 = np.array([F_b_for_k_from(tail_b_hi, k) for k in Ks])  # P_bright_hi(N >= k)
# F_d already computed: P_dark(N < k)

# Average fidelity (same priors)
F_avg_035_allk = pri_b * F_b_035 + pri_d * F_d
k_opt_035 = int(Ks[np.argmax(F_avg_035_allk)])

# Fidelities at that threshold
Fb_035_at_own_k = F_b_for_k_from(tail_b_hi, k_opt_035)
Fd_at_own_k     = F_d_for_k_from(cdf_d,     k_opt_035)
Favg_035_at_own = pri_b * Fb_035_at_own_k + pri_d * Fd_at_own_k

# Also report the overlap (misclassification probability) at that threshold
# overlap = pri_b * P_bright_hi(N < k) + pri_d * P_dark(N >= k)
overlap_035 = pri_b * (1.0 - Fb_035_at_own_k) + pri_d * (1.0 - Fd_at_own_k)

print("\n=== Optimal threshold using NA=0.35 ===")
print(f"k*_opt = {k_opt_035}")
print(
    f"F_b(NA=0.35, k*) = {Fb_035_at_own_k:.6f}   "
    f"F_d(k*) = {Fd_at_own_k:.6f}   "
    f"F_avg(k*) = {Favg_035_at_own:.6f}"
)

# =========================
# Plot bars + optimal threshold
# =========================
bar_w = 0.28  # a bit narrower so three series can overlap visibly

plt.figure(figsize=(10,5))
plt.bar(N, pmf_dark,            width=bar_w, color='none', edgecolor='crimson', linewidth=1.6, label='Dark Ion')
plt.bar(N, pmf_bright,          width=bar_w, alpha=0.55, label='Bright Ion (NA=0.28)')
plt.bar(N, pmf_bright_hiNA,     width=bar_w, alpha=0.55, label='Bright Ion (NA=0.35)')
plt.axvline(k_opt, linestyle='--', linewidth=2, color='k', label='Threshold')

plt.xlabel('Number of photons')
plt.ylabel('Probability  P(N)')
plt.xlim(-1, 40)
plt.legend()
plt.tight_layout()
plt.show()




