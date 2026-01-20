# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 11:41:47 2025

@author: iontrap
"""

# poisson_bright_dark_min_bars_threshold.py
# Bars for Poisson PMFs (bright per NA + dark), a vertical threshold line,
# and fidelity numbers at that threshold.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import poisson

# =========================
# EDIT THESE PARAMETERS
# =========================
lambda_baseline = 20.0    # mean bright counts at NA_baseline
NA_baseline = 0.28
eta_baseline = 0.020      # collection efficiency at NA_baseline (e.g., 2.0%)

cases = [
    (0.28, 0.020),        # baseline
    (0.35, 0.032),        # higher NA example
]

lambda_dark = 1.44        # mean dark counts (assumed NA-independent)
lambda_bg   = 0.0         # optional constant background (not scaled)

THRESHOLD_K = 9          # <<<----- decision rule: "bright" if N >= THRESHOLD_K

# =========================
# helpers
# =========================
def scale_lambda(lam_base, eta_base, eta_new, lam_bg=0.0):
    """Scale bright mean with collection efficiency; keep background unscaled."""
    signal = max(lam_base - lam_bg, 0.0)
    return signal * (eta_new / eta_base) + lam_bg

def F_bright(k_star, lam_b):
    """Bright fidelity: P(N >= k* | bright)."""
    return float(1 - poisson.cdf(k_star - 1, lam_b))

def F_dark(k_star, lam_d):
    """Dark fidelity: P(N < k* | dark)."""
    return float(poisson.cdf(k_star - 1, lam_d))

# =========================
# compute means
# =========================
rows = []
for NA, eta in cases:
    lam_b = scale_lambda(lambda_baseline, eta_baseline, eta, lambda_bg)
    rows.append({"NA": NA, "eta": eta, "lambda_b": lam_b})
df = pd.DataFrame(rows)

# x-axis (counts)
max_lambda = max(df["lambda_b"].max(), lambda_dark)
N = np.arange(0, int(np.ceil(max_lambda + 6*np.sqrt(max_lambda))) + 1)

# =========================
# plot (bars) + threshold
# =========================
num_series = len(cases) + 1  # +1 for dark
group_width = 0.9
bar_w = group_width / num_series
offsets = np.linspace(-group_width/2 + bar_w/2, group_width/2 - bar_w/2, num_series)

plt.figure(figsize=(10, 5))

# dark bars
pmf_dark = poisson.pmf(N, lambda_dark)
plt.bar(N + offsets[0], pmf_dark, width=bar_w, color='none',
        edgecolor='crimson', linewidth=1.6, label=f'Dark (λ_d={lambda_dark:.2f})')

# bright bars
for i, (_, row) in enumerate(df.iterrows(), start=1):
    pmf_b = poisson.pmf(N, row["lambda_b"])
    plt.bar(N + offsets[i], pmf_b, width=bar_w, alpha=0.6,
            label=f'Bright NA {row.NA:.2f} (λ_b={row.lambda_b:.1f})')

# vertical threshold line
plt.axvline(THRESHOLD_K, linestyle='--', linewidth=2, label=f'Threshold k*={THRESHOLD_K}', color='k')

plt.xlabel('Detected counts N')
plt.ylabel('Poisson PMF  P(N|λ)')
plt.title('Poisson distributions (bars) with decision threshold')
plt.legend()
plt.tight_layout()
plt.show()

# =========================
# fidelities at threshold
# =========================
print("\n=== Fidelities at threshold k* =", THRESHOLD_K, "===")
Fd = F_dark(THRESHOLD_K, lambda_dark)
print(f"Dark fidelity F_d = P(N < k* | dark, λ_d={lambda_dark:.2f}) = {Fd:.6f}")

for _, r in df.iterrows():
    Fb = F_bright(THRESHOLD_K, r.lambda_b)
    Favg = 0.5 * (Fb + Fd)  # equal priors
    print(f"Bright (NA {r.NA:.2f}, λ_b={r.lambda_b:.1f}): "
          f"F_b = P(N >= k* | bright) = {Fb:.6f}  |  F_avg = {Favg:.6f}")

