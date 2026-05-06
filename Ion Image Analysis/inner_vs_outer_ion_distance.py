# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 13:00:51 2026

@author: iontrap
"""

# -*- coding: utf-8 -*-
"""
Plot d1 and d2 vs V_center from a set of ion-chain images.

Assumes you already have a separate script/module containing:
    - find_ion_centers(...)
    - nearest_neighbor_distances_projected(...)

For an odd-number ion chain (e.g. 11 ions):
    d1 = average of the two spacings adjacent to the center ion
    d2 = average of the two outermost spacings
"""

import numpy as np
import matplotlib.pyplot as plt
import tifffile

# Import from your previous script
from projected_axial_ion_distance import find_ion_centers, nearest_neighbor_distances_projected


# --------------------------------------------------
# User input: file paths and corresponding V_center values
# --------------------------------------------------
image_paths = [
    r"Z:\Lab Data\EMCCD\Images\raw_img_Tue Mar 31 2026_09.01.35_1463.tif",
    r"Z:\Lab Data\EMCCD\Images\raw_img_Tue Mar 31 2026_09.04.07_1465.tif",
    r"Z:\Lab Data\EMCCD\Images\raw_img_Tue Mar 31 2026_09.04.31_1466.tif",
    r"Z:\Lab Data\EMCCD\Images\raw_img_Tue Mar 31 2026_09.05.06_1467.tif",
    r"Z:\Lab Data\EMCCD\Images\raw_img_Tue Mar 31 2026_09.05.58_1468.tif",
]

V_center_values = np.array([
    0,
    -0.2,
    -0.4,
    -0.6,
    -0.75,
    
], dtype=float)


# --------------------------------------------------
# Helper: load tif and average if stack
# --------------------------------------------------
def load_avg_image(tif_path):
    data = tifffile.imread(tif_path)
    return np.mean(data, axis=0) if data.ndim == 3 else data


# --------------------------------------------------
# Helper: compute d1 and d2 from spacing list
# --------------------------------------------------
def extract_d1_d2(dists):
    """
    Automatically extract d1 and d2 for either odd or even ion chains.

    Definitions
    -----------
    d2:
        Average of the two edge spacings:
            d2 = (dists[0] + dists[-1]) / 2

    d1:
        A symmetric 'center spacing' quantity:
        - If number of ions is odd:
            average of the two spacings around the center ion
        - If number of ions is even:
            average of the two middle spacings

    Parameters
    ----------
    dists : 1D ndarray
        Nearest-neighbor spacings. If there are N ions, len(dists) = N - 1.

    Returns
    -------
    d1 : float
        Center spacing metric
    d2 : float
        Edge spacing metric
    info : dict
        Extra information about which indices were used
    """
    n_spacings = len(dists)
    n_ions = n_spacings + 1

    if n_spacings < 2:
        raise ValueError("Need at least 3 ions to define d1 and d2.")
        
    d2 = 0.5*(dists[0] + dists[-1])    

    # Edge spacing metric
    d2L = dists[0] 
    d2R = dists[-1]

    # Center spacing metric
    if n_ions % 2 == 1:
        # Odd number of ions: average spacings around the center ion
        # Example: 9 ions -> 8 spacings -> use d3 and d4
        left_idx = n_spacings // 2 - 1
        right_idx = n_spacings // 2
        d1 = 0.5*(dists[left_idx] + dists[right_idx])
        definition = "odd ion chain: spacings around center ion"
    else:
        # Even number of ions: use the single middle spacing
        center_idx = n_spacings // 2

        d1 = dists[center_idx]
        left_idx = center_idx
        right_idx = center_idx

        definition = "even ion chain: single middle spacing"

    info = {
        "n_ions": n_ions,
        "n_spacings": n_spacings,
        "d1_indices": (left_idx, right_idx),
        "d2_indices": (0, n_spacings - 1),
        "definition": definition,
    }

    return d1, d2, d2L, d2R, info


# --------------------------------------------------
# Main analysis
# --------------------------------------------------
if len(image_paths) != len(V_center_values):
    raise ValueError("image_paths and V_center_values must have the same length.")

d1_values = []
d2_values = []
d2L_values = []
d2R_values = []

for tif_path, Vc in zip(image_paths, V_center_values):
    avg_img = load_avg_image(tif_path)

    centers = find_ion_centers(
        avg_img,
        sigma=0.8,
        threshold_rel=0.5,
        min_separation=3,
        roi_half_width=4
    )

    dists, projected_coords, axis_unit, fit_params = nearest_neighbor_distances_projected(centers)

    if len(centers) == 0:
        raise RuntimeError(f"No ions found in image:\n{tif_path}")

    print(f"\nFile: {tif_path}")
    print(f"V_center = {Vc}")
    print(f"Found {len(centers)} ions")
    print("Projected nearest-neighbor spacings (px):")
    print(np.array2string(dists, precision=3))

    d1, d2, d2L, d2R, info = extract_d1_d2(dists)

    d1_values.append(d1)
    d2_values.append(d2)
    d2L_values.append(d2L)
    d2R_values.append(d2R)

    print(f"d1 = {d1:.3f} px")
    print(f"d2 = {d2:.3f} px")
    print(f"Detected {info['n_ions']} ions")
    print(f"d1 used spacing indices {info['d1_indices']} ({info['definition']})")
    print(f"d2 used spacing indices {info['d2_indices']} (edge spacings)")

d1_values = np.array(d1_values)
d2_values = np.array(d2_values)
d2L_values = np.array(d2L_values)
d2R_values = np.array(d2R_values)

# --------------------------------------------------
# Best-fit lines
# --------------------------------------------------
m_d1, b_d1 = np.polyfit(V_center_values, d1_values, 1)
m_d2, b_d2 = np.polyfit(V_center_values, d2_values, 1)
m_d2L, b_d2L = np.polyfit(V_center_values, d2L_values, 1)
m_d2R, b_d2R = np.polyfit(V_center_values, d2R_values, 1)

xfit = np.linspace(np.min(V_center_values), np.max(V_center_values), 300)

yfit_d1  = m_d1  * xfit + b_d1
yfit_d2  = m_d2  * xfit + b_d2
yfit_d2L = m_d2L * xfit + b_d2L
yfit_d2R = m_d2R * xfit + b_d2R

print("\nBest-fit line slopes:")
print(f"d1 slope  = {m_d1:.6f} px / V")
print(f"d2 slope  = {m_d2:.6f} px / V")
print(f"d2L slope = {m_d2L:.6f} px / V")
print(f"d2R slope = {m_d2R:.6f} px / V")


# --------------------------------------------------
# Plot 1: d1 and d2
# --------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(V_center_values, d1_values, 'o-', label='d1 (center-adjacent spacing)')
plt.plot(V_center_values, d2_values, 's-', label='d2 (outermost spacing)')
plt.xlabel("V_center")
plt.ylabel("Spacing (pixels)")
plt.title("d1 and d2 vs V_center")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


# --------------------------------------------------
# Plot 2: d1, d2L, d2R
# --------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(V_center_values, d1_values, 'o-', label='d1 (center-adjacent spacing)')
plt.plot(V_center_values, d2L_values, 'o-', label='d2L (left outer spacing)')
plt.plot(V_center_values, d2R_values, 'o-', label='d2R (right outer spacing)')
plt.xlabel("V_center")
plt.ylabel("Spacing (pixels)")
plt.title("d1 and d2 vs V_center")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


# --------------------------------------------------
# Plot 3: all data + best-fit lines
# --------------------------------------------------
plt.figure(figsize=(7, 5))

# raw data
plt.plot(V_center_values, d1_values, 'o', label='d1 data')
plt.plot(V_center_values, d2_values, 'o', label='d2 data')
plt.plot(V_center_values, d2L_values, 'o', label='d2L data')
plt.plot(V_center_values, d2R_values, 'o', label='d2R data')

# best-fit lines
plt.plot(xfit, yfit_d1,  '-', label=f'd1 fit (slope={m_d1:.3f})')
plt.plot(xfit, yfit_d2,  '-', label=f'd2 fit (slope={m_d2:.3f})')
plt.plot(xfit, yfit_d2L, '-', label=f'd2L fit (slope={m_d2L:.3f})')
plt.plot(xfit, yfit_d2R, '-', label=f'd2R fit (slope={m_d2R:.3f})')

plt.xlabel("V_center")
plt.ylabel("Spacing (pixels)")
plt.title("d1 and d2 vs V_center")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()