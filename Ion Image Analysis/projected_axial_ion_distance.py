# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 12:21:12 2026

@author: iontrap
"""

import numpy as np
import matplotlib.pyplot as plt
import tifffile
from scipy.ndimage import gaussian_filter, maximum_filter, label, center_of_mass


# --------------------------------------------------
# Find ion centers
# --------------------------------------------------
def find_ion_centers(img, sigma=0.65, threshold_rel=0.5, min_separation=3, roi_half_width=4):
    """
    Find ion centers using local maxima + local intensity-weighted centroiding.

    Parameters
    ----------
    img : 2D ndarray
        Image
    sigma : float
        Gaussian smoothing width
    threshold_rel : float
        Relative threshold as fraction of max(smoothed image)
    min_separation : int
        Minimum allowed peak separation in pixels
    roi_half_width : int
        Half-width of local box used for centroiding

    Returns
    -------
    centers : ndarray, shape (N, 2)
        Ion centers in (y, x) format
    """
    smoothed = gaussian_filter(img, sigma=sigma)
    threshold = threshold_rel * np.max(smoothed)

    # Candidate local maxima
    mf = maximum_filter(smoothed, size=min_separation)
    peak_mask = (smoothed == mf) & (smoothed > threshold)

    peak_y, peak_x = np.where(peak_mask)

    if len(peak_x) == 0:
        return np.empty((0, 2))

    # Sort peaks by brightness, brightest first
    peak_vals = smoothed[peak_y, peak_x]
    order = np.argsort(peak_vals)[::-1]
    peak_y = peak_y[order]
    peak_x = peak_x[order]

    # Non-maximum suppression by distance
    accepted = []
    for y0, x0 in zip(peak_y, peak_x):
        keep = True
        for ya, xa in accepted:
            if np.hypot(x0 - xa, y0 - ya) < min_separation:
                keep = False
                break
        if keep:
            accepted.append((y0, x0))

    centers = []

    for y0, x0 in accepted:
        y_min = max(0, y0 - roi_half_width)
        y_max = min(img.shape[0], y0 + roi_half_width + 1)
        x_min = max(0, x0 - roi_half_width)
        x_max = min(img.shape[1], x0 + roi_half_width + 1)

        roi = img[y_min:y_max, x_min:x_max].astype(float)

        # Optional local background subtraction
        bg = np.median(roi)
        roi_bgsub = roi - bg
        roi_bgsub[roi_bgsub < 0] = 0

        if roi_bgsub.sum() <= 0:
            continue

        yy, xx = np.indices(roi_bgsub.shape)
        yc = (yy * roi_bgsub).sum() / roi_bgsub.sum() + y_min
        xc = (xx * roi_bgsub).sum() / roi_bgsub.sum() + x_min

        centers.append((yc, xc))

    centers = np.array(centers)

    if len(centers) == 0:
        return np.empty((0, 2))

    # Sort left to right
    centers = centers[np.argsort(centers[:, 1])]

    return centers


# --------------------------------------------------
# Fit chain axis and project ion centers onto it
# --------------------------------------------------
def project_onto_chain_axis(centers):
    """
    Fit a straight line to the ion centers and project each center
    onto that best-fit axis.

    Parameters
    ----------
    centers : ndarray of shape (N, 2)
        Ion centers in (y, x) format

    Returns
    -------
    projected_coords : ndarray
        Scalar coordinate of each ion along the fitted chain axis
    axis_unit : ndarray of shape (2,)
        Unit vector along the fitted chain axis in (x, y) form
    fit_params : tuple
        (slope, intercept) for y = slope*x + intercept
    """
    if len(centers) < 2:
        return np.array([]), np.array([1.0, 0.0]), (0.0, 0.0)

    y = centers[:, 0]
    x = centers[:, 1]

    # Fit line y = m*x + b
    m, b = np.polyfit(x, y, 1)

    # Unit vector along fitted line, in (x, y) coordinates
    axis_vec = np.array([1.0, m])
    axis_unit = axis_vec / np.linalg.norm(axis_vec)

    # Put coordinates into (x, y) form for projection
    coords_xy = np.column_stack((x, y))

    # Project each point onto the axis
    projected_coords = coords_xy @ axis_unit

    # Sort projected coordinates from left to right
    projected_coords = np.sort(projected_coords)

    return projected_coords, axis_unit, (m, b)


# --------------------------------------------------
# Compute nearest-neighbor distances along chain axis
# --------------------------------------------------
def nearest_neighbor_distances_projected(centers):
    if len(centers) < 2:
        return np.array([]), np.array([]), np.array([1.0, 0.0]), (0.0, 0.0)

    projected_coords, axis_unit, fit_params = project_onto_chain_axis(centers)

    dists = np.diff(projected_coords)

    return dists, projected_coords, axis_unit, fit_params


# --------------------------------------------------
# Compute spacing statistics
# --------------------------------------------------
def spacing_stats(dists, ddof=0):
    """
    Returns mean spacing, standard deviation, and std/mean ratio.

    Parameters
    ----------
    dists : array-like
        Nearest-neighbor distances
    ddof : int
        0 for population std, 1 for sample std
    """
    if len(dists) == 0:
        return np.nan, np.nan, np.nan

    mean_spacing = np.mean(dists)
    std_spacing = np.std(dists, ddof=ddof)

    if mean_spacing == 0:
        ratio = np.nan
    else:
        ratio = std_spacing / mean_spacing

    return mean_spacing, std_spacing, ratio


# --------------------------------------------------
# Load image
# --------------------------------------------------
tif_path = r"Z:\Lab Data\EMCCD\Images\raw_img_Thu Apr  2 2026_14.29.06_1506.tif"
data = tifffile.imread(tif_path)

# Average if stack
avg_img = np.mean(data, axis=0) if data.ndim == 3 else data


# --------------------------------------------------
# Find ions + projected distances
# --------------------------------------------------
centers = find_ion_centers(avg_img)

dists, projected_coords, axis_unit, fit_params = nearest_neighbor_distances_projected(centers)

mean_spacing, std_spacing, spacing_ratio = spacing_stats(dists, ddof=0)

# Inner-only spacing stats (exclude spacings involving outermost ions)
if len(dists) >= 3:
    inner_dists = dists[1:-1]
else:
    inner_dists = np.array([])

mean_spacing_inner, std_spacing_inner, spacing_ratio_inner = spacing_stats(inner_dists, ddof=0)

m, b = fit_params


# --------------------------------------------------
# Print results
# --------------------------------------------------
print("\nIon positions (x, y):")
for i, (y, x) in enumerate(centers):
    print(f"Ion {i}: x={x:.2f}, y={y:.2f}")

print("\nProjected ion coordinates along fitted chain axis:")
for i, s in enumerate(projected_coords):
    print(f"Ion {i}: s = {s:.3f} pixels")

print("\nNearest-neighbor distances along chain axis (pixels):")
for i, d in enumerate(dists):
    print(f"{i} → {i+1}: {d:.3f} pixels")

print("\nBest-fit chain axis:")
print(f"slope m = {m:.6f}")
print(f"intercept b = {b:.6f}")
print(f"axis unit vector (x, y) = ({axis_unit[0]:.6f}, {axis_unit[1]:.6f})")

print("\nSpacing statistics (all ions):")
print(f"Mean spacing             = {mean_spacing:.6f} pixels")
print(f"Std. dev. of spacing     = {std_spacing:.6f} pixels")
print(f"Uniformity               = {(1 - spacing_ratio)*100:.6f} %")

print("\nSpacing statistics (excluding outermost ions):")
print(f"Mean spacing             = {mean_spacing_inner:.6f} pixels")
print(f"Std. dev. of spacing     = {std_spacing_inner:.6f} pixels")
print(f"Uniformity               = {(1 - spacing_ratio_inner)*100:.6f} %")


# --------------------------------------------------
# Plot result: all ions included
# --------------------------------------------------
plt.figure(figsize=(8, 4))
plt.imshow(avg_img, cmap="gray", origin="lower")

if len(centers) > 0:
    x = centers[:, 1]
    y = centers[:, 0]

    plt.scatter(x, y, color="cyan", s=20, label="Ion centers")

    for i, (xi, yi) in enumerate(zip(x, y)):
        plt.text(xi + 1, yi + 1, str(i), color="red", fontsize=10)

    # Plot best-fit chain axis
    x_line = np.linspace(np.min(x) - 5, np.max(x) + 5, 200)
    y_line = m * x_line + b
    plt.plot(x_line, y_line, 'y--', lw=1.5, label="Best-fit axis")

stats_text = (
    f"Mean spacing = {mean_spacing:.3f} px\n"
    f"Std spacing  = {std_spacing:.3f} px\n"
    f"Uniformity   = {(1 - spacing_ratio)*100:.5f} %"
)

plt.text(
    -0.02, 1.05, stats_text,
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)

plt.title("Ion Spacing Along Best-Fit Chain Axis (All Ions)")
plt.legend(loc="lower right")
plt.tight_layout()
plt.show()


# --------------------------------------------------
# Plot result: excluding outermost ions in uniformity
# --------------------------------------------------
plt.figure(figsize=(8, 4))
plt.imshow(avg_img, cmap="gray", origin="lower")

if len(centers) > 0:
    x = centers[:, 1]
    y = centers[:, 0]

    plt.scatter(x, y, color="cyan", s=20, label="Ion centers")

    for i, (xi, yi) in enumerate(zip(x, y)):
        plt.text(xi + 1, yi + 1, str(i), color="red", fontsize=10)

    # Plot best-fit chain axis
    x_line = np.linspace(np.min(x) - 5, np.max(x) + 5, 200)
    y_line = m * x_line + b
    plt.plot(x_line, y_line, 'y--', lw=1.5, label="Best-fit axis")

stats_text_inner = (
    f"Mean spacing = {mean_spacing_inner:.3f} px\n"
    f"Std spacing  = {std_spacing_inner:.3f} px\n"
    f"Uniformity   = {(1 - spacing_ratio_inner)*100:.5f} %\n"
    f"(outer spacings excluded)"
)

plt.text(
    0.02, 0.98, stats_text_inner,
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)

plt.title("Ion Spacing Along Best-Fit Chain Axis (Inner Spacings Only)")
plt.legend(loc="lower right")
plt.tight_layout()
plt.show()


# --------------------------------------------------
# Optional: plot spacing vs ion-pair index
# --------------------------------------------------
if len(dists) > 0:
    plt.figure(figsize=(6, 4))
    plt.plot(np.arange(len(dists)), dists, 'o-')
    plt.xlabel("Ion pair index")
    plt.ylabel("Projected spacing (pixels)")
    plt.title("Nearest-Neighbor Spacing Along Chain Axis")
    plt.tight_layout()
    plt.show()