import numpy as np
import matplotlib.pyplot as plt
import tifffile
from scipy.ndimage import gaussian_filter, maximum_filter, label, center_of_mass


# --------------------------------------------------
# Find ion centers
# --------------------------------------------------
def find_ion_centers(img, sigma=1.0, threshold_rel=0.2, min_separation=4):

    smoothed = gaussian_filter(img, sigma=sigma)

    threshold = threshold_rel * np.max(smoothed)

    local_max = smoothed == maximum_filter(smoothed, size=min_separation)
    detected = local_max & (smoothed > threshold)

    labels, num = label(detected)

    if num == 0:
        return np.empty((0, 2))

    centers = np.array(center_of_mass(smoothed, labels, range(1, num + 1)))

    # sort left → right (by x)
    centers = centers[np.argsort(centers[:, 1])]

    return centers


# --------------------------------------------------
# Compute nearest-neighbor distances
# --------------------------------------------------
def nearest_neighbor_distances(centers):
    if len(centers) < 2:
        return np.array([])

    dx = np.diff(centers[:, 1])
    dy = np.diff(centers[:, 0])

    return np.sqrt(dx**2 + dy**2)


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
tif_path = r"Z:\Lab Data\EMCCD\Images\raw_img_Fri Mar 27 2026_12.00.25_1418.tif"
data = tifffile.imread(tif_path)

# Average if stack
avg_img = np.mean(data, axis=0) if data.ndim == 3 else data

# --------------------------------------------------
# Find ions + distances
# --------------------------------------------------
centers = find_ion_centers(avg_img)

dists = nearest_neighbor_distances(centers)

mean_spacing, std_spacing, spacing_ratio = spacing_stats(dists, ddof=0)

# --------------------------------------------------
# Print results
# --------------------------------------------------
print("\nIon positions (x, y):")
for i, (y, x) in enumerate(centers):
    print(f"Ion {i}: x={x:.2f}, y={y:.2f}")

print("\nNearest-neighbor distances (pixels):")
for i, d in enumerate(dists):
    print(f"{i} → {i+1}: {d:.3f} pixels")

print("\nSpacing statistics:")
print(f"Mean spacing             = {mean_spacing:.6f} pixels")
print(f"Std. dev. of spacing     = {std_spacing:.6f} pixels")
print(f"Uniformity  = {(1 - spacing_ratio)*100:.6f} %")

# --------------------------------------------------
# Plot result
# --------------------------------------------------
plt.figure(figsize=(8, 4))
plt.imshow(avg_img, cmap="gray", origin="lower")

if len(centers) > 0:
    x = centers[:, 1]
    y = centers[:, 0]
    plt.scatter(x, y, color="cyan", s=40)

    for i, (xi, yi) in enumerate(zip(x, y)):
        plt.text(xi + 1, yi + 1, str(i), color="red", fontsize=10)

# Text box with spacing stats
stats_text = (
    f"Mean spacing = {mean_spacing:.3f} px\n"
    f"Std spacing  = {std_spacing:.3f} px\n"
    f"Uniformity = {(1 - spacing_ratio)*100:.5f} %"
)

plt.text(
    0.02, 0.98, stats_text,
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)

plt.title("Ion Spacing")
plt.tight_layout()
plt.show()