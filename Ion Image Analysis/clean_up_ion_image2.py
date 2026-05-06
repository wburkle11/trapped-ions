# -*- coding: utf-8 -*-
"""
Ion crystal cleanup with automatic ion-focused ROI.

Pipeline:
1. Load full raw image
2. Automatically find compact ion-like peaks
3. Crop ROI around ion crystal only
4. Run tilt_image on ROI image
5. Build averaged PSF
6. Run Richardson-Lucy deconvolution
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.optimize import curve_fit
from scipy.signal import fftconvolve
from scipy.ndimage import gaussian_filter, maximum_filter
from imageio.v2 import imread, imwrite
from pathlib import Path

from tilt_image import tilt_image


# ============================================================
# USER SETTINGS
# ============================================================

file_path = r"Z:\Lab Data\EMCCD\Images\raw_img_Tue Apr 28 2026_14.09.34_1462.tif"

# ----------------------------
# Auto-ROI settings
# ----------------------------
USE_AUTO_ROI = True
EXPECTED_N_IONS = 7

AUTO_ROI_PADDING_X = 80
AUTO_ROI_PADDING_Y = 55

BANDPASS_SMALL_SIGMA = 1.2
BANDPASS_LARGE_SIGMA = 12.0

PEAK_FOOTPRINT = 9
PEAK_PERCENTILE = 97.5
MAX_PEAKS_TO_TEST = 120
CLUSTER_RADIUS = 300

# ----------------------------
# Tilt-image settings
# ----------------------------
THRESHOLD_FRAC = 0.35
MIN_SEPARATION = 5
FIT_HALF_WIDTH_TILT = 4
HX_DEFAULT = 7
HY_DEFAULT = 7
SAVE_ROTATED = True

# ----------------------------
# Display settings
# ----------------------------
DISPLAY_GAMMA = 1.0

# Plotting-only crop for the final rotated-vs-cleaned figure.
# This does NOT affect saved data.
USE_CENTER_CROP_DISPLAY = True
CENTER_CROP_FRACTION = 0.65     # smaller = more zoomed in; try 0.55-0.85
CENTER_CROP_SIZE = None         # set to an integer pixel size to override fraction

# ----------------------------
# Gaussian-fit settings
# ----------------------------
FIT_SIGMA_INIT = 2.0
FIT_BG_MODE = "median"

# ----------------------------
# RL deconvolution settings
# ----------------------------
RL_ITERATIONS = 13
EPS = 1e-12

# ----------------------------
# Averaged-PSF settings
# ----------------------------
USE_PSF_BACKGROUND_SUBTRACTION = True
PSF_SELECTION_MODE = "all"       # "all", "middle", or "manual"
MIDDLE_SKIP = 2
PSF_SELECTION_IONS = [1, 2, 3]   # 1-based ion numbers

# ----------------------------
# PSF crosstalk suppression
# ----------------------------
USE_PSF_WEIGHT_MASK = True
PSF_MASK_SIGMA_X_SCALE = 0.8
PSF_MASK_SIGMA_Y_SCALE = 1.0


# ============================================================
# HELPERS
# ============================================================

def display_image(img, gamma=1.0):
    disp = img.astype(float).copy()
    disp -= disp.min()
    if disp.max() > 0:
        disp /= disp.max()
    disp = disp ** gamma
    return disp


def center_crop_display(img, crop_size=None, crop_fraction=0.65):
    """
    Plotting-only center crop.

    This does NOT modify the image data or the saved RL output.
    It only returns a centered crop for nicer display.

    Parameters
    ----------
    img : 2D array
        Image to crop.
    crop_size : int or None
        If given, use this square crop size in pixels.
    crop_fraction : float
        If crop_size is None, crop to this fraction of the smaller image dimension.
        Example: 0.65 keeps the central 65% of the image.
    """
    ny, nx = img.shape

    if crop_size is None:
        crop_size = int(crop_fraction * min(ny, nx))

    crop_size = min(crop_size, ny, nx)

    cx = nx // 2
    cy = ny // 2

    x1 = max(0, cx - crop_size // 2)
    x2 = min(nx, x1 + crop_size)
    y1 = max(0, cy - crop_size // 2)
    y2 = min(ny, y1 + crop_size)

    # Adjust if we hit an edge
    x1 = max(0, x2 - crop_size)
    y1 = max(0, y2 - crop_size)

    return img[y1:y2, x1:x2]


def load_image_average(path):
    img = imread(path).astype(float)

    if img.ndim == 3:
        return np.mean(img, axis=0)

    return img.copy()


def gaussian_2d(coords, A, x0, y0, sx, sy, B):
    x, y = coords
    model = A * np.exp(
        -((x - x0) ** 2 / (2 * sx ** 2) +
          (y - y0) ** 2 / (2 * sy ** 2))
    ) + B
    return model.ravel()


def fit_local_gaussian(image, x_guess, y_guess, hx, hy):
    ny, nx = image.shape

    x1 = max(0, int(round(x_guess)) - hx)
    x2 = min(nx, int(round(x_guess)) + hx + 1)
    y1 = max(0, int(round(y_guess)) - hy)
    y2 = min(ny, int(round(y_guess)) + hy + 1)

    patch = image[y1:y2, x1:x2].copy()
    yy, xx = np.mgrid[y1:y2, x1:x2]

    if FIT_BG_MODE == "median":
        B0 = np.median(patch)
    else:
        B0 = np.min(patch)

    A0 = np.max(patch) - B0
    if A0 <= 0:
        A0 = 1.0

    p0 = [A0, x_guess, y_guess, FIT_SIGMA_INIT, FIT_SIGMA_INIT, B0]

    lower_bounds = [0, x1 - 0.5, y1 - 0.5, 0.5, 0.5, -np.inf]
    upper_bounds = [np.inf, x2 - 0.5, y2 - 0.5, 10.0, 10.0, np.inf]

    popt, _ = curve_fit(
        gaussian_2d,
        (xx, yy),
        patch.ravel(),
        p0=p0,
        bounds=(lower_bounds, upper_bounds),
        maxfev=10000
    )

    A_fit, x0_fit, y0_fit, sx_fit, sy_fit, B_fit = popt

    return {
        "patch": patch,
        "bounds": (x1, x2, y1, y2),
        "fit": {
            "A": A_fit,
            "x0": x0_fit,
            "y0": y0_fit,
            "sx": sx_fit,
            "sy": sy_fit,
            "B": B_fit
        },
        "local_center": (x0_fit - x1, y0_fit - y1)
    }


def extract_patch(image, x0, y0, hx, hy):
    ny, nx = image.shape

    x0 = int(round(x0))
    y0 = int(round(y0))

    x1 = max(0, x0 - hx)
    x2 = min(nx, x0 + hx + 1)
    y1 = max(0, y0 - hy)
    y2 = min(ny, y0 + hy + 1)

    patch = image[y1:y2, x1:x2].copy()
    return patch, (x1, x2, y1, y2)


def richardson_lucy(image, psf, iterations=RL_ITERATIONS, epsilon=EPS):
    image = image.astype(float)
    psf = psf.astype(float)

    psf_sum = np.sum(psf)
    if psf_sum <= 0:
        raise ValueError("PSF sum must be positive.")

    psf = psf / psf_sum

    estimate = np.full_like(image, 0.5)
    psf_mirror = psf[::-1, ::-1]

    for _ in range(iterations):
        blurred = fftconvolve(estimate, psf, mode="same")
        ratio = image / (blurred + epsilon)
        correction = fftconvolve(ratio, psf_mirror, mode="same")
        estimate *= correction
        estimate[estimate < 0] = 0

    return estimate


def make_auto_roi_image(
    file_path,
    expected_n_ions=EXPECTED_N_IONS,
    padding_x=AUTO_ROI_PADDING_X,
    padding_y=AUTO_ROI_PADDING_Y,
    save_roi=True
):
    raw = load_image_average(file_path)
    ny, nx = raw.shape

    work = raw.astype(float).copy()
    bg = np.percentile(work, 20)
    work -= bg
    work[work < 0] = 0

    if work.max() <= 0:
        raise RuntimeError("Auto ROI failed: no positive signal after background subtraction.")

    # Bandpass filtering:
    # small blur keeps ion-sized spots
    # large blur estimates broad scatter/background
    small = gaussian_filter(work, sigma=BANDPASS_SMALL_SIGMA)
    large = gaussian_filter(work, sigma=BANDPASS_LARGE_SIGMA)

    band = small - large
    band[band < 0] = 0

    if band.max() <= 0:
        raise RuntimeError("Auto ROI failed: bandpass image has no signal.")

    # Local maxima
    local_max = band == maximum_filter(band, size=PEAK_FOOTPRINT)

    positive_band = band[band > 0]
    peak_threshold = np.percentile(positive_band, PEAK_PERCENTILE)

    peak_mask = local_max & (band > peak_threshold)

    ys, xs = np.where(peak_mask)

    if len(xs) == 0:
        raise RuntimeError(
            "Auto ROI failed: no peaks found. "
            "Try lowering PEAK_PERCENTILE."
        )

    peak_vals = band[ys, xs]

    # Keep strongest candidate peaks
    order = np.argsort(peak_vals)[::-1]
    order = order[:MAX_PEAKS_TO_TEST]

    xs = xs[order]
    ys = ys[order]
    peak_vals = peak_vals[order]

    candidates = np.column_stack([xs, ys])

    # Pick densest compact bright cluster
    best_score = -np.inf
    best_cluster = None

    for i, p in enumerate(candidates):
        dx = candidates[:, 0] - p[0]
        dy = candidates[:, 1] - p[1]
        dist = np.sqrt(dx**2 + dy**2)

        cluster_idx = np.where(dist < CLUSTER_RADIUS)[0]

        if expected_n_ions is not None and len(cluster_idx) >= expected_n_ions:
            sub_vals = peak_vals[cluster_idx]
            sub_order = np.argsort(sub_vals)[::-1][:expected_n_ions]
            cluster_idx = cluster_idx[sub_order]

        cluster = candidates[cluster_idx]

        if len(cluster) < 1:
            continue

        x_span = cluster[:, 0].max() - cluster[:, 0].min()
        y_span = cluster[:, 1].max() - cluster[:, 1].min()
        total_brightness = np.sum(peak_vals[cluster_idx])

        compactness_penalty = 1 + 0.02 * x_span + 0.05 * y_span
        score = total_brightness * len(cluster_idx) / compactness_penalty

        if score > best_score:
            best_score = score
            best_cluster = cluster

    if best_cluster is None:
        raise RuntimeError("Auto ROI failed: no valid ion-like cluster found.")

    x_min = max(0, int(best_cluster[:, 0].min()) - padding_x)
    x_max = min(nx, int(best_cluster[:, 0].max()) + padding_x + 1)

    y_min = max(0, int(best_cluster[:, 1].min()) - padding_y)
    y_max = min(ny, int(best_cluster[:, 1].max()) + padding_y + 1)

    roi = raw[y_min:y_max, x_min:x_max].copy()

    if save_roi:
        in_path = Path(file_path)
        roi_path = in_path.with_name(in_path.stem + "_ionCrystal_autoROI.tif")
        imwrite(roi_path, roi.astype(np.float32))
    else:
        roi_path = None

    print("\n================ ION-CRYSTAL AUTO ROI ================")
    print(f"Original image shape: {raw.shape}")
    print("Chosen ion-like peak positions:")
    for x, y in best_cluster:
        print(f"    x={x}, y={y}")
    print(f"ROI bounds: x=[{x_min}:{x_max}], y=[{y_min}:{y_max}]")
    print(f"ROI shape: {roi.shape}")
    if roi_path is not None:
        print(f"Saved ROI image to: {roi_path}")

    # Diagnostic plot
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), dpi=150)

    axes[0].imshow(display_image(raw, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
    axes[0].add_patch(
        Rectangle(
            (x_min - 0.5, y_min - 0.5),
            x_max - x_min,
            y_max - y_min,
            edgecolor="cyan",
            facecolor="none",
            linewidth=1.5
        )
    )
    axes[0].scatter(best_cluster[:, 0], best_cluster[:, 1], c="cyan", marker="x")
    axes[0].set_title("Full Raw Image with Ion-Crystal Auto ROI")
    axes[0].axis("off")

    axes[1].imshow(display_image(band, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
    axes[1].scatter(best_cluster[:, 0], best_cluster[:, 1], c="cyan", marker="x")
    axes[1].set_title("Bandpassed Image Used for Ion Detection")
    axes[1].axis("off")

    axes[2].imshow(display_image(roi, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
    axes[2].set_title("Final Auto-Cropped Ion ROI")
    axes[2].axis("off")

    plt.tight_layout()
    plt.show()

    return roi_path, roi, (x_min, x_max, y_min, y_max), raw


# ============================================================
# AUTO ROI FIRST
# ============================================================

if USE_AUTO_ROI:
    roi_path, roi_image, roi_bounds, raw_full_image = make_auto_roi_image(
        file_path=file_path,
        expected_n_ions=EXPECTED_N_IONS,
        padding_x=AUTO_ROI_PADDING_X,
        padding_y=AUTO_ROI_PADDING_Y,
        save_roi=True
    )

    analysis_file_path = str(roi_path)

else:
    analysis_file_path = file_path
    roi_bounds = None


# ============================================================
# TILT PREPROCESSING ON ROI IMAGE
# ============================================================

tilt_result = tilt_image(
    file_path=analysis_file_path,
    threshold_frac=THRESHOLD_FRAC,
    min_separation=MIN_SEPARATION,
    fit_half_width=FIT_HALF_WIDTH_TILT,
    hx_default=HX_DEFAULT,
    hy_default=HY_DEFAULT,
    save_rotated=SAVE_ROTATED
)

raw_image = tilt_result["raw_image"]
rotated_image = tilt_result["rotated_image"]
psf_settings = tilt_result["psf_settings"]

print(f"\nRotation angle: {tilt_result['angle_deg']:.6f} degrees")
if tilt_result["rotated_path"] is not None:
    print(f"Saved rotated image to: {tilt_result['rotated_path']}")

print("\nPSF settings used:")
print("psf_settings = [")
for s in psf_settings:
    print(f'    {{"x0": {s["x0"]}, "y0": {s["y0"]}, "hx": {s["hx"]}, "hy": {s["hy"]}}},')
print("]")


# ============================================================
# FIGURE 1: ROI RAW VS ROTATED
# ============================================================

fig1, axes1 = plt.subplots(2, 1, figsize=(16, 10), dpi=150)

axes1[0].imshow(display_image(raw_image, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
axes1[0].set_title("Auto-ROI Raw Image", fontsize=18)
axes1[0].axis("off")

axes1[1].imshow(display_image(rotated_image, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
axes1[1].set_title("Auto-ROI Rotated Image", fontsize=18)
axes1[1].axis("off")

plt.tight_layout()
plt.show()


# ============================================================
# FIGURE 2: ROTATED IMAGE WITH PSF BOXES
# ============================================================

fig2, ax2 = plt.subplots(figsize=(16, 5), dpi=150)
ax2.imshow(display_image(rotated_image, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
ax2.set_title("Rotated ROI Image with PSF Crops", fontsize=18)
ax2.axis("off")

for i, s in enumerate(psf_settings, start=1):
    x0, y0, hx, hy = s["x0"], s["y0"], s["hx"], s["hy"]

    rect = Rectangle(
        (x0 - hx - 0.5, y0 - hy - 0.5),
        2 * hx + 1,
        2 * hy + 1,
        linewidth=1.2,
        edgecolor="cyan",
        facecolor="none"
    )
    ax2.add_patch(rect)
    ax2.scatter([x0], [y0], c="cyan", marker="x", s=35)
    ax2.text(
        x0,
        y0 - hy - 1,
        str(i),
        color="white",
        fontsize=8,
        ha="center",
        va="bottom"
    )

plt.tight_layout()
plt.show()


# ============================================================
# FIGURE 3: CROPS WITH GAUSSIAN-FIT CENTERS
# ============================================================

fit_results = []

n_ions = len(psf_settings)
ncols = 4
nrows = int(np.ceil(n_ions / ncols))

fig3, axes3 = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), dpi=150)
axes3 = np.atleast_2d(axes3)

for i in range(nrows * ncols):
    r = i // ncols
    c = i % ncols
    ax = axes3[r, c]

    if i >= n_ions:
        ax.axis("off")
        continue

    s = psf_settings[i]
    x0, y0, hx, hy = s["x0"], s["y0"], s["hx"], s["hy"]

    result = fit_local_gaussian(rotated_image, x0, y0, hx, hy)
    fit_results.append(result)

    patch = result["patch"]
    x_fit_local, y_fit_local = result["local_center"]

    ax.imshow(display_image(patch, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
    ax.scatter([hx], [hy], c="cyan", marker="x", s=35)
    ax.scatter([x_fit_local], [y_fit_local], c="lime", marker="+", s=50)

    ax.set_title(
        f"Ion {i+1}\n"
        f"manual=({x0},{y0})\n"
        f"fit=({result['fit']['x0']:.2f},{result['fit']['y0']:.2f})",
        fontsize=8
    )
    ax.axis("off")

plt.tight_layout()
plt.show()


# ============================================================
# BUILD GLOBAL AVERAGED PSF
# ============================================================

psf_crops = []
used_indices = []

for i, s in enumerate(psf_settings):
    x0, y0, hx, hy = s["x0"], s["y0"], s["hx"], s["hy"]
    patch, _ = extract_patch(rotated_image, x0, y0, hx, hy)

    psf = patch.copy()

    if USE_PSF_BACKGROUND_SUBTRACTION:
        psf = psf - np.min(psf)
        psf[psf < 0] = 0

    if USE_PSF_WEIGHT_MASK:
        yy, xx = np.indices(psf.shape)

        cx = (psf.shape[1] - 1) / 2
        cy = (psf.shape[0] - 1) / 2

        mask_sigma_x = PSF_MASK_SIGMA_X_SCALE * hx
        mask_sigma_y = PSF_MASK_SIGMA_Y_SCALE * hy

        weight_mask = np.exp(
            -((xx - cx) ** 2 / (2 * mask_sigma_x ** 2) +
              (yy - cy) ** 2 / (2 * mask_sigma_y ** 2))
        )

        psf *= weight_mask

    if np.sum(psf) <= 0:
        continue

    psf = psf / np.sum(psf)
    psf_crops.append(psf)
    used_indices.append(i)

if len(psf_crops) == 0:
    raise RuntimeError("No valid PSF crops found.")

if PSF_SELECTION_MODE == "middle" and len(psf_crops) > 2 * MIDDLE_SKIP:
    psf_crops = psf_crops[MIDDLE_SKIP:-MIDDLE_SKIP]
    used_indices = used_indices[MIDDLE_SKIP:-MIDDLE_SKIP]

elif PSF_SELECTION_MODE == "manual":
    selected_zero_based = [i - 1 for i in PSF_SELECTION_IONS]

    filtered_psfs = []
    filtered_indices = []

    for psf, idx in zip(psf_crops, used_indices):
        if idx in selected_zero_based:
            filtered_psfs.append(psf)
            filtered_indices.append(idx)

    psf_crops = filtered_psfs
    used_indices = filtered_indices

if len(psf_crops) == 0:
    raise RuntimeError("PSF selection removed all PSFs. Check PSF_SELECTION_MODE.")

avg_psf = np.mean(np.stack(psf_crops, axis=0), axis=0)
avg_psf = avg_psf / np.sum(avg_psf)

print(f"\nUsing {len(psf_crops)} ions to build averaged PSF.")
print("Ion indices used for averaged PSF:", [i + 1 for i in used_indices])


# ============================================================
# GLOBAL RL WITH AVERAGED PSF
# ============================================================

global_rl = richardson_lucy(
    rotated_image,
    avg_psf,
    iterations=RL_ITERATIONS,
    epsilon=EPS
)

input_path = Path(file_path)

if USE_AUTO_ROI:
    rl_path = input_path.with_name(input_path.stem + "_autoROI_rotated_rl.tif")
else:
    rl_path = input_path.with_name(input_path.stem + "_rl.tif")

imwrite(rl_path, global_rl.astype(np.float32))
print(f"Saved RL image to: {rl_path}")


# ============================================================
# FIGURE 4: ROTATED VS GLOBAL RL
# ============================================================

fig4, axes4 = plt.subplots(2, 1, figsize=(16, 10), dpi=150)

if USE_CENTER_CROP_DISPLAY:
    rotated_for_display = center_crop_display(
        rotated_image,
        crop_size=CENTER_CROP_SIZE,
        crop_fraction=CENTER_CROP_FRACTION
    )
    rl_for_display = center_crop_display(
        global_rl,
        crop_size=CENTER_CROP_SIZE,
        crop_fraction=CENTER_CROP_FRACTION
    )
else:
    rotated_for_display = rotated_image
    rl_for_display = global_rl

axes4[0].imshow(display_image(rotated_for_display, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
axes4[0].set_title("rotated raw image with ROI, centered display crop", fontsize=18)
axes4[0].axis("off")

axes4[1].imshow(display_image(rl_for_display, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
axes4[1].set_title("cleaned image, centered display crop", fontsize=18)
axes4[1].axis("off")

plt.tight_layout()
plt.show()


# ============================================================
# FIGURE 5: INDIVIDUAL PSFS + AVERAGED PSF
# ============================================================

n_psfs = len(psf_crops)

fig5, axes5 = plt.subplots(
    1,
    n_psfs + 1,
    figsize=(2.4 * (n_psfs + 1), 3),
    dpi=150
)

if n_psfs == 1:
    axes5 = np.array([axes5])

for i, psf in enumerate(psf_crops):
    axes5[i].imshow(display_image(psf, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
    axes5[i].set_title(f"Ion {used_indices[i] + 1}", fontsize=9)
    axes5[i].axis("off")

axes5[-1].imshow(display_image(avg_psf, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
axes5[-1].set_title("Averaged PSF", fontsize=9)
axes5[-1].axis("off")

plt.tight_layout()
plt.show()