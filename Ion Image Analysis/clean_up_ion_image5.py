# -*- coding: utf-8 -*-
"""
Ion crystal cleanup with MANUALLY defined ROI.

Pipeline:
1. Load full raw image
2. Crop a manually defined ROI around the ion crystal
3. Save that ROI
4. Run tilt_image on the ROI
5. Build averaged empirical PSF from the detected ions
6. Run Richardson-Lucy deconvolution

The main change from the auto-ROI version is that the ROI is set by hand
near the top of this file.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.optimize import curve_fit
from scipy.signal import fftconvolve
from scipy.ndimage import gaussian_filter
from imageio.v2 import imread, imwrite
from pathlib import Path

from tilt_image import tilt_image


# ============================================================
# USER SETTINGS
# ============================================================

file_path = r"Z:\Lab Data\EMCCD\Images\raw_img_Wed Apr 29 2026_15.07.13_1545.tif"


# ============================================================
# MANUAL ROI SETTINGS -- SET THESE FIRST
# ============================================================

USE_MANUAL_ROI = True

# Option A: define ROI by explicit pixel bounds in the full raw image.
#
# Coordinates are full-image pixel coordinates:
#   x increases left -> right
#   y increases top  -> bottom
#
# Crop used by the code:
#   raw[MANUAL_Y_MIN:MANUAL_Y_MAX, MANUAL_X_MIN:MANUAL_X_MAX]
#
# Make this large enough to include:
#   - the full ion crystal
#   - the visible coma/tails
#   - a little surrounding background
#
# But avoid including bright background scatter far away from the crystal.
MANUAL_X_MIN = 180
MANUAL_X_MAX = 400
MANUAL_Y_MIN = 220
MANUAL_Y_MAX = 300

# Option B: define ROI by center and size.
# If USE_CENTER_SIZE_ROI=True, the bounds above are ignored.
USE_CENTER_SIZE_ROI = False
MANUAL_CENTER_X = 250
MANUAL_CENTER_Y = 190
MANUAL_ROI_WIDTH = 240
MANUAL_ROI_HEIGHT = 150

SAVE_MANUAL_ROI = True


# ============================================================
# TILT-IMAGE SETTINGS
# ============================================================

THRESHOLD_FRAC = 0.35
MIN_SEPARATION = 5
FIT_HALF_WIDTH_TILT = 4
HX_DEFAULT = 5
HY_DEFAULT = 5
SAVE_ROTATED = True


# ============================================================
# DISPLAY SETTINGS
# ============================================================

DISPLAY_GAMMA = 1.0
DISPLAY_CMAP = "inferno"       # try "inferno", "magma", "cividis", "turbo"
DISPLAY_INTERPOLATION = "none" # try "none", "bilinear", or "bicubic"


# ============================================================
# GAUSSIAN-FIT SETTINGS
# ============================================================

FIT_SIGMA_INIT = 2.0
FIT_BG_MODE = "median"


# ============================================================
# RL DECONVOLUTION SETTINGS
# ============================================================

RL_ITERATIONS = 14
EPS = 1e-12


# ============================================================
# AVERAGED-PSF SETTINGS
# ============================================================

USE_PSF_BACKGROUND_SUBTRACTION = True
PSF_SELECTION_MODE = "all"       # "all", "middle", or "manual"
MIDDLE_SKIP = 2
PSF_SELECTION_IONS = [1, 2, 3]   # 1-based ion numbers


# ============================================================
# PSF CROSSTALK SUPPRESSION
# ============================================================

USE_PSF_WEIGHT_MASK = True
PSF_MASK_SIGMA_X_SCALE = 0.8
PSF_MASK_SIGMA_Y_SCALE = 1.0


# ============================================================
# HELPERS
# ============================================================

def display_image(img, gamma=1.0):
    disp = img.astype(float).copy()
    disp -= np.nanmin(disp)
    if np.nanmax(disp) > 0:
        disp /= np.nanmax(disp)
    return disp ** gamma


def load_image_average(path):
    img = imread(path).astype(float)

    if img.ndim == 3:
        print(f"Image stack detected: {img.shape}. Using average over frames.")
        return np.mean(img, axis=0)

    return img.copy()


def get_manual_roi_bounds(raw):
    """
    Convert manual ROI settings into safe image bounds.
    """
    ny, nx = raw.shape

    if USE_CENTER_SIZE_ROI:
        x_min = int(round(MANUAL_CENTER_X - MANUAL_ROI_WIDTH / 2))
        x_max = int(round(MANUAL_CENTER_X + MANUAL_ROI_WIDTH / 2))
        y_min = int(round(MANUAL_CENTER_Y - MANUAL_ROI_HEIGHT / 2))
        y_max = int(round(MANUAL_CENTER_Y + MANUAL_ROI_HEIGHT / 2))
    else:
        x_min = int(round(MANUAL_X_MIN))
        x_max = int(round(MANUAL_X_MAX))
        y_min = int(round(MANUAL_Y_MIN))
        y_max = int(round(MANUAL_Y_MAX))

    x_min = max(0, min(nx - 1, x_min))
    x_max = max(x_min + 1, min(nx, x_max))
    y_min = max(0, min(ny - 1, y_min))
    y_max = max(y_min + 1, min(ny, y_max))

    return x_min, x_max, y_min, y_max


def make_manual_roi_image(file_path, save_roi=True):
    """
    Load the full image, crop the manually defined ROI, save it,
    and show a diagnostic plot.
    """
    raw = load_image_average(file_path)
    x_min, x_max, y_min, y_max = get_manual_roi_bounds(raw)

    roi = raw[y_min:y_max, x_min:x_max].copy()

    if save_roi:
        in_path = Path(file_path)
        roi_path = in_path.with_name(in_path.stem + "_manualROI.tif")
        imwrite(roi_path, roi.astype(np.float32))
    else:
        roi_path = None

    print("\n================ MANUAL ROI ================")
    print(f"Original image shape: {raw.shape}")
    print(f"Manual ROI bounds: x=[{x_min}:{x_max}], y=[{y_min}:{y_max}]")
    print(f"Manual ROI shape: {roi.shape}")
    if roi_path is not None:
        print(f"Saved manual ROI image to: {roi_path}")

    fig, axes = plt.subplots(2, 1, figsize=(8, 9), dpi=150)

    axes[0].imshow(
        display_image(raw, DISPLAY_GAMMA),
        cmap=DISPLAY_CMAP,
        interpolation=DISPLAY_INTERPOLATION
    )
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
    axes[0].set_title("Full raw image")
    axes[0].axis("off")

    axes[1].imshow(
        display_image(roi, DISPLAY_GAMMA),
        cmap=DISPLAY_CMAP,
        interpolation=DISPLAY_INTERPOLATION
    )
    axes[1].set_title("manually chosen ROI")
    axes[1].axis("off")

    plt.tight_layout()
    plt.show()

    return roi_path, roi, (x_min, x_max, y_min, y_max), raw


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

    B0 = np.median(patch) if FIT_BG_MODE == "median" else np.min(patch)

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


# ============================================================
# MANUAL ROI FIRST
# ============================================================

if USE_MANUAL_ROI:
    roi_path, roi_image, roi_bounds, raw_full_image = make_manual_roi_image(
        file_path=file_path,
        save_roi=SAVE_MANUAL_ROI
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

axes1[0].imshow(display_image(raw_image, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
axes1[0].set_title("Manual-ROI Raw Image", fontsize=18)
axes1[0].axis("off")

axes1[1].imshow(display_image(rotated_image, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
axes1[1].set_title("Manual-ROI Rotated Image", fontsize=18)
axes1[1].axis("off")

plt.tight_layout()
plt.show()


# ============================================================
# FIGURE 2: ROTATED IMAGE WITH PSF BOXES
# ============================================================

fig2, ax2 = plt.subplots(figsize=(16, 5), dpi=150)
ax2.imshow(display_image(rotated_image, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
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

    ax.imshow(display_image(patch, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
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
rl_path = input_path.with_name(input_path.stem + "_manualROI_rotated_rl.tif")

imwrite(rl_path, global_rl.astype(np.float32))
print(f"Saved RL image to: {rl_path}")


# ============================================================
# FIGURE 4A: PIXEL COORDINATES (DENSE TICKS)
# ============================================================

fig4a, axes4a = plt.subplots(2, 1, figsize=(10, 10), dpi=150)

for ax, img, title in zip(
    axes4a,
    [rotated_image, global_rl],
    ["raw image (px)", "cleaned image (px)"]
):
    ax.imshow(
        display_image(img, DISPLAY_GAMMA),
        cmap=DISPLAY_CMAP,
        interpolation=DISPLAY_INTERPOLATION
    )

    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Pixel X")
    ax.set_ylabel("Pixel Y")

    ny, nx = img.shape

    # Dense ticks
    ax.set_xticks(np.arange(0, nx, 20))
    ax.set_yticks(np.arange(0, ny, 20))

    ax.tick_params(labelsize=8)
    ax.set_aspect("equal")

plt.tight_layout()
plt.show()

# ============================================================
# FIGURE 4B: PHYSICAL UNITS (µm)
# ============================================================

PIXEL_SIZE_UM = 0.65  # <-- your calibration

fig4b, axes4b = plt.subplots(2, 1, figsize=(10, 10), dpi=150)

for ax, img, title in zip(
    axes4b,
    [rotated_image, global_rl],
    ["raw image (um)", "cleaned image (um)"]
):
    ax.imshow(
        display_image(img, DISPLAY_GAMMA),
        cmap=DISPLAY_CMAP,
        interpolation=DISPLAY_INTERPOLATION
    )

    ny, nx = img.shape

    # Choose reasonable spacing in microns
    tick_um = 10  # every 10 µm

    tick_px = int(tick_um / PIXEL_SIZE_UM)

    xticks = np.arange(0, nx, tick_px)
    yticks = np.arange(0, ny, tick_px)

    ax.set_xticks(xticks)
    ax.set_yticks(yticks)

    ax.set_xticklabels(np.round(xticks * PIXEL_SIZE_UM, 1))
    ax.set_yticklabels(np.round(yticks * PIXEL_SIZE_UM, 1))

    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_title(title + " (physical units)", fontsize=14)

    ax.tick_params(labelsize=8)
    ax.set_aspect("equal")

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
    axes5[i].imshow(display_image(psf, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
    axes5[i].set_title(f"Ion {used_indices[i] + 1}", fontsize=9)
    axes5[i].axis("off")

axes5[-1].imshow(display_image(avg_psf, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
axes5[-1].set_title("Averaged PSF", fontsize=9)
axes5[-1].axis("off")

plt.tight_layout()
plt.show()
