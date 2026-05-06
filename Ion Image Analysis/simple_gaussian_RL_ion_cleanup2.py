# -*- coding: utf-8 -*-
"""
Simple ion image cleanup using Gaussian PSF + Richardson-Lucy deconvolution.

Idea:
1. Load TIFF, optionally average frames
2. Crop around ions automatically
3. Estimate and subtract smooth/background scatter
4. Build a simple 2D Gaussian PSF
5. Run Richardson-Lucy deconvolution
6. Suppress residual low-level background for display/saving

This is intentionally simpler than the older version:
- no PSF extraction from each ion
- no fitted ion-by-ion PSF averaging
- fewer tuning knobs
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.signal import fftconvolve
import tifffile
from imageio.v2 import imwrite
from pathlib import Path


# ============================================================
# USER SETTINGS
# ============================================================

file_path = r"Z:\Lab Data\EMCCD\Images\raw_img_Tue Apr 28 2026_14.30.29_1467.tif"

# Use "first" for a single-shot image, "average" for cleaner multi-frame stacks.
FRAME_MODE = "first"       # "first" or "average"

# Auto-crop settings
USE_AUTO_ROI = True
ROI_PADDING_X = 80
ROI_PADDING_Y = 50

# Background subtraction
BACKGROUND_PERCENTILE = 15
BACKGROUND_SMOOTH_SIGMA = 15     # larger = smoother estimated background

# Gaussian PSF model in pixels.
# Increase SIGMA_X / SIGMA_Y if the cleaned image has ringing or over-sharp artifacts.
# Decrease them if the cleaned image is still too blurry.
PSF_SIGMA_X = 2.0
PSF_SIGMA_Y = 3.0
PSF_SIZE = 25                    # odd number; try 21, 25, or 31

# Richardson-Lucy settings
RL_ITERATIONS = 4             # try 5-12; too high amplifies scatter
EPS = 1e-12

# Final display/background cleanup
FINAL_BG_PERCENTILE = 5
FINAL_SOFT_DENOISE = True
FINAL_DENOISE_SIGMA = 0.7
DISPLAY_GAMMA = 1.0

# Plotting only: makes zoomed image look less blocky.
# This does NOT change the saved TIFF data.
PLOT_INTERPOLATION = "bicubic"   # try "bilinear", "bicubic", or "none" 

SAVE_OUTPUT = True

# Final cleaned-image zoom plot
MAKE_FINAL_ZOOM_PLOT = True
ZOOM_PADDING_X = 60
ZOOM_PADDING_Y = 50
ZOOM_THRESHOLD_PERCENTILE = 99.5


# ============================================================
# BASIC HELPERS
# ============================================================

def normalize_for_display(img, gamma=1.0):
    img = img.astype(float).copy()
    img -= np.nanmin(img)
    if np.nanmax(img) > 0:
        img /= np.nanmax(img)
    return img ** gamma


def load_tiff(path, frame_mode="first"):
    data = tifffile.imread(path).astype(float)

    print("\n================ TIFF LOAD ================")
    print("Path:", path)
    print("Shape:", data.shape)

    if data.ndim == 2:
        return data

    if data.ndim == 3:
        if frame_mode.lower() == "average":
            print("Using average over all frames.")
            return np.mean(data, axis=0)
        else:
            print("Using first frame only.")
            return data[0]

    raise RuntimeError(f"Unexpected TIFF shape: {data.shape}")


def subtract_background(img, percentile=20, smooth_sigma=20):
    """
    Removes camera offset + broad background scatter.

    First subtract a constant percentile offset, then subtract a heavily
    smoothed version of the image to remove broad glow/scatter.
    """
    work = img.astype(float).copy()

    offset = np.percentile(work, percentile)
    work -= offset
    work[work < 0] = 0

    broad_bg = gaussian_filter(work, sigma=smooth_sigma)
    work = work - broad_bg
    work[work < 0] = 0

    return work, broad_bg, offset


def auto_roi_from_signal(img, padding_x=80, padding_y=50):
    """
    Simple robust auto-ROI:
    threshold the cleaned signal and crop around the bright pixels.
    """
    if img.max() <= 0:
        raise RuntimeError("No positive signal found for ROI.")

    thresh = np.percentile(img[img > 0], 98.0)
    mask = img > thresh

    ys, xs = np.where(mask)
    if len(xs) == 0:
        raise RuntimeError("Auto ROI failed. Try lowering the ROI threshold.")

    ny, nx = img.shape

    x1 = max(0, xs.min() - padding_x)
    x2 = min(nx, xs.max() + padding_x + 1)
    y1 = max(0, ys.min() - padding_y)
    y2 = min(ny, ys.max() + padding_y + 1)

    return int(x1), int(x2), int(y1), int(y2)


def make_gaussian_psf(size=21, sigma_x=2.0, sigma_y=2.0):
    """
    Make an elliptical Gaussian PSF.

    This is the cleanest starting model if the goal is to reverse a roughly
    Gaussian blur/aberration without building a noisy empirical PSF from ions.
    """
    if size % 2 == 0:
        raise ValueError("PSF_SIZE must be odd.")

    ax = np.arange(size) - size // 2
    xx, yy = np.meshgrid(ax, ax)

    psf = np.exp(
        -0.5 * ((xx / sigma_x) ** 2 + (yy / sigma_y) ** 2)
    )

    psf /= psf.sum()
    return psf


def richardson_lucy(image, psf, iterations=8, eps=1e-12):
    """
    Basic Richardson-Lucy deconvolution.

    More iterations sharpen the image but also amplify noise/scatter.
    """
    image = image.astype(float)
    image[image < 0] = 0

    psf = psf.astype(float)
    psf /= psf.sum()

    estimate = np.full_like(image, np.mean(image[image > 0]) if np.any(image > 0) else 1.0)
    psf_mirror = psf[::-1, ::-1]

    for _ in range(iterations):
        blurred = fftconvolve(estimate, psf, mode="same")
        ratio = image / (blurred + eps)
        correction = fftconvolve(ratio, psf_mirror, mode="same")
        estimate *= correction
        estimate[estimate < 0] = 0

    return estimate


def final_cleanup(img, bg_percentile=5, soft_denoise=True, denoise_sigma=0.4):
    """
    Remove residual low-level haze after RL.
    """
    out = img.astype(float).copy()

    bg = np.percentile(out[out > 0], bg_percentile) if np.any(out > 0) else 0
    out -= bg
    out[out < 0] = 0

    if soft_denoise:
        out = gaussian_filter(out, sigma=denoise_sigma)

    return out


def get_zoom_bounds(img, padding_x=30, padding_y=20, n_pixels=100):
    work = img.copy()

    # Ignore weak background
    threshold = np.percentile(work[work > 0], 95)
    work = np.where(work > threshold, work, 0)

    flat = work.ravel()
    idx = np.argpartition(flat, -n_pixels)[-n_pixels:]
    ys, xs = np.unravel_index(idx, work.shape)

    ny, nx = work.shape

    x1 = max(0, xs.min() - padding_x)
    x2 = min(nx, xs.max() + padding_x + 1)
    y1 = max(0, ys.min() - padding_y)
    y2 = min(ny, ys.max() + padding_y + 1)

    return x1, x2, y1, y2


# ============================================================
# MAIN
# ============================================================

raw = load_tiff(file_path, FRAME_MODE)

# First pass background subtraction on full image for ROI finding
clean_for_roi, broad_bg_full, offset_full = subtract_background(
    raw,
    percentile=BACKGROUND_PERCENTILE,
    smooth_sigma=BACKGROUND_SMOOTH_SIGMA
)

if USE_AUTO_ROI:
    x1, x2, y1, y2 = auto_roi_from_signal(
        clean_for_roi,
        padding_x=ROI_PADDING_X,
        padding_y=ROI_PADDING_Y
    )
    raw_roi = raw[y1:y2, x1:x2]
    print(f"Auto ROI: x=[{x1}:{x2}], y=[{y1}:{y2}], shape={raw_roi.shape}")
else:
    raw_roi = raw

# Background subtract inside the ROI
analysis_img, broad_bg_roi, offset_roi = subtract_background(
    raw_roi,
    percentile=BACKGROUND_PERCENTILE,
    smooth_sigma=BACKGROUND_SMOOTH_SIGMA
)

# Build simple Gaussian PSF
psf = make_gaussian_psf(
    size=PSF_SIZE,
    sigma_x=PSF_SIGMA_X,
    sigma_y=PSF_SIGMA_Y
)

# Deconvolve
rl = richardson_lucy(
    analysis_img,
    psf,
    iterations=RL_ITERATIONS,
    eps=EPS
)

# Final cleanup
rl_clean = final_cleanup(
    rl,
    bg_percentile=FINAL_BG_PERCENTILE,
    soft_denoise=FINAL_SOFT_DENOISE,
    denoise_sigma=FINAL_DENOISE_SIGMA
)

# Save
input_path = Path(file_path)
if SAVE_OUTPUT:
    out_path = input_path.with_name(input_path.stem + "_simple_gaussian_RL_cleaned.tif")
    imwrite(out_path, rl_clean.astype(np.float32))
    print("Saved cleaned image to:", out_path)

# ============================================================
# DIAGNOSTIC PLOTS
# ============================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=150)

axes[0].imshow(normalize_for_display(raw_roi, DISPLAY_GAMMA), cmap="inferno", interpolation=PLOT_INTERPOLATION)
axes[0].set_title("Raw ROI")
axes[0].axis("off")

axes[1].imshow(normalize_for_display(analysis_img, DISPLAY_GAMMA), cmap="inferno", interpolation=PLOT_INTERPOLATION)
axes[1].set_title("Background subtracted")
axes[1].axis("off")

axes[2].imshow(normalize_for_display(psf, 0.7), cmap="inferno", interpolation=PLOT_INTERPOLATION)
axes[2].set_title("Gaussian PSF")
axes[2].axis("off")

plt.tight_layout()
plt.show()


# ============================================================
# FINAL CLEANED IMAGE: OWN PLOT
# ============================================================

plt.figure(figsize=(9, 7), dpi=180)
plt.imshow(normalize_for_display(rl_clean, DISPLAY_GAMMA), cmap="inferno", interpolation=PLOT_INTERPOLATION)
plt.title("Final cleaned image", fontsize=16)
plt.axis("off")
plt.tight_layout()
plt.show()


# ============================================================
# FINAL CLEANED IMAGE: ZOOMED ON ION CRYSTAL
# ============================================================

if MAKE_FINAL_ZOOM_PLOT:
    zx1, zx2, zy1, zy2 = get_zoom_bounds(
        rl_clean,
        padding_x=ZOOM_PADDING_X,
        padding_y=ZOOM_PADDING_Y,
        n_pixels=400
    )

    rl_zoom = rl_clean[zy1:zy2, zx1:zx2]

    plt.figure(figsize=(9, 6), dpi=180)
    plt.imshow(
        normalize_for_display(rl_zoom, DISPLAY_GAMMA),
        cmap="inferno",
        interpolation=PLOT_INTERPOLATION
    )
    plt.title("Final cleaned image, zoomed on ion crystal", fontsize=16)
    plt.axis("off")
    plt.tight_layout()
    plt.show()

    print(f"Zoom bounds in ROI coordinates: x=[{zx1}:{zx2}], y=[{zy1}:{zy2}]")
