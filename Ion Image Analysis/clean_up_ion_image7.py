# -*- coding: utf-8 -*-
"""
Robust ion-crystal cleanup with manual ROI and improved PSF selection.

Main idea:
1. Load full raw image
2. Manually crop ROI around crystal
3. Detect ion-like peaks inside the ROI
4. Keep many detected ions for diagnostics
5. Select only GOOD PSF candidates:
   - not too close to ROI boundary
   - sufficiently bright
   - sufficiently isolated from neighbors
   - not badly off-center in the local crop
6. Build averaged empirical PSF only from selected good candidates
7. Run Richardson-Lucy deconvolution

This version is designed to work better for:
- small crystals, where most ions can be used
- large dense crystals, where using every detected ion gives terrible PSFs
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.signal import fftconvolve
from imageio.v2 import imread, imwrite
from pathlib import Path


# ============================================================
# USER SETTINGS
# ============================================================

file_path = r"Z:\Lab Data\EMCCD\Images\raw_img_Tue Apr 28 2026_16.31.48_1491.tif"


# ============================================================
# MANUAL ROI SETTINGS -- SET THESE FIRST
# ============================================================

USE_MANUAL_ROI = True

# Full-image pixel coordinates.
# Crop used:
#     raw[MANUAL_Y_MIN:MANUAL_Y_MAX, MANUAL_X_MIN:MANUAL_X_MAX]
#
# Make this large enough to include the whole crystal and coma/tails.
MANUAL_X_MIN = 120
MANUAL_X_MAX = 480
MANUAL_Y_MIN = 200
MANUAL_Y_MAX = 320

# Alternative center-size ROI. If True, explicit bounds above are ignored.
USE_CENTER_SIZE_ROI = False
MANUAL_CENTER_X = 300
MANUAL_CENTER_Y = 260
MANUAL_ROI_WIDTH = 360
MANUAL_ROI_HEIGHT = 120

SAVE_MANUAL_ROI = True


# ============================================================
# IMAGE / DISPLAY SETTINGS
# ============================================================

FRAME_MODE = "average"          # "first" or "average" for TIFF stacks
DISPLAY_GAMMA = 1.0
DISPLAY_CMAP = "inferno"
DISPLAY_INTERPOLATION = "none"

PIXEL_SIZE_UM = 0.4


# ============================================================
# ION DETECTION SETTINGS
# ============================================================

# Bandpass: small blur keeps ion-sized spots, large blur removes broad background.
BANDPASS_SMALL_SIGMA = 1.0
BANDPASS_LARGE_SIGMA = 8.0

# Peak threshold. Lower for dim images, higher if background spots are being detected.
PEAK_PERCENTILE = 98.0

# Non-maximum suppression radius. Should be comparable to nearest-neighbor spacing / 2.
# For small crystals: 4-6 is usually fine.
# For dense large crystals: 4-5 may be needed.
MIN_PEAK_SEPARATION = 5

# Limit total detected peaks for diagnostics.
MAX_DETECTED_PEAKS = 300


# ============================================================
# PSF SELECTION SETTINGS
# ============================================================

# Half-width of each PSF crop.
# Crop size will be (2*PSF_HALF_WIDTH_Y+1, 2*PSF_HALF_WIDTH_X+1).
PSF_HALF_WIDTH_X = 7
PSF_HALF_WIDTH_Y = 7

# Number of high-quality PSFs to average.
# For small crystals, this may use most/all ions.
# For large crystals, this prevents using hundreds of bad overlapping crops.
MAX_PSFS_TO_USE = 50

# Minimum number to proceed. If fewer pass strict isolation, the code relaxes automatically.
MIN_PSFS_TO_USE = 5

# Preferred isolation distance from nearest neighbor, in pixels.
# For dense crystals, this may be impossible; code will relax if needed.
PREFERRED_MIN_NEIGHBOR_DISTANCE = 12

# Brightness filter. Candidate peak amplitude must be above this fraction of the brightest peak.
MIN_PEAK_FRACTION_OF_MAX = 0.20

# Reject candidates whose local center-of-mass is too far from the crop center.
# This catches cases where the detected x is on a tail rather than the ion core.
MAX_CENTER_OFFSET_PIXELS = 3.0

# Use a soft mask on each PSF crop before averaging.
USE_PSF_WEIGHT_MASK = True
PSF_MASK_SIGMA_X_SCALE = 0.75
PSF_MASK_SIGMA_Y_SCALE = 0.85

# Background subtraction for each PSF crop.
PSF_BACKGROUND_PERCENTILE = 15


# ============================================================
# RICHARDSON-LUCY SETTINGS
# ============================================================

RL_ITERATIONS = 10
EPS = 1e-12


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def display_image(img, gamma=1.0):
    disp = img.astype(float).copy()
    disp -= np.nanmin(disp)
    if np.nanmax(disp) > 0:
        disp /= np.nanmax(disp)
    return disp ** gamma


def load_image(path, frame_mode="average"):
    data = imread(path).astype(float)

    if data.ndim == 2:
        print(f"Loaded single image: shape={data.shape}")
        return data.copy()

    if data.ndim == 3:
        print(f"Loaded stack: shape={data.shape}")
        if frame_mode.lower() == "first":
            print("Using first frame.")
            return data[0].copy()
        else:
            print("Using average over frames.")
            return np.mean(data, axis=0)

    raise RuntimeError(f"Unexpected image shape: {data.shape}")


def get_manual_roi_bounds(raw):
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


def crop_manual_roi(file_path):
    raw = load_image(file_path, FRAME_MODE)
    x_min, x_max, y_min, y_max = get_manual_roi_bounds(raw)

    roi = raw[y_min:y_max, x_min:x_max].copy()

    input_path = Path(file_path)
    roi_path = input_path.with_name(input_path.stem + "_manualROI.tif")

    if SAVE_MANUAL_ROI:
        imwrite(roi_path, roi.astype(np.float32))

    print("\n================ MANUAL ROI ================")
    print(f"Full image shape: {raw.shape}")
    print(f"ROI bounds: x=[{x_min}:{x_max}], y=[{y_min}:{y_max}]")
    print(f"ROI shape: {roi.shape}")
    print(f"Saved ROI: {roi_path}")

    fig, axes = plt.subplots(2, 1, figsize=(9, 9), dpi=150)

    axes[0].imshow(display_image(raw, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
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
    axes[0].set_title("Full raw image with manual ROI")
    axes[0].axis("off")

    axes[1].imshow(display_image(roi, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
    axes[1].set_title("Manual ROI used for cleanup")
    axes[1].axis("off")

    plt.tight_layout()
    plt.show()

    return raw, roi, roi_path, (x_min, x_max, y_min, y_max)


def make_bandpass_image(img):
    work = img.astype(float).copy()
    bg = np.percentile(work, 10)
    work -= bg
    work[work < 0] = 0

    small = gaussian_filter(work, sigma=BANDPASS_SMALL_SIGMA)
    large = gaussian_filter(work, sigma=BANDPASS_LARGE_SIGMA)

    band = small - large
    band[band < 0] = 0

    return band


def detect_peaks(img):
    """
    Detect ion-like peaks by bandpass + local maxima + greedy nonmax suppression.
    """
    band = make_bandpass_image(img)

    positive = band[band > 0]
    if positive.size == 0:
        raise RuntimeError("Bandpass image has no positive signal.")

    threshold = np.percentile(positive, PEAK_PERCENTILE)

    local_max = band == maximum_filter(band, size=2 * MIN_PEAK_SEPARATION + 1)
    mask = local_max & (band >= threshold)

    ys, xs = np.where(mask)
    vals = band[ys, xs]

    if len(xs) == 0:
        raise RuntimeError("No peaks detected. Lower PEAK_PERCENTILE.")

    # Sort by strength.
    order = np.argsort(vals)[::-1]
    xs = xs[order]
    ys = ys[order]
    vals = vals[order]

    # Greedy non-maximum suppression.
    kept = []
    for x, y, v in zip(xs, ys, vals):
        if len(kept) >= MAX_DETECTED_PEAKS:
            break

        if not kept:
            kept.append((x, y, v))
            continue

        prev = np.array([[p[0], p[1]] for p in kept])
        d = np.sqrt((prev[:, 0] - x)**2 + (prev[:, 1] - y)**2)

        if np.all(d >= MIN_PEAK_SEPARATION):
            kept.append((x, y, v))

    peaks = np.array(kept, dtype=float)

    print("\n================ PEAK DETECTION ================")
    print(f"Detected peaks: {len(peaks)}")
    print(f"Peak threshold percentile: {PEAK_PERCENTILE}")

    return peaks, band


def nearest_neighbor_distances(peaks):
    if len(peaks) <= 1:
        return np.full(len(peaks), np.inf)

    xy = peaks[:, :2]
    dists = np.full(len(peaks), np.inf)

    for i in range(len(peaks)):
        dx = xy[:, 0] - xy[i, 0]
        dy = xy[:, 1] - xy[i, 1]
        d = np.sqrt(dx**2 + dy**2)
        d[i] = np.inf
        dists[i] = np.min(d)

    return dists


def extract_patch(img, x, y, hx, hy):
    ny, nx = img.shape

    xi = int(round(x))
    yi = int(round(y))

    x1 = xi - hx
    x2 = xi + hx + 1
    y1 = yi - hy
    y2 = yi + hy + 1

    if x1 < 0 or y1 < 0 or x2 > nx or y2 > ny:
        return None, None

    return img[y1:y2, x1:x2].copy(), (x1, x2, y1, y2)


def crop_center_of_mass(patch):
    """
    Estimate subpixel center of the bright spot inside a patch using a
    background-subtracted center of mass.
    """
    p = patch.astype(float).copy()
    p -= np.percentile(p, PSF_BACKGROUND_PERCENTILE)
    p[p < 0] = 0

    if p.sum() <= 0:
        cy = (patch.shape[0] - 1) / 2
        cx = (patch.shape[1] - 1) / 2
        return cx, cy, 0.0

    yy, xx = np.indices(p.shape)
    cx = np.sum(xx * p) / np.sum(p)
    cy = np.sum(yy * p) / np.sum(p)

    peak_signal = np.max(p)
    return cx, cy, peak_signal


def select_good_psf_candidates(img, peaks):
    """
    Select good PSF candidates from detected peaks.

    Uses a score that prefers:
    - bright peaks
    - isolated peaks
    - centered local crops
    - not near boundaries

    This avoids using every ion in large dense crystals.
    """
    nn_dist = nearest_neighbor_distances(peaks)
    max_val = np.max(peaks[:, 2])

    candidates = []

    for i, (x, y, val) in enumerate(peaks):
        if val < MIN_PEAK_FRACTION_OF_MAX * max_val:
            continue

        patch, bounds = extract_patch(img, x, y, PSF_HALF_WIDTH_X, PSF_HALF_WIDTH_Y)
        if patch is None:
            continue

        cx_local, cy_local, peak_signal = crop_center_of_mass(patch)

        crop_cx = PSF_HALF_WIDTH_X
        crop_cy = PSF_HALF_WIDTH_Y
        center_offset = np.sqrt((cx_local - crop_cx)**2 + (cy_local - crop_cy)**2)

        if center_offset > MAX_CENTER_OFFSET_PIXELS:
            continue

        # Score: bright, isolated, well-centered.
        isolation_score = min(nn_dist[i] / PREFERRED_MIN_NEIGHBOR_DISTANCE, 1.5)
        centering_score = max(0.1, 1.0 - center_offset / MAX_CENTER_OFFSET_PIXELS)

        score = val * isolation_score * centering_score

        candidates.append({
            "index": i,
            "x": x,
            "y": y,
            "val": val,
            "nn_dist": nn_dist[i],
            "center_offset": center_offset,
            "score": score,
            "patch": patch,
            "bounds": bounds,
        })

    # Prefer sufficiently isolated candidates.
    strict = [c for c in candidates if c["nn_dist"] >= PREFERRED_MIN_NEIGHBOR_DISTANCE]

    if len(strict) >= MIN_PSFS_TO_USE:
        pool = strict
        mode = "strict isolated"
    else:
        pool = candidates
        mode = "relaxed"

    if len(pool) == 0:
        raise RuntimeError(
            "No usable PSF candidates. Try lowering PEAK_PERCENTILE, "
            "lowering MIN_PEAK_FRACTION_OF_MAX, or increasing MAX_CENTER_OFFSET_PIXELS."
        )

    pool = sorted(pool, key=lambda c: c["score"], reverse=True)
    selected = pool[:MAX_PSFS_TO_USE]

    print("\n================ PSF SELECTION ================")
    print(f"Candidate PSFs after quality cuts: {len(candidates)}")
    print(f"Selection mode: {mode}")
    print(f"Selected PSFs: {len(selected)}")
    print("Selected peak numbers:", [int(c["index"] + 1) for c in selected])

    return selected, candidates


def build_average_psf(selected):
    psfs = []

    for c in selected:
        psf = c["patch"].astype(float).copy()

        # Robust local background subtraction.
        psf -= np.percentile(psf, PSF_BACKGROUND_PERCENTILE)
        psf[psf < 0] = 0

        if USE_PSF_WEIGHT_MASK:
            yy, xx = np.indices(psf.shape)
            cx = (psf.shape[1] - 1) / 2
            cy = (psf.shape[0] - 1) / 2

            sig_x = PSF_MASK_SIGMA_X_SCALE * PSF_HALF_WIDTH_X
            sig_y = PSF_MASK_SIGMA_Y_SCALE * PSF_HALF_WIDTH_Y

            mask = np.exp(
                -((xx - cx)**2 / (2 * sig_x**2) +
                  (yy - cy)**2 / (2 * sig_y**2))
            )
            psf *= mask

        if psf.sum() <= 0:
            continue

        psf /= psf.sum()
        psfs.append(psf)

    if len(psfs) == 0:
        raise RuntimeError("All selected PSF crops were empty after background subtraction.")

    avg = np.mean(np.stack(psfs, axis=0), axis=0)
    avg /= avg.sum()

    return avg, psfs


def richardson_lucy(image, psf, iterations=RL_ITERATIONS, eps=EPS):
    image = image.astype(float).copy()
    image[image < 0] = 0

    psf = psf.astype(float).copy()
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


def plot_detected_and_selected(img, peaks, selected):
    fig, ax = plt.subplots(figsize=(16, 5), dpi=150)
    ax.imshow(display_image(img, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
    ax.set_title("Detected peaks and selected PSFs")
    ax.axis("off")

    # All detected peaks: small white dots / labels.
    ax.scatter(peaks[:, 0], peaks[:, 1], c="white", marker=".", s=10, alpha=0.5)

    selected_indices = set(c["index"] for c in selected)

    for i, (x, y, val) in enumerate(peaks):
        if i in selected_indices:
            rect = Rectangle(
                (x - PSF_HALF_WIDTH_X - 0.5, y - PSF_HALF_WIDTH_Y - 0.5),
                2 * PSF_HALF_WIDTH_X + 1,
                2 * PSF_HALF_WIDTH_Y + 1,
                linewidth=1.5,
                edgecolor="cyan",
                facecolor="none"
            )
            ax.add_patch(rect)
            ax.scatter([x], [y], c="cyan", marker="x", s=35)
            ax.text(x, y - PSF_HALF_WIDTH_Y - 2, str(i + 1), color="white",
                    fontsize=8, ha="center", va="bottom")
        else:
            # label fewer unselected peaks to avoid clutter
            pass

    plt.tight_layout()
    plt.show()


def plot_psf_gallery(selected, avg_psf):
    n = len(selected)
    ncols = min(6, n)
    nrows = int(np.ceil((n + 1) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(2.4 * ncols, 2.4 * nrows), dpi=150)
    axes = np.atleast_1d(axes).ravel()

    for ax in axes:
        ax.axis("off")

    for k, c in enumerate(selected):
        ax = axes[k]
        ax.imshow(display_image(c["patch"], DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
        ax.set_title(
            f"peak {int(c['index']+1)}\n"
            f"NN={c['nn_dist']:.1f}, off={c['center_offset']:.1f}",
            fontsize=8
        )
        ax.axis("off")

    axes[n].imshow(display_image(avg_psf, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
    axes[n].set_title("Averaged PSF", fontsize=9)
    axes[n].axis("off")

    plt.tight_layout()
    plt.show()


def plot_final_images(raw_roi, cleaned):
    # Pixel coordinates
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), dpi=150)

    for ax, img, title in zip(
        axes,
        [raw_roi, cleaned],
        ["raw ROI image (pixels)", "cleaned image (pixels)"]
    ):
        ax.imshow(display_image(img, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
        ax.set_title(title)
        ax.set_xlabel("Pixel X")
        ax.set_ylabel("Pixel Y")

        ny, nx = img.shape
        ax.set_xticks(np.arange(0, nx, 20))
        ax.set_yticks(np.arange(0, ny, 20))
        ax.tick_params(labelsize=8)
        ax.set_aspect("equal")

    plt.tight_layout()
    plt.show()

    # Physical units
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), dpi=150)

    tick_um = 10
    tick_px = max(1, int(round(tick_um / PIXEL_SIZE_UM)))

    for ax, img, title in zip(
        axes,
        [raw_roi, cleaned],
        ["raw ROI image (µm)", "cleaned image (µm)"]
    ):
        ax.imshow(display_image(img, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
        ax.set_title(title)
        ax.set_xlabel("X (µm)")
        ax.set_ylabel("Y (µm)")

        ny, nx = img.shape
        xticks = np.arange(0, nx, tick_px)
        yticks = np.arange(0, ny, tick_px)

        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.set_xticklabels(np.round(xticks * PIXEL_SIZE_UM, 1))
        ax.set_yticklabels(np.round(yticks * PIXEL_SIZE_UM, 1))

        ax.tick_params(labelsize=8)
        ax.set_aspect("equal")

    plt.tight_layout()
    plt.show()


# ============================================================
# MAIN
# ============================================================

raw_full, roi, roi_path, roi_bounds = crop_manual_roi(file_path)

peaks, band = detect_peaks(roi)

fig, axes = plt.subplots(2, 1, figsize=(12, 7), dpi=150)
axes[0].imshow(display_image(roi, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
axes[0].scatter(peaks[:, 0], peaks[:, 1], c="cyan", marker="x", s=20)
axes[0].set_title("ROI with detected peaks")
axes[0].axis("off")

axes[1].imshow(display_image(band, DISPLAY_GAMMA), cmap=DISPLAY_CMAP, interpolation=DISPLAY_INTERPOLATION)
axes[1].scatter(peaks[:, 0], peaks[:, 1], c="cyan", marker="x", s=20)
axes[1].set_title("Bandpassed ROI used for peak detection")
axes[1].axis("off")
plt.tight_layout()
plt.show()

selected, all_candidates = select_good_psf_candidates(roi, peaks)

plot_detected_and_selected(roi, peaks, selected)

avg_psf, psf_list = build_average_psf(selected)

plot_psf_gallery(selected, avg_psf)

global_rl = richardson_lucy(
    roi,
    avg_psf,
    iterations=RL_ITERATIONS,
    eps=EPS
)

input_path = Path(file_path)
rl_path = input_path.with_name(input_path.stem + "_manualROI_robustPSF_RL.tif")
imwrite(rl_path, global_rl.astype(np.float32))
print(f"\nSaved cleaned image to: {rl_path}")

plot_final_images(roi, global_rl)
