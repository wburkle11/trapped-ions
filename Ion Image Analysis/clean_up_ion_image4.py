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

file_path = r"Z:\Lab Data\EMCCD\Images\raw_img_Wed Apr 29 2026_11.16.20_1513.tif"

# ----------------------------
# Auto-ROI settings
# ----------------------------
USE_AUTO_ROI = True
EXPECTED_N_IONS = 10

AUTO_ROI_PADDING_X = 70
AUTO_ROI_PADDING_Y = 70

# Robust compact-cluster ROI settings.
# This prevents isolated bright background/scatter peaks from pulling the ROI away.
ION_CLUSTER_RADIUS = 45          # pixels; increase to 60-80 for very spread-out crystals
MIN_IONS_IN_CLUSTER = 3
CENTER_ROI_ON_IONS = True
CENTERED_ROI_WIDTH = None        # None -> cluster span + padding
CENTERED_ROI_HEIGHT = None       # None -> cluster span + padding

BANDPASS_SMALL_SIGMA = 1.2
BANDPASS_LARGE_SIGMA = 12.0

PEAK_FOOTPRINT = 9
PEAK_PERCENTILE = 98.5
MAX_PEAKS_TO_TEST = 120
CLUSTER_RADIUS = ION_CLUSTER_RADIUS

# ----------------------------
# Tilt-image settings
# ----------------------------
THRESHOLD_FRAC = 0.35
MIN_SEPARATION = 5
FIT_HALF_WIDTH_TILT = 4
HX_DEFAULT = 5
HY_DEFAULT = 5
SAVE_ROTATED = True

# ----------------------------
# Display settings
# ----------------------------
DISPLAY_GAMMA = 1.0

# ----------------------------
# Gaussian-fit settings
# ----------------------------
FIT_SIGMA_INIT = 2.0
FIT_BG_MODE = "median"

# ----------------------------
# RL deconvolution settings
# ----------------------------
RL_ITERATIONS = 12
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


def load_image_average(path):
    img = imread(path).astype(float)

    if img.ndim == 3:
        return np.mean(img, axis=0)

    return img.copy()


def centered_bounds(cx, cy, width, height, nx, ny):
    """
    Return safe image bounds for a crop centered on (cx, cy).
    """
    width = int(round(min(width, nx)))
    height = int(round(min(height, ny)))

    x_min = int(round(cx - width / 2))
    y_min = int(round(cy - height / 2))
    x_max = x_min + width
    y_max = y_min + height

    if x_min < 0:
        x_max -= x_min
        x_min = 0
    if x_max > nx:
        x_min -= (x_max - nx)
        x_max = nx

    if y_min < 0:
        y_max -= y_min
        y_min = 0
    if y_max > ny:
        y_min -= (y_max - ny)
        y_max = ny

    return int(max(0, x_min)), int(x_max), int(max(0, y_min)), int(y_max)


def choose_compact_ion_cluster(candidates, peak_vals,
                               cluster_radius=45,
                               min_ions=3,
                               expected_n_ions=None):
    """
    Choose the densest compact group of ion-like peaks.

    This rejects isolated bright background spots by requiring peaks to be
    close to each other. The old logic could include far-away bright scatter
    when EXPECTED_N_IONS was larger than the number of compact ion peaks found.
    """
    best_score = -np.inf
    best_cluster = None
    best_vals = None

    for p in candidates:
        dx = candidates[:, 0] - p[0]
        dy = candidates[:, 1] - p[1]
        dist = np.sqrt(dx**2 + dy**2)

        cluster_idx = np.where(dist <= cluster_radius)[0]

        if len(cluster_idx) < min_ions:
            continue

        # If many local maxima exist inside the compact ion patch, keep the
        # expected number of strongest peaks, but only from this compact patch.
        if expected_n_ions is not None and len(cluster_idx) > expected_n_ions:
            local_vals = peak_vals[cluster_idx]
            keep = np.argsort(local_vals)[::-1][:expected_n_ions]
            cluster_idx = cluster_idx[keep]

        cluster = candidates[cluster_idx]
        vals = peak_vals[cluster_idx]

        x_span = cluster[:, 0].max() - cluster[:, 0].min() + 1
        y_span = cluster[:, 1].max() - cluster[:, 1].min() + 1
        area = max(x_span * y_span, 1)

        total_brightness = np.sum(vals)

        # Dense, bright, compact groups win.
        score = total_brightness * len(cluster_idx) / area

        if score > best_score:
            best_score = score
            best_cluster = cluster
            best_vals = vals

    # Fallback: if no cluster passes min_ions, use the single brightest peak.
    if best_cluster is None:
        i = np.argmax(peak_vals)
        best_cluster = candidates[[i]]
        best_vals = peak_vals[[i]]

    return best_cluster, best_vals


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
    ion_cluster_radius=ION_CLUSTER_RADIUS,
    min_ions_in_cluster=MIN_IONS_IN_CLUSTER,
    center_roi_on_ions=CENTER_ROI_ON_IONS,
    centered_roi_width=CENTERED_ROI_WIDTH,
    centered_roi_height=CENTERED_ROI_HEIGHT,
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

    # Pick densest compact bright cluster.
    # This intentionally rejects isolated bright background/scatter peaks.
    best_cluster, best_cluster_vals = choose_compact_ion_cluster(
        candidates=candidates,
        peak_vals=peak_vals,
        cluster_radius=ion_cluster_radius,
        min_ions=min_ions_in_cluster,
        expected_n_ions=expected_n_ions
    )

    if best_cluster is None:
        raise RuntimeError("Auto ROI failed: no valid ion-like cluster found.")

    cluster_x_min = int(best_cluster[:, 0].min())
    cluster_x_max = int(best_cluster[:, 0].max())
    cluster_y_min = int(best_cluster[:, 1].min())
    cluster_y_max = int(best_cluster[:, 1].max())

    cluster_cx = float(np.mean(best_cluster[:, 0]))
    cluster_cy = float(np.mean(best_cluster[:, 1]))

    if center_roi_on_ions:
        cluster_width = cluster_x_max - cluster_x_min + 1
        cluster_height = cluster_y_max - cluster_y_min + 1

        roi_width = centered_roi_width
        roi_height = centered_roi_height

        if roi_width is None:
            roi_width = cluster_width + 2 * padding_x
        if roi_height is None:
            roi_height = cluster_height + 2 * padding_y

        x_min, x_max, y_min, y_max = centered_bounds(
            cluster_cx,
            cluster_cy,
            roi_width,
            roi_height,
            nx,
            ny
        )
    else:
        x_min = max(0, cluster_x_min - padding_x)
        x_max = min(nx, cluster_x_max + padding_x + 1)

        y_min = max(0, cluster_y_min - padding_y)
        y_max = min(ny, cluster_y_max + padding_y + 1)

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
    print(f"Compact cluster radius used: {ion_cluster_radius} px")
    print(f"Detected compact ion-cluster peaks: {len(best_cluster)}")
    print(f"Ion-cluster center: x={cluster_cx:.2f}, y={cluster_cy:.2f}")
    print(f"Centered ROI enabled: {center_roi_on_ions}")
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
        ion_cluster_radius=ION_CLUSTER_RADIUS,
        min_ions_in_cluster=MIN_IONS_IN_CLUSTER,
        center_roi_on_ions=CENTER_ROI_ON_IONS,
        centered_roi_width=CENTERED_ROI_WIDTH,
        centered_roi_height=CENTERED_ROI_HEIGHT,
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

axes4[0].imshow(display_image(rotated_image, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
axes4[0].set_title("rotated raw image with ROI", fontsize=18)
axes4[0].axis("off")

axes4[1].imshow(display_image(global_rl, DISPLAY_GAMMA), cmap="inferno", interpolation="none")
axes4[1].set_title("cleaned image", fontsize=18)
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