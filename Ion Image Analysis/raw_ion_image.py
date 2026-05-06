
# -*- coding: utf-8 -*-
"""
Simple TIFF viewer for trapped-ion image data

What it does:
1) Loads a .tif file
2) Prints the data shape and dtype
3) Displays the image

If the TIFF contains a stack of frames, it will show:
- the first frame
- the average over all frames
"""

import numpy as np
import matplotlib.pyplot as plt
import tifffile


# ============================================================
# USER INPUT: path to your TIFF file
# ============================================================
tif_path = r"Z:\Lab Data\EMCCD\Images\raw_img_Wed Apr 29 2026_15.07.13_1545.tif"

# ============================================================
# LOAD TIFF
# ============================================================
data = tifffile.imread(tif_path)

print("Loaded TIFF file:")
print("Shape:", data.shape)
print("Dtype:", data.dtype)


# ============================================================
# DISPLAY IMAGE
# ============================================================
if data.ndim == 2:
    # Single image: shape = (height, width)
    plt.figure(figsize=(8, 6))
    plt.imshow(data, cmap="inferno", origin="lower")
    plt.tight_layout()
    plt.show()

elif data.ndim == 3:
    # Stack of images: shape = (n_frames, height, width)
    first_frame = data[0]
    avg_frame = np.mean(data, axis=0)

    plt.figure(figsize=(8, 6))
    plt.imshow(first_frame, cmap="inferno", origin="lower")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.imshow(avg_frame, cmap="inferno", origin="lower")
    plt.tight_layout()
    plt.show()

else:
    print("Unexpected TIFF shape. Could not automatically display.")