# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 13:24:51 2026

@author: iontrap
"""

import numpy as np
import matplotlib.pyplot as plt
import tifffile

# --------------------------------------------------
# Load raw image
# --------------------------------------------------
tif_path = r"Z:\Lab Data\EMCCD\Images\raw_img_Mon Apr 27 2026_17.20.53_1442.tif"


data = tifffile.imread(tif_path)

print("Raw shape:", data.shape)
print("dtype:", data.dtype)

# If this is an image stack, average over frames
if data.ndim == 3:
    img = np.mean(data, axis=0)
else:
    img = data.astype(float)

# --------------------------------------------------
# Display full image
# --------------------------------------------------
plt.figure(figsize=(7, 6))
plt.imshow(img, cmap="magma", origin="lower")
plt.colorbar(label="Counts")
plt.title("Full raw image")
plt.xlabel("x pixel")
plt.ylabel("y pixel")
plt.show()

# --------------------------------------------------
# ROI around the 3 ions
# Adjust these numbers if needed
# --------------------------------------------------
x_min, x_max = 180, 415
y_min, y_max = 160, 340

roi = img[y_min:y_max, x_min:x_max]

# --------------------------------------------------
# Display ROI
# --------------------------------------------------
plt.figure(figsize=(6, 5))
plt.imshow(roi, cmap="magma", origin="lower",
           extent=[x_min, x_max, y_min, y_max])
#plt.colorbar(label="Counts")
#plt.title("ROI containing 3 ions")
plt.xlabel("x pixel")
plt.ylabel("y pixel")
plt.show()

print("ROI shape:", roi.shape)
print("ROI x range:", x_min, x_max)
print("ROI y range:", y_min, y_max)