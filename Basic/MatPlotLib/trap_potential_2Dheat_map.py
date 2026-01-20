import numpy as np
import matplotlib.pyplot as plt

# Seed for reproducibility
np.random.seed(2)

# Data
data = [[0.131, 0.178, 0.224, 0.28, 0.32, 0.37, 0.42, 0.46, 0.51, 0.56],
        [0.144, 0.184, 0.24, 0.29, 0.34, 0.39, 0.44, 0.49, 0.54, 0.58],
        [0.15, 0.2, 0.25, 0.31, 0.35, 0.41, 0.47, 0.51, 0.57, 0.6],
        [0.16, 0.211, 0.264, 0.32, 0.38, 0.43, 0.49, 0.54, 0.58, 0.63],
        [0.171, 0.228, 0.284, 0.35, 0.4, 0.45, 0.51, 0.56, 0.61, 0.65]]

# Labels
xlabs = [200, 225, 250, 275, 300, 325, 350, 375, 400, 425]
ylabs = [500, 475, 450, 425, 400]

# Heat map
fig, ax = plt.subplots(figsize = (8, 8))
im = ax.imshow(data, cmap='Blues')

# Add the color bar
cbar = ax.figure.colorbar(im, ax = ax)
cbar.ax.set_ylabel("Trap Depth in Axial Direction (eV)", rotation = -90, va = "bottom")

# Add the labels
ax.set_xticks(np.arange(len(xlabs)), labels=xlabs)
ax.set_yticks(np.arange(len(ylabs)), labels=ylabs)

# naming the x axis
plt.xlabel('Distance Between DC Fiber Tips (µm)')
# naming the y axis
plt.ylabel('Distance Between RF Rods (µm)')

# Rotate the labels of the X-axis
plt.setp(ax.get_xticklabels(), rotation=40,
         ha="right", rotation_mode="anchor")
plt.show()