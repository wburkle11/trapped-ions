import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def custom_colormap():
    magma_cmap = plt.get_cmap('gist_rainbow')
    colors = magma_cmap(np.linspace(0, 1, 256))

    # Change color of the first index to white
    colors[0] = [0, 0, 0, 1]  # RGB values with alpha = 1

    return LinearSegmentedColormap.from_list('custom_colormap', colors)


# Data for Radial Trap Depth
data = [[0, 0, 0, 0, 0.04],
        [0, 0, 0, 0.03, 0.09],
        [0, 0, 0, 0.08, 0.16],
        [0, 0, 0.06, 0.16, 0.31],
        [0, 0.03, 0.17, 0.36, 0.61],
        [0, 0.17, 0.48, 0.93, 1.48]]
        #[0.16, 0.93, 2.18, 3.94, 6.15]]

# Data for Axial Trap Depth
data1 = [[0, 0, 0, 0, 0.13],
        [0, 0, 0, 0.13, 0.13],
        [0, 0, 0, 0.14, 0.14],
        [0, 0, 0.13, 0.14, 0.14],
        [0, 0.13, 0.14, 0.15, 0.17],
        [0, 0.14, 0.16, 0.18, 0.21]]
        #[0.14, 0.19, 0.26, 0.37, 0.5]]

# Data for Secular Frequencies
data2 = [[0, 0, 0, 0, 0.24],
        [0, 0, 0, 0.14, 0.40],
        [0, 0, 0, 0.36, 0.57],
        [0, 0, 0.29, 0.57, 0.80],
        [0, 0.14, 0.57, 0.87, 1.15],
        [0, 0.57, 1.01, 1.41, 1.80]]
        #[0.57, 1.41, 2.18, 2.94, 3.70]]

# Labels
xlabs = [50, 100, 150, 200, 250]
ylabs = [35, 30, 25, 20, 15, 10]

# Create custom colormap
cmap = custom_colormap()

# Heat map
fig, ax = plt.subplots(figsize = (8, 8))
im = ax.imshow(data2, cmap=cmap)

# Add the color bar
cbar = ax.figure.colorbar(im, ticks=np.arange(0, np.max(data2)+1, 0.10), ax = ax)
cbar.ax.set_ylabel("Radial Secular Frequency (MHz)", rotation = -90, va = "bottom")

# Add the labels
ax.set_xticks(np.arange(len(xlabs)), labels=xlabs)
ax.set_yticks(np.arange(len(ylabs)), labels=ylabs)

# naming the x axis
plt.xlabel('RF Voltage Amplitude (V)')
# naming the y axis
plt.ylabel('Resonant Frequency (MHz)')

# Rotate the labels of the X-axis
plt.setp(ax.get_xticklabels(), rotation=40,
         ha="right", rotation_mode="anchor")
plt.show()