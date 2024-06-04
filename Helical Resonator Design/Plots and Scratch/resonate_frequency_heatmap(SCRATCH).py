import numpy as np
import matplotlib.pyplot as plt

# Seed for reproducibility
np.random.seed(2)

# Data
data = [[0, 0, 0, 0, 0.49],
        [0, 0, 0, 0.42, 0.63],
        [0, 0, 0, 0.60, 0.82],
        [0, 0, 0.53, 0.83, 1.09],
        [0, 0.42, 0.83, 1.18, 1.52],
        [0, 0.83, 1.35, 1.85, 2.34],
        [0.83, 1.8, 2.8, 3.79, 4.75]]

# Labels
xlabs = [50, 100, 150, 200, 250]
ylabs = [35, 30, 25, 20, 15, 10, 5]

# Heat map
fig, ax = plt.subplots(figsize = (8, 8))
im = ax.imshow(data, cmap='magma')

# Add the color bar
cbar = ax.figure.colorbar(im, ax = ax)
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