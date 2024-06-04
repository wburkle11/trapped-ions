import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

# Seed for reproducibility
np.random.seed(2)

# Data
data = [[0, 0, 0, 0, 0.04],
        [0, 0, 0, 0.03, 0.09],
        [0, 0, 0, 0.08, 0.16],
        [0, 0, 0.06, 0.16, 0.31],
        [0, 0.03, 0.17, 0.36, 0.61],
        [0, 0.17, 0.48, 0.93, 1.48],
        [0.16, 0.93, 2.18, 3.94, 6.15]]



# Labels
xlabs = [50, 100, 150, 200, 250]
ylabs = [35, 30, 25, 20, 15, 10, 5]


# Heat map
fig, ax = plt.subplots(figsize = (8, 8))
#cmap = mpl.colors.ListedColormap(['royalblue', 'cyan', 'yellow', 'orange', 'red', 'green', 'black'])
#cmap.set_under('black')
#cmap.set_over('green')
cmap = ListedColormap(["darkorange", "gold", "lawngreen", "lightseagreen"])

magma = mpl.colormaps['magma']
newcolors = magma(np.linspace(0, 1, 256))
pink = np.array([248/256, 24/256, 148/256, 1])
newcolors[:25, :] = pink
newcmp = ListedColormap(newcolors)

im = ax.imshow(data, cmap=newcmp)


# Add the color bar
cbar = ax.figure.colorbar(im, ax = ax)
cbar.ax.set_ylabel("Radial Trap Depth (eV)", rotation = -90, va = "bottom")


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