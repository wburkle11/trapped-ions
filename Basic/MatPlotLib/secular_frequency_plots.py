import numpy as np
import matplotlib.pyplot as plt

# List Data
x_data = [200, 225, 250, 275, 300, 325, 350, 375, 400, 425]
data1_rad = [0.72, 0.69, 0.69, 0.7, 0.73, 0.76, 0.81, 0.84, 0.87, 0.9]
data2_ax = [637.7, 661.6, 672.9, 673.9, 667.9, 656.4, 641.3, 623.7, 604.6, 584.7]

#Data for Chaning Optical Fiber Diameter (plot 2)
#x_data = [280, 265, 250, 235, 220]
#data1_rad = [0.72, 0.65, 0.55, 0, 0]
#data2_ax = [637.7, 664.08, 694.6, 726.8, 760.8]

#Data for Changing Rod Diameter (plot 3)
#x_data = [100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300]
#data1_rad = [0.58, 0.82, 1.02, 1.22, 1.42, 1.62, 1.84, 2.07, 2.34, 2.63, 2.95]
#data2_ax = [661.17, 673.26, 684.90, 697.49, 710.03, 723.8, 738.27, 753.78, 770.8, 789.61, 810.94]

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Distance between Optical Fiber Faces (µm)')
ax1.set_ylabel('Radial Secular Frequency (MHz)', color=color)
ax1.plot(x_data, data1_rad, color='red', linestyle='dotted', marker='o', markerfacecolor='red', markersize=5)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(0.68, 0.91)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Axial Secular Frequency (kHz)', color=color)  # we already handled the x-label with ax1
ax2.plot(x_data, data2_ax, color='blue', linestyle='dashdot', marker='o', markerfacecolor='blue', markersize=5)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

plt.show()