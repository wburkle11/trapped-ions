# Point Plot for QFactor vs. Trap Capacitance
import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# x-axis data
x = [0, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]

# y-axis data
y = [1000, 380.9, 116.3, 56.8, 40.0, 32.0, 27.2, 24.0, 21.6, 19.9, 18.4, 17.3, 16.3, 15.5, 14.8, 14.1, 13.6, 13.1, 12.7, 12.3, 11.9, 11.6]

z = [1000, 1000, 946, 498.2, 359.3, 290.7, 249.2, 221.0, 200.3, 184.5, 171.8, 161.4, 152.6, 145.1, 138.7, 133.0, 128.0, 123.4, 119.4, 115.7, 112.3, 109.3]

#z1 = [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 974.8, 930.9, 893.2, 860.2, 831.0, 804.9, 781.4, 760.1, 740.7, 722.8, 706.3, 691.0]

fig, ax = plt.subplots(figsize=(10, 8))
ax.plot(x, y)
ax.plot(x, z, ls='-')
#ax.plot(x, z1, ls='--')
# find line of best fit
# a, b = np.polyfit(x, y, 1)
# setting x and y axis range
#ax.secondary_yaxis('right')
ax.xaxis.set_minor_locator(AutoMinorLocator())
# ax.xaxis.set_major_formatter('{x:.0f}')
ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.secondary_yaxis.set_minor_locator(AutoMinorLocator())

ax.set_yticks(np.arange(0, 1200, step=200))
ax.set_xticks(np.arange(0, 120, step=20))
# plt.plot(x, math.exp(-b))
# naming the x axis
plt.xlabel('Trap Capacitance [pF]')
# naming the y axis
plt.ylabel('Q')

# Adding gridlines
plt.grid()
# function to show the plot
plt.show()