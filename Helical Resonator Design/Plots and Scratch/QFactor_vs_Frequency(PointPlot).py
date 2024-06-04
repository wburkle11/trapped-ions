# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 11:36:32 2024

@author: iontrap
"""

#Point Plot for Resonant Frequency vs. Trap Capacitance
import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#x-axis data
x = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]


#y-axis data
y = [72.2, 48.3, 38.8, 33.3, 29.6, 26.9, 24.9, 23.2, 21.9, 20.7, 19.8, 18.9, 18.1, 17.5, 16.9, 16.3, 15.8, 15.4, 15.0, 14.6, 14.2]


fig, ax = plt.subplots(figsize = (10, 8))
ax.plot(x, y)
#find line of best fit
#a, b = np.polyfit(x, y, 1)
# setting x and y axis range
ax.secondary_yaxis('right')
ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.xaxis.set_major_formatter('{x:.0f}')
ax.yaxis.set_minor_locator(AutoMinorLocator())
#ax.secondary_yaxis.set_minor_locator(AutoMinorLocator())

ax.set_yticks(np.arange(0, 120, step=20))
ax.set_xticks(np.arange(0, 120, step=20))
#plt.plot(x, math.exp(-b))
# naming the x axis
plt.xlabel('Trap Capacitance [pF]')
# naming the y axis
plt.ylabel('Frequency [MHz]')
#Adding gridlines 
plt.grid()
# function to show the plot
plt.show()
