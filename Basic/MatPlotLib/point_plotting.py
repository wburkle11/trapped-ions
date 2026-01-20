# Code for simple point plots 

import numpy as np
import matplotlib.pyplot as plt 

#x-axis data
x = [100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300]
             
#y-axis data
y = [0.140, 0.145, 0.148, 0.155, 0.160, 0.168, 0.176, 0.182, 0.192, 0.202, 0.215]

#Note: len(x) = len(y)

# plotting the points 
plt.plot(x, y, color='black', linestyle='None', linewidth = 2,
         marker='o', markerfacecolor='black', markersize=5)             

# setting x and y axis range
plt.yticks(np.arange(0.100, 0.240, step=0.01))
plt.xticks(np.arange(100, 320, step=20))


# naming the x axis
plt.xlabel('Diameter of RF Rods (µm)')
# naming the y axis
plt.ylabel('Trap Depth (eV)')
 
#Adding gridlines 
plt.grid()
# function to show the plot
plt.show()
