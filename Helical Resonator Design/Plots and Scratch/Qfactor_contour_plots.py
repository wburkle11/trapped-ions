# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 16:27:43 2024

@author: iontrap
"""


# Implementation of matplotlib function 
import matplotlib.pyplot as plt 
import numpy as np 
import math
x = [0.44, 0.47, 0.5, 0.53, 0.56, 0.59, 0.62, 0.66, 0.71, 0.83, 0.01, 0.14, 0.16, 0.20] 
d = [0.045, 0.048, 0.051, 0.054, 0.057, 0.06, 0.063, 0.067, 0.072, 0.085, 0.01, 0.014, 0.017, 0.020] 
 
  
# Creating 2-D grid of features 
[X, Y] = np.meshgrid(x, d) 
  
fig, ax = plt.subplots(1, 1) 
print([X, Y]) 
#Z = X**2-Y**2
#Z = math.cos(math.sqrt((X**2+Y**2)))

# plots filled contour plot 
ax.contourf(X, Y, Z) 
  
ax.set_title('Q Factor') 
ax.set_xlabel('d/D') 
ax.set_ylabel('d (m)') 
  
plt.show()