import numpy as np
import matplotlib.pyplot as plt

# List Data
x_data_z0 = np.arange(200, 450, 25)
#print(x_data_z0)
y_data_z0 = [0.123, 0.175, 0.223, 0.272, 0.317, 0.369, 0.415, 0.462, 0.508, 0.549]

x_data_rodthick = np.arange(100, 320, 20)
#print(x_data_rod_d0)
y_data_rodthick = [0.140, 0.145, 0.148, 0.155, 0.16, 0.168, 0.176, 0.182, 0.192, 0.202, 0.215]

x_data_roddist = np.arange(275, 525, 25)
#print(x_data_roddist)
y_data_roddist = [0.333, 0.275, 0.235, 0.207, 0.189, 0.172, 0.159, 0.15, 0.140, 0.134]

# plotting the points
plt.plot(x_data_z0, y_data_z0, label='DC Electrode Gap', color='blue', linestyle='solid', linewidth=1,
         marker='o', markerfacecolor='blue', markersize=3)

plt.plot(x_data_rodthick, y_data_rodthick, label='Rod Diameter', color='red', linestyle='solid', linewidth=1,
         marker='o', markerfacecolor='red', markersize=3)

plt.plot(x_data_roddist, y_data_roddist, label='Distance Between Rods', color='green', linestyle='solid', linewidth=1,
         marker='o', markerfacecolor='green', markersize=3)

# setting x and y axis range
plt.yticks(np.arange(0, 0.65, step=0.05))
plt.xticks(np.arange(100, 600, step=50))

# naming the x axis
plt.xlabel('Distance (µm)')
# naming the y axis
plt.ylabel('Axial Trap Depth (eV)')

plt.legend()

# Adding gridlines
plt.grid()
# function to show the plot
plt.show()

