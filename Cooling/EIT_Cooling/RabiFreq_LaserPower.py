import numpy as np
import matplotlib.pyplot as plt
from sympy import S, sqrt as sy_sqrt
import math

# -----------------------------
# Physical constants
# -----------------------------
omega = 751e12 * 2 * math.pi  # Laser frequency (rad/s)
gamma = 19.6e6 * 2 * math.pi  # Linewidth (rad/s)
hbar = 1.06e-34               # Reduced Planck constant (J·s)
c = 3.0e8                     # Speed of light (m/s)
eps = 8.8e-12                 # Vacuum permittivity (F/m)

# -----------------------------
# Clebsch-Gordan coefficients
# -----------------------------
threedrive = S(1) / sy_sqrt(3)
threeprobe = -S(1) / sy_sqrt(3)
sixdrive = S(1) / sy_sqrt(6)
sixprobe = S(1) / sy_sqrt(6)

# -----------------------------
# Dipole matrix element
# -----------------------------
sqrt_factor = sy_sqrt((2 * 0 + 1) * (2 * 1 + 1) * (2 * S(1)/2 + 1))
dipole_prefactor = sy_sqrt(3 * math.pi * eps * hbar * c**3 * gamma / omega**3)
phase = (-1)**(1/2 + 1/2 + 1 + 1)

drive_element = phase * sqrt_factor * threedrive * sixdrive * dipole_prefactor
probe_element = phase * sqrt_factor * threeprobe * sixprobe * dipole_prefactor

# Evaluate to float
drive_element_val = float(drive_element.evalf())
probe_element_val = float(probe_element.evalf())

# -----------------------------
# Laser beam waists to compare
# -----------------------------
beam_waists = [45e-6, 60e-6, 80e-6]  # in meters  # <<< NEW
colors = ['blue', 'orange', 'red']  # one color per beam waist  # <<< NEW

# -----------------------------
# Rabi frequency ranges
# -----------------------------
rabidrive_array = np.linspace(1e6 * 2 * np.pi, 50e6 * 2 * np.pi, 200)
rabiprobe_array = np.linspace(1e6 * 2 * np.pi, 50e6 * 2 * np.pi, 200)

# -----------------------------
# Electric fields
# -----------------------------
E_drive = (-1 * rabidrive_array * hbar) / drive_element_val
E_probe = (-1 * rabiprobe_array * hbar) / probe_element_val

# -----------------------------
# Plotting
# -----------------------------
fig, (ax1) = plt.subplots(1, 1, figsize=(8, 6), sharex=True)

# Loop over beam waists
for radius, color in zip(beam_waists, colors):  # <<< NEW
    spot = math.pi * radius**2

    # Intensities and powers
    Idrive = E_drive**2 * eps * c / 2
    Iprobe = E_probe**2 * eps * c / 2

    Pdrive = Idrive * spot * 1e6  # µW
    Pprobe = Iprobe * spot * 1e6  # µW

    label = f"Beam waist = {int(radius*1e6)} µm"

    ax1.plot(Pdrive, rabidrive_array / (2 * np.pi * 1e6), color=color, label=label)
    #ax2.plot(Pprobe, rabiprobe_array / (2 * np.pi * 1e6), color=color, label=label)

# Labels and styling
ax1.set_ylabel('Rabi Frequency (MHz)')
ax1.set_title('Drive Laser Power vs Rabi Frequency')
ax1.grid(True)
ax1.legend()

ax1.set_xlabel('Laser Power (µW)')
#ax2.set_ylabel('Rabi Frequency (MHz)')
#ax2.set_title('Probe Laser Power vs Rabi Frequency')
#ax2.grid(True)
#ax2.legend()

# Set custom x-ticks
max_power = max(Pdrive.max(), Pprobe.max())
xticks = np.linspace(0, max_power, 15)  # 8 evenly spaced ticks, adjust as needed

ax1.set_xticks(xticks)
#ax2.set_xticks(xticks)

plt.tight_layout()
plt.show()
