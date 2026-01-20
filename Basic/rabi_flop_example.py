import numpy as np
import matplotlib.pyplot as plt

# Time array (0 to 2 μs)
t = np.linspace(0, 5e-5, 1000)  # seconds

# Rabi frequency (in Hz)
Omega = 2 * np.pi * 46e3  # 1 MHz

# State populations
P_0 = np.cos(Omega * t / 2)**2
P_1 = np.sin(Omega * t / 2)**2

# π-pulse time
t_pi = np.pi / Omega  # seconds

# Plot
plt.figure(figsize=(8, 5))
plt.plot(t * 1e6, P_1, color='red')

# Vertical line at π/Omega
plt.axvline(t_pi * 1e6, color='black', linestyle='--',
            label=r"$\pi - time$")

# Aesthetics
plt.title(r"Rabi Flopping at $\Omega_{n_i,\,n_i-1} = 46kHz$", fontsize=16)
plt.xlabel("Time (μs)", fontsize=12)
plt.ylabel("Photon Counts", fontsize=12)
plt.ylim(0, 1.05)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

