import numpy as np
import matplotlib.pyplot as plt

tau_g = 20e-6
tau_d = 100e-6
t_c   = 100e-6

delta_max = 300e3 * 2*np.pi   # rad/s
delta_min = 25e3 * 2*np.pi    # rad/s

V0 = 1.0
j = 3

T = 2*tau_g + 2*tau_d + t_c
dt = 1e-8
t = np.arange(0, T, dt)

def delta_down(t):
    b = delta_max**(-j)
    c = (2/tau_d)*(delta_min**(-j) - delta_max**(-j))
    g = 0.5*t - (tau_d/(4*np.pi))*np.sin(2*np.pi*t/tau_d)
    return (b + c*g)**(-1/j)

def delta_up(t):
    return delta_down(tau_d - t)

delta = np.zeros_like(t)
V = np.zeros_like(t)

for i, ti in enumerate(t):
    if ti < tau_g:
        delta[i] = delta_max
    elif ti < tau_g + tau_d:
        delta[i] = delta_down(ti - tau_g)
    elif ti < tau_g + tau_d + t_c:
        delta[i] = delta_min
    elif ti < tau_g + tau_d + t_c + tau_d:
        delta[i] = delta_up(ti - (tau_g + tau_d + t_c))
    else:
        delta[i] = delta_max

    if ti < tau_g:
        V[i] = V0 * np.sin(0.5*np.pi*ti/tau_g)**2
    elif ti < tau_g + tau_d + t_c + tau_d:
        V[i] = V0
    else:
        x = ti - (tau_g + tau_d + t_c + tau_d)
        V[i] = V0 * np.sin(0.5*np.pi*(1 - x/tau_g))**2

# Correct phase: phi(t) = integral delta(t) dt
phi = np.cumsum(delta) * dt

signal = V * np.sin(phi)

fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

ax[0].plot(t*1e6, delta/(2*np.pi*1e3))
ax[0].set_ylabel("Detuning (kHz)")
ax[0].set_title("Smooth Detuning Ramp")

ax[1].plot(t*1e6, V)
ax[1].set_ylabel("Amplitude")

ax[2].plot(t*1e6, signal)
ax[2].set_ylabel("V(t) sin(phi(t))")
ax[2].set_xlabel("Time (µs)")

plt.tight_layout()
plt.show()