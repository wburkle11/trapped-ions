# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Parameters copied from full_ms_smoothgate_AWGdrive.py
# ============================================================
dt_awg = 2e-6       # DDS update interval
tau_g  = 20e-6      # amplitude ramp on/off
tau_d  = 100e-6     # detuning ramp down/up
t_c    = 100e-6     # hold at delta_min
j      = 3

delta_max_hz = 300e3
delta_min_hz = 25e3
amp_max_frac = 0.40

f_carrier_hz  = 222e6
f_mode_hz     = 1e6
f_blue_ref_hz = f_carrier_hz + f_mode_hz   # 223 MHz
f_red_ref_hz  = f_carrier_hz - f_mode_hz   # 221 MHz

# "Scope" simulation time step
dt_scope = 0.2e-9   # 0.2 ns = 5 GS/s equivalent
T = 2*tau_g + 2*tau_d + t_c
t = np.arange(0, T, dt_scope)

# ============================================================
# Smooth functions
# ============================================================
def smooth_delta_down(t, tau_d, delta_min_hz, delta_max_hz, j=3):
    b = delta_max_hz**(-j)
    c = (2.0 / tau_d) * (delta_min_hz**(-j) - delta_max_hz**(-j))
    g = 0.5*t - (tau_d/(4.0*np.pi))*np.sin(2.0*np.pi*t/tau_d)
    return (b + c*g)**(-1.0/j)

def smooth_delta_up(t, tau_d, delta_min_hz, delta_max_hz, j=3):
    return smooth_delta_down(tau_d - t, tau_d, delta_min_hz, delta_max_hz, j=j)

def sin2_amp_envelope(t, tau_g, amp_max_frac):
    x = np.clip(t / tau_g, 0.0, 1.0)
    return amp_max_frac * np.sin(0.5 * np.pi * x)**2

# ============================================================
# Build AWG command sequence
# ============================================================
n_steps = int(np.ceil(T / dt_awg))
t_steps = np.arange(n_steps) * dt_awg

delta_cmd = np.zeros_like(t_steps)
amp_cmd   = np.zeros_like(t_steps)

for i, ti in enumerate(t_steps):
    if ti < tau_g:
        delta_cmd[i] = delta_max_hz
    elif ti < tau_g + tau_d:
        delta_cmd[i] = smooth_delta_down(
            ti - tau_g, tau_d, delta_min_hz, delta_max_hz, j=j
        )
    elif ti < tau_g + tau_d + t_c:
        delta_cmd[i] = delta_min_hz
    elif ti < tau_g + tau_d + t_c + tau_d:
        delta_cmd[i] = smooth_delta_up(
            ti - (tau_g + tau_d + t_c), tau_d, delta_min_hz, delta_max_hz, j=j
        )
    else:
        delta_cmd[i] = delta_max_hz

    if ti < tau_g:
        amp_cmd[i] = sin2_amp_envelope(ti, tau_g, amp_max_frac)
    elif ti < tau_g + tau_d + t_c + tau_d:
        amp_cmd[i] = amp_max_frac
    else:
        x = ti - (tau_g + tau_d + t_c + tau_d)
        amp_cmd[i] = sin2_amp_envelope(tau_g - x, tau_g, amp_max_frac)

# ============================================================
# Zero-order hold from AWG step grid to scope grid
# ============================================================
idx = np.minimum((t / dt_awg).astype(int), len(t_steps) - 1)
delta_t = delta_cmd[idx]
amp_t   = amp_cmd[idx]

f_blue_t = f_blue_ref_hz + delta_t
f_red_t  = f_red_ref_hz  - delta_t
f_beat_t = 0.5 * (f_blue_t - f_red_t)   # = 1 MHz + delta(t)

# ============================================================
# Phase accumulation
# ============================================================
phi_blue = 2*np.pi * np.cumsum(f_blue_t * dt_scope)
phi_red  = 2*np.pi * np.cumsum(f_red_t  * dt_scope)

phi_blue -= phi_blue[0]
phi_red  -= phi_red[0]

# ============================================================
# AWG output
# ============================================================
v_out = amp_t * np.sin(phi_blue) + amp_t * np.sin(phi_red)

# Slow envelope proxy
envelope_proxy = 2 * amp_t * np.cos(0.5 * (phi_blue - phi_red))

# ============================================================
# Scope-style display helper:
# compress into horizontal pixel bins and plot min/max
# ============================================================
def make_scope_envelope(t, y, n_bins=1000):
    edges = np.linspace(t[0], t[-1], n_bins + 1)
    t_mid = 0.5 * (edges[:-1] + edges[1:])
    y_min = np.empty(n_bins)
    y_max = np.empty(n_bins)

    start = 0
    N = len(t)

    for i in range(n_bins):
        left = edges[i]
        right = edges[i + 1]

        while start < N and t[start] < left:
            start += 1
        stop = start
        while stop < N and t[stop] < right:
            stop += 1

        if stop > start:
            ys = y[start:stop]
            y_min[i] = ys.min()
            y_max[i] = ys.max()
        else:
            y_min[i] = np.nan
            y_max[i] = np.nan

    return t_mid, y_min, y_max

# ============================================================
# Choose a scope window like your photo
# Example: a few microseconds total span
# ============================================================
t_scope_start = 5e-6
t_scope_stop  = 8e-6
m_scope = (t >= t_scope_start) & (t <= t_scope_stop)

t_scope_view = t[m_scope]
v_scope_view = v_out[m_scope]

t_bin, v_min, v_max = make_scope_envelope(t_scope_view, v_scope_view, n_bins=1200)

# ============================================================
# Plotting
# ============================================================
fig, ax = plt.subplots(5, 1, figsize=(11, 13), sharex=False)

# 1) Commanded detuning
ax[0].step(t_steps * 1e6, delta_cmd / 1e3, where='post')
ax[0].set_ylabel("Detuning (kHz)")
ax[0].set_title("Commanded detuning")

# 2) Beat frequency
ax[1].step(t * 1e6, f_beat_t / 1e6, where='mid')
ax[1].set_ylabel("Beat freq (MHz)")
ax[1].set_title("Expected slow modulation frequency = 1 MHz + delta(t)")

# 3) Commanded amplitude
ax[2].step(t_steps * 1e6, amp_cmd, where='post')
ax[2].set_ylabel("Amp (frac)")
ax[2].set_title("Commanded amplitude envelope")

# 4) Full raw waveform over full gate
decim_full = max(1, len(t) // 30000)
ax[3].plot(t[::decim_full] * 1e6, v_out[::decim_full], lw=0.6)
ax[3].set_ylabel("V_out")
ax[3].set_title("Full simulated AWG output")

# 5) Scope-like display
ax[4].fill_between(t_bin * 1e6, v_min, v_max, alpha=0.9)
ax[4].plot(t_bin * 1e6, 0.5*(v_min + v_max), lw=0.5)
ax[4].set_xlabel("Time (us)")
ax[4].set_ylabel("Displayed V")
ax[4].set_title("Scope-like min/max display over selected window")

plt.tight_layout()
plt.show()