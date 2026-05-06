# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 18:44:18 2026

@author: iontrap
"""
# -*- coding: utf-8 -*-
"""
Simple scope-waveform + Hilbert detuning extraction

What this script does
---------------------
1) Loads the oscilloscope CSV
2) Plots the measured waveform directly (like the scope trace)
3) Uses the Hilbert transform to extract instantaneous frequency
4) Converts that into detuning magnitude relative to a known 1 MHz reference
5) Plots detuning vs time

Expected physical behavior
--------------------------
The extracted detuning plot should show approximately:
300 kHz -> 25 kHz -> 300 kHz

@author: iontrap
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import hilbert, savgol_filter


# ============================================================
# LOAD CSV
# ============================================================
def load_scope_csv(csv_path):
    """
    Load a scope CSV that may have header text before numeric data.

    Returns
    -------
    t : ndarray
        Time array [s]
    v : ndarray
        Voltage array [V]
    df : DataFrame
        Parsed numeric data
    """
    with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    data_start = None
    sep = ","

    for i, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue

        for sep_try in [",", ";", "\t"]:
            parts = [p.strip() for p in s.split(sep_try)]
            if len(parts) >= 2:
                try:
                    float(parts[0])
                    float(parts[1])
                    data_start = i
                    sep = sep_try
                    break
                except ValueError:
                    pass
        if data_start is not None:
            break

    if data_start is None:
        raise ValueError("Could not find numeric data in CSV.")

    df = pd.read_csv(
        csv_path,
        skiprows=data_start,
        sep=sep,
        header=None,
        engine="python",
        encoding="utf-8",
        on_bad_lines="skip",
    )

    df = df.iloc[:, :2].copy()
    df.columns = ["time", "voltage"]

    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df["voltage"] = pd.to_numeric(df["voltage"], errors="coerce")
    df = df.dropna()

    t = df["time"].to_numpy()
    v = df["voltage"].to_numpy()

    if len(t) < 20:
        raise ValueError("Too few valid data points found in scope file.")

    return t, v, df


# ============================================================
# HILBERT ANALYSIS
# ============================================================
def extract_instantaneous_frequency(t, v, smooth_window=301, polyorder=3):
    """
    Extract instantaneous frequency from a real-valued waveform using
    the analytic signal from the Hilbert transform.

    Parameters
    ----------
    t : ndarray
        Time [s]
    v : ndarray
        Voltage [V]
    smooth_window : int
        Smoothing window for f_inst
    polyorder : int
        Polynomial order for Savitzky-Golay smoothing

    Returns
    -------
    env : ndarray
        Envelope
    phi : ndarray
        Unwrapped phase [rad]
    f_inst : ndarray
        Instantaneous frequency [Hz]
    """
    # Remove DC offset
    v0 = v - np.mean(v)

    # Analytic signal
    z = hilbert(v0)
    env = np.abs(z)
    phi = np.unwrap(np.angle(z))

    # Instantaneous frequency
    omega_inst = np.gradient(phi, t)
    f_inst = omega_inst / (2 * np.pi)

    # Smooth the extracted frequency a bit
    n = len(f_inst)
    w = int(smooth_window)
    if w >= n:
        w = n - 1
    if w % 2 == 0:
        w -= 1
    if w >= 5:
        f_inst = savgol_filter(f_inst, window_length=w, polyorder=polyorder)

    return env, phi, f_inst


# ============================================================
# PLOTTING
# ============================================================
def plot_scope_waveform(t, v, title="Measured Scope Waveform"):
    """
    Plot the raw scope waveform.
    """
    t_us = 1e6 * t

    plt.figure(figsize=(10, 4.5))
    plt.plot(t_us, v, lw=1.0)
    plt.xlabel("time (µs)")
    plt.ylabel("Voltage (V)")
    plt.title(title)
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.show()


def plot_detuning_ramp(t, detuning_hz, title="Extracted Detuning Ramp"):
    """
    Plot detuning magnitude vs time.
    """
    t_us = 1e6 * t
    detuning_khz = detuning_hz / 1e3

    plt.figure(figsize=(10, 4.5))
    plt.plot(t_us, detuning_khz, lw=2.5)
    plt.xlabel("time (µs)")
    plt.ylabel(r"$\delta(t)$ (kHz)")
    plt.title(title)
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.show()


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    # --------------------------------------------------------
    # USER INPUTS
    # --------------------------------------------------------
    csv_path = r"Z:\Users\Wes\AWG Data\smooth gate data\smooth_gate_ramp_data_4.csv"

    # Known reference frequency for the programmed chirp
    reference_hz = 1e6

    # Smoothing for the instantaneous-frequency trace
    smooth_window = 301
    polyorder = 3

    # --------------------------------------------------------
    # LOAD DATA
    # --------------------------------------------------------
    t, v, df = load_scope_csv(csv_path)

    dt = np.median(np.diff(t))
    fs = 1.0 / dt

    print(f"Loaded {len(t)} points")
    print(f"Sample spacing dt = {dt:.3e} s")
    print(f"Sample rate fs = {fs/1e6:.3f} MHz")

    # --------------------------------------------------------
    # PLOT RAW WAVEFORM
    # --------------------------------------------------------
    plot_scope_waveform(t, v, title="Measured Scope Waveform")

    # --------------------------------------------------------
    # HILBERT EXTRACTION
    # --------------------------------------------------------
    env, phi, f_inst = extract_instantaneous_frequency(
        t, v,
        smooth_window=smooth_window,
        polyorder=polyorder,
    )

    print(f"Instantaneous frequency range: "
          f"{np.min(f_inst)/1e6:.6f} to {np.max(f_inst)/1e6:.6f} MHz")

    # --------------------------------------------------------
    # DETUNING RELATIVE TO 1 MHz
    # --------------------------------------------------------
    # Use absolute value because the measured chirp may sit either
    # above or below the 1 MHz reference, while the physical detuning
    # magnitude is positive.
    detuning_hz = np.abs(f_inst - reference_hz)

    print(f"Extracted detuning range: "
          f"{np.min(detuning_hz)/1e3:.2f} to {np.max(detuning_hz)/1e3:.2f} kHz")

    # --------------------------------------------------------
    # PLOT DETUNING RAMP
    # --------------------------------------------------------
    plot_detuning_ramp(
        t,
        detuning_hz,
        title="Extracted Detuning Ramp Relative to 1 MHz"
    )