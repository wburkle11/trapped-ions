# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 09:46:28 2026

@author: iontrap
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import hilbert, savgol_filter
from pathlib import Path


def load_scope_csv(csv_path):
    csv_path = r"Z:\Users\Wes\AWG Data\smooth gate data\smooth_gate_ramp_data_1.csv"

    # Read raw lines first
    with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # Find the first line that looks like actual numeric data
    data_start = None
    for i, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue

        # Try comma-separated numeric pair
        parts = [p.strip() for p in s.split(",")]
        if len(parts) >= 2:
            try:
                float(parts[0])
                float(parts[1])
                data_start = i
                break
            except ValueError:
                pass

        # Try semicolon-separated numeric pair
        parts = [p.strip() for p in s.split(";")]
        if len(parts) >= 2:
            try:
                float(parts[0])
                float(parts[1])
                data_start = i
                break
            except ValueError:
                pass

    if data_start is None:
        raise ValueError("Could not find where numeric scope data begins in the CSV.")

    # Detect delimiter from the first numeric line
    first_data_line = lines[data_start].strip()
    if ";" in first_data_line:
        sep = ";"
    elif "\t" in first_data_line:
        sep = "\t"
    else:
        sep = ","

    df = pd.read_csv(
        csv_path,
        skiprows=data_start,
        sep=sep,
        header=None,
        engine="python",
        encoding="utf-8",
        on_bad_lines="skip",
    )

    # Keep first two columns only
    df = df.iloc[:, :2].copy()
    df.columns = ["time", "voltage"]

    # Convert to numeric
    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df["voltage"] = pd.to_numeric(df["voltage"], errors="coerce")
    df = df.dropna()

    t = df["time"].to_numpy()
    v = df["voltage"].to_numpy()

    if len(t) < 10:
        raise ValueError("Found too few numeric data points after parsing.")

    return t, v, df


def estimate_instantaneous_frequency(t, v, smooth_window=501, polyorder=3,
                                     amp_threshold_frac=0.20):
    """
    Estimate instantaneous frequency from a real-valued waveform using
    the analytic signal (Hilbert transform).

    Parameters
    ----------
    t : ndarray
        Time array [s]
    v : ndarray
        Voltage array [V]
    smooth_window : int
        Window length for Savitzky-Golay smoothing; must be odd
    polyorder : int
        Polynomial order for Savitzky-Golay filter
    amp_threshold_frac : float
        Fraction of max envelope below which frequency estimate is masked

    Returns
    -------
    env : ndarray
        Signal envelope
    phi : ndarray
        Unwrapped instantaneous phase [rad]
    f_inst : ndarray
        Instantaneous frequency [Hz]
    good_mask : ndarray
        Boolean mask where amplitude is large enough to trust f_inst
    """
    # Remove DC offset
    v0 = v - np.mean(v)

    # Analytic signal
    z = hilbert(v0)
    env = np.abs(z)
    phi = np.unwrap(np.angle(z))

    # dphi/dt -> instantaneous angular frequency
    omega_inst = np.gradient(phi, t)
    f_inst = omega_inst / (2 * np.pi)

    # Smooth the noisy derivative if possible
    if smooth_window >= len(f_inst):
        smooth_window = len(f_inst) - 1
    if smooth_window % 2 == 0:
        smooth_window -= 1
    if smooth_window >= 5:
        f_inst = savgol_filter(f_inst, window_length=smooth_window, polyorder=polyorder)

    # Mask out regions where the amplitude is too low
    env_thresh = amp_threshold_frac * np.max(env)
    good_mask = env > env_thresh

    return env, phi, f_inst, good_mask


def choose_reference_frequency(t, f_inst, good_mask, method="median"):
    """
    Choose a reference frequency to define delta(t) = f_inst - f_ref.

    method:
        'median'  -> median over trusted region
        'start'   -> average of first 10% trusted points
        'middle'  -> average of middle 20% trusted points
    """
    f_good = f_inst[good_mask]
    t_good = t[good_mask]

    if len(f_good) < 10:
        raise ValueError("Not enough trusted points to determine a reference frequency.")

    if method == "median":
        return np.median(f_good)

    elif method == "start":
        n = max(5, len(f_good) // 10)
        return np.mean(f_good[:n])

    elif method == "middle":
        n = len(f_good)
        i1 = n * 4 // 10
        i2 = n * 6 // 10
        return np.mean(f_good[i1:i2])

    else:
        raise ValueError(f"Unknown method '{method}'")


def plot_detuning_analysis(t, v, env, f_inst, delta, good_mask,
                           title="Detuning Ramp Extracted from Scope Data"):
    """
    Make a 4-panel diagnostic plot:
      1) raw waveform
      2) envelope
      3) instantaneous frequency
      4) detuning
    """
    t_us = 1e6 * t
    f_khz = 1e-3 * f_inst
    delta_khz = 1e-3 * delta

    fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(t_us, v, lw=1)
    axes[0].set_ylabel("Voltage (V)")
    axes[0].set_title(title)

    axes[1].plot(t_us, env, lw=1)
    axes[1].set_ylabel("Envelope")

    axes[2].plot(t_us[good_mask], f_khz[good_mask], lw=1.5)
    axes[2].set_ylabel("Inst. freq (kHz)")

    axes[3].plot(t_us[good_mask], delta_khz[good_mask], lw=2)
    axes[3].set_ylabel("Detuning (kHz)")
    axes[3].set_xlabel("Time (µs)")

    for ax in axes:
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def plot_stylized_detuning(t, delta, good_mask,
                           title="Extracted Detuning Ramp"):
    t_us = 1e6 * t
    delta_khz = 1e-3 * delta

    # Restrict to trusted region
    tt = t_us[good_mask]
    dd = delta_khz[good_mask]

    plt.figure(figsize=(8, 4.5))
    plt.plot(tt, dd, lw=3)
    plt.xlabel("time (µs)")
    plt.ylabel(r"$\delta(t)$ (kHz)")
    plt.title(title)
    plt.grid(True, alpha=0.25)

    # ---- minimal fix ----
    ymin, ymax = plt.ylim()
    plt.ylim(bottom=min(0, ymin))
    # ---------------------

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # ===========================
    # USER INPUTS
    # ===========================
    csv_path = "scope_data.csv"   # <-- replace with your actual CSV filename

    # How to define the reference frequency:
    # 'median' is the safest generic choice.
    reference_method = "median"

    # If you already know the intended center/reference frequency, set it here.
    # Example: known_reference_hz = 222e6
    # If None, the script estimates one automatically.
    known_reference_hz = 100e3

    # ===========================
    # LOAD DATA
    # ===========================
    t, v, df = load_scope_csv(csv_path)

    # Basic sanity check for uniform sampling
    dt = np.median(np.diff(t))
    fs = 1.0 / dt
    #print(f"Loaded {len(t)} points")
    #print(f"Estimated sample spacing dt = {dt:.3e} s")
    #print(f"Estimated sample rate fs = {fs/1e6:.3f} MHz")

    # ===========================
    # EXTRACT INSTANTANEOUS FREQUENCY
    # ===========================
    env, phi, f_inst, good_mask = estimate_instantaneous_frequency(
        t, v,
        smooth_window=301,
        polyorder=3,
        amp_threshold_frac=0.10,
    )

    # ===========================
    # DEFINE DETUNING
    # ===========================
    if known_reference_hz is not None:
        f_ref = known_reference_hz
    else:
        f_ref = choose_reference_frequency(t, f_inst, good_mask,
                                           method=reference_method)

    delta = f_inst - f_ref

    print(f"Reference frequency used: {f_ref/1e6:.6f} MHz")
    print(f"Detuning range over trusted region: "
          f"{np.min(delta[good_mask])/1e3:.2f} to {np.max(delta[good_mask])/1e3:.2f} kHz")

    # ===========================
    # PLOTS
    # ===========================
    plot_detuning_analysis(t, v, env, f_inst, delta, good_mask)

    # Cleaner figure
    plot_stylized_detuning(t, delta, good_mask)