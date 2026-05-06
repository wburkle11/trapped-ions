# -*- coding: utf-8 -*-
"""
Extract detuning ramp from oscilloscope waveform using Hilbert transform.

This version:
1) Detects RF burst start
2) Defines the true detuning-ramp window
3) Crops a SLIGHTLY LARGER analysis window around that ramp
4) Runs Hilbert on that larger window
5) Plots the detuning ramp with a little extra real edge data

@author: iontrap
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import hilbert, savgol_filter


# ============================================================
# CSV LOADING
# ============================================================
def load_scope_csv(csv_path):
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
# ENVELOPE / START DETECTION
# ============================================================
def analytic_envelope(v):
    v0 = v - np.mean(v)
    z = hilbert(v0)
    return np.abs(z)


def smooth_array(y, window_pts=101, polyorder=2):
    n = len(y)
    w = int(window_pts)

    if w >= n:
        w = n - 1
    if w % 2 == 0:
        w -= 1
    if w < 5:
        return y.copy()

    return savgol_filter(y, window_length=w, polyorder=polyorder)


def find_rf_start_time(t, v, threshold_frac=0.40, env_smooth_us=2.0):
    env = analytic_envelope(v)

    dt = np.median(np.diff(t))
    smooth_pts = max(5, int(round(env_smooth_us * 1e-6 / dt)))
    if smooth_pts % 2 == 0:
        smooth_pts += 1

    env_smooth = smooth_array(env, window_pts=smooth_pts, polyorder=2)

    thresh = threshold_frac * np.max(env_smooth)
    above = np.where(env_smooth > thresh)[0]

    if len(above) == 0:
        raise ValueError("Could not detect RF start from envelope threshold.")

    i0 = above[0]
    t_rf_start = t[i0]

    return t_rf_start, env, env_smooth, thresh


# ============================================================
# WINDOWING
# ============================================================
def crop_to_time_window(t, v, t_start, duration):
    t_stop = t_start + duration
    mask = (t >= t_start) & (t <= t_stop)

    if np.count_nonzero(mask) < 20:
        raise ValueError("Too few points in cropped window.")

    t_win = t[mask]
    v_win = v[mask]

    return t_win, v_win, mask


# ============================================================
# INSTANTANEOUS FREQUENCY
# ============================================================
def instantaneous_frequency_hilbert(t, v, smooth_window=301, polyorder=3):
    v0 = v - np.mean(v)

    z = hilbert(v0)
    env = np.abs(z)
    phi = np.unwrap(np.angle(z))

    omega_inst = np.gradient(phi, t)
    f_inst = omega_inst / (2 * np.pi)

    n = len(f_inst)
    w = int(smooth_window)
    if w >= n:
        w = n - 1
    if w % 2 == 0:
        w -= 1
    if w >= 5:
        f_inst = savgol_filter(f_inst, window_length=w, polyorder=polyorder)

    return env, phi, f_inst


def choose_reference_frequency(f_inst, method="middle"):
    n = len(f_inst)
    if n < 10:
        raise ValueError("Too few points to define reference frequency.")

    if method == "median":
        return np.median(f_inst)

    if method == "start":
        m = max(5, n // 10)
        return np.mean(f_inst[:m])

    if method == "middle":
        i1 = n * 4 // 10
        i2 = n * 6 // 10
        return np.mean(f_inst[i1:i2])

    raise ValueError(f"Unknown method '{method}'")
    
def find_left_peak_time_us(t_rel_us, delta_win, search_window_us=25.0, smooth_pts=21):
    """
    Find the left-hand local maximum of the detuning curve near the start.

    Parameters
    ----------
    t_rel_us : ndarray
        Time axis in microseconds, referenced to the current detuning start.
    delta_win : ndarray
        Detuning in Hz.
    search_window_us : float
        Only search for the left peak in 0 <= t <= search_window_us.
    smooth_pts : int
        Small smoothing window to reduce noise before peak finding.

    Returns
    -------
    t_peak_us : float
        Time of the left local maximum in microseconds.
    """
    mask = (t_rel_us >= 0.0) & (t_rel_us <= search_window_us)
    if np.count_nonzero(mask) < 5:
        raise ValueError("Not enough points in left-peak search window.")

    t_search = t_rel_us[mask]
    d_search = delta_win[mask]

    # smooth slightly before finding the maximum
    n = len(d_search)
    w = int(smooth_pts)
    if w >= n:
        w = n - 1
    if w % 2 == 0:
        w -= 1
    if w >= 5:
        d_search_smooth = savgol_filter(d_search, window_length=w, polyorder=2)
    else:
        d_search_smooth = d_search.copy()

    i_peak = np.argmax(d_search_smooth)
    t_peak_us = t_search[i_peak]

    return t_peak_us

def find_first_target_crossing_time_us(t_rel_us, delta_win,
                                       target_khz=300.0,
                                       search_window_us=30.0,
                                       smooth_pts=21):
    """
    Find the first time the detuning reaches the target value.

    Parameters
    ----------
    t_rel_us : ndarray
        Time axis in microseconds, referenced to current detuning start.
    delta_win : ndarray
        Detuning in Hz.
    target_khz : float
        Target detuning value in kHz.
    search_window_us : float
        Only search for the crossing near the left edge.
    smooth_pts : int
        Small smoothing window to reduce noise before crossing detection.

    Returns
    -------
    t_cross_us : float
        First time (in us) that the smoothed detuning reaches target_khz.
    """
    mask = (t_rel_us >= 0.0) & (t_rel_us <= search_window_us)
    if np.count_nonzero(mask) < 5:
        raise ValueError("Not enough points in target-crossing search window.")

    t_search = t_rel_us[mask]
    d_search_khz = delta_win[mask] / 1e3

    # Smooth slightly before finding crossing
    n = len(d_search_khz)
    w = int(smooth_pts)
    if w >= n:
        w = n - 1
    if w % 2 == 0:
        w -= 1
    if w >= 5:
        d_smooth = savgol_filter(d_search_khz, window_length=w, polyorder=2)
    else:
        d_smooth = d_search_khz.copy()

    target = float(target_khz)

    # Find first point where detuning is at or above target
    idx = np.where(d_smooth >= target)[0]
    if len(idx) == 0:
        raise ValueError(f"Did not find crossing of {target_khz:.3f} kHz in search window.")

    i = idx[0]

    # If crossing is at very first point, just return that time
    if i == 0:
        return t_search[0]

    # Linear interpolation between neighboring points
    t1, t2 = t_search[i - 1], t_search[i]
    y1, y2 = d_smooth[i - 1], d_smooth[i]

    if y2 == y1:
        return t2

    frac = (target - y1) / (y2 - y1)
    t_cross_us = t1 + frac * (t2 - t1)

    return t_cross_us


# ============================================================
# PLOTTING
# ============================================================
def plot_rf_detection(t, v, env, env_smooth, thresh, t_rf_start, t_det_start, t_det_stop,
                      t_analysis_start, t_analysis_stop):
    t_us = 1e6 * t

    plt.figure(figsize=(10, 6))
    plt.plot(t_us, v, lw=1.0, label="Waveform")
    plt.plot(t_us, env_smooth, lw=2.0, label="Smoothed envelope")
    plt.axhline(thresh, color="k", ls=":", label="Envelope threshold")
    plt.axvline(1e6 * t_rf_start, color="r", ls="--", label="RF start")
    plt.axvline(1e6 * t_det_start, color="g", ls="--", label="Detuning start")
    plt.axvline(1e6 * t_det_stop, color="m", ls="--", label="Detuning end")
    plt.axvline(1e6 * t_analysis_start, color="c", ls=":", label="Analysis start")
    plt.axvline(1e6 * t_analysis_stop, color="c", ls=":", label="Analysis end")
    plt.xlabel("Scope time (µs)")
    plt.ylabel("Signal / Envelope")
    plt.title("Detected RF Burst and Analysis Window")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_detuning_diagnostics(t_rel_us, v_win, env_win, f_inst_win, delta_win):
    f_khz = 1e-3 * f_inst_win
    delta_khz = 1e-3 * delta_win

    fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(t_rel_us, v_win, lw=1)
    axes[0].set_ylabel("Voltage (V)")
    axes[0].set_title("Detuning Analysis Window Diagnostics")

    axes[1].plot(t_rel_us, env_win, lw=1)
    axes[1].set_ylabel("Envelope")

    axes[2].plot(t_rel_us, f_khz, lw=1.5)
    axes[2].set_ylabel("Inst. freq (kHz)")

    axes[3].plot(t_rel_us, delta_khz, lw=2.2)
    axes[3].set_ylabel("Detuning (kHz)")
    axes[3].set_xlabel("time relative to detuning-ramp start (µs)")

    for ax in axes:
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def plot_detuning_ramp_only(t_rel_us, delta_win, plot_left_us=-5.0, plot_right_us=305.0):
    delta_khz = 1e-3 * delta_win

    plt.figure(figsize=(8, 4.5))
    plt.plot(t_rel_us, delta_khz, lw=3)

    plt.xlim(plot_left_us, plot_right_us)

    plt.xlabel("time (µs)")
    plt.ylabel(r"$\delta(t)$ (kHz)")
    plt.title("Extracted Detuning Ramp")
    plt.grid(True, alpha=0.25)

    ymin, ymax = plt.ylim()
    plt.ylim(bottom=min(0, ymin))

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

    tau_g = 20e-6
    tau_d = 100e-6
    t_c   = 100e-6

    # true detuning-ramp duration
    detuning_duration = 2 * tau_d + t_c   # 300 us

    # extra data on BOTH sides for Hilbert analysis
    analysis_pad = 20e-6   # try 10 us; can increase to 15 us if needed

    # final displayed range relative to detuning start
    plot_left_us = -20.0
    plot_right_us = 315.0

    start_threshold_frac = 0.20
    env_smooth_us = 2.0

    smooth_window = 301
    polyorder = 3

    known_reference_hz = 1e6
    reference_method = "middle"

    # --------------------------------------------------------
    # LOAD SCOPE DATA
    # --------------------------------------------------------
    t, v, df = load_scope_csv(csv_path)

    dt = np.median(np.diff(t))
    fs = 1.0 / dt

    print(f"Loaded {len(t)} points")
    print(f"Sample spacing dt = {dt:.3e} s")
    print(f"Sample rate fs = {fs/1e6:.3f} MHz")

    # --------------------------------------------------------
    # FIND RF BURST START
    # --------------------------------------------------------
    t_rf_start, env_full, env_smooth_full, env_thresh = find_rf_start_time(
        t, v,
        threshold_frac=start_threshold_frac,
        env_smooth_us=env_smooth_us,
    )

    # true detuning ramp
    t_det_start = t_rf_start + tau_g
    t_det_stop = t_det_start + detuning_duration

    # widened analysis window
    t_analysis_start = t_det_start - analysis_pad
    t_analysis_stop = t_det_stop + analysis_pad
    analysis_duration = t_analysis_stop - t_analysis_start

    print(f"Detected RF start         : {t_rf_start * 1e6:.3f} us")
    print(f"Detuning ramp start used  : {t_det_start * 1e6:.3f} us")
    print(f"Detuning ramp end used    : {t_det_stop * 1e6:.3f} us")
    print(f"Analysis start used       : {t_analysis_start * 1e6:.3f} us")
    print(f"Analysis end used         : {t_analysis_stop * 1e6:.3f} us")

    plot_rf_detection(
        t, v, env_full, env_smooth_full, env_thresh,
        t_rf_start, t_det_start, t_det_stop,
        t_analysis_start, t_analysis_stop
    )

    # --------------------------------------------------------
    # CROP TO WIDER ANALYSIS WINDOW
    # --------------------------------------------------------
    t_win, v_win, win_mask = crop_to_time_window(
        t, v,
        t_start=t_analysis_start,
        duration=analysis_duration,
    )

    # --------------------------------------------------------
    # HILBERT ANALYSIS ON WIDER WINDOW
    # --------------------------------------------------------
    env_win, phi_win, f_inst_win = instantaneous_frequency_hilbert(
        t_win, v_win,
        smooth_window=smooth_window,
        polyorder=polyorder,
    )

    # --------------------------------------------------------
    # DEFINE DETUNING
    # --------------------------------------------------------
    if known_reference_hz is not None:
        f_ref = known_reference_hz
    else:
        f_ref = choose_reference_frequency(f_inst_win, method=reference_method)

    # Use this sign convention if the waveform is programmed below the reference
    delta_win = f_inst_win - f_ref

    # --------------------------------------------------------
    # REDEFINE TIME RELATIVE TO TRUE DETUNING START
    # --------------------------------------------------------
    t_rel = t_win - t_det_start
    t_rel_us = 1e6 * t_rel
    
    # --------------------------------------------------------
    # ANCHOR TIME AXIS TO FIRST TIME DETUNING REACHES TARGET
    # --------------------------------------------------------
    target_detuning_khz = 298.0

    t_target_us = find_first_target_crossing_time_us(
        t_rel_us,
        delta_win,
        target_khz=target_detuning_khz,
        search_window_us=30.0,
        smooth_pts=21
    )

    t_anchor_us = t_rel_us - t_target_us

    print(f"First crossing of {target_detuning_khz:.3f} kHz occurs at: {t_target_us:.3f} us")
    
    
    # PLOTS
    # --------------------------------------------------------
    plot_detuning_diagnostics(t_anchor_us, v_win, env_win, f_inst_win, delta_win)

    plot_detuning_ramp_only(
        t_anchor_us,
        delta_win,
        plot_left_us=0.0,
        plot_right_us=305.0
    )