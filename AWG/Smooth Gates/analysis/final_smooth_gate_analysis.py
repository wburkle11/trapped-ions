
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import hilbert, butter, sosfiltfilt, savgol_filter
from scipy.signal import find_peaks


# ============================================================
# USER SETTINGS
# ============================================================
csv_path = r"Z:\Users\Wes\AWG Data\smooth gate data\full_smooth_gate_data_2.csv"

carrier_center_hz = 222e6          # center frequency
expected_base_beat_hz = 1e6        # 1 MHz
expected_delta_max_hz = 0.35e6      # adjust if needed
plot_only_active_region = False    # set True if you want to crop to pulse region

# bandpass for the squared-signal beat extraction
# Squaring produces the difference frequency:
#   f_diff(t) = 2 * (1 MHz + delta(t))
# so expected band is around ~2 MHz, broadened by delta(t)
band_margin_hz = 0.25e6

# smoothing for extracted instantaneous frequency
sg_window = 301    # must be odd; auto-fixed below if needed
sg_poly = 3


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
        raise ValueError("Too few valid data points found.")

    return t, v


# ============================================================
# HELPERS
# ============================================================
def make_odd(n, minimum=5):
    n = int(n)
    if n < minimum:
        n = minimum
    if n % 2 == 0:
        n += 1
    return n


def bandpass_filter(y, fs, f_low, f_high, order=4):
    nyq = 0.5 * fs
    if f_low <= 0:
        raise ValueError("Bandpass lower edge must be > 0.")
    if f_high >= nyq:
        raise ValueError(
            f"Bandpass upper edge {f_high/1e6:.3f} MHz exceeds Nyquist {nyq/1e6:.3f} MHz."
        )
    sos = butter(order, [f_low / nyq, f_high / nyq], btype="bandpass", output="sos")
    return sosfiltfilt(sos, y)


def smooth(y, window, poly=3):
    window = make_odd(min(window, len(y) - 1))
    if window >= len(y):
        window = len(y) - 1
        if window % 2 == 0:
            window -= 1
    if window < poly + 2:
        return y.copy()
    return savgol_filter(y, window_length=window, polyorder=poly)


def estimate_active_region(v, frac=0.10):
    """
    Use Hilbert envelope of raw trace to identify the main pulse region.
    """
    v0 = v - np.mean(v)
    env = np.abs(hilbert(v0))
    env = smooth(env, 1001, poly=3)

    if np.max(env) <= 0:
        return 0, len(v)

    env_norm = env / np.max(env)
    mask = env_norm > frac

    if not np.any(mask):
        return 0, len(v)

    idx = np.where(mask)[0]
    gaps = np.where(np.diff(idx) > 1)[0]

    starts = [idx[0]] + [idx[g + 1] for g in gaps]
    stops = [idx[g] + 1 for g in gaps] + [idx[-1] + 1]

    lengths = [b - a for a, b in zip(starts, stops)]
    k = int(np.argmax(lengths))
    return starts[k], stops[k]


# ============================================================
# CORE EXTRACTION
# ============================================================
def extract_beat_frequency(t, v, expected_base_beat_hz, expected_delta_max_hz, band_margin_hz):
    """
    For a symmetric bichromatic signal with tones:
        f_blue = 222 MHz + (1 MHz + delta)
        f_red  = 222 MHz - (1 MHz + delta)

    Their separation is:
        f_diff = f_blue - f_red = 2*(1 MHz + delta)

    If we square the measured signal, this difference-frequency component
    appears explicitly. We bandpass around that component, Hilbert-transform
    it, then differentiate the analytic phase to get instantaneous f_diff(t).

    Finally:
        f_beat(t) = f_diff(t)/2 = 1 MHz + delta(t)
    """
    dt = np.median(np.diff(t))
    fs = 1.0 / dt

    v0 = v - np.mean(v)
    vsq = v0**2
    vsq -= np.mean(vsq)

    f_diff_center = 2.0 * expected_base_beat_hz
    f_diff_halfspan = 2.0 * expected_delta_max_hz + band_margin_hz

    f_low = max(1e3, f_diff_center - f_diff_halfspan)
    f_high = f_diff_center + f_diff_halfspan

    vsq_bp = bandpass_filter(vsq, fs, f_low, f_high, order=4)

    z = hilbert(vsq_bp)
    phase = np.unwrap(np.angle(z))
    f_diff_inst = np.gradient(phase, t) / (2.0 * np.pi)

    f_diff_smooth = smooth(f_diff_inst, sg_window, poly=sg_poly)
    f_beat = 0.5 * f_diff_smooth
    delta = f_beat - expected_base_beat_hz

    return {
        "fs": fs,
        "vsq": vsq,
        "vsq_bp": vsq_bp,
        "f_diff_inst": f_diff_inst,
        "f_diff_smooth": f_diff_smooth,
        "f_beat": f_beat,
        "delta": delta,
        "f_low": f_low,
        "f_high": f_high,
    }

def analyze_squared_spectrum(t, vsq, n_peaks=10):
    """
    Compute FFT of squared signal and print dominant frequency components.

    Parameters
    ----------
    t : ndarray
        Time array [s]
    vsq : ndarray
        Squared signal (mean-subtracted)
    n_peaks : int
        Number of strongest peaks to report
    """
    dt = np.median(np.diff(t))
    fs = 1.0 / dt
    N = len(vsq)

    # FFT
    freqs = np.fft.fftfreq(N, d=dt)
    spectrum = np.fft.fft(vsq)

    # Only keep positive frequencies
    mask = freqs > 0
    freqs = freqs[mask]
    power = np.abs(spectrum[mask])

    # Normalize for readability
    power /= np.max(power)

    # Find peaks
    peaks, _ = find_peaks(power, height=0.05)  # adjust threshold if needed

    peak_freqs = freqs[peaks]
    peak_powers = power[peaks]

    # Sort by strength
    idx = np.argsort(peak_powers)[::-1]

    print("\n=== Dominant frequency components in squared signal ===")
    for i in idx[:n_peaks]:
        print(f"{peak_freqs[i]/1e6:10.6f} MHz   (relative amplitude: {peak_powers[i]:.3f})")

    return freqs, power

# ============================================================
# PLOTTING
# ============================================================
def plot_waveform(t, v, title="Scope Data"):
    plt.figure(figsize=(10, 4))
    plt.plot(1e6 * t, v, lw=1.0)
    plt.xlabel("Time (µs)")
    plt.ylabel("Voltage (V)")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_squared_bandpassed(t, vsq_bp, title="Scope Data (squared, filtered)"):
    plt.figure(figsize=(10, 4))
    plt.plot(1e6 * t, vsq_bp, lw=1.0)
    plt.xlabel("Time (µs)")
    plt.ylabel("Beatnote Tone Amplitude")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_beat_and_delta(t, f_beat, delta):
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    ax[0].plot(1e6 * t, f_beat / 1e6, lw=2)
    ax[0].set_ylabel("MHz")
    ax[0].set_title(r"Extracted Beatnote")
    ax[0].grid(True, alpha=0.3)

    ax[1].plot(1e6 * t, delta / 1e3, lw=2)
    ax[1].set_xlabel("Time (µs)")
    ax[1].set_ylabel("kHz")
    ax[1].set_title(r"Detuning from Mode: $\delta(t)$")
    ax[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def plot_inferred_tones(t, f_blue, f_red):
    plt.figure(figsize=(10, 4.5))
    plt.plot(1e6 * t, f_blue / 1e6, lw=2, label="222 MHz + (1 MHz + δ(t))")
    plt.plot(1e6 * t, f_red / 1e6, lw=2, label="222 MHz - (1 MHz + δ(t))")
    plt.xlabel("Time (µs)")
    plt.ylabel("Frequency (MHz)")
    plt.title("Smoothgate RF to AOM")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t, v = load_scope_csv(csv_path)

    dt = np.median(np.diff(t))
    fs = 1.0 / dt
    nyq = 0.5 * fs


    # Optional: isolate main pulse region for cleaner plots
    if plot_only_active_region:
        i0, i1 = estimate_active_region(v, frac=0.10)
        t = t[i0:i1]
        v = v[i0:i1]
        print(f"Using active region only: indices {i0}:{i1}")

    # 1) measured waveform
    plot_waveform(t, v)

    # 2) extract beat frequency and delta(t)
    out = extract_beat_frequency(
        t=t,
        v=v,
        expected_base_beat_hz=expected_base_beat_hz,
        expected_delta_max_hz=expected_delta_max_hz,
        band_margin_hz=band_margin_hz,
    )
    
    # Analyze spectrum of squared signal
    #freqs, power = analyze_squared_spectrum(t, out["vsq"])


    f_beat = out["f_beat"]     # = 1 MHz + delta(t)
    delta = out["delta"]

    # 3) infer the two actual tones
    f_blue = carrier_center_hz + f_beat
    f_red = carrier_center_hz - f_beat

    # plots
    plot_squared_bandpassed(t, out["vsq_bp"])
    plot_beat_and_delta(t, f_beat, delta)
    plot_inferred_tones(t, f_blue, f_red)