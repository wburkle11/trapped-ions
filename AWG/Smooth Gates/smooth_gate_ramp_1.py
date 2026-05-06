"""
Paper-inspired smooth detuning ramp for a Spectrum DDS AWG.

This approximates the smooth-gate detuning ramp from
arXiv:2510.17286 by updating the DDS carrier frequency
in many small steps.

Important:
- This is a piecewise approximation, not a truly continuous chirp.
- For a real MS gate, you will likely need TWO carriers
  (red and blue sidebands), not just one.
"""

import numpy as np
import spcm
from spcm import units


def smooth_delta_down(t, tau_d, delta_min_hz, delta_max_hz, j=3):
    """
    Smooth detuning ramp from delta_max -> delta_min over 0 <= t <= tau_d
    based on the paper's Eq. (18)-(19) form with j=3 recommended.
    All deltas in Hz.
    """
    b = delta_max_hz ** (-j)
    c = (2.0 / tau_d) * (delta_min_hz ** (-j) - delta_max_hz ** (-j))
    g = 0.5 * t - (tau_d / (4.0 * np.pi)) * np.sin(2.0 * np.pi * t / tau_d)
    return (b + c * g) ** (-1.0 / j)


def smooth_delta_up(t, tau_d, delta_min_hz, delta_max_hz, j=3):
    """
    Reverse ramp from delta_min -> delta_max over 0 <= t <= tau_d.
    """
    return smooth_delta_down(tau_d - t, tau_d, delta_min_hz, delta_max_hz, j=j)


def sin2_amp_envelope(t, tau_g, amp_max_frac):
    """
    Smooth amplitude turn-on/off envelope using sin^2.
    Returns amplitude in fractional units [0, 1].
    """
    x = np.clip(t / tau_g, 0.0, 1.0)
    return amp_max_frac * np.sin(0.5 * np.pi * x) ** 2


# -----------------------------
# User parameters
# -----------------------------

# DDS update interval
dt = 2e-6  # 2 us per command update

# Gate shape parameters
tau_g = 20e-6      # amplitude ramp-on / ramp-off time
tau_d = 100e-6     # detuning ramp-down and ramp-up time
t_c  = 100e-6         # optional hold time at delta_min; paper says tc=0 often best

# Paper-inspired choice
j = 3

# Carrier mapping:
# f_out(t) = f_center + delta(t)
# Replace f_center with whatever your AWG frequency should be referenced to.
f_center_hz = 1.0e6

delta_max_hz = 300e3
delta_min_hz = 25e3

amp_max_frac = 0.40   # 40% amplitude

# -----------------------------
# Open Spectrum card
# -----------------------------
card: spcm.Card
with spcm.Card(card_type=spcm.SPCM_TYPE_AO) as card:

    card.card_mode(spcm.SPC_REP_STD_DDS)

    channels = spcm.Channels(card)
    channels.enable(True)
    channels.output_load(50 * units.ohm)
    channels.amp(500 * units.mV)

    # Required for DDS clocking
    card.write_setup()

    dds = spcm.DDS(card, channels=channels)
    dds.reset()

    # Timer trigger source: one queued command executed every dt
    dds.trg_src(spcm.SPCM_DDS_TRG_SRC_TIMER)
    dds.trg_timer(dt * units.s)

    # Keep DDS ramp step sizes small / benign.
    # We are mostly doing explicit point-by-point updates anyway.
    dds.freq_ramp_stepsize(1)
    dds.amp_ramp_stepsize(1)

    # Use carrier 0
    carr = dds[0]

    # Start from zero output at large detuning
    carr.amp(0 * units.percent)
    carr.freq((f_center_hz + delta_max_hz) * units.Hz)
    dds.exec_at_trg()

    # -----------------------------
    # Step 1: amplitude ramp on at delta_max
    # -----------------------------
    n_g = max(2, int(np.ceil(tau_g / dt)))
    tvals = np.linspace(0.0, tau_g, n_g, endpoint=True)

    for t in tvals:
        amp_frac = sin2_amp_envelope(t, tau_g, amp_max_frac)
        carr.freq((f_center_hz + delta_max_hz) * units.Hz)
        carr.amp((100.0 * amp_frac) * units.percent)
        dds.exec_at_trg()

    # -----------------------------
    # Step 2: ramp detuning down, amplitude constant
    # -----------------------------
    n_d = max(2, int(np.ceil(tau_d / dt)))
    tvals = np.linspace(0.0, tau_d, n_d, endpoint=True)

    for t in tvals:
        delta_t = smooth_delta_down(
            t=t,
            tau_d=tau_d,
            delta_min_hz=delta_min_hz,
            delta_max_hz=delta_max_hz,
            j=j,
        )
        carr.freq((f_center_hz + delta_t) * units.Hz)
        carr.amp((100.0 * amp_max_frac) * units.percent)
        dds.exec_at_trg()

    # -----------------------------
    # Step 3: optional hold at delta_min
    # -----------------------------
    if t_c > 0:
        n_c = max(1, int(np.ceil(t_c / dt)))
        for _ in range(n_c):
            carr.freq((f_center_hz + delta_min_hz) * units.Hz)
            carr.amp((100.0 * amp_max_frac) * units.percent)
            dds.exec_at_trg()

    # -----------------------------
    # Step 4: ramp detuning back up
    # -----------------------------
    for t in tvals:
        delta_t = smooth_delta_up(
            t=t,
            tau_d=tau_d,
            delta_min_hz=delta_min_hz,
            delta_max_hz=delta_max_hz,
            j=j,
        )
        carr.freq((f_center_hz + delta_t) * units.Hz)
        carr.amp((100.0 * amp_max_frac) * units.percent)
        dds.exec_at_trg()

    # -----------------------------
    # Step 5: amplitude ramp off at delta_max
    # -----------------------------
    for t in tvals:
        amp_frac = sin2_amp_envelope(tau_g - t, tau_g, amp_max_frac)
        carr.freq((f_center_hz + delta_max_hz) * units.Hz)
        carr.amp((100.0 * amp_frac) * units.percent)
        dds.exec_at_trg()

    # Write queued command list to card
    dds.write_to_card()

    # Start DDS command execution
    card.start(spcm.M2CMD_CARD_ENABLETRIGGER, spcm.M2CMD_CARD_FORCETRIGGER)

    input("Press Enter to Exit")