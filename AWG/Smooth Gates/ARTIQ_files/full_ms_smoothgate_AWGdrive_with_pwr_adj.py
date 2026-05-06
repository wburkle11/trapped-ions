"""
Symmetric bichromatic smooth detuning ramp for a Spectrum DDS AWG.

This module is intended to be CALLED by an ARTIQ experiment.
It programs two DDS tones:

    blue tone: f_blue(t) = f_blue_ref + delta(t)
    red  tone: f_red(t)  = f_red_ref  - delta(t)

so the tones stay symmetric about the carrier frequency.
"""

import numpy as np
import spcm
from spcm import units


def smooth_delta_down(t, tau_d, delta_min_hz, delta_max_hz, j=3):
    """
    Smooth detuning ramp from delta_max -> delta_min over 0 <= t <= tau_d.
    All deltas are in Hz.
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
    Smooth amplitude envelope using sin^2.
    Returns amplitude in fractional units [0, amp_max_frac].
    """
    x = np.clip(t / tau_g, 0.0, 1.0)
    return amp_max_frac * np.sin(0.5 * np.pi * x) ** 2


def program_smooth_gate_awg(
    dt,
    tau_g,
    tau_d,
    t_c,
    j,
    delta_max_hz,
    delta_min_hz,
    amp_max_frac,
    f_red_ref_hz,
    f_carrier_hz,
    f_blue_ref_hz,
    red_amp_scale=1.0,
    blue_amp_scale=1.0,
    output_amp_mV=500.0,
    trigger_mode="ext",
):
    """
    Program the Spectrum DDS AWG for one smooth-gate waveform.

    Parameters
    ----------
    dt : float
        DDS update interval [s]
    tau_g : float
        amplitude ramp time [s]
    tau_d : float
        detuning ramp time [s]
    t_c : float
        center hold time at delta_min [s]
    j : int
        smoothness exponent
    delta_max_hz : float
        max detuning magnitude [Hz]
    delta_min_hz : float
        min detuning magnitude [Hz]
    amp_max_frac : float
        overall max fractional DDS amplitude [0, 1]. The larger of the red/blue
        balanced tones reaches this value.
    red_amp_scale : float
        dimensionless red-tone scale factor, usually in [0, 1]
    blue_amp_scale : float
        dimensionless blue-tone scale factor, usually in [0, 1]
    f_red_ref_hz : float
        red-side reference frequency [Hz]
    f_carrier_hz : float
        carrier frequency [Hz]
    f_blue_ref_hz : float
        blue-side reference frequency [Hz]
    output_amp_mV : float
        analog output range in mV
    trigger_mode : str
        "ext" for external trigger from ARTIQ,
        "timer" for self-clocked testing

    Returns
    -------
    card, dds, total_gate_time_s
    """
    red_amp_scale = float(red_amp_scale)
    blue_amp_scale = float(blue_amp_scale)

    def blue_amp_percent(amp_frac):
        return 100.0 * amp_frac * blue_amp_scale

    def red_amp_percent(amp_frac):
        return 100.0 * amp_frac * red_amp_scale

    def f_blue(delta_hz):
        return f_blue_ref_hz + delta_hz

    def f_red(delta_hz):
        return f_red_ref_hz - delta_hz

    total_gate_time_s = 2.0 * tau_g + 2.0 * tau_d + t_c

    card = spcm.Card(card_type=spcm.SPCM_TYPE_AO)
    card.__enter__()

    card.card_mode(spcm.SPC_REP_STD_DDS)

    channels = spcm.Channels(card)
    channels.enable(True)
    channels.output_load(50 * units.ohm)
    channels.amp(output_amp_mV * units.mV)

    card.write_setup()

    dds = spcm.DDS(card, channels=channels)
    dds.reset()

    dds.freq_ramp_stepsize(1)
    dds.amp_ramp_stepsize(1)

    blue = dds[0]
    red = dds[1]

    blue.amp(0 * units.percent)
    red.amp(0 * units.percent)
    blue.freq(f_blue(delta_max_hz) * units.Hz)
    red.freq(f_red(delta_max_hz) * units.Hz)

    if trigger_mode == "ext":
        dds.trg_timer(dt * units.s)
        dds.trg_src(spcm.SPCM_DDS_TRG_SRC_TIMER)
        dds.exec_at_trg()

    elif trigger_mode == "timer":
        dds.trg_src(spcm.SPCM_DDS_TRG_SRC_TIMER)
        dds.trg_timer(dt * units.s)
        dds.exec_at_trg()

    else:
        raise ValueError(f"Unsupported trigger_mode: {trigger_mode}")

    # Step 1: amplitude ramp on at delta_max
    n_g = max(2, int(np.ceil(tau_g / dt)))
    tvals_g = np.linspace(0.0, tau_g, n_g, endpoint=True)

    for t in tvals_g:
        amp_frac = sin2_amp_envelope(t, tau_g, amp_max_frac)

        blue.freq(f_blue(delta_max_hz) * units.Hz)
        red.freq(f_red(delta_max_hz) * units.Hz)

        blue.amp(blue_amp_percent(amp_frac) * units.percent)
        red.amp(red_amp_percent(amp_frac) * units.percent)
        dds.exec_at_trg()

    # Step 2: ramp inward
    n_d = max(2, int(np.ceil(tau_d / dt)))
    tvals_d = np.linspace(0.0, tau_d, n_d, endpoint=True)

    for t in tvals_d:
        delta_t = smooth_delta_down(
            t=t,
            tau_d=tau_d,
            delta_min_hz=delta_min_hz,
            delta_max_hz=delta_max_hz,
            j=j,
        )

        blue.freq(f_blue(delta_t) * units.Hz)
        red.freq(f_red(delta_t) * units.Hz)

        blue.amp(blue_amp_percent(amp_max_frac) * units.percent)
        red.amp(red_amp_percent(amp_max_frac) * units.percent)
        dds.exec_at_trg()

    # Step 3: hold at delta_min
    if t_c > 0:
        n_c = max(1, int(np.ceil(t_c / dt)))
        for _ in range(n_c):
            blue.freq(f_blue(delta_min_hz) * units.Hz)
            red.freq(f_red(delta_min_hz) * units.Hz)

            blue.amp((100.0 * amp_max_frac) * units.percent)
            red.amp((100.0 * amp_max_frac) * units.percent)
            dds.exec_at_trg()

    # Step 4: ramp outward
    for t in tvals_d:
        delta_t = smooth_delta_up(
            t=t,
            tau_d=tau_d,
            delta_min_hz=delta_min_hz,
            delta_max_hz=delta_max_hz,
            j=j,
        )

        blue.freq(f_blue(delta_t) * units.Hz)
        red.freq(f_red(delta_t) * units.Hz)

        blue.amp(blue_amp_percent(amp_max_frac) * units.percent)
        red.amp(red_amp_percent(amp_max_frac) * units.percent)
        dds.exec_at_trg()

    # Step 5: amplitude ramp off at delta_max
    for t in tvals_g:
        amp_frac = sin2_amp_envelope(tau_g - t, tau_g, amp_max_frac)

        blue.freq(f_blue(delta_max_hz) * units.Hz)
        red.freq(f_red(delta_max_hz) * units.Hz)

        blue.amp(blue_amp_percent(amp_frac) * units.percent)
        red.amp(red_amp_percent(amp_frac) * units.percent)
        dds.exec_at_trg()

    dds.write_to_card()

    if trigger_mode == "ext":
        card.start(spcm.M2CMD_CARD_START | spcm.M2CMD_CARD_ENABLETRIGGER)
    elif trigger_mode == "timer":
        card.start(spcm.M2CMD_CARD_START | spcm.M2CMD_CARD_ENABLETRIGGER)

    return card, dds, total_gate_time_s


def close_awg(card):
    """
    Cleanly close the AWG card handle if it exists.
    """
    if card is not None:
        try:
            card.__exit__(None, None, None)
        except Exception:
            pass
