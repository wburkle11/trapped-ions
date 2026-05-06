# smooth_ms_awg.py

import numpy as np
import spcm
from spcm import units


def smooth_delta_down(t, tau_d, delta_min_hz, delta_max_hz, j=3):
    b = delta_max_hz ** (-j)
    c = (2.0 / tau_d) * (delta_min_hz ** (-j) - delta_max_hz ** (-j))
    g = 0.5 * t - (tau_d / (4.0 * np.pi)) * np.sin(2.0 * np.pi * t / tau_d)
    return (b + c * g) ** (-1.0 / j)


def smooth_delta_up(t, tau_d, delta_min_hz, delta_max_hz, j=3):
    return smooth_delta_down(tau_d - t, tau_d, delta_min_hz, delta_max_hz, j=j)


def sin2_amp_envelope(t, tau_g, amp_max_frac):
    x = np.clip(t / tau_g, 0.0, 1.0)
    return amp_max_frac * np.sin(0.5 * np.pi * x) ** 2


def program_smooth_ms_awg(
    *,
    blue_sideband_ref_mhz,
    red_sideband_ref_mhz,
    delta_max_khz,
    delta_min_khz,
    tau_d_us,
    tau_g_us,
    hold_us,
    blue_amp_frac,
    red_amp_frac,
    dt_us=2.0,
    j=3,
    output_mv=500,
    load_ohm=50,
):
    f_blue_ref_hz = float(blue_sideband_ref_mhz) * 1e6
    f_red_ref_hz  = float(red_sideband_ref_mhz) * 1e6

    delta_max_hz = float(delta_max_khz) * 1e3
    delta_min_hz = float(delta_min_khz) * 1e3

    tau_d = float(tau_d_us) * 1e-6
    tau_g = float(tau_g_us) * 1e-6
    t_c   = float(hold_us) * 1e-6
    dt    = float(dt_us) * 1e-6

    blue_amp_frac = float(blue_amp_frac)
    red_amp_frac  = float(red_amp_frac)

    def f_blue(delta_hz):
        return f_blue_ref_hz + delta_hz

    def f_red(delta_hz):
        return f_red_ref_hz - delta_hz

    n_g = max(2, int(np.ceil(tau_g / dt)))
    n_d = max(2, int(np.ceil(tau_d / dt)))
    n_c = max(1, int(np.ceil(t_c / dt))) if t_c > 0 else 0

    tvals_g = np.linspace(0.0, tau_g, n_g, endpoint=True)
    tvals_d = np.linspace(0.0, tau_d, n_d, endpoint=True)

    card = spcm.Card(card_type=spcm.SPCM_TYPE_AO)
    try:
        card.card_mode(spcm.SPC_REP_STD_DDS)

        channels = spcm.Channels(card)
        channels.enable(True)
        channels.output_load(load_ohm * units.ohm)
        channels.amp(output_mv * units.mV)
        card.write_setup()

        dds = spcm.DDS(card, channels=channels)
        dds.reset()

        dds.trg_src(spcm.SPCM_DDS_TRG_SRC_TIMER)
        dds.trg_timer(dt * units.s)

        dds.freq_ramp_stepsize(1)
        dds.amp_ramp_stepsize(1)

        blue = dds[0]
        red = dds[1]

        blue.amp(0 * units.percent)
        red.amp(0 * units.percent)
        blue.freq(f_blue(delta_max_hz) * units.Hz)
        red.freq(f_red(delta_max_hz) * units.Hz)
        dds.exec_at_trg()

        # 1) amplitude ramp on
        for t in tvals_g:
            blue_a = sin2_amp_envelope(t, tau_g, blue_amp_frac)
            red_a  = sin2_amp_envelope(t, tau_g, red_amp_frac)

            blue.freq(f_blue(delta_max_hz) * units.Hz)
            red.freq(f_red(delta_max_hz) * units.Hz)
            blue.amp((100.0 * blue_a) * units.percent)
            red.amp((100.0 * red_a) * units.percent)
            dds.exec_at_trg()

        # 2) inward ramp
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
            blue.amp((100.0 * blue_amp_frac) * units.percent)
            red.amp((100.0 * red_amp_frac) * units.percent)
            dds.exec_at_trg()

        # 3) hold
        for _ in range(n_c):
            blue.freq(f_blue(delta_min_hz) * units.Hz)
            red.freq(f_red(delta_min_hz) * units.Hz)
            blue.amp((100.0 * blue_amp_frac) * units.percent)
            red.amp((100.0 * red_amp_frac) * units.percent)
            dds.exec_at_trg()

        # 4) outward ramp
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
            blue.amp((100.0 * blue_amp_frac) * units.percent)
            red.amp((100.0 * red_amp_frac) * units.percent)
            dds.exec_at_trg()

        # 5) amplitude ramp off
        for t in tvals_g:
            blue_a = sin2_amp_envelope(tau_g - t, tau_g, blue_amp_frac)
            red_a  = sin2_amp_envelope(tau_g - t, tau_g, red_amp_frac)

            blue.freq(f_blue(delta_max_hz) * units.Hz)
            red.freq(f_red(delta_max_hz) * units.Hz)
            blue.amp((100.0 * blue_a) * units.percent)
            red.amp((100.0 * red_a) * units.percent)
            dds.exec_at_trg()

        dds.write_to_card()
        card.start(spcm.M2CMD_CARD_ENABLETRIGGER)

        return {
            "card": card,
            "blue_amp_frac": blue_amp_frac,
            "red_amp_frac": red_amp_frac,
        }

    except Exception:
        card.close()
        raise


def close_awg(awg_info):
    card = awg_info.get("card", None)
    if card is not None:
        card.close()