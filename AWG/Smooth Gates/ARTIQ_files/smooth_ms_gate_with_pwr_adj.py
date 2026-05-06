from artiq.experiment import *
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

# Import lab code
from subartiq_lib.experiment.utilities import scan_to_array, us_to_ns
from subartiq_lib.analysis.analysis import compute_rsb_omegas
from subartiq_lib.analysis.fits import estimate_rabi_fit_params
from subartiq_lib.analysis.functions import rabi_cosine_exp_decay, blackman_window, tukey_window
from subartiq_lib.experiment.calibrations import rf_power_mw_to_phaser_amp_single_tone, phaser_amp_single_tone_to_rf_power_mw
from subartiq_lib.experiment.prepare import prepare_red_blue_phaser_amps, calculate_sbc_time_mu
from full_ms_smoothgate_AWGdrive import program_smooth_gate_awg, close_awg


class MolmerSorensen(EnvExperiment):
    """Smooth MS Gate"""

    # Import heavily repeated code
    from subartiq_lib.experiment.build import \
        build_sbc_defaults, build_meas_basis_defaults, import_hardware
    from subartiq_lib.experiment.experiment import \
        (import_exp_parameters, setup_369_aoms_defaults,
         doppler_cool, detect_cool, doppler_cool_leave_on, phaser_setup, phaser_multi_ch_setup, import_phaser_duc_freqs,
         set_phaser_ch_frequency, clear_phase_accum, raman_pulse,
         raman_pulse_outer_arm, raman_pulse_inner_arm,
         set_phaser_osc_frequency,
         setup_camera_grabber_rois, camera_detect_cool,
         align_timeline_to_frame, phaser_sbc_sched_multi_mode,
         set_phaser_oscs_amp_phase_clr, raman_pulse)

    def build(self):
        # Hardware setup
        self.import_hardware()
        self.import_exp_parameters()

        # --------------------------------------------------
        # Required TTL devices
        # Replace "ttlX" with your actual device-db names
        # --------------------------------------------------
        #self.setattr_device("AWG_trig")
        #self.setattr_device("outer_arm_rf_switch")

        # --------------------------------------------------
        # General experiment options
        # --------------------------------------------------
        self.setattr_argument("use_camera", BooleanValue(False))
        self.setattr_argument("use_smooth_gate", BooleanValue(True))

        # RF-switch
        self.setattr_argument("use_rf_switch", BooleanValue(True))
        self.setattr_argument(
            "rf_switch_settle_us",
            NumberValue(1.0, precision=3, step=0.1, min=0.0)
        )

        self.setattr_argument(
            "phaser_duc_freq",
            NumberValue(
                self.phaser_duc_freq_default[0] / 1e6,
                precision=6,
                step=1,
                min=self.phaser_trf_freq / 1e6 - 200,
                max=self.phaser_trf_freq / 1e6 + 200
            )
        )

        # --------------------------------------------------
        # Smooth-gate scan parameter
        # ms_pulse_times is now interpreted as t_c
        # --------------------------------------------------
        self.setattr_argument(
            "ms_pulse_times",
            Scannable(
                default=RangeScan(1, 500, 20, randomize=False),
                global_min=0.01,
                precision=2
            )
        )
        self.setattr_argument("n_reps", NumberValue(100, precision=0, step=10))

        # --------------------------------------------------
        # Frequencies shared by regular MS and smooth gate
        # --------------------------------------------------
        self.setattr_argument("rsb_aom_freq",
                              NumberValue(220.95, precision=6, step=0.1, max=250, min=180))
        self.setattr_argument("carrier_aom_freq",
                              NumberValue(221.99562, precision=6, step=10, max=260, min=200))
        self.setattr_argument("bsb_aom_freq",
                              NumberValue(223.05, precision=6, step=0.1, max=250, min=190))

        # --------------------------------------------------
        # Static MS/phaser parameters
        # --------------------------------------------------
        self.setattr_argument("phaser_amp",
                              NumberValue(self.phaser_amp_max,
                                          max=self.phaser_amp_max,
                                          min=0.0, precision=5, step=1))

        # --------------------------------------------------
        # Smooth-gate waveform parameters
        # --------------------------------------------------
        self.setattr_argument("smooth_dt_us",
                              NumberValue(2.0, precision=3, step=0.1, min=0.01))
        self.setattr_argument("smooth_tau_g_us",
                              NumberValue(20.0, precision=3, step=1.0, min=0.01))
        self.setattr_argument("smooth_tau_d_us",
                              NumberValue(100.0, precision=3, step=1.0, min=0.01))
        self.setattr_argument("smooth_j",
                              NumberValue(3, precision=0, step=1, min=1))
        self.setattr_argument("smooth_delta_max_khz",
                              NumberValue(300.0, precision=3, step=1.0, min=0.001))
        self.setattr_argument("smooth_delta_min_khz",
                              NumberValue(25.0, precision=3, step=1.0, min=0.001))
        self.setattr_argument("smooth_amp_max_frac",
                              NumberValue(0.40, precision=4, step=0.01, min=0.0, max=1.0))

        # Same red/blue tone-balancing structure used in the normal MS gate.
        # pwr_adj splits the total RF power between the red and blue tones.
        self.setattr_argument("smooth_pwr_adj",
                              NumberValue(0.49, precision=6, step=0.00001, max=1, min=0))

        # Measurement basis
        self.setattr_argument("measurement_basis",
                              EnumerationValue(['x', 'y', 'z'], 'z'),
                              group='measurement-basis')
        self.setattr_argument("carrier_pi_half_pulse_time",
                              NumberValue(1.2, precision=4, step=1),
                              group='measurement-basis')

        self.build_sbc_defaults()

    def prepare(self):
        # --------------------------------------------------
        # Shared setup / conversions
        # --------------------------------------------------
        self.phaser_duc_freq_default = np.array(self.import_phaser_duc_freqs())
        print("phaser_duc_freq_default =", self.phaser_duc_freq_default)

        # ms_pulse_times is now interpreted as the scanned hold time t_c [us]
        self.ms_pulse_times = scan_to_array(self.ms_pulse_times)
        self.smooth_t_c_us = np.array(self.ms_pulse_times, dtype=float)

        self.n_steps = len(self.smooth_t_c_us)
        self.carrier_pi_half_pulse_time_mu = us_to_ns(self.carrier_pi_half_pulse_time)
        self.rf_switch_settle_mu = us_to_ns(self.rf_switch_settle_us)

        # --------------------------------------------------
        # Smooth-gate timing
        # --------------------------------------------------
        self.smooth_tau_g_us = float(self.smooth_tau_g_us)
        self.smooth_tau_d_us = float(self.smooth_tau_d_us)

        self.smooth_gate_total_time_us = (
            2.0 * self.smooth_tau_g_us
            + 2.0 * self.smooth_tau_d_us
            + self.smooth_t_c_us
        )
        self.smooth_gate_total_time_mu = np.array(us_to_ns(self.smooth_gate_total_time_us), dtype=np.int64)

        # --------------------------------------------------
        # Sideband-cooling prep
        # --------------------------------------------------
        self.sbc_rsb_aom_freqs = scan_to_array(self.sbc_rsb_aom_freqs)
        self.sbc_cycles = scan_to_array(self.sbc_cycles).astype(int)
        self.sbc_pulse_times = scan_to_array(self.sbc_pulse_times)

        self.rsb_aom_detunings = self.sbc_rsb_aom_freqs * MHz - self.phaser_duc_freq * MHz

        self.sbc_pulse_schedule_mu = [[] for _ in range(len(self.sbc_rsb_aom_freqs))]
        for i in range(len(self.sbc_rsb_aom_freqs)):
            self.sbc_pulse_schedule_mu[i] = [us_to_ns(self.sbc_pulse_times[i])] * self.sbc_cycles[i]

        if len(self.sbc_rsb_aom_freqs) != len(self.sbc_cycles):
            print("Error! The number of RSB frequencies does not match the number of SBC cycles")
        if len(self.sbc_rsb_aom_freqs) != len(self.sbc_pulse_times):
            print("Error! The number of RSB frequencies does not match the number of SBC pulse times")

        # --------------------------------------------------
        # Smooth-gate red/blue amplitude balancing
        # --------------------------------------------------
        # This mirrors the normal MS-gate structure:
        #   phaser_amp-like total scale -> total RF power -> pwr_adj split -> red/blue amps
        #
        # IMPORTANT: smooth_amp_max_frac is still the overall AWG DDS amplitude scale.
        # smooth_red_amp_scale and smooth_blue_amp_scale are dimensionless relative
        # weights normalized so the larger tone remains at smooth_amp_max_frac.
        self.smooth_total_power_mw = phaser_amp_single_tone_to_rf_power_mw(self.smooth_amp_max_frac)
        self.smooth_red_amp_raw, self.smooth_blue_amp_raw = prepare_red_blue_phaser_amps(
            self.smooth_pwr_adj,
            total_power_mw=self.smooth_total_power_mw
        )

        self.smooth_amp_norm = max(self.smooth_red_amp_raw, self.smooth_blue_amp_raw)
        if self.smooth_amp_norm > 0.0:
            self.smooth_red_amp_scale = self.smooth_red_amp_raw / self.smooth_amp_norm
            self.smooth_blue_amp_scale = self.smooth_blue_amp_raw / self.smooth_amp_norm
        else:
            self.smooth_red_amp_scale = 0.0
            self.smooth_blue_amp_scale = 0.0

        print("smooth red raw amp", self.smooth_red_amp_raw)
        print("smooth blue raw amp", self.smooth_blue_amp_raw)
        print("smooth red amp scale", self.smooth_red_amp_scale)
        print("smooth blue amp scale", self.smooth_blue_amp_scale)

        self.set_dataset("smooth_ms.pwr_adj", float(self.smooth_pwr_adj), broadcast=True)
        self.set_dataset("smooth_ms.red_amp_raw", float(self.smooth_red_amp_raw), broadcast=True)
        self.set_dataset("smooth_ms.blue_amp_raw", float(self.smooth_blue_amp_raw), broadcast=True)
        self.set_dataset("smooth_ms.red_amp_scale", float(self.smooth_red_amp_scale), broadcast=True)
        self.set_dataset("smooth_ms.blue_amp_scale", float(self.smooth_blue_amp_scale), broadcast=True)

        # --------------------------------------------------
        # Common AWG waveform parameters
        # --------------------------------------------------
        self.awg_waveform_common = {
            "dt": self.smooth_dt_us * 1e-6,
            "tau_g": self.smooth_tau_g_us * 1e-6,
            "tau_d": self.smooth_tau_d_us * 1e-6,
            "j": int(self.smooth_j),
            "delta_max_hz": self.smooth_delta_max_khz * 1e3,
            "delta_min_hz": self.smooth_delta_min_khz * 1e3,
            "amp_max_frac": float(self.smooth_amp_max_frac),
            "red_amp_scale": float(self.smooth_red_amp_scale),
            "blue_amp_scale": float(self.smooth_blue_amp_scale),

            "f_red_ref_hz": self.rsb_aom_freq * 1e6,
            "f_carrier_hz": self.carrier_aom_freq * 1e6,
            "f_blue_ref_hz": self.bsb_aom_freq * 1e6,

            "output_amp_mV": 500.0,
            "trigger_mode": "ext",
        }

        # --------------------------------------------------
        # Save datasets for debugging / analysis
        # --------------------------------------------------
        self.set_dataset("smooth_ms.t_c_us", self.smooth_t_c_us, broadcast=True)
        self.set_dataset("smooth_ms.total_gate_time_us", self.smooth_gate_total_time_us, broadcast=True)

        self.set_dataset(
            "smooth_ms.awg_inputs",
            [
                float(self.smooth_dt_us),
                float(self.smooth_tau_g_us),
                float(self.smooth_tau_d_us),
                float(self.smooth_j),
                float(self.smooth_delta_max_khz),
                float(self.smooth_delta_min_khz),
                float(self.smooth_amp_max_frac),
                float(self.smooth_pwr_adj),
                float(self.smooth_red_amp_scale),
                float(self.smooth_blue_amp_scale),
                float(self.rsb_aom_freq),
                float(self.carrier_aom_freq),
                float(self.bsb_aom_freq),
            ],
            broadcast=True
        )

        self.set_dataset("smooth_ms.rf_switch_settle_us", float(self.rf_switch_settle_us), broadcast=True)
        self.set_dataset("smooth_ms.use_rf_switch", int(self.use_rf_switch), broadcast=True)

        self.set_dataset('smooth_ms.sbc_rsb_aom_freqs', self.sbc_rsb_aom_freqs, broadcast=True)
        self.set_dataset('smooth_ms.sbc_cycles', self.sbc_cycles, broadcast=True)
        self.set_dataset('smooth_ms.sbc_pulse_times', self.sbc_pulse_times, broadcast=True)

        # --------------------------------------------------
        # Initialize data arrays
        # --------------------------------------------------
        self.mean_pmt_counts = np.int64(0)
        self.set_dataset("smooth_ms.pmt_counts", [np.nan] * self.n_steps * self.n_reps, broadcast=True)
        self.set_dataset("smooth_ms.mean_counts", [np.nan] * self.n_steps, broadcast=True)

        # --------------------------------------------------
        # Timing estimates
        # --------------------------------------------------
        self.sbc_time_mu = calculate_sbc_time_mu(
            self.rsb_aom_detunings,
            self.sbc_cycles,
            self.sbc_pulse_schedule_mu,
            self.op_time
        )

        total_time_us = self.n_steps * self.n_reps * (
            self.dc_time
            + self.op_time
            + self.det_time
            + 2.0 * self.rf_switch_settle_us
            + 0.5 * self.smooth_gate_total_time_us[-1]
            + self.sbc_time_mu * 1e-3
        )
        print("Estimated execution time (s):", total_time_us / 1e6)

        self.timeline_correction = np.int64(320) - np.int64(
            ((self.dc_time + self.op_time) * 1e3 + self.sbc_time_mu) % 320
        )

        # --------------------------------------------------
        # Host-side AWG handles
        # --------------------------------------------------
        self._awg_card = None
        self._awg_dds = None

    @kernel
    def set_rf_switch_to_awg(self):
        
            self.outer_arm_rf_switch.off()

    @kernel
    def set_rf_switch_to_phaser(self):

            self.outer_arm_rf_switch.on()

    @kernel
    def exp_setup(self):
        self.core.break_realtime()
        delay(self.delay * us)
        self.core.reset()
        delay(self.delay * us)

        self.setup_369_aoms_defaults()

        self.phaser_multi_ch_setup(self.phaser_duc_freq_default)

        # keep only carrier oscillator for analysis pulses
        self.set_phaser_osc_frequency(
            0, 0,
            self.carrier_aom_freq * MHz - self.phaser_duc_freq_default[0]
        )
        delay_mu(self.phaser_param_change_delay_mu)

        delay_mu(self.phaser_sample_window_mu - self.phaser_param_change_delay_mu)

        # Default RF chain state at experiment start: PHASER
        if self.use_rf_switch:
            self.set_rf_switch_to_phaser()
            delay_mu(self.rf_switch_settle_mu)

    @rpc
    def program_awg_for_point(self, i):
        params = dict(self.awg_waveform_common)
        params["t_c"] = self.smooth_t_c_us[i] * 1e-6

        if self._awg_card is not None:
            close_awg(self._awg_card)
            self._awg_card = None
            self._awg_dds = None

        self._awg_card, self._awg_dds, total_gate_time_s = program_smooth_gate_awg(**params)

        self.current_smooth_gate_time_mu = np.int64(total_gate_time_s * 1e9)

    @kernel
    def run(self):
        # One-time experiment setup
        self.exp_setup()

        # Record SBC pulse schedule to DMA
        with self.core_dma.record(self.sbc_dma_handle):
            self.phaser_sbc_sched_multi_mode(0, 1, self.rsb_aom_detunings, self.sbc_pulse_schedule_mu)

        sbc_pulse_handle = self.core_dma.get_handle(self.sbc_dma_handle)

        delay(1 * s)
        print("running smooth_ms_gate")

        for i in range(self.n_steps):

            # Host-side: program AWG waveform for this scan point
            self.core.break_realtime()
            self.program_awg_for_point(i)
            self.core.break_realtime()

            for j in range(self.n_reps):

                # Align to experimental frame
                self.align_timeline_to_frame()

                # Doppler cool
                self.doppler_cool()

                # Optical pumping + clear phaser phase
                with parallel:
                    self.optical_pumping_aom.sw.pulse(self.op_time * us)
                    with sequential:
                        self.clear_phase_accum(0, [0, 1, 2, 3],
                                               [0.0, 0.0, 0.0, 0.0],
                                               [0.0, 0.0, 0.0, 0.0])
                        self.clear_phase_accum(1, [0, 1, 2, 3],
                                               [0.0, 0.0, 0.0, 0.0],
                                               [0.0, 0.0, 0.0, 0.0])

                # Pulsed sideband cooling on PHASER path
                self.core_dma.playback_handle(sbc_pulse_handle)

                # Timeline correction
                delay_mu(self.timeline_correction + 320)

                # --------------------------------------------------
                # Switch outer-arm RF chain: PHASER -> AWG
                # --------------------------------------------------
                if self.use_rf_switch:
                    self.set_rf_switch_to_awg()
                    delay_mu(self.rf_switch_settle_mu)

                # --------------------------------------------------
                # Smooth gate: AWG-triggered
                # --------------------------------------------------
                with parallel:
                    with sequential:
                        delay_mu(self.inner355_arm_delay_mu)
                        self.inner355_aom_sw.on()
                    with sequential:
                        delay_mu(self.inner355_arm_delay_mu)
                        self.AWG_trig.pulse(200 * ns)

                delay_mu(self.smooth_gate_total_time_mu[i])
                self.inner355_aom_sw.off()

                # --------------------------------------------------
                # Switch outer-arm RF chain back: AWG -> PHASER
                # --------------------------------------------------
                if self.use_rf_switch:
                    self.set_rf_switch_to_phaser()
                    delay_mu(self.rf_switch_settle_mu)

                # --------------------------------------------------
                # Measurement basis rotation
                # --------------------------------------------------
                if self.measurement_basis == 'x':
                    self.raman_pulse(0, [0], [self.phaser_amp_max],
                                     [0.75], self.carrier_pi_half_pulse_time_mu)
                elif self.measurement_basis == 'y':
                    self.raman_pulse(0, [0], [self.phaser_amp_max],
                                     [0.0], self.carrier_pi_half_pulse_time_mu)

                # --------------------------------------------------
                # Detection
                # --------------------------------------------------
                if self.use_camera:
                    self.camera_detect_cool(self.det_time)
                else:
                    pmt_counts = self.detect_cool(self.det_time)
                    self.mean_pmt_counts += pmt_counts
                    self.mutate_dataset("smooth_ms.pmt_counts", i * self.n_reps + j, pmt_counts)

            self.mutate_dataset("smooth_ms.mean_counts", i, self.mean_pmt_counts / self.n_reps)
            self.mean_pmt_counts = np.int64(0)

        # Leave system in safe/default state
        if self.use_rf_switch:
            self.set_rf_switch_to_phaser()
            delay_mu(self.rf_switch_settle_mu)

        # Leave cooling on at the end
        self.doppler_cool_leave_on()

    def analyze(self):
        pulse_times = np.array(self.get_dataset("smooth_ms.t_c_us"))
        reps = int(self.n_reps)

        if self.use_camera:
            camera_counts = np.array(self.get_dataset("smooth_ms.camera_counts"))

            for i in range(len(self.rois)):
                mean_counts = np.mean(
                    np.reshape(camera_counts[:, i], (len(pulse_times), -1)),
                    axis=1
                )
                plt.plot(pulse_times, mean_counts, label=f"roi {i}")

            plt.xlabel(r"$t_c$ ($\mu$s)")
            plt.ylabel("camera counts")
            plt.grid()
            plt.legend()
            plt.show()

            self.set_dataset("smooth_ms.mean_counts", mean_counts, broadcast=True)
            return

        pmt_counts = np.array(self.get_dataset("smooth_ms.pmt_counts"))
        mean_counts = np.mean(np.reshape(pmt_counts, (-1, reps)), axis=1)
        self.set_dataset("smooth_ms.mean_counts", mean_counts, broadcast=True)

        amp_est, bg_est, freq_est, phase_est = estimate_rabi_fit_params(pulse_times, mean_counts)
        decay_est = pulse_times[-1]

        guess = [amp_est, bg_est, freq_est, decay_est, phase_est]
        popt, pcov = curve_fit(
            rabi_cosine_exp_decay,
            pulse_times,
            mean_counts,
            p0=guess,
            absolute_sigma=True
        )

        fit_cos_decay = rabi_cosine_exp_decay(pulse_times, *popt)

        self.set_dataset("smooth_ms.fit", fit_cos_decay, broadcast=True)
        self.set_dataset("smooth_ms.fit_params", popt, broadcast=True)
        self.set_dataset("smooth_ms.fit_cov", pcov, broadcast=True)

        print("fit parameters", popt)
        print("fit uncertainties", np.sqrt(np.diag(pcov)))
        print("Oscillation frequency (kHz)",
              np.round(popt[2] / 2 / np.pi * 1000, 6), "+/-",
              np.round(np.sqrt(pcov[2, 2]) / 2 / np.pi * 1000, 6))
        print("Characteristic decay time (us)",
              np.round(popt[3], 5), "+/-",
              np.round(np.sqrt(pcov[3, 3]), 6))
