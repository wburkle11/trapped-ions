# SmoothMSGate.py

from artiq.experiment import *
from smooth_ms_awg import program_smooth_ms_awg, close_awg


class SmoothMSGate(EnvExperiment):
    """Smooth MS Gate"""

    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl_awg_trigger")

        # Sideband references from scan
        self.setattr_argument(
            "blue_sideband_ref_mhz",
            NumberValue(default=223.000000, unit="MHz", ndecimals=6)
        )
        self.setattr_argument(
            "red_sideband_ref_mhz",
            NumberValue(default=221.000000, unit="MHz", ndecimals=6)
        )

        # Smooth gate parameters
        self.setattr_argument(
            "delta_max_khz",
            NumberValue(default=300.0, unit="kHz", ndecimals=3)
        )
        self.setattr_argument(
            "delta_min_khz",
            NumberValue(default=25.0, unit="kHz", ndecimals=3)
        )
        self.setattr_argument(
            "tau_d_us",
            NumberValue(default=100.0, unit="us", ndecimals=3)
        )
        self.setattr_argument(
            "tau_g_us",
            NumberValue(default=20.0, unit="us", ndecimals=3)
        )
        self.setattr_argument(
            "hold_us",
            NumberValue(default=100.0, unit="us", ndecimals=3)
        )
        self.setattr_argument(
            "blue_amp_frac",
            NumberValue(default=0.40, min=0.0, max=1.0, ndecimals=4)
        )
        self.setattr_argument(
            "red_amp_frac",
            NumberValue(default=0.40, min=0.0, max=1.0, ndecimals=4)
        )

        # DDS step timing
        self.setattr_argument(
            "dt_us",
            NumberValue(default=2.0, unit="us", ndecimals=3)
        )

        # TTL trigger pulse width
        self.setattr_argument(
            "awg_trigger_ns",
            NumberValue(default=200.0, unit="ns", ndecimals=1)
        )

    def prepare(self):
        self.awg_info = program_smooth_ms_awg(
            blue_sideband_ref_mhz=self.blue_sideband_ref_mhz,
            red_sideband_ref_mhz=self.red_sideband_ref_mhz,
            delta_max_khz=self.delta_max_khz,
            delta_min_khz=self.delta_min_khz,
            tau_d_us=self.tau_d_us,
            tau_g_us=self.tau_g_us,
            hold_us=self.hold_us,
            blue_amp_frac=self.blue_amp_frac,
            red_amp_frac=self.red_amp_frac,
            dt_us=self.dt_us,
    )

        print("AWG programmed successfully")
        print(f"Blue sideband ref (MHz): {self.awg_info['blue_sideband_ref_mhz']}")
        print(f"Red sideband ref (MHz):  {self.awg_info['red_sideband_ref_mhz']}")
        print(f"Total time (us):         {self.awg_info['total_time_us']}")
        print(f"Total steps:             {self.awg_info['total_steps']}")
        print(f"Blue amp frac: {self.blue_amp_frac}")
        print(f"Red amp frac:  {self.red_amp_frac}")

    @kernel
    def run(self):
        self.core.reset()
        self.core.break_realtime()

        # Cooling / pumping / state prep goes here
        delay(1 * us)

        # Launch preloaded AWG waveform
        self.AWG_trigger.pulse(self.awg_trigger_ns * ns)

        # Continue experiment
        delay(1 * us)

    def analyze(self):
        close_awg(self.awg_info)