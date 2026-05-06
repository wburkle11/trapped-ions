import subprocess
import re
from dataclasses import dataclass

#_NUM and _parse_float handle the SNMP return values to make them user
#friendly for Python and readible
_NUM = re.compile(r"([-+]?\d+(\.\d+)?([eE][-+]?\d+)?)")

def _parse_float(s: str) -> float:
    m = _NUM.search(s)
    if not m:
        raise ValueError(f"Could not parse number from: {s}")
    return float(m.group(1))


@dataclass
class ChannelStatus:
    slot: int
    ch: int
    idx: int
    on: int
    vset: float
    v: float
    i: float

    def __str__(self) -> str:
        state = "ON" if self.on else "OFF"
        return (
            f"MPOD slot {self.slot}, channel {self.ch} "
            f"(idx {self.idx})\n"
            f"  State   : {state}\n"
            f"  Vset    : {self.vset:.6f} V\n"
            f"  Vmeas   : {self.v:.6f} V\n"
            f"  Imeas   : {self.i:.6e} A"
        )


class MPOD:
    def __init__(self,
                 #ip address for MPOD crate;
                 #if using different crate, change the ip address
                 #which can be found in the properties tab of the Easy LV HV software
                 #after booting up with USB
                 #see lab writeup

                 #If on Windows, you will need to make sure and manually set the ip address
                 #to the corresponding ethernet port and create a routing priority table for 
                 #multiple ethernet ports so that Windows knows where to send packets to
                 #See lab writeup
                 ip: str = "140.182.217.76",
                 community: str = "guru"):
        self.ip = ip
        self.community = community
        self.BASE = "1.3.6.1.4.1.19947"

    # GUI channel → SNMP index which takes care of appropriate mapping of
    #channel to index
    def idx(self, slot: int, ch: int) -> int:
        return slot * 100 + ch + 1

    def _run(self, cmd: list[str]) -> str:
        return subprocess.check_output(
            cmd, stderr=subprocess.STDOUT, text=True
        ).strip()

    def _get(self, oid: str) -> str:
        return self._run([
            "snmpget", "-Ovq", "-v2c", "-c",
            self.community, self.ip, oid
        ])

    def _set(self, oid: str, typecode: str, value: str) -> None:
        self._run([
            "snmpset", "-Ovq", "-v2c", "-c",
            self.community, self.ip, oid, typecode, value
        ])

    # ---------- Public API ----------

    #prints out relevant channel information; can add more configurations if needed
    #see mpod-snmp-guide.pdf
    def read(self, slot: int, ch: int) -> ChannelStatus:
        i = self.idx(slot, ch)
        return ChannelStatus(
            slot=slot,
            ch=ch,
            idx=i,
            on=int(_parse_float(self._get(f"{self.BASE}.1.3.2.1.9.{i}"))),
            vset=_parse_float(self._get(f"{self.BASE}.1.3.2.1.10.{i}")),
            v=_parse_float(self._get(f"{self.BASE}.1.3.2.1.6.{i}")),
            i=_parse_float(self._get(f"{self.BASE}.1.3.2.1.7.{i}")),
        )

    #set voltage - format (slot, channel) 
    #slot corresponds to the module's slot on the crate starting at 0
    #channel corresponds to physical channel port on module that is listed
    def set_voltage(self, slot: int, ch: int, volts: float) -> None:
        self._set(
            f"{self.BASE}.1.3.2.1.10.{self.idx(slot, ch)}",
            "F", str(volts)
        )

    #A helper function I added as the measured voltage takes time to ramp up
    #def wait_for_vmeas(self, slot: int, ch: int, target: float | None = None,
    #              tol: float = 0.02, timeout: float, poll: float = 0.2):
    #    #Wait until measured voltage is close to target (or close to Vset if target=None),
    #    #then return ChannelStatus (so you can print it normally).
        
    #    t0 = time.time()
    #    last = None

    #    while time.time() - t0 < timeout:
    #        s = self.read(slot, ch)
    #        last = s

    #        tv = target if target is not None else s.vset

    #        if abs(s.v - tv) <= tol:
    #            return s

    #        time.sleep(poll)

    #    #If it never reaches tolerance, return the last read so you can see what happened
    #    return last
    
    #set current limit
    def set_current_limit(self, slot: int, ch: int, amps: float) -> None:
        self._set(
            f"{self.BASE}.1.3.2.1.12.{self.idx(slot, ch)}",
            "F", str(amps)
        )

    #turns channel on
    def on(self, slot: int, ch: int) -> None:
        self._set(
            f"{self.BASE}.1.3.2.1.9.{self.idx(slot, ch)}",
            "i", "1"
        )
    
    #turns channel off
    def off(self, slot: int, ch: int) -> None:
        self._set(
            f"{self.BASE}.1.3.2.1.9.{self.idx(slot, ch)}",
            "i", "0"
        )
    #Sustain channel voltage for timehold seconds (call this function instead of
    # wait_for_vmeas if you want to sustain your target voltage for a prolonged time)
    #def sustain_voltage(self, slot: int, ch: int, target: float | None = None,
    #               tol: float = 0.02, timehold: float, timeout: float, poll: float = 0.2):
        
    #    wait_for_vmeas(self,slot,ch,target,tol,timeout,poll)
    #    while time.time() - timeout < timehold:
    #        s = self.read(slot, ch)
    #        last = s
    #        time.sleep(poll)

