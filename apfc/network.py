"""The plant: an 11 kV supply, a 1600 kVA transformer and a 415 V main bus
carrying direct-on-line motors and six-pulse drives.

Everything is per phase, referred to 415 V, at harmonic order h.
"""

import math
from dataclasses import dataclass

V_LL = 415.0
F1 = 50.0
W1 = 2 * math.pi * F1

# ---- supply ------------------------------------------------------------------
GRID_SK_MVA = 350.0
GRID_R_OVER_X = 0.1
TR_KVA = 1600.0
TR_UK = 0.0625
TR_PK_KW = 15.5
MV_KV = 11.0


def source_impedance():
    """Grid + transformer at 415 V, fundamental (ohm)."""
    zg = V_LL ** 2 / (GRID_SK_MVA * 1e6)
    xg = zg / math.sqrt(1 + GRID_R_OVER_X ** 2)
    zb = V_LL ** 2 / (TR_KVA * 1e3)
    rt = TR_PK_KW / TR_KVA * zb
    xt = math.sqrt((TR_UK * zb) ** 2 - rt ** 2)
    return complex(GRID_R_OVER_X * xg + rt, xg + xt)


def source_at(h, z1=None):
    """Resistance rises roughly with sqrt(h) from eddy-current and skin effect."""
    z1 = z1 or source_impedance()
    return complex(z1.real * math.sqrt(h), z1.imag * h)


def short_circuit_mva():
    return V_LL ** 2 / abs(source_impedance()) / 1e6


# ---- loads ---------------------------------------------------------------------
@dataclass
class Plant:
    motor_kw: float = 550.0      # DOL induction motors at full plant load
    motor_q_mag: float = 440.0   # magnetising kvar, roughly constant with load
    motor_q_leak: float = 293.0  # leakage kvar at full load, grows with load^2
    drive_kw: float = 450.0      # six-pulse VFDs
    drive_dpf: float = 0.95      # displacement pf of a diode front end
    motor_kva_rated: float = 890.0
    motor_xd2_pu: float = 0.17   # locked-rotor (subtransient) reactance
    motor_x_over_r: float = 8.0

    def pq(self, motor_frac=1.0, drive_frac=1.0):
        pm = self.motor_kw * motor_frac
        qm = self.motor_q_mag + self.motor_q_leak * motor_frac ** 2
        pd = self.drive_kw * drive_frac
        qd = pd * math.tan(math.acos(self.drive_dpf))
        return pm + pd, qm + qd, pd

    def motor_at(self, h):
        x = self.motor_xd2_pu * V_LL ** 2 / (self.motor_kva_rated * 1e3)
        return complex(x / self.motor_x_over_r * math.sqrt(h), x * h)


# Harmonic current of a six-pulse drive with a 3 % AC line reactor, percent of
# the drive's fundamental current (typical maker application data). The sign is
# the phase of a six-pulse bridge, used only to draw waveforms.
DRIVE_SPECTRUM = {5: (35.0, -1), 7: (12.0, 1), 11: (7.5, -1), 13: (5.0, 1),
                  17: (3.5, -1), 19: (3.0, 1), 23: (2.0, -1), 25: (2.0, 1)}
HARMONICS = sorted(DRIVE_SPECTRUM)


def drive_harmonic_currents(drive_kw, dpf=0.95):
    i1 = drive_kw * 1e3 / (math.sqrt(3) * V_LL * dpf)
    return i1, {h: i1 * pct / 100 for h, (pct, _) in DRIVE_SPECTRUM.items()}
