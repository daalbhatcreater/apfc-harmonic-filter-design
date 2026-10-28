"""Frequency-domain harmonic penetration at the 415 V bus.

The drives are harmonic current sources. At each order h the bus voltage is the
injected current times the parallel impedance of supply, motors and
compensation branches; the share that flows back into the supply is what the
utility sees.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .network import (HARMONICS, MV_KV, V_LL, Plant, drive_harmonic_currents,
                      source_at, source_impedance)

V_PH = V_LL / math.sqrt(3)

# IEEE 519-2022 Table 2, 120 V to 69 kV: Isc/IL -> (h<11, 11-17, 17-23, 23-35, TDD) %
IEEE519_CURRENT = [
    (20, (4.0, 2.0, 1.5, 0.6, 5.0)),
    (50, (7.0, 3.5, 2.5, 1.0, 8.0)),
    (100, (10.0, 4.5, 4.0, 1.5, 12.0)),
    (1000, (12.0, 5.5, 5.0, 2.0, 15.0)),
    (math.inf, (15.0, 7.0, 6.0, 2.5, 20.0)),
]
# IEEE 519-2022 Table 1: bus voltage -> (individual %, THD %)
IEEE519_VOLTAGE = {"lv": (5.0, 8.0), "mv": (3.0, 5.0)}


def current_limits(isc_over_il):
    for bound, limits in IEEE519_CURRENT:
        if isc_over_il < bound:
            return limits


def band_limit(h, limits):
    if h < 11:
        return limits[0]
    if h < 17:
        return limits[1]
    if h < 23:
        return limits[2]
    return limits[3]


@dataclass
class Result:
    name: str
    p_kw: float
    q_kvar: float
    q_comp_kvar: float
    v_h: dict              # bus voltage, % of fundamental
    i_supply_h: dict       # amps into the supply
    i_branch_h: dict = field(default_factory=dict)
    i1_supply: float = 0.0
    thd_v: float = 0.0
    thd_i: float = 0.0
    tdd: float = 0.0
    dpf: float = 0.0
    true_pf: float = 0.0
    isc_over_il: float = 0.0
    tdd_limit: float = 0.0
    individual_ok: bool = True
    duty: list = field(default_factory=list)

    @property
    def tdd_ok(self):
        return self.tdd <= self.tdd_limit

    @property
    def thd_v_ok(self):
        return self.thd_v <= IEEE519_VOLTAGE["lv"][1] and max(self.v_h.values()) <= \
            IEEE519_VOLTAGE["lv"][0]


def solve(name, branches, plant=Plant(), motor_frac=1.0, drive_frac=1.0):
    p, q, p_drive = plant.pq(motor_frac, drive_frac)
    q_comp = sum(b.q_eff_kvar for b in branches)
    _, inj = drive_harmonic_currents(p_drive, plant.drive_dpf)
    z1 = source_impedance()

    v_h, i_s, i_b = {}, {}, {b.name: {} for b in branches}
    for h in HARMONICS:
        zs = source_at(h, z1)
        y = 1 / zs + 1 / plant.motor_at(h) + sum(1 / b.z(h) for b in branches)
        v = inj[h] / y
        v_h[h] = 100 * abs(v) / V_PH
        i_s[h] = abs(v / zs)
        for b in branches:
            i_b[b.name][h] = abs(v / b.z(h))

    q_net = q - q_comp
    s1 = math.hypot(p, q_net)
    i1 = s1 * 1e3 / (math.sqrt(3) * V_LL)
    ih = math.sqrt(sum(x * x for x in i_s.values()))
    r = Result(name, p, q, q_comp, v_h, i_s, i_b, i1)
    r.thd_v = math.sqrt(sum(x * x for x in v_h.values()))
    r.thd_i = 100 * ih / i1
    r.dpf = p / s1
    r.true_pf = p * 1e3 / (math.sqrt(3) * V_LL * math.hypot(i1, ih))
    return r


def evaluate_ieee519(r, i_l, isc_a):
    """TDD against the maximum demand current i_l at this bus."""
    r.isc_over_il = isc_a / i_l
    limits = current_limits(r.isc_over_il)
    r.tdd_limit = limits[4]
    r.tdd = 100 * math.sqrt(sum(x * x for x in r.i_supply_h.values())) / i_l
    r.individual_ok = all(100 * i / i_l <= band_limit(h, limits) for h, i in r.i_supply_h.items())
    return limits


def capacitor_duty(r, branches):
    """IEEE 18-2012 limits: 110 % rms voltage, 120 % peak voltage,
    135 % rms current, 135 % kvar, each relative to the nameplate."""
    rows = []
    for b in branches:
        i1 = V_PH / abs(b.z(1))
        vc = {1: i1 * b.x_c}
        ic = {1: i1}
        for h, i in r.i_branch_h[b.name].items():
            vc[h] = i * b.x_c / h
            ic[h] = i
        u_r = b.cap_rated_v / math.sqrt(3)
        i_r = b.rated_kvar_at_cap_voltage * 1e3 / (3 * u_r)
        v_rms = math.sqrt(sum(v * v for v in vc.values()))
        v_pk = math.sqrt(2) * sum(vc.values())
        i_rms = math.sqrt(sum(i * i for i in ic.values()))
        kvar = 3 * sum(vc[h] * ic[h] for h in vc) / 1e3
        rows.append({
            "branch": b.name, "v_rms_pct": 100 * v_rms / u_r,
            "v_peak_pct": 100 * v_pk / (math.sqrt(2) * u_r),
            "i_rms_pct": 100 * i_rms / i_r,
            "kvar_pct": 100 * kvar / b.rated_kvar_at_cap_voltage,
        })
    for row in rows:
        row["ok"] = (row["v_rms_pct"] <= 110 and row["v_peak_pct"] <= 120
                     and row["i_rms_pct"] <= 135 and row["kvar_pct"] <= 135)
    r.duty = rows
    return rows


def frequency_scan(branches, plant=Plant(), orders=None):
    """|Z| seen by a harmonic source at the 415 V bus."""
    orders = np.linspace(1, 25, 1200) if orders is None else orders
    z1 = source_impedance()
    out = []
    for h in orders:
        y = 1 / source_at(h, z1) + 1 / plant.motor_at(h) + sum(1 / b.z(h) for b in branches)
        out.append(abs(1 / y))
    return np.array(orders), np.array(out)


def waveform(r, i1_peak_scale=True, cycles=2, n=2000, plant=Plant()):
    """Supply current waveform rebuilt from its harmonic phasors (amps)."""
    from .network import DRIVE_SPECTRUM
    t = np.linspace(0, cycles / 50, n)
    i = math.sqrt(2) * r.i1_supply * np.sin(2 * math.pi * 50 * t)
    for h, amp in r.i_supply_h.items():
        sign = DRIVE_SPECTRUM[h][1]
        i += sign * math.sqrt(2) * amp * np.sin(2 * math.pi * 50 * h * t)
    return t, i


def mv_view(r):
    """The same harmonic currents seen at 11 kV through the Dyn transformer
    (5th, 7th, 11th... all pass; only the magnitude scales)."""
    return {h: i * V_LL / (MV_KV * 1e3) for h, i in r.i_supply_h.items()}
