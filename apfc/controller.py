"""24-hour simulation of the APFC relay switching capacitor steps against a
varying plant load.

The relay measures P and Q every `dt` seconds and works to keep the power
factor inside a band around the target:
  - switch a step in when the reactive deficit exceeds the C/k threshold
    (a fraction of the smallest step), choosing the largest step that fits
  - switch a step out when the plant goes leading by more than the threshold
  - a step that has just been switched out is locked out for the discharge
    time (IEC 60831 discharge to 75 V), and consecutive operations are spaced
    by the switching delay
  - equal steps are rotated so they wear evenly
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .network import Plant


@dataclass
class RelaySettings:
    target_pf: float = 0.99
    ck_fraction: float = 0.65
    switching_delay_s: float = 20.0
    discharge_s: float = 60.0
    dt_s: float = 10.0


def load_profile(seed=7, dt_s=10.0):
    """Three-shift plant: motor and drive loading as fractions of full load."""
    rng = np.random.default_rng(seed)
    t = np.arange(0, 24 * 3600, dt_s)
    hour = t / 3600
    # base shift pattern (night shift lighter, lunch and shift-change dips)
    base = np.interp(hour, [0, 6, 7, 9, 12.5, 13.5, 14, 17.5, 18.5, 22, 23, 24],
                     [0.45, 0.45, 0.75, 0.95, 0.95, 0.60, 0.95, 1.00, 0.70, 0.70, 0.45, 0.45])
    # minute-scale process variation, smoothed
    noise = rng.normal(0, 0.06, len(t))
    kernel = np.ones(12) / 12
    noise = np.convolve(noise, kernel, mode="same")
    motor = np.clip(base + noise, 0.25, 1.0)
    drive = np.clip(base * 0.95 + 0.6 * noise + rng.normal(0, 0.02, len(t)), 0.25, 1.0)
    return t, motor, drive


@dataclass
class Simulation:
    t: np.ndarray
    p: np.ndarray
    q: np.ndarray
    q_comp: np.ndarray
    steps_on: np.ndarray
    operations: list = field(default_factory=list)

    @property
    def pf_before(self):
        return self.p / np.hypot(self.p, self.q)

    @property
    def pf_after(self):
        qn = self.q - self.q_comp
        return np.sign(qn + 1e-9) * self.p / np.hypot(self.p, qn)

    def energy(self, dt):
        h = dt / 3600
        kwh = np.sum(self.p) * h
        kvah_before = np.sum(np.hypot(self.p, self.q)) * h
        kvah_after = np.sum(np.hypot(self.p, self.q - self.q_comp)) * h
        return {"kwh": kwh, "kvah_before": kvah_before, "kvah_after": kvah_after,
                "avg_pf_before": kwh / kvah_before, "avg_pf_after": kwh / kvah_after}


def simulate(steps, fixed_kvar=0.0, plant=Plant(), s: RelaySettings = None, seed=7):
    """`steps`: switchable effective kvar per step; `fixed_kvar`: permanently
    connected compensation (tuned filters)."""
    s = s or RelaySettings()
    t, motor, drive = load_profile(seed, s.dt_s)
    n = len(t)
    on = [False] * len(steps)
    off_since = [-1e9] * len(steps)
    last_op = -1e9
    use_count = [0] * len(steps)
    tan_t = math.tan(math.acos(s.target_pf))
    threshold = s.ck_fraction * min(steps)

    p_arr, q_arr, qc_arr, steps_arr = (np.zeros(n) for _ in range(4))
    ops = []
    for k in range(n):
        p, q, _ = plant.pq(motor[k], drive[k])
        q_comp = fixed_kvar + sum(st for st, o in zip(steps, on) if o)
        deficit = (q - q_comp) - p * tan_t      # >0: still lagging below target
        now = t[k]
        if now - last_op >= s.switching_delay_s:
            if deficit > threshold:
                ready = [i for i, o in enumerate(on) if not o
                         and now - off_since[i] >= s.discharge_s
                         and steps[i] <= deficit + threshold]
                if ready:
                    # largest step that fits, least used among equals
                    i = max(ready, key=lambda j: (steps[j], -use_count[j]))
                    on[i] = True
                    use_count[i] += 1
                    last_op = now
                    ops.append((now, "in", i))
            elif deficit < -threshold:
                live = [i for i, o in enumerate(on) if o and steps[i] <= -deficit + threshold]
                if live:
                    i = max(live, key=lambda j: (steps[j], use_count[j]))
                    on[i] = False
                    off_since[i] = now
                    last_op = now
                    ops.append((now, "out", i))
        q_comp = fixed_kvar + sum(st for st, o in zip(steps, on) if o)
        p_arr[k], q_arr[k], qc_arr[k] = p, q, q_comp
        steps_arr[k] = sum(on)
    return Simulation(t, p_arr, q_arr, qc_arr, steps_arr, ops)
