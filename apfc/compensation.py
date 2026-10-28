"""Shunt compensation branches: plain capacitors, detuned (reactor-capacitor)
steps and single-tuned harmonic filters.

All branches are star-equivalent per phase at 415 V. `q_eff_kvar` is what the
branch delivers at 415 V and 50 Hz, which is what the APFC controller and the
power factor see.
"""

import math
from dataclasses import dataclass

from .network import V_LL


@dataclass
class Branch:
    name: str
    q_eff_kvar: float
    tuning: float = None      # series resonance order; None for a plain capacitor
    quality: float = 40.0     # X/R of the reactor at the tuning frequency
    cap_rated_v: float = 415.0

    @property
    def x_c(self):
        """Capacitive reactance at 50 Hz giving q_eff at 415 V."""
        q = self.q_eff_kvar * 1e3
        if self.tuning is None:
            return V_LL ** 2 / q
        # X_c - X_L = V^2 / Q  with  X_L = X_c / tuning^2
        return V_LL ** 2 / q / (1 - 1 / self.tuning ** 2)

    @property
    def x_l(self):
        return 0.0 if self.tuning is None else self.x_c / self.tuning ** 2

    @property
    def r(self):
        return 0.0 if self.tuning is None else self.tuning * self.x_l / self.quality

    def z(self, h):
        return complex(self.r, h * self.x_l - self.x_c / h)

    @property
    def p_detuning(self):
        return self.x_l / self.x_c

    @property
    def cap_voltage_fundamental(self):
        """Capacitor terminal voltage at 50 Hz, raised by the series reactor."""
        return V_LL * self.x_c / (self.x_c - self.x_l)

    @property
    def rated_kvar_at_cap_voltage(self):
        """Nameplate kvar of the capacitor bank at its rated voltage."""
        return self.cap_rated_v ** 2 / self.x_c / 1e3


def plain(name, q, cap_rated_v=440.0):
    return Branch(name, q, cap_rated_v=cap_rated_v)


def detuned(name, q, p=0.07, cap_rated_v=480.0):
    """Reactor-capacitor step tuned below the 5th: p = 7 % gives 189 Hz."""
    return Branch(name, q, tuning=1 / math.sqrt(p), quality=40.0, cap_rated_v=cap_rated_v)


def tuned(name, q, order, quality=40.0, cap_rated_v=525.0):
    """Single-tuned filter, tuned slightly below the target harmonic so that
    capacitor ageing (C falls, f rises) does not push it onto the harmonic."""
    return Branch(name, q, tuning=order, quality=quality, cap_rated_v=cap_rated_v)


def steps_for(q_needed, step=100.0, small=(25.0, 50.0)):
    """Step table: one each of the small steps plus enough large ones."""
    n = math.ceil((q_needed - sum(small)) / step)
    return list(small) + [step] * max(n, 0)
