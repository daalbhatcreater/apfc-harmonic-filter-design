import math

import numpy as np
import pytest

from apfc.compensation import detuned, plain, steps_for, tuned
from apfc.controller import simulate
from apfc.harmonics import current_limits, evaluate_ieee519, solve
from apfc.network import V_LL, Plant, short_circuit_mva


def test_short_circuit_level():
    # 350 MVA grid in series with 1600 kVA at 6.25 % -> about 23.9 MVA at 415 V
    assert short_circuit_mva() == pytest.approx(23.86, abs=0.05)


@pytest.mark.parametrize("branch", [plain("p", 100), detuned("d", 100), tuned("t", 100, 4.7)])
def test_branches_deliver_their_kvar_at_50hz(branch):
    q = V_LL ** 2 / -branch.z(1).imag / 1e3
    assert q == pytest.approx(100, rel=1e-6)


def test_tuned_branch_resonates_at_its_order():
    f = tuned("t", 300, 4.7)
    assert abs(f.z(4.7).imag) < 1e-9
    assert detuned("d", 100).tuning == pytest.approx(3.78, abs=0.01)


def test_detuned_capacitor_voltage_rise():
    # 7 % reactor: capacitor sees 415 / 0.93
    assert detuned("d", 100).cap_voltage_fundamental == pytest.approx(415 / 0.93, rel=1e-6)


def test_ieee519_table():
    assert current_limits(10)[4] == 5.0
    assert current_limits(23.7)[4] == 8.0
    assert current_limits(347)[4] == 15.0


def _run(branches):
    r = solve("x", branches, Plant())
    isc = short_circuit_mva() * 1e6 / (math.sqrt(3) * V_LL)
    evaluate_ieee519(r, r.i1_supply, isc)
    return r


def test_plain_capacitors_amplify_and_filters_reduce_distortion():
    a = _run([])
    b = _run([plain("p", 750)])
    c = _run([detuned("d", 750)])
    d = _run([tuned("f5", 300, 4.7), tuned("f7", 150, 6.7), detuned("d", 300)])
    assert b.thd_v > 2 * a.thd_v and b.tdd > 3 * a.tdd
    assert d.thd_v < c.thd_v < a.thd_v
    assert d.tdd < c.tdd < a.tdd
    assert d.tdd_ok and d.individual_ok and d.thd_v_ok
    assert not a.tdd_ok and not b.tdd_ok


def test_step_table_covers_the_need():
    steps = steps_for(738)
    assert sum(steps) >= 738 and steps[:2] == [25.0, 50.0]


def test_apfc_holds_target_through_the_day():
    sim = simulate(steps_for(738))
    e = sim.energy(10.0)
    assert e["avg_pf_before"] < 0.75
    assert e["avg_pf_after"] > 0.98
    assert np.mean(np.abs(sim.pf_after) >= 0.95) > 0.99
