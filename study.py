"""APFC sizing, harmonic penetration and filter design for the plant.

    python study.py
"""

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from apfc.compensation import detuned, plain, steps_for, tuned
from apfc.controller import simulate
from apfc.harmonics import (capacitor_duty, current_limits, evaluate_ieee519,
                            frequency_scan, mv_view, solve, waveform)
from apfc.network import (GRID_SK_MVA, HARMONICS, MV_KV, V_LL, Plant,
                          short_circuit_mva)

OUT = Path(__file__).parent / "results"
PLANT = Plant()
FILTER_5_KVAR, FILTER_7_KVAR = 300.0, 150.0


def configurations():
    p, q, _ = PLANT.pq()
    q_need = q - p * math.tan(math.acos(0.99))
    steps = steps_for(q_need)
    # at full load the relay settles on the step combination closest to the need
    q_full = min((sum(c) for c in _subsets(steps)), key=lambda x: abs(x - q_need))
    fixed = FILTER_5_KVAR + FILTER_7_KVAR
    hybrid_steps = steps_for(q_need - fixed)
    q_hybrid = min((sum(c) for c in _subsets(hybrid_steps)),
                   key=lambda x: abs(x + fixed - q_need))
    # name: (branches at full load, switched steps, fixed kvar)
    return {
        "A. Uncompensated": ([], None, 0.0),
        "B. Plain capacitor APFC": ([plain("Plain capacitors", q_full)], steps, 0.0),
        "C. Detuned APFC (7 %, 189 Hz)": ([detuned("Detuned steps", q_full)], steps, 0.0),
        "D. 5th + 7th tuned filters + detuned APFC": (
            [tuned("5th filter (4.7)", FILTER_5_KVAR, 4.7),
             tuned("7th filter (6.7)", FILTER_7_KVAR, 6.7),
             detuned("Detuned steps", q_hybrid)], hybrid_steps, fixed),
    }, q_need


def _subsets(values):
    from itertools import combinations
    for n in range(len(values) + 1):
        yield from combinations(values, n)


def harmonic_study(configs):
    isc_lv = short_circuit_mva() * 1e6 / (math.sqrt(3) * V_LL)
    isc_mv = GRID_SK_MVA * 1e6 / (math.sqrt(3) * MV_KV * 1e3)
    x_grid_mv = (MV_KV * 1e3) ** 2 / (GRID_SK_MVA * 1e6) * 0.995
    results, rows, duty_rows = {}, [], []
    for name, (branches, _, _) in configs.items():
        r = solve(name, branches, PLANT)
        evaluate_ieee519(r, r.i1_supply, isc_lv)
        duty = capacitor_duty(r, branches)
        # the same currents at the 11 kV point of common coupling
        i_mv = mv_view(r)
        il_mv = r.i1_supply * V_LL / (MV_KV * 1e3)
        tdd_mv = 100 * math.sqrt(sum(i * i for i in i_mv.values())) / il_mv
        v_ph_mv = MV_KV * 1e3 / math.sqrt(3)
        thdv_mv = 100 * math.sqrt(sum((i * h * x_grid_mv) ** 2 for h, i in i_mv.items())) / v_ph_mv
        mv_limits = current_limits(isc_mv / il_mv)
        results[name] = r
        rows.append({
            "Configuration": name,
            "Compensation (kvar)": round(r.q_comp_kvar),
            "Displacement pf": round(r.dpf, 3),
            "True pf": round(r.true_pf, 3),
            "THDv 415 V (%)": round(r.thd_v, 2),
            "Worst Vh 415 V (%)": round(max(r.v_h.values()), 2),
            "TDD 415 V (%)": round(r.tdd, 2),
            "Isc/IL 415 V": round(r.isc_over_il, 1),
            "TDD limit 415 V (%)": r.tdd_limit,
            "415 V current OK": r.tdd_ok and r.individual_ok,
            "415 V voltage OK": r.thd_v_ok,
            "THDv 11 kV (%)": round(thdv_mv, 2),
            "TDD 11 kV (%)": round(tdd_mv, 2),
            "Isc/IL 11 kV": round(isc_mv / il_mv),
            "TDD limit 11 kV (%)": mv_limits[4],
            "Caps within IEEE 18": all(d["ok"] for d in duty) if duty else "",
        })
        for d in duty:
            duty_rows.append({"Configuration": name, "Branch": d["branch"],
                              "V rms (%)": round(d["v_rms_pct"], 1),
                              "V peak (%)": round(d["v_peak_pct"], 1),
                              "I rms (%)": round(d["i_rms_pct"], 1),
                              "kvar (%)": round(d["kvar_pct"], 1), "Within limits": d["ok"]})
    return results, pd.DataFrame(rows), pd.DataFrame(duty_rows)


def filter_table():
    rows = []
    for b in (detuned("Detuned step, 100 kvar", 100.0),
              tuned("5th filter (4.7)", FILTER_5_KVAR, 4.7),
              tuned("7th filter (6.7)", FILTER_7_KVAR, 6.7)):
        rows.append({
                "Branch": b.name, "Effective kvar at 415 V": round(b.q_eff_kvar),
                "Tuned to": f"{b.tuning:.2f} ({b.tuning * 50:.0f} Hz)",
                "p = XL/XC (%)": round(100 * b.p_detuning, 2),
                "C per phase, star (uF)": round(1e6 / (2 * math.pi * 50 * b.x_c), 1),
                "L per phase (mH)": round(1e3 * b.x_l / (2 * math.pi * 50), 4),
                "R per phase (mOhm)": round(1e3 * b.r, 3),
                "Capacitor voltage at 50 Hz (V)": round(b.cap_voltage_fundamental),
                "Capacitor rating (V)": b.cap_rated_v,
                "Nameplate kvar at rated V": round(b.rated_kvar_at_cap_voltage),
        })
    return pd.DataFrame(rows)


def plot_scan(path):
    fig, ax = plt.subplots(figsize=(9, 5))
    for label, branches, style in [
        ("Uncompensated", [], dict(color="black")),
        ("Plain capacitors 200 kvar", [plain("x", 200)], dict(color="tab:red", ls=":")),
        ("Plain capacitors 400 kvar", [plain("x", 400)], dict(color="tab:red", ls="--")),
        ("Plain capacitors 725 kvar", [plain("x", 725)], dict(color="tab:red")),
        ("Detuned 7 %, 725 kvar", [detuned("x", 725)], dict(color="tab:blue")),
        ("5th + 7th filters + detuned", [tuned("a", 300, 4.7), tuned("b", 150, 6.7),
                                         detuned("c", 275)], dict(color="tab:green")),
    ]:
        h, z = frequency_scan(branches, PLANT)
        ax.plot(h, 1e3 * z, label=label, **style)
    for h in HARMONICS:
        ax.axvline(h, color="grey", lw=0.5, ls=":")
    ax.set_yscale("log")
    ax.set_xlim(1, 25)
    ax.set_xlabel("Harmonic order")
    ax.set_ylabel("|Z| at 415 V bus (mOhm)")
    ax.set_title("Frequency scan: impedance seen by the drives")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", lw=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_spectrum(results, path):
    names = list(results)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    w = 0.2
    for k, name in enumerate(names):
        r = results[name]
        vals = [100 * r.i_supply_h[h] / r.i1_supply for h in HARMONICS]
        ax.bar(np.arange(len(HARMONICS)) + (k - 1.5) * w, vals, w, label=name)
    ax.set_xticks(range(len(HARMONICS)), [str(h) for h in HARMONICS])
    ax.set_xlabel("Harmonic order")
    ax.set_ylabel("Supply current (% of fundamental)")
    ax.set_title("Harmonic current drawn from the transformer at full load")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", lw=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_waveforms(results, path):
    fig, ax = plt.subplots(figsize=(9, 4))
    for name, colour in [("A. Uncompensated", "black"), ("B. Plain capacitor APFC", "tab:red"),
                         ("D. 5th + 7th tuned filters + detuned APFC", "tab:green")]:
        t, i = waveform(results[name])
        ax.plot(1e3 * t, i, color=colour, label=f"{name} (THDi {results[name].thd_i:.1f} %)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Phase current (A)")
    ax.set_title("Transformer LV current rebuilt from the harmonic solution")
    ax.legend(fontsize=8)
    ax.grid(True, lw=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_day(sim, path, title):
    h = sim.t / 3600
    fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    ax[0].plot(h, sim.p, label="kW")
    ax[0].plot(h, sim.q, label="kvar, plant")
    ax[0].plot(h, sim.q_comp, label="kvar, compensation")
    ax[0].set_ylabel("kW / kvar")
    ax[0].legend(fontsize=8, loc="lower right")
    ax[1].plot(h, sim.pf_before, color="tab:red", label="pf without APFC")
    ax[1].plot(h, np.abs(sim.pf_after), color="tab:green", label="pf with APFC")
    ax[1].axhline(0.99, color="grey", ls=":", lw=0.8)
    ax[1].set_ylim(0.5, 1.01)
    ax[1].set_ylabel("Power factor")
    ax[1].legend(fontsize=8, loc="lower right")
    ax[2].step(h, sim.steps_on, where="post", color="tab:purple")
    ax[2].set_ylabel("Steps in service")
    ax[2].set_xlabel("Hour of day")
    for a in ax:
        a.grid(True, lw=0.3)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    configs, q_need = configurations()
    p, q, _ = PLANT.pq()

    results, harm, duty = harmonic_study(configs)
    filters = filter_table()

    day_rows = []
    for name, (_, steps, fixed) in configs.items():
        if steps is None or name.startswith("B"):
            continue
        sim = simulate(steps, fixed, PLANT)
        e = sim.energy(10.0)
        key = "detuned" if name.startswith("C") else "hybrid"
        plot_day(sim, OUT / f"apfc_day_{key}.png", name)
        day_rows.append({
            "Scheme": name, "Switched steps (kvar)": " + ".join(f"{s:g}" for s in steps),
            "Fixed filters (kvar)": fixed, "kWh/day": round(e["kwh"]),
            "kVAh/day without": round(e["kvah_before"]), "kVAh/day with": round(e["kvah_after"]),
            "kVAh saved (%)": round(100 * (1 - e["kvah_after"] / e["kvah_before"]), 1),
            "Average pf without": round(e["avg_pf_before"], 3),
            "Average pf with": round(e["avg_pf_after"], 3),
            "Switching operations/day": len(sim.operations),
            "Time below 0.98 (%)": round(100 * np.mean(np.abs(sim.pf_after) < 0.98), 1),
        })
    day = pd.DataFrame(day_rows)

    plot_scan(OUT / "frequency_scan.png")
    plot_spectrum(results, OUT / "harmonic_spectrum.png")
    plot_waveforms(results, OUT / "current_waveforms.png")

    basis = pd.DataFrame([
        ("Plant at full load", f"{p:.0f} kW, {q:.0f} kvar, pf {p / math.hypot(p, q):.3f}"),
        ("Drives (six-pulse, 3 % reactor)", f"{PLANT.drive_kw:.0f} kW"),
        ("Short-circuit level at 415 V", f"{short_circuit_mva():.1f} MVA"),
        ("kvar to reach pf 0.99 at full load", f"{q_need:.0f} kvar"),
        ("Parallel resonance with plain capacitors", "h = sqrt(Ssc / Qc): " + ", ".join(
            f"{qc} kvar -> {math.sqrt(short_circuit_mva() * 1e3 / qc):.2f}" for qc in (200, 400, 725))),
    ], columns=["Item", "Value"])

    with pd.ExcelWriter(OUT / "apfc_harmonics.xlsx") as xl:
        basis.to_excel(xl, sheet_name="Basis", index=False)
        harm.to_excel(xl, sheet_name="Harmonics", index=False)
        duty.to_excel(xl, sheet_name="Capacitor duty", index=False)
        filters.to_excel(xl, sheet_name="Filter design", index=False)
        day.to_excel(xl, sheet_name="24 h APFC", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    for title, df in [("Basis", basis), ("Harmonics at full load", harm),
                      ("Capacitor duty (IEEE 18)", duty), ("Filter and detuned step design", filters),
                      ("24 h APFC simulation", day)]:
        print(title)
        print(df.to_string(index=False))
        print()


if __name__ == "__main__":
    main()
