"""Generate a phase/neutral teaching figure only from audited saved evidence."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "studies/australian_pv"


def load_evidence(folder):
    env = json.loads((folder / "environment.json").read_text())
    audit = json.loads((folder / "audit.json").read_text())
    if not audit["all_planned_accepted"]:
        raise ValueError("all predeclared points must pass the saved-point audit")
    for path, digest in env["source_sha256"].items():
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != digest:
            raise ValueError(f"stale physical/setup source: {path}")
    with (folder / "cases.csv").open() as stream:
        cases = list(csv.DictReader(stream))
    with (folder / "devices.csv").open() as stream:
        devices = list(csv.DictReader(stream))
    if len(cases) != audit["planned"] or len({c["case_id"] for c in cases}) != len(cases):
        raise ValueError("case coverage mismatch")
    if {c["case_id"] for c in cases} != set(audit["cases"]):
        raise ValueError("unplanned or missing case")
    seen = set()
    for c in cases:
        truth = audit["cases"][c["case_id"]]
        if c["accepted"] != "true" or not truth["accepted"]:
            raise ValueError("unaccepted case")
        for key in ("p_poc_W", "q_poc_var", "vmax_sensed_V", "vmin_sensed_V", "max_vuf_percent", "max_neutral_current_A"):
            if not np.isfinite(float(c[key])) or abs(float(c[key]) - truth[key]) > 1e-5:
                raise ValueError(f"case metric mismatch: {key}")
    for d in devices:
        key = (d["case_id"], d["device_id"])
        if key in seen:
            raise ValueError("duplicate device")
        seen.add(key)
        if d["accepted"] != "true":
            raise ValueError("unaccepted device")
        truth = audit["cases"][d["case_id"]]["devices"]
        if d["model_present"] == "false":
            if truth or float(d["p_poc_W"]) != 0 or float(d["neutral_total_rms_current_A"]) != 0:
                raise ValueError("disconnected record mismatch")
            continue
        t = truth[d["device_id"]]
        for column, canonical in (("p_poc_W", "p_poc_W"), ("q_poc_var", "q_poc_var"),
                                  ("neutral_total_rms_current_A", "neutral_current_A"),
                                  ("vmax_sensed_V", "vmax_sensed_V"), ("vuf_percent", "vuf_percent")):
            if not np.isfinite(float(d[column])) or abs(float(d[column]) - t[canonical]) > 1e-5:
                raise ValueError(f"device metric mismatch: {column}")
    expected = {(c["case_id"], device) for c in cases for device in
                (audit["cases"][c["case_id"]]["devices"] or
                 audit["cases"]["matched_s2_p3"]["devices"])}
    if seen != expected:
        raise ValueError("device coverage mismatch")
    return cases, devices, audit


def generate(folder, destination):
    cases, devices, audit = load_evidence(folder)
    point = json.loads((folder / "points/matched_s2_p3.json").read_text())
    selected = [d for d in devices if d["case_id"] == "matched_s2_p3"]
    customer = max(selected, key=lambda d: float(d["vmax_sensed_V"]))
    ident = customer["device_id"]
    def ph(t):
        return complex(t["vr"], t["vi"])
    bus = point["bus"][customer["bus"]]
    voltage = np.array([ph(bus[p]) - ph(bus["n"]) for p in "abc"])
    record = point["devices"][ident]
    current = np.array([complex(*p) for p in record["grid_current"]])
    power = voltage * np.conj(current)
    rawq = 2000 * np.interp(abs(voltage), [207, 220, 240, 258], [.44, 0, 0, -.60])
    colors = ["#c54a3e", "#2d8b57", "#396bb0"]
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "svg.hashsalt": "four-leg-phase-control"})
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.5), gridspec_kw={"width_ratios": [1.1, 1, 1.15]})
    ax = axes[0]
    order = sorted(selected, key=lambda d: float(d["vmax_sensed_V"]))
    for k, phase in enumerate("abc"):
        volts = []
        for d in order:
            b = point["bus"][d["bus"]]
            volts.append(abs(ph(b[phase]) - ph(b["n"])))
        ax.plot(range(1, len(order)+1), volts, "o-", color=colors[k], label=phase.upper(), ms=3)
    ax.axhline(240, color="#777", ls=":", lw=1)
    ax.axhline(253, color="#777", ls="--", lw=1)
    ax.text(1, 240.4, "VV absorption starts", color="#666", fontsize=8)
    ax.text(1, 253.4, "VW starts / illustrative upper criterion", color="#666", fontsize=8)
    ax.set(title="a  Each phase senses its own voltage", xlabel="Customer (ordered by maximum voltage)", ylabel="Phase-neutral voltage (V)", ylim=(239.4, 255))
    ax.legend(ncol=3, loc="center left", bbox_to_anchor=(.01, .32), frameon=False)
    ax = axes[1]
    x = np.arange(3)
    ax.bar(x-.17, rawq/1000, width=.34, color=colors, alpha=.35, label="Requested")
    ax.bar(x+.17, power.imag/1000, width=.34, color=colors, label="Delivered")
    ax.set(xticks=x, xticklabels=list("ABC"), ylabel="Phase Q (kvar)", title="b  Q priority preserves the phase vector", ylim=(-1.38, .08))
    ax.axhline(0, color="#888", lw=.7)
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="lower left", bbox_to_anchor=(.01, .24))
    ax.text(.03, .035, f"Highest-voltage customer:\n{ident}\nQ service = {record['reactive_scale']:.3f}; equal phase P = {power.real.mean():.0f} W",
            transform=ax.transAxes, va="bottom", fontsize=8)
    ax = axes[2]
    for k, z in enumerate(current):
        ax.annotate("", xy=(z.real, z.imag), xytext=(0, 0), arrowprops={"arrowstyle": "->", "color": colors[k], "lw": 2})
        ax.text(z.real*1.10, z.imag*1.10, f"I{list('ABC')[k]}", color=colors[k])
    neutral = -sum(current)
    ax.annotate("", xy=(neutral.real, neutral.imag), xytext=(0, 0), arrowprops={"arrowstyle": "->", "color": "#222", "lw": 2})
    ax.annotate("IN", xy=(neutral.real, neutral.imag), xytext=(-3, 2), fontsize=9,
                arrowprops={"arrowstyle": "-", "color": "#444", "lw": .7})
    ax.text(.02, .02, f"IN = −(IA + IB + IC)\n|IN| = {abs(neutral):.2f} A; rating 10 A RMS", transform=ax.transAxes, fontsize=9)
    ax.axhline(0, color="#ddd", lw=.7); ax.axvline(0, color="#ddd", lw=.7)
    ax.set(xlim=(-11, 11), ylim=(-10, 10), aspect="equal", xlabel="Real current (A RMS)", ylabel="Imaginary current (A RMS)",
           title="c  The neutral carries the vector sum")
    fig.suptitle("Four-leg per-phase VV+VW: equal active powers, unequal currents", fontsize=14, y=.98)
    fig.text(.5, .025, "Full PV · source 1.03 pu · 1% VUF · fixed phase/DC sizing · explicit neutral · quasi-static counterfactual", ha="center", fontsize=10)
    fig.tight_layout(rect=(0, .07, 1, .92))
    destination.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(destination / f"four_leg_phase_control.{ext}", dpi=180,
                    metadata={"Date": None} if ext == "svg" else {})
    plt.close(fig)
    lines = ["# Australian four-leg research-control results", "",
             "Generated from eleven accepted saved-point audits. This is an explicit four-wire topology and policy counterfactual; phase and DC sizing are unchanged. Neutral rating is 10 A RMS and neutral filter R/X equals the phase arm. These assumptions do not establish Australian connection conformity or dynamic stability.", "",
             "| Case | POC P (kW) | POC Q (kvar) | Maximum PN voltage (V) | Maximum neutral (A RMS) |",
             "|---|---:|---:|---:|---:|"]
    for c in cases:
        lines.append(f"| {c['case_id']} | {float(c['p_poc_W'])/1000:.3f} | {float(c['q_poc_var'])/1000:.3f} | {float(c['vmax_sensed_V']):.3f} | {float(c['max_neutral_current_A']):.3f} |")
    headline = next(c for c in cases if c["case_id"] == "matched_s2_p3")
    lines += ["", f"The full-PV, 1.03 pu, 1% VUF point exports {float(headline['p_poc_W'])/1000:.3f} kW with maximum PN voltage {float(headline['vmax_sensed_V']):.3f} V. All phase Q requests are delivered, while a shared source couples active dispatch. The maximum voltage still exceeds the illustrative 253 V planning criterion: per-phase VV is not a voltage-constrained OPF or a guarantee of compliance.", "",
              f"The independent saved-phasor audit recomputes {audit['physical_checks_recomputed']:,} hardware/policy checks for {audit['physical_devices']} energized inverter records. Fixed-PQ replay checks network injection consistency in the same engine; it is not an external circuit or EMT validation.", "",
              "![Phase voltages, reactive requests and neutral-current vector sum](figures/four_leg_phase_control.png)", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path, nargs="?", default=STUDY / "results/four_leg")
    parser.add_argument("--figures", type=Path, default=STUDY / "figures")
    parser.add_argument("--report", type=Path, default=STUDY / "four_leg_results.md")
    args = parser.parse_args()
    args.report.write_text(generate(args.folder, args.figures))
