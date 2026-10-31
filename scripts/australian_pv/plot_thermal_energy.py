"""Audit mixed-controller CSVs and publish the fixed-sizing model comparison."""
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "studies/australian_pv"


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def audit_records(folder):
    env = json.loads((folder / "environment.json").read_text())
    truth = json.loads((folder / "audit.json").read_text())
    if not truth["all_planned_accepted"]:
        raise ValueError("mixed campaign incomplete or rejected")
    for file, digest in env["source_sha256"].items():
        if hashlib.sha256((REPO / file).read_bytes()).hexdigest() != digest:
            raise ValueError(f"stale source: {file}")
    cases = read_csv(folder / "cases.csv")
    devices = read_csv(folder / "devices.csv")
    if len(cases) != truth["planned"] or len({r["case_id"] for r in cases}) != len(cases) or {r["case_id"] for r in cases} != set(truth["cases"]):
        raise ValueError("case coverage mismatch")
    def close(left, right):
        if not math.isfinite(float(left)) or abs(float(left)-right) > 1e-5:
            raise ValueError("CSV metric differs from saved-phasor evidence")
    for row in cases:
        expected = truth["cases"][row["case_id"]]
        if row["accepted"] != "true" or not expected["accepted"]:
            raise ValueError("unaccepted case")
        for key in ("p_poc_W", "q_poc_var", "vmax_sensed_V", "maximum_capacitor_model_A"):
            close(row[key], expected[key])
    ids = sorted(truth["cases"]["fundamental_s1"]["devices"])
    labels = {ident: env["controllers"][i % 4] for i, ident in enumerate(ids)}
    seen = set()
    for row in devices:
        key = (row["case_id"], row["device_id"])
        if key in seen or row["device_id"] not in labels or row["controller"] != labels[row["device_id"]]:
            raise ValueError("device/controller coverage mismatch")
        seen.add(key)
        if row["accepted"] != "true":
            raise ValueError("unaccepted device")
        if row["dc_energy_model"] != ("nonlinear" if row["case_id"].startswith("energy_") else "linearized"):
            raise ValueError("DC model label mismatch")
        if row["pwm_strategy"] != ("NONE" if row["case_id"].startswith("fundamental_") else "CENTERED"):
            raise ValueError("PWM label mismatch")
        active = truth["cases"][row["case_id"]]["devices"]
        if row["model_present"] == "false":
            if row["device_id"] in active or float(row["p_poc_W"]) != 0 or float(row["q_poc_var"]) != 0:
                raise ValueError("disconnected device mismatch")
            continue
        expected = active[row["device_id"]]
        for key in ("p_poc_W", "q_poc_var", "p_dc_W", "p_available_W", "vmax_sensed_V", "capacitor_thermal_current_A"):
            close(row[key], expected[key])
        key = "capacitor_carrier_audit_current_A"
        if expected[key] is None:
            if row[key] != "":
                raise ValueError("unexpected PWM audit value")
        else:
            close(row[key], expected[key])
    if seen != {(c["case_id"], i) for c in cases for i in ids}:
        raise ValueError("missing device")
    return cases, devices, truth


def publish():
    folder = ROOT / "results/thermal_energy"
    cases, devices, truth = audit_records(folder)
    emt_folder = ROOT / "results/thermal_energy_emt"
    emt = json.loads((emt_folder / "results.json").read_text())
    receipt = json.loads((emt_folder / "audit.json").read_text())
    if not receipt["accepted"] or receipt["results_sha256"] != hashlib.sha256((emt_folder / "results.json").read_bytes()).hexdigest():
        raise ValueError("EMT audit required")
    allocated = [c for ident, c in emt["cases"].items() if ident != "nonlinear_boundary" and not ident.startswith("calibration_")]
    if len(allocated) != 27 or not all(c["modes"]["switched"]["finest"]["all_engineering_checks_passed"] for c in allocated):
        raise ValueError("all 27 PWM-allocated EMT points must satisfy engineering limits")
    baseline = json.loads((ROOT / "results/emt/results.json").read_text())
    old_receipt = json.loads((ROOT / "results/emt/audit.json").read_text())
    if not old_receipt["accepted"] or old_receipt["results_sha256"] != hashlib.sha256((ROOT / "results/emt/results.json").read_bytes()).hexdigest():
        raise ValueError("matched baseline EMT audit required")
    data = {r["case_id"]: r for r in cases}
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "svg.hashsalt": "thermal-energy-mixed"})
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), layout="constrained")
    colors = ["#777777", "#c36b2d", "#207b8d"]
    models = ["fundamental", "pwm", "energy"]
    labels = ["Fundamental only", "PWM thermal", "PWM + nonlinear energy"]
    ax = axes[0, 0]
    for j, (model, label, color) in enumerate(zip(models, labels, colors)):
        ax.bar(np.arange(2)+(j-1)*.24, [float(data[f"{model}_s{s}"]["p_poc_W"])/1000 for s in (1, 2)], .24, color=color, label=label)
    ax.axhline(37, color="#888", ls=":", label="Available DC: 37 kW")
    ax.axhline(0, color="#aaa", lw=.6)
    ax.set(xticks=[0, 1], xticklabels=["1.00 pu source", "1.03 pu source"], ylabel="Fleet POC net active injection (kW)", ylim=(-2, 45), title="a  Thermal feasibility changes the dispatch")
    ax.annotate("≈0 kW under reactive priority", xy=(1, 0), xytext=(1, 5),
                ha="center", fontsize=8, arrowprops={"arrowstyle": "->", "color": "#555"})
    ax.legend(frameon=False, fontsize=8)
    ax = axes[0, 1]
    groups = [[], [], [], []]
    for ident, case in emt["cases"].items():
        if ident.startswith("calibration_"):
            continue
        value = case["modes"]["switched"]["finest"]["capacitor_total_rms_A"]
        group = 1 if ident.startswith("mixed_energy_s1") else 2 if ident.startswith("mixed_energy_s2") else 3 if ident == "nonlinear_boundary" else 0
        groups[group].append(value)
    ax.bar(range(4), [max(g) for g in groups], color=[colors[2], colors[2], colors[2], colors[0]])
    ax.axhline(1, color="#b2182b", ls="--", label="Fixed 1 A RMS budget")
    ax.set(xticks=range(4), xticklabels=["Allocated fixtures", "Mixed 1.00 pu", "Mixed 1.03 pu", "LF-only boundary"], ylabel="Largest switched capacitor current (A RMS)", title="b  Allocated points are checked in EMT", ylim=(0, 5.2))
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1, 0]
    for model, label, color in zip(models, labels, colors):
        ax.plot([1, 1.03], [float(data[f"{model}_s{s}"]["vmax_sensed_V"]) for s in (1, 2)], "o-", color=color, label=label)
    ax.axhline(253, color="#888", ls="--", label="Illustrative 253 V criterion")
    ax.set(xlabel="Source positive-sequence voltage (pu)", ylabel="Highest customer PN voltage (V)", title="c  All four controller types share one feeder")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1, 1]
    old = abs(complex(*baseline["cases"]["binding_dc_capacitor"]["modes"]["averaged"]["finest"]["dc_2omega_peak"]))
    new = abs(complex(*emt["cases"]["nonlinear_boundary"]["modes"]["averaged"]["finest"]["dc_2omega_peak"]))
    ax.bar([0, 1], [(old-7.5)*1e6, (new-7.5)*1e6], color=[colors[0], colors[2]])
    ax.axhline(0, color="#b2182b", ls="--", label="7.5 V limit")
    ax.set(xticks=[0, 1], xticklabels=["Linearized allocation", "Nonlinear-energy allocation"], ylabel="Averaged EMT 2ω excess above 7.5 V (µV)", title="d  Energy bounds correct the small boundary error", ylim=(-10, 120))
    for ax in axes.flat:
        ax.grid(axis="y", alpha=.15)
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle("Mixed Australian controllers with fixed phase/DC sizing\nMean · worst phase · sequence · four-leg phase control; 12 sites with matched PV availability groups", fontsize=13)
    for ext in ("png", "svg"):
        path = ROOT / f"figures/thermal_energy_mixed.{ext}"
        fig.savefig(path, dpi=180, metadata={"Date": None} if ext == "svg" else {})
        if ext == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    plt.close(fig)
    lines = ["# Mixed-controller PWM and nonlinear-energy results", "",
             "All eight predeclared coupled cases are accepted. Each feeder contains mean, worst-phase, sequence and four-leg controllers simultaneously. Phase/DC sizing is unchanged. Each controller group receives the same 25/60/100% availability pattern; locations differ, so these mixed cases do not establish a causal ranking of controller laws.", "",
             "| Case | POC net P (kW) | POC Q (kvar) | Maximum PN voltage (V) | Maximum capacitor model (A RMS) |",
             "|---|---:|---:|---:|---:|"]
    for r in cases:
        lines.append(f"| {r['case_id']} | {float(r['p_poc_W'])/1000:.3f} | {float(r['q_poc_var'])/1000:.3f} | {float(r['vmax_sensed_V']):.3f} | {float(r['maximum_capacitor_model_A']):.6f} |")
    allocation_values = [c["modes"]["switched"]["finest"]["capacitor_total_rms_A"] for ident, c in emt["cases"].items() if ident != "nonlinear_boundary" and not ident.startswith("calibration_")]
    lines += ["", f"Saved-point reconstruction passes {truth['conditions']} conditions and {truth['physical_checks']} SI physical checks, plus same-engine fixed-PQ replay. The EMT matrix covers all 24 energized inverter points from both nonlinear-energy daytime mixed cases and four fixed-boundary references, at three meshes in averaged and switched modes. All 27 PWM-allocated references meet the declared switched engineering checks; their largest capacitor RMS is {max(allocation_values):.6f} A. The separate LF-only boundary is retained as a switching-thermal negative control.", "",
              f"At the matched ripple boundary, averaged time-domain 2ω voltage is {old:.9f} V for linearized allocation and {new:.9f} V for nonlinear-energy allocation. The specified 750 V mean is preserved. The nonlinear correction is small here; the large dispatch change is caused by enforcing the previously omitted PWM thermal component.", "",
              "The initial 5% allowance passed four balanced probes but exceeded 1 A at three four-leg mixed points (up to 1.0059 A). Those earlier points remain explicit negative controls. The default model allowance is now 10%, with sizing unchanged; it is not a universal bound on AC ripple or device dynamics. An averaged three-wire common-mode projection error exposed by the mixed cases was corrected and receives an independent KCL/energy regression.", "",
              "The mixed night case retains reactive service and imports loss power at zero PV. Partial disconnection removes four complete AC models. No 253 V constraint is imposed by these local policies; the plotted criterion is illustrative. EMT uses fixed POC voltages, ideal centered PWM, constant DC input power and fixed mean converter loss; feeder/service dynamics, LCL/split-link dynamics and semiconductor losses remain outside this validation.", "",
              "![Thermal dispatch, EMT stress, feeder voltages and nonlinear boundary correction](figures/thermal_energy_mixed.png)", "",
              "See [setup and reproduction](thermal_energy.md).", ""]
    (ROOT / "thermal_energy_results.md").write_text("\n".join(lines))


if __name__ == "__main__":
    publish()
