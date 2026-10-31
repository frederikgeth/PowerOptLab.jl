"""Publish the audited matched PF/OPF benchmark: figures and a report whose claims are computed from data.

Every quantitative sentence in the generated report is derived from the audited
records here, so prose cannot drift from the CSVs (roadmap R7).
"""
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from audit_engineering_benchmark_emt import audit as audit_matched_emt
from network_si_audit import audit_network

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO/"studies/australian_pv"
CONDITIONS = ["nominal", "high_source", "rotated_unbalance", "night", "oversized_dc"]
LABELS = ["Nominal", "High source", "Rotated 2% U₂", "Night", "Oversized DC"]
MODES = ["local", "opf_equal_service", "opf_vv_released"]
MODE_LABELS = {"local": "Local control", "opf_equal_service": "OPF, same volt-var service",
               "opf_vv_released": "OPF, volt-var released"}
COLORS = {"local": "#167d9a", "opf_equal_service": "#6a9f58", "opf_vv_released": "#c47832"}
PLANNING_V, SUSTAINED_V = 253., 258.
# Diagnostic columns: timing and model size are machine-specific and not audited.
DIAGNOSTIC = {"elapsed_s", "iterations", "nlp_variables", "nlp_constraints"}


def equal_cell(cell, expected):
    if expected is None:
        return cell == ""
    if isinstance(expected, bool):
        return cell == str(expected).lower()
    if isinstance(expected, (int, float)):
        return math.isfinite(float(cell)) and abs(float(cell)-expected) <= 1e-5
    return cell == str(expected)


def binding_counts(devices):
    """Same rule as engineering_binding_counts in run_engineering_benchmark.jl."""
    counts = {}
    for row in devices:
        labels = []
        for name in filter(None, row["binding_constraints"].split("|")):
            label = name[:-2] if name[-2:] in ("_1", "_2", "_3") and "_current" in name else name
            if label not in labels:
                labels.append(label)
        for label in labels or ["none"]:
            counts[label] = counts.get(label, 0)+1
    return counts


def lv_voltage_max(point, network):
    sources = {s["bus"] for s in network["voltage_source"].values()}
    return max(math.hypot(ts[p]["vr"]-ts["n"]["vr"], ts[p]["vi"]-ts["n"]["vi"])
               for bus, ts in point["bus"].items() if bus not in sources and "n" in ts
               for p in "abc" if p in ts)


def audited_records(folder):
    env = json.loads((folder/"environment.json").read_text())
    truth = json.loads((folder/"audit.json").read_text())
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    planned = len(CONDITIONS)*len(MODES)
    if not truth["all_planned_accepted"] or truth["accepted"] != planned:
        raise ValueError("complete accepted matched benchmark required")
    if truth["environment_sha256"] != digest(folder/"environment.json"):
        raise ValueError("benchmark environment changed")
    if truth["audit_code_sha256"] != digest(REPO/"scripts/australian_pv/audit_engineering_benchmark.jl"):
        raise ValueError("stale saved-point audit")
    if env["network_sha256"] != digest(ROOT/"networks/LV3_55bus.bmopf.json"):
        raise ValueError("network input changed")
    if env["Project_sha256"] != digest(REPO/"Project.toml"):
        raise ValueError("Project changed")
    if env["Manifest_sha256"] != digest(folder/"dependency_manifest.toml"):
        raise ValueError("recorded Manifest changed")
    for file, sha in env["source_sha256"].items():
        if digest(REPO/file) != sha:
            raise ValueError(f"stale benchmark source: {file}")
    for file, sha in truth["record_sha256"].items():
        if digest(folder/file) != sha:
            raise ValueError(f"benchmark record changed: {file}")
    if not all(a["accepted"] for a in truth["agreement"].values()) or set(truth["agreement"]) != set(CONDITIONS):
        raise ValueError("PF/OPF agreement record incomplete or rejected")
    network = json.loads((ROOT/"networks/LV3_55bus.bmopf.json").read_text())
    expected = {(condition, mode) for condition in CONDITIONS for mode in MODES}
    with (folder/"cases.csv").open() as stream:
        cases = list(csv.DictReader(stream))
    with (folder/"devices.csv").open() as stream:
        devices = list(csv.DictReader(stream))
    if {(r["case_id"], r["mode"]) for r in cases} != expected or len(cases) != planned:
        raise ValueError("case CSV coverage changed")
    for row in cases:
        key = row["case_id"]+"_"+row["mode"]
        evidence = truth["cases"][key]
        point = json.loads((folder/f"points/{key}.json").read_text())
        vmax = lv_voltage_max(point, network)
        sensed = max(d["vmax_sensed_V"] for d in evidence["devices"])
        derived = dict(maximum_lv_voltage_V=vmax, planning_253V_exceeded=vmax > PLANNING_V+1e-6,
                       sustained_258V_exceeded=sensed > SUSTAINED_V,
                       binding_counts=";".join(sorted(f"{k}={v}" for k, v in binding_counts(evidence["devices"]).items())))
        for column, cell in row.items():
            if column in ("case_id", "mode") or column in DIAGNOSTIC:
                continue
            if column in derived:
                value = derived[column]
            elif column in evidence:
                value = evidence[column]
            else:
                value = sum(d[column] for d in evidence["devices"])
            if not equal_cell(cell, value):
                raise ValueError(f"case CSV differs from audited phasors: {column}")
    seen = set()
    for row in devices:
        key = row["case_id"]+"_"+row["mode"]
        evidence = truth["cases"][key]
        saved = next(d for d in evidence["devices"] if d["device_id"] == row["device_id"])
        identity = (key, row["device_id"])
        if identity in seen:
            raise ValueError("duplicate device CSV row")
        seen.add(identity)
        for field, cell in row.items():
            if not equal_cell(cell, evidence["accepted"] if field == "accepted" else saved[field]):
                raise ValueError(f"device CSV differs from audited phasors: {field}")
    if seen != {(key, d["device_id"]) for key, c in truth["cases"].items() for d in c["devices"]}:
        raise ValueError("missing device CSV row")
    return env, truth, {(r["case_id"], r["mode"]): r for r in cases}


def save_figure(fig, stem):
    for suffix in ("png", "svg"):
        path = ROOT/f"figures/{stem}.{suffix}"
        fig.savefig(path, dpi=180, metadata={"Date": None} if suffix == "svg" else {})
        if suffix == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    plt.close(fig)


def fmt_counts(text):
    """'dc_source_power_max=12;planning_voltage=3' -> 'PV availability 12, planning voltage 3'."""
    names = {"dc_source_power_max": "PV availability", "dc_source_power_min": "zero DC source",
             "poc_apparent_power": "apparent power", "converter_current": "phase current",
             "grid_current": "phase current", "planning_voltage": "253 V planning",
             "capacitor_current": "capacitor current", "dc_voltage_ripple": "DC 2ω ripple",
             "modulation": "modulation", "none": "none"}
    parts = []
    for item in filter(None, text.split(";")):
        name, count = item.split("=")
        parts.append(f"{names.get(name, name)} {count}")
    return ", ".join(parts)


def legacy_rows():
    path = ROOT/"results/engineering_benchmark_legacy_fixture/cases.csv"
    if not path.exists():
        return None
    with path.open() as stream:
        return {(r["case_id"], r["mode"]): r for r in csv.DictReader(stream)}


def publish():
    folder = ROOT/"results/engineering_benchmark"
    env, truth, cases = audited_records(folder)
    emt_folder = ROOT/"results/engineering_benchmark_emt"
    inputs = ROOT/"results/engineering_benchmark_emt_inputs.json"
    emt = audit_matched_emt(inputs, emt_folder, folder)
    simulation = json.loads((emt_folder/"results.json").read_text())
    references = json.loads(inputs.read_text())
    data = truth["cases"]
    network = json.loads((ROOT/"networks/LV3_55bus.bmopf.json").read_text())
    utilization = 0.
    for condition in CONDITIONS:
        for mode in MODES:
            point = json.loads((folder/f"points/{condition}_{mode}.json").read_text())
            result = audit_network(network, point, point["case"], voltage_max=point["planning_voltage_max_V"])
            if not result["accepted"]:
                raise ValueError("publication-time feeder check rejected a saved point")
            utilization = max(utilization, result["maximum_conductor_utilization"])
    P = {(c, m): data[f"{c}_{m}"]["p_poc_W"] for c in CONDITIONS for m in MODES}
    Q = {(c, m): data[f"{c}_{m}"]["q_poc_var"] for c in CONDITIONS for m in MODES}
    V = {(c, m): float(cases[(c, m)]["maximum_lv_voltage_V"]) for c in CONDITIONS for m in MODES}
    available = {c: sum(d["p_available_W"] for d in data[f"{c}_local"]["devices"]) for c in CONDITIONS}
    limits = env["operating_limits"]

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "svg.hashsalt": "engineering-benchmark"})
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5), layout="constrained")
    x = np.arange(len(CONDITIONS));w = .26
    for k, mode in enumerate(MODES):
        offset = (k-1)*w
        axes[0, 0].bar(x+offset, [P[(c, mode)]/1000 for c in CONDITIONS], w, color=COLORS[mode], label=MODE_LABELS[mode])
        axes[0, 1].bar(x+offset, [-Q[(c, mode)]/1000 for c in CONDITIONS], w, color=COLORS[mode])
        axes[1, 0].plot(x, [V[(c, mode)] for c in CONDITIONS], "o-", color=COLORS[mode], label=MODE_LABELS[mode])
    axes[0, 0].scatter(x, [available[c]/1000 for c in CONDITIONS], marker="_", s=900, color="k", label="PV available", zorder=3)
    axes[1, 0].axhline(PLANNING_V, color="#b2182b", ls="--", lw=1, label="253 V planning limit")
    axes[1, 0].axhline(SUSTAINED_V, color="#7f0000", ls=":", lw=1, label="258 V sustained-overvoltage setting")
    for ax in (axes[0, 0], axes[0, 1], axes[1, 0]):
        ax.set_xticks(x, LABELS)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].axhline(0, color="#777", lw=.7);axes[0, 1].axhline(0, color="#777", lw=.7)
    axes[0, 0].set(title="a  Export: nearly identical unless the feeder limit binds", ylabel="Total PV export at POC (kW)")
    axes[0, 0].legend(frameon=False, fontsize=8, loc="upper left")
    axes[0, 1].set(title="b  Reactive absorption (volt-var service)", ylabel="Total PV absorption (kvar)")
    axes[1, 0].set(title="c  Highest LV phase-neutral voltage", ylabel="Voltage (V)")
    axes[1, 0].legend(frameon=False, fontsize=8, loc="lower left")
    names = ["dc_source_power_max", "poc_apparent_power", "converter_current", "planning_voltage",
             "capacitor_current", "dc_voltage_ripple"]
    short = ["PV\navailable", "Apparent\npower", "Phase\ncurrent", "253 V\nplanning", "Capacitor\ncurrent", "DC 2ω\nripple"]
    rows, ylabels = [], []
    for condition, label in zip(CONDITIONS, LABELS):
        if condition == "night":
            continue
        for mode in MODES:
            counts = dict(item.split("=") for item in filter(None, cases[(condition, mode)]["binding_counts"].split(";")))
            rows.append([int(counts.get(n, 0)) for n in names])
            ylabels.append(f"{label} · {['local', 'OPF same', 'OPF released'][MODES.index(mode)]}")
    image = axes[1, 1].imshow(rows, vmin=0, vmax=12, cmap="YlOrBr", aspect="auto")
    for j, row in enumerate(rows):
        for k, value in enumerate(row):
            if value:
                axes[1, 1].text(k, j, str(value), ha="center", va="center", fontsize=8,
                                color="white" if value > 7 else "black")
    axes[1, 1].set(xticks=range(len(names)), xticklabels=short, yticks=range(len(rows)), yticklabels=ylabels,
                   title="d  Sites at each limit (daytime cases, 12 sites)")
    axes[1, 1].tick_params(axis="y", labelsize=8)
    fig.colorbar(image, ax=axes[1, 1], label="Number of sites", shrink=.8)
    fig.suptitle("Matched hardware and network: local control, same-service OPF and volt-var-released OPF", fontsize=13)
    save_figure(fig, "engineering_coordination")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), layout="constrained")
    for mode in MODES:
        group = [c for c in references["cases"] if c["parent_mode"] == mode]
        axes[0].scatter([c["parent_capacitor_model_A"] for c in group],
            [simulation["cases"][c["case_id"]]["modes"]["switched"]["finest"]["capacitor_total_rms_A"] for c in group],
            label=MODE_LABELS[mode], color=COLORS[mode], s=24)
    top = limits["capacitor_rms_A"]*1.05
    axes[0].plot([0, top], [0, top], ":", color="#777", label="Model = switched EMT")
    axes[0].axhline(limits["capacitor_rms_A"], color="#b2182b", ls="--", label=f"{limits['capacitor_rms_A']:g} A rating")
    axes[0].set(xlabel="Modelled capacitor current bound (A RMS)", ylabel="Switched EMT capacitor current (A RMS)",
                title="a  Model bound versus switched circuit", xlim=(0, top), ylim=(0, top))
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    for j, condition in enumerate(CONDITIONS):
        for k, mode in enumerate(MODES):
            selected = [simulation["cases"][c["case_id"]]["modes"]["switched"]["finest"] for c in references["cases"]
                        if c["parent_case_id"] == condition and c["parent_mode"] == mode]
            axes[1].plot([min(r["v_dc_min_V"] for r in selected), max(r["v_dc_max_V"] for r in selected)],
                         [j+(k-1)*.2]*2, "|-", color=COLORS[mode], lw=2, markersize=8)
    for value in (limits["dc_min_V"], limits["dc_max_V"]):
        axes[1].axvline(value, color="#b2182b", ls="--", lw=1)
    axes[1].set(yticks=range(len(CONDITIONS)), yticklabels=LABELS, ylim=(len(CONDITIONS)-.5, -.5),
                xlim=(limits["dc_min_V"]-2, limits["dc_max_V"]+2),
                xlabel="Switched DC-link voltage range (V)", title="b  DC-link voltage within its band")
    agreement = truth["agreement"]
    dv = [max(agreement[c]["max_voltage_difference_V"], 1e-12) for c in CONDITIONS]
    dq = [max(agreement[c]["max_site_Q_difference_var"], 1e-12) for c in CONDITIONS]
    axes[2].bar(x-.18, dv, .36, color="#167d9a", label="max |ΔV| (V)")
    axes[2].bar(x+.18, dq, .36, color="#c47832", label="max site |ΔQ| (var)")
    axes[2].axhline(env["agreement_tolerance"]["voltage_V"], color="#167d9a", ls="--", lw=1, label="tolerance, V")
    axes[2].axhline(env["agreement_tolerance"]["reactive_var"], color="#c47832", ls="--", lw=1, label="tolerance, var")
    axes[2].set(yscale="log", xticks=x, title="c  PF/OPF agreement at equal service",
                ylabel="Local fixed point vs OPF difference")
    axes[2].set_xticks(x, LABELS, rotation=20)
    axes[2].legend(frameon=False, fontsize=8)
    fig.suptitle(f"Checks: {emt['references']} fixed-POC EMT references ({emt['integrations']} integrations) and the "
                 "simulation/optimisation agreement", fontsize=12)
    save_figure(fig, "engineering_coordination_emt")

    # ---- Report: every number and limit statement below is computed from the audited records.
    def kw(v): return f"{v/1000:.3f}"
    lines = ["# Matched benchmark: local control, same-service OPF and volt-var-released OPF", "",
        "*Generated by `scripts/australian_pv/plot_engineering_benchmark.py` from audited records; do not edit by hand.* "
        "Overview, fixture and limitations: [study README](README.md).", "",
        "Twelve 6 kVA inverters (mean, worst-phase, sequence and four-leg VV/VW controllers, three of each) on the CSIRO "
        "`LV3_55bus` feeder. Hardware, PV availability, loads and network limits are identical across modes. "
        f"Fixture: capacitor rating {limits['capacitor_rms_A']:g} A RMS, phase/neutral peak {limits['phase_peak_A']:g} A, "
        f"DC band {limits['dc_min_V']:g}–{limits['dc_max_V']:g} V; rated-capability gate "
        f"{'passed' if env['rated_capability_accepted'] else 'FAILED'}.", "",
        "![Export, reactive service, voltage and limiting constraints](figures/engineering_coordination.png)", "",
        "| Condition | PV available (kW) | Local P (kW) | Same-service OPF P (kW) | Released OPF P (kW) | "
        "Released − local (kW) | Local Q (kvar) | Released Q (kvar) | Max LV voltage local / same / released (V) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    for c, label in zip(CONDITIONS, LABELS):
        lines.append(f"| {label} | {kw(available[c])} | {kw(P[(c, 'local')])} | {kw(P[(c, 'opf_equal_service')])} | "
            f"{kw(P[(c, 'opf_vv_released')])} | {kw(P[(c, 'opf_vv_released')]-P[(c, 'local')])} | "
            f"{Q[(c, 'local')]/1000:.3f} | {Q[(c, 'opf_vv_released')]/1000:.3f} | "
            f"{V[(c, 'local')]:.2f} / {V[(c, 'opf_equal_service')]:.2f} / {V[(c, 'opf_vv_released')]:.2f} |")
    lines += ["", "Positive P/Q are injection; negative Q is absorption. Night export is negative (converter losses imported).", "",
        "**What limits each site** (number of sites at each limit; one site can be at several):", ""]
    for c, label in zip(CONDITIONS, LABELS):
        parts = [f"{['local', 'same-service OPF', 'released OPF'][k]}: {fmt_counts(cases[(c, m)]['binding_counts'])}"
                 for k, m in enumerate(MODES)]
        lines.append(f"- {label}: " + "; ".join(parts) + ".")
    daytime = [c for c in CONDITIONS if c != "night"]
    gains = [(P[(c, 'opf_vv_released')]-P[(c, 'local')])/P[(c, 'local')] for c in daytime]
    exceed = [c for c in CONDITIONS if cases[(c, "local")]["planning_253V_exceeded"] == "true"]
    sustained = [(c, m) for c in CONDITIONS for m in MODES if cases[(c, m)]["sustained_258V_exceeded"] == "true"]
    equal_gap = max(abs(P[(c, 'opf_equal_service')]-P[(c, 'local')]) for c in CONDITIONS if c not in exceed)
    lines += ["", "**Findings computed from these records:**", "",
        f"- Releasing volt-var changes daytime export by {100*min(gains):.2f}% to {100*max(gains):.2f}% of local export.",
        f"- Where local control meets 253 V, same-service OPF export differs from local export by at most {equal_gap:.3f} W.",
        ("- Local control exceeds the 253 V planning limit in: " +
         ", ".join(f"{LABELS[CONDITIONS.index(c)]} ({V[(c, 'local')]:.2f} V)" for c in exceed) +
         ". In these cases same-service OPF must curtail to meet the limit "
         f"({', '.join(kw(P[(c, 'opf_equal_service')]-P[(c, 'local')])+' kW' for c in exceed)}), while released OPF "
         f"changes export by {', '.join(kw(P[(c, 'opf_vv_released')]-P[(c, 'local')])+' kW' for c in exceed)} "
         "using less absorption.") if exceed else "- Local control meets the 253 V planning limit in every condition.",
        "- No site exceeds the 258 V sustained-overvoltage setting (snapshot treated as its 10-minute average)." if not sustained
        else "- Sites exceed the 258 V sustained-overvoltage setting in: " + ", ".join(f"{c}/{m}" for c, m in sustained) + ".",
        f"- Largest native conductor utilisation across all points: {100*utilization:.1f}%.", ""]
    worst = {k: max(agreement[c][k] for c in CONDITIONS) for k in
             ("max_voltage_difference_V", "max_site_P_difference_W", "max_site_Q_difference_var")}
    lines += ["**Simulation/optimisation agreement.** Without the 253 V row, the same-service OPF reproduces the local "
        f"fixed point in all {len(CONDITIONS)} conditions: largest bus-voltage difference {worst['max_voltage_difference_V']:.1e} V, "
        f"site P {worst['max_site_P_difference_W']:.1e} W, site Q {worst['max_site_Q_difference_var']:.1e} var "
        f"(tolerances {env['agreement_tolerance']['voltage_V']:g} V, {env['agreement_tolerance']['power_W']:g} W, "
        f"{env['agreement_tolerance']['reactive_var']:g} var).", ""]
    ratios, counts = [], []
    for c in CONDITIONS:
        history = json.loads((folder/f"{c}_local_status.json").read_text())["history"]
        v = [r["voltage_residual_inf_V"] for r in history if r.get("voltage_residual_inf_V") is not None]
        counts.append(len(history))
        ratios += [b/a for a, b in zip(v, v[1:]) if a > 0 and b > 0]
    lines += [f"**Local fixed point.** Converged in {min(counts)}–{max(counts)} iterations without relaxation. The largest "
        f"ratio between successive voltage residuals is {max(ratios):.3f} (empirical contraction indicator on this stiff feeder; "
        "not a uniqueness or stability proof).", ""]
    legacy = legacy_rows()
    if legacy:
        lines += ["**Sensitivity: legacy fixture** (1 A capacitor budget, 14 A peak; fails the rated-capability gate). "
            "Shown only to explain the earlier headline result; it does not represent nameplate-capable hardware.", "",
            "| Condition | Local P (kW) | Released OPF P (kW) | Released − local (kW) |", "|---|---:|---:|---:|"]
        for c, label in zip(CONDITIONS, LABELS):
            l, r = legacy.get((c, "local")), legacy.get((c, "opf_vv_released"))
            if l and r:
                lines.append(f"| {label} | {kw(float(l['p_poc_W']))} | {kw(float(r['p_poc_W']))} | "
                             f"{kw(float(r['p_poc_W'])-float(l['p_poc_W']))} |")
        lines.append("")
    cov = emt["thermal_prediction_coverage"]
    margin = min(v["predicted_A"]-v["switched_A"] for v in cov.values())
    peak = max(v["switched_A"] for v in cov.values())
    lines += ["![EMT enclosure, DC-link voltage and PF/OPF agreement](figures/engineering_coordination_emt.png)", "",
        f"**EMT check.** {emt['references']} references (the highest modelled capacitor stress per condition, mode and controller "
        f"group), each integrated averaged and switched at three step sizes ({emt['integrations']} integrations). "
        f"{'All' if emt['all_thermal_predictions_enclose'] else 'Not all'} modelled capacitor-current bounds enclose the switched "
        f"measurement (largest switched value {peak:.3f} A against the {limits['capacitor_rms_A']:g} A rating; smallest margin "
        f"{margin:.3f} A). The EMT circuit is separately coded in this repository and holds the POC voltage fixed; it is not a "
        "PSCAD/PLECS comparison or a hardware test, and it omits PLL, controller and feeder dynamics. The thermal bound was "
        "revised after earlier night-time misses and is checked here with the same selection rule, so this is not an "
        "out-of-sample validation.", "",
        "**Solve times** (diagnostic, this machine): " + "; ".join(
            f"{MODE_LABELS[m]} median {np.median([float(cases[(c, m)]['elapsed_s']) for c in CONDITIONS]):.1f} s"
            for m in MODES) +
        f". OPF size: {cases[('nominal', 'opf_vv_released')]['nlp_variables']} variables, "
        f"{cases[('nominal', 'opf_vv_released')]['nlp_constraints']} constraints for 12 sites. "
        "Scaling beyond one feeder is not yet measured.", "",
        "Per-site records: [devices.csv](results/engineering_benchmark/devices.csv) "
        "(requested and delivered Q, curtailment, losses, named limits). Reproduce: see the [study README](README.md).", ""]
    (ROOT/"engineering_benchmark.md").write_text("\n".join(lines))
    artifacts = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in
        [ROOT/"engineering_benchmark.md", *(ROOT/f"figures/{stem}.{suffix}" for stem in
         ("engineering_coordination", "engineering_coordination_emt") for suffix in ("png", "svg"))]}
    (folder/"publication.json").write_text(json.dumps(dict(artifacts=artifacts,
        plot_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        audit_sha256=hashlib.sha256((folder/"audit.json").read_bytes()).hexdigest(),
        matched_emt_audit_sha256=hashlib.sha256((emt_folder/"matched_audit.json").read_bytes()).hexdigest()), indent=2)+"\n")
    print("Published matched PF/OPF benchmark report and figures")


if __name__ == "__main__":
    publish()
