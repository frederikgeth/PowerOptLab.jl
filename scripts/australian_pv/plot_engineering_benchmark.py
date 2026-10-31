"""Publish audited matched export/service tradeoffs, binding limits and EMT gaps."""
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
CONDITIONS = ["nominal", "high_source", "rotated_unbalance", "night"]
LABELS = ["Nominal", "High source", "Rotated 2% U₂", "Night"]
SHORT = dict(mean_vv_vw="mean", worst_vv_vw="worst", sequence_vv_vw="sequence", phase_vv_vw="four-leg")


def equal_cell(cell, expected):
    if expected is None:
        return cell == ""
    if isinstance(expected, bool):
        return cell == str(expected).lower()
    if isinstance(expected, (int, float)):
        return math.isfinite(float(cell)) and abs(float(cell)-expected) <= 1e-5
    return cell == str(expected)


def audited_records(folder):
    env = json.loads((folder/"environment.json").read_text())
    truth = json.loads((folder/"audit.json").read_text())
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    if not truth["all_planned_accepted"] or truth["accepted"] != 8:
        raise ValueError("complete accepted matched benchmark required")
    if truth["environment_sha256"] != digest(folder/"environment.json"):
        raise ValueError("benchmark environment changed")
    if truth["audit_code_sha256"] != digest(REPO/"scripts/australian_pv/audit_engineering_benchmark.jl"):
        raise ValueError("stale saved-point audit")
    if env["network_sha256"] != digest(ROOT/"networks/LV3_55bus.bmopf.json"):
        raise ValueError("network input changed")
    for name in ("Project", "Manifest"):
        if env[name+"_sha256"] != digest(REPO/(name+".toml")):
            raise ValueError(f"{name} changed")
    for file, sha in env["source_sha256"].items():
        if digest(REPO/file) != sha:
            raise ValueError(f"stale benchmark source: {file}")
    for file, sha in truth["record_sha256"].items():
        if digest(folder/file) != sha:
            raise ValueError(f"benchmark record changed: {file}")
    expected = {(condition, mode) for condition in CONDITIONS for mode in ("local", "coordinated")}
    with (folder/"cases.csv").open() as stream:
        cases = list(csv.DictReader(stream))
    with (folder/"devices.csv").open() as stream:
        devices = list(csv.DictReader(stream))
    if {(r["case_id"], r["mode"]) for r in cases} != expected or len(cases) != 8:
        raise ValueError("case CSV coverage changed")
    for row in cases:
        evidence = truth["cases"][row["case_id"]+"_"+row["mode"]]
        for key, cell in row.items():
            if key in ("case_id", "mode"):
                continue
            value = evidence[key] if key in evidence else sum(d[key] for d in evidence["devices"])
            if not equal_cell(cell, value):
                raise ValueError(f"case CSV differs from audited phasors: {key}")
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
    return env, truth


def save_figure(fig, stem):
    for suffix in ("png", "svg"):
        path = ROOT/f"figures/{stem}.{suffix}"
        fig.savefig(path, dpi=180, metadata={"Date": None} if suffix == "svg" else {})
        if suffix == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    plt.close(fig)


def publish():
    folder = ROOT/"results/engineering_benchmark"
    env, truth = audited_records(folder)
    emt_folder = ROOT/"results/engineering_benchmark_emt"
    emt = audit_matched_emt(ROOT/"results/engineering_benchmark_emt_inputs.json", emt_folder, folder)
    simulation = json.loads((emt_folder/"results.json").read_text())
    references = json.loads((ROOT/"results/engineering_benchmark_emt_inputs.json").read_text())
    data = truth["cases"]
    network = json.loads((ROOT/"networks/LV3_55bus.bmopf.json").read_text())
    network_checks = []
    for condition in CONDITIONS:
        for mode in ("local", "coordinated"):
            point = json.loads((folder/f"points/{condition}_{mode}.json").read_text())
            result = audit_network(network, point, point["case"])
            if not result["accepted"]:
                raise ValueError("publication-time independent feeder rejection")
            network_checks.append(result)
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "svg.hashsalt": "engineering-benchmark"})
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), layout="constrained")
    colors = {"local": "#167d9a", "coordinated": "#c47832"}
    x = np.arange(4)
    for mode, offset in (("local", -.18), ("coordinated", .18)):
        axes[0, 0].bar(x+offset, [data[c+"_"+mode]["p_poc_W"]/1000 for c in CONDITIONS], .35,
                       color=colors[mode], label="Local Q priority" if mode == "local" else "Coordinated, Q free")
        axes[0, 1].bar(x+offset, [-data[c+"_"+mode]["q_poc_var"]/1000 for c in CONDITIONS], .35, color=colors[mode])
    for ax in axes[0]:
        ax.set_xticks(x, LABELS)
        ax.axhline(0, color="#777", lw=.7)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].set(title="a  Coordination recovers active injection", ylabel="Total PV POC injection (kW)")
    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[0, 1].set(title="b  Reactive absorption is given up", ylabel="Total PV reactive absorption (kvar)")
    baseline = {r["device_id"]: r for r in data["high_source_local"]["devices"]}
    coordinated = {r["device_id"]: r for r in data["high_source_coordinated"]["devices"]}
    ids = sorted(baseline)
    y = np.arange(12)
    recovery = [(coordinated[i]["p_poc_W"]-baseline[i]["p_poc_W"])/1000 for i in ids]
    axes[1, 0].barh(y, recovery, color=colors["coordinated"])
    labels = [f"{i.split('_')[1][2:]} · {SHORT[baseline[i]['local_controller']]} · {baseline[i]['p_available_W']/1000:g}" for i in ids]
    axes[1, 0].set(yticks=y, yticklabels=labels, xlabel="Additional POC injection (kW)",
                   title="c  High-source recovery at every site\nConnection · local controller · available DC kW")
    axes[1, 0].invert_yaxis()
    axes[1, 0].spines[["top", "right"]].set_visible(False)
    utilization = []
    for i in ids:
        l, c = baseline[i], coordinated[i]
        utilization.append([l["p_dc_W"]/l["p_available_W"], l["capacitor_thermal_A"],
            c["p_dc_W"]/c["p_available_W"], c["capacitor_thermal_A"],
            max((750-c["dc_minimum_bound_V"])/15, (c["dc_maximum_bound_V"]-750)/15),
            max(c["phase_peak_bound_A"], c["neutral_peak_bound_A"])/14])
    image = axes[1, 1].imshow(utilization, vmin=0, vmax=1, cmap="YlOrBr", aspect="auto")
    axes[1, 1].set(yticks=y, yticklabels=[""]*12, xticks=range(6),
        xticklabels=["Local\nDC source", "Local\ncapacitor", "OPF\nDC source", "OPF\ncapacitor", "OPF\nDC band", "OPF\nAC peak"],
        title="d  Same high-source sites: which limit is tight?")
    for j, row in enumerate(utilization):
        for k, value in enumerate(row):
            if value >= .995:
                axes[1, 1].text(k, j, "*", color="white", ha="center", va="center", fontsize=14)
    fig.colorbar(image, ax=axes[1, 1], label="Modeled fraction of limit (* ≥99.5%)", shrink=.85)
    fig.suptitle("Matched hardware and network limits · different service obligations\nExport recovery must be read alongside the released VV/VW support", fontsize=14)
    save_figure(fig, "engineering_coordination")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), layout="constrained")
    for mode in ("local", "coordinated"):
        group = [c for c in references["cases"] if c["parent_mode"] == mode]
        axes[0].scatter([c["parent_capacitor_model_A"] for c in group],
            [simulation["cases"][c["case_id"]]["modes"]["switched"]["finest"]["capacitor_total_rms_A"] for c in group],
            label="Local" if mode == "local" else "Coordinated", color=colors[mode], s=28)
    axes[0].plot([0, 1.05], [0, 1.05], ":", color="#777", label="Prediction = measurement")
    axes[0].axhline(1, color="#b2182b", ls="--", label="Fixed 1 A hardware budget")
    axes[0].annotate("Night OPF: prediction misses\nAC-induced DC ripple", xy=(.012, .30), xytext=(.23, .43),
                     color="#a52e27", fontsize=9, arrowprops={"arrowstyle": "->", "color": "#a52e27"})
    axes[0].set(xlabel="Calibrated thermal model (A RMS)", ylabel="Switched EMT capacitor current (A RMS)",
                title="a  Hardware passes; prediction is not a bound", xlim=(-.03, 1.07), ylim=(0, 1.08))
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    for j, condition in enumerate(CONDITIONS):
        for mode, offset in (("local", -.13), ("coordinated", .13)):
            selected = [simulation["cases"][c["case_id"]]["modes"]["switched"]["finest"] for c in references["cases"]
                        if c["parent_case_id"] == condition and c["parent_mode"] == mode]
            axes[1].plot([min(r["v_dc_min_V"] for r in selected), max(r["v_dc_max_V"] for r in selected)],
                         [j+offset]*2, "|-", color=colors[mode], lw=2, markersize=9)
    axes[1].axvspan(735, 765, color="#e3efe5", alpha=.5)
    for value in (735, 765):
        axes[1].axvline(value, color="#b2182b", ls="--", lw=1)
    axes[1].set(yticks=range(4), yticklabels=LABELS, ylim=(3.5, -.5), xlim=(733, 767),
                xlabel="Switched DC-voltage range across selected sites (V)", title="b  Absolute DC limits have margin")
    fig.suptitle("Independent EMT retains both hardware outcomes and model discrepancies\n32 fixed-POC references · 192 averaged/switched integrations", fontsize=13)
    save_figure(fig, "engineering_coordination_emt")

    lines = ["# Matched mixed-controller PF versus coordinated OPF", "",
        "The same 12 sites, physical hardware, PV availability, loads, source conditions, PWM/energy models and network ceilings "
        "are used in both modes. Local mean, worst-phase, sequence and four-leg controllers operate together with reactive priority. "
        "The coordinated objective maximizes total POC active injection with Q free and the original VV/VW commands released. "
        "This is a comparison of service choices under matched physical capability; it is not an equal-reactive-service dominance result or a global nonconvex bound.", "",
        "The common 735–765 V DC band and 14 A phase/neutral peak ceilings are declared study operating assumptions. "
        "No manufacturer peak rating is inferred from the unchanged 10 A RMS rating.", "",
        "![Export, reactive-service tradeoff, per-site recovery and limiting constraints](figures/engineering_coordination.png)", "",
        "Positive P/Q mean injection; negative values mean import/absorption. Available DC is 37 kW in each daytime case and zero at night. "
        "The 253 V planning ceiling is checked at all 55 LV buses, alongside native conductor/transformer limits. It is an illustrative study criterion, not a regional conformity claim.", "",
        "| Condition | Local P (kW) | OPF P (kW) | Recovery (kW) | Local Q (kvar) | OPF Q (kvar) | Local curtailed DC (kW) | OPF curtailed DC (kW) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for condition, label in zip(CONDITIONS, LABELS):
        l, c = data[condition+"_local"], data[condition+"_coordinated"]
        lines.append(f"| {label} | {l['p_poc_W']/1000:.3f} | {c['p_poc_W']/1000:.3f} | {(c['p_poc_W']-l['p_poc_W'])/1000:.3f} | "
            f"{l['q_poc_var']/1000:.3f} | {c['q_poc_var']/1000:.3f} | {l['curtailment_W']/1000:.3f} | {c['curtailment_W']/1000:.3f} |")
    max_cap = max(c["modes"]["switched"]["finest"]["capacitor_total_rms_A"] for c in simulation["cases"].values())
    lines += ["", "In the high-source case, the local reactive-priority allocations use almost all their capacitor thermal budget for reactive support. "
        "Releasing those commands permits active injection. Higher-availability sites then reach the capacitor budget, while low-availability sites reach their DC source ceiling. "
        "DC-voltage and AC peak bounds retain margin. Panel d normalizes DC excursions by the 15 V headroom around 750 V; its stars mark near-limit utilization, not a dual sensitivity certificate.", "",
        f"These cases are inverter-limited: the largest native conductor utilization is {100*max(r['maximum_conductor_utilization'] for r in network_checks):.2f}%, "
        f"and the largest LV phase-neutral voltage is {max(c['value'] for r in network_checks for c in r['checks'] if c['name'].startswith('planning_voltage/')):.3f} V. "
        "The comparison therefore does not establish a coordination benefit under active feeder congestion.", "",
        "The [complete per-site CSV](results/engineering_benchmark/devices.csv) includes requested and achieved Q, reactive-service error, DC curtailment, "
        "converter/filter losses and named binding limits. Night recovery is reduced loss import, with no PV curtailment.", "",
        "![EMT capacitor-model discrepancy and absolute DC-voltage margin](figures/engineering_coordination_emt.png)", "",
        f"All **{emt['references']} representative references and {emt['integrations']} integrations** pass the declared averaged/switched mechanism, mesh and hardware checks. "
        f"The largest switched capacitor RMS is **{max_cap:.6f} A** against the unchanged 1 A budget. The parent/circuit audit passes "
        f"{emt['reference_checks']} selection/availability checks and {emt['circuit_checks']} circuit-record checks.", "",
        "**A retained thermal-prediction failure:** the four coordinated night references predict approximately 0.012 A capacitor RMS, "
        "while switched EMT measures 0.276–0.316 A. The frozen-fundamental-current predictor omits the contribution generated by AC switching ripple. "
        "Its 10% calibrated allowance does not make it a conservative bound in this low-current regime. Those points pass the fixed 1 A hardware budget; "
        "their prediction failures remain explicit in `matched_audit.json`. The analytical AC-current and combined DC-voltage envelopes are separate from this calibrated thermal estimate.", "",
        "EMT holds each selected POC phasor boundary fixed, uses constant-power DC input and fixed mean converter loss, and includes finite-window DC drift. "
        "It does not simulate PLL/current/service-controller or feeder dynamics, establish periodic steady state, or certify all-time voltage excursions. "
        "The independent SI feeder audit reconstructs fundamental-frequency currents; carrier harmonics are not propagated through feeder ampacity constraints. "
        "The carrier approximation and its declared P/Q/sequence/modulation domain are documented in the [domain validation](engineering_domain.md); "
        "the retained LF-only DC-voltage failure is in the [boundary experiment](engineering_boundary.md).", "",
        "The solver initializes four-leg neutral channels at zero, then sequentially removes those temporary rows and retains accepted primal/dual starts. "
        "The final solve has full matched capability and the pure export objective. A separate [initialization ablation](results/engineering_benchmark_numerical_control/status.json) "
        "retains an `ITERATION_LIMIT` outcome under the same final constraints; successful SI feeder checks do not promote it.", "",
        f"Saved-point audits pass {truth['conditions']} record conditions, {truth['physical_checks']} SI inverter checks and {truth['network_conditions']} "
        "independent SI network checks. The campaign retains every raw allocation/OPF stage and point/source fingerprints. "
        "The authoritative PF/OPF inputs were generated from clean committed source; sizing is unchanged.", "",
        "Reproduce into new output paths:", "", "```sh",
        "julia --project=. scripts/australian_pv/run_engineering_benchmark.jl /tmp/engineering_benchmark python3",
        "julia --project=. scripts/australian_pv/audit_engineering_benchmark.jl /tmp/engineering_benchmark python3",
        "julia --project=. scripts/australian_pv/prepare_engineering_benchmark_emt.jl /tmp/engineering_benchmark /tmp/engineering_emt_inputs.json",
        "julia --project=. scripts/australian_pv/emt_circuit.jl /tmp/engineering_emt_inputs.json /tmp/engineering_emt",
        "python3 scripts/australian_pv/audit_engineering_benchmark_emt.py /tmp/engineering_emt_inputs.json /tmp/engineering_emt /tmp/engineering_benchmark",
        "python3 scripts/australian_pv/plot_engineering_benchmark.py", "```", "",
        "The plotting command checks the committed CSVs against saved-point evidence and reruns the independent matched EMT audit before publication. "
        "Python requires NumPy and Matplotlib. Generated figures and this report use committed inputs; new runs retain their own fingerprints.", ""]
    (ROOT/"engineering_benchmark.md").write_text("\n".join(lines))
    artifacts = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in
        [ROOT/"engineering_benchmark.md", *(ROOT/f"figures/{stem}.{suffix}" for stem in
         ("engineering_coordination", "engineering_coordination_emt") for suffix in ("png", "svg"))]}
    (folder/"publication.json").write_text(json.dumps(dict(artifacts=artifacts,
        plot_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        audit_sha256=hashlib.sha256((folder/"audit.json").read_bytes()).hexdigest(),
        matched_emt_audit_sha256=hashlib.sha256((emt_folder/"matched_audit.json").read_bytes()).hexdigest()), indent=2)+"\n")
    print("Published matched PF/OPF tradeoffs, per-inverter limits and explicit EMT prediction failures")


if __name__ == "__main__":
    publish()
