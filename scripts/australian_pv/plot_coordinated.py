"""Audit coordinated CSV metrics against saved-phasor evidence and plot comparisons."""
import csv
import hashlib
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "studies/australian_pv"
VOLTAGE_TOLERANCE = 1e-5+253*1e-8


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def close(left, right, tolerance=.0101):
    assert math.isfinite(left) and math.isfinite(right)
    assert abs(left-right) <= tolerance, (left, right)


def audit_records():
    path = ROOT / "results/coordinated"
    (path / "record_audit.json").write_text(json.dumps({"all_planned_accepted": False, "reason": "AUDIT_INCOMPLETE"})+"\n")
    env = json.loads((path / "environment.json").read_text())
    audit = json.loads((path / "audit.json").read_text())
    summary = json.loads((path / "summary.json").read_text())
    assert hashlib.sha256((REPO / "scripts/australian_pv/audit_coordinated.jl").read_bytes()).hexdigest() == audit["audit_code_sha256"]
    for file, digest in env["source_sha256"].items():
        assert hashlib.sha256((REPO / file).read_bytes()).hexdigest() == digest, file
    for name in ("Project", "Manifest"):
        assert hashlib.sha256((REPO / f"{name}.toml").read_bytes()).hexdigest() == env[f"{name}_sha256"]
    assert hashlib.sha256((ROOT / "networks/LV3_55bus.bmopf.json").read_bytes()).hexdigest() == env["network_sha256"]
    cases, devices = rows(path / "cases.csv"), rows(path / "devices.csv")
    assert len(cases) == summary["completed"]
    assert len({r["case_id"] for r in cases}) == len(cases)
    assert len(devices) == 12*len(cases)
    assert len({(r["case_id"], r["device_id"]) for r in devices}) == len(devices)
    for case in cases:
        cid = case["case_id"]
        assert (case["accepted"] == "true") == summary["cases"][cid]["accepted"]
        selected = [r for r in devices if r["case_id"] == cid]
        assert all(r["accepted"] == case["accepted"] for r in selected)
        if case["accepted"] != "true":
            continue
        evidence = audit["cases"][cid]
        assert evidence["accepted"] and case["pq_replay_accepted"] == "true"
        assert float(case["pq_replay_difference_V"]) <= .001
        for key in ("p_poc_W", "q_poc_var", "vmax_sensed_V"):
            close(float(case[key]), evidence[key], 1e-5)
        for row in selected:
            point = evidence["devices"][row["device_id"]]
            assert row["physical_checked"] == row["physical_feasible"] == "true"
            for key, value in point.items():
                close(float(row[key]), value)
            close(float(row["p_available_W"]), 5000*float(case["pv_fraction"]))
            close(float(row["p_dc_W"])-float(row["p_poc_W"]), sum(float(row[k]) for k in
                  ("converter_loss_W", "capacitor_loss_W", "filter_loss_W")))
        for key in ("p_poc_W", "q_poc_var"):
            close(float(case[key]), sum(float(r[key]) for r in selected))
    accepted_ids = {r["case_id"] for r in cases if r["accepted"] == "true"}
    for name in ("physical_checks.csv", "network_checks.csv"):
        for check in rows(path / name):
            if check["case_id"] not in accepted_ids:
                continue
            value, lo, hi, tolerance = [float(check[k]) for k in ("value", "lower", "upper", "tolerance")]
            assert math.isfinite(value) and math.isfinite(tolerance) and tolerance > 0
            assert lo-tolerance <= value <= hi+tolerance and check["passed"] == "true"
    complete = summary["all_planned_accepted"] and audit["all_planned_accepted"]
    assert len(accepted_ids) == summary["accepted"] == audit["accepted_saved_points"]
    record = {"all_planned_accepted": complete, "cases_checked": len(cases), "accepted": len(accepted_ids),
              "device_rows_checked": len(devices), "source_hashes_verified": len(env["source_sha256"]),
              "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scope": "CSV completeness, source freshness, physical/network inequalities and metrics compared with independent saved-phasor audit"}
    (path / "record_audit.json").write_text(json.dumps(record, indent=2)+"\n")
    assert complete, "Matched figures withheld: incomplete or unresolved benchmark cases"
    return env, cases, devices


def publish():
    env, cases, devices = audit_records()
    auto = ROOT / "results/voltage_aware_snapshots"
    autoenv = json.loads((auto / "environment.json").read_text())
    assert json.loads((auto / "summary.json").read_text())["all_planned_accepted"]
    for file, digest in autoenv["source_sha256"].items():
        assert env["source_sha256"][file] == digest
    assert autoenv["hardware"] == env["hardware"] and autoenv["control_settings"] == env["control_settings"]
    automatic = rows(auto / "cases.csv")
    selections = []
    for mode, label in (("mean_vv_vw", "Mean VV+VW"), ("worst_vv_vw", "Worst-phase VV+VW"), ("sequence_vv_vw", "Sequence VV+VW")):
        case = next(r for r in automatic if r["case_id"] == f"matched_s2_p3_{mode}")
        selections.append((label, float(case["p_poc_W"])/1000, float(case["q_poc_var"])/1000, float(case["vmax_sensed_V"])))
    for allocation, label in (("balanced_converter", "Coordinated balanced"), ("free", "Coordinated free I₂")):
        for suffix, cap in (("unrestricted", "No extra cap"), ("253V", "253 V cap")):
            case = next(r for r in cases if r["case_id"] == f"s2_p3_{allocation}_{suffix}")
            selections.append((f"{label}; {cap}", float(case["p_poc_W"])/1000, float(case["q_poc_var"])/1000, float(case["vmax_sensed_V"])))
    lines = ["# Generated coordinated OPF comparison", "", f"Executed benchmark source revision: `{env['git_revision']}`.", "",
             "All 24 predeclared cases pass the independent physical/policy, native SI equation, additional customer-voltage and fixed-PQ replay gates. "
             "Saved phasors are independently re-audited. These are feasible local nonconvex results, not global upper bounds.", "",
             "## Full PV, 1.03 pu source, 1% aligned upstream VUF", "",
             "Hardware is identical. The 253 V upper phase-to-neutral limit is an illustrative additional planning criterion. "
             "Autonomous VV+VW uses reactive priority; the benchmarks release those Q laws. Free I₂ additionally releases sequence allocation.", "",
             "| Operation | Fleet P (kW) | Fleet Q (kvar) | Maximum customer phase (V) | Meets 253 V criterion? |",
             "| --- | ---: | ---: | ---: | --- |"]
    for label, p, q, v in selections:
        lines.append(f"| {label} | {p:.4f} | {q:.4f} | {v:.4f} | {'Yes' if v <= 253+VOLTAGE_TOLERANCE else 'No'} |")
    lines += ["", "An autonomous point above 253 V cannot serve as a feasible export comparator under that added criterion. "
              "The balanced-current reference isolates coordination/free-Q capability; the free-I₂ reference includes extra actuation. "
              "The latter is research capability and does not establish regional permission for deliberate unbalanced injection.", "",
              "## Active-power accounting at the 253 V criterion", "",
              "| Allocation | Available DC (kW) | Unused DC (kW) | Converter/capacitor/filter loss (kW) | POC P (kW) |",
              "| --- | ---: | ---: | ---: | ---: |"]
    for allocation, label in (("balanced_converter", "Balanced"), ("free", "Free I₂")):
        selected = [r for r in devices if r["case_id"] == f"s2_p3_{allocation}_253V"]
        total = lambda key: sum(float(r[key]) for r in selected)/1000
        available, unused = total("p_available_W"), total("p_available_W")-total("p_dc_W")
        loss = sum(total(k) for k in ("converter_loss_W", "capacitor_loss_W", "filter_loss_W"))
        close(available-unused-loss, total("p_poc_W"), .0002)
        lines.append(f"| {label} | {available:.4f} | {unused:.4f} | {loss:.4f} | {total('p_poc_W'):.4f} |")
    lines += ["", "DC-link availability excludes upstream source/DC–DC losses. All cases are energized; these snapshots are not chronological energy or EMT results.", "",
              "![Coordinated capability and the additional voltage criterion](figures/coordinated_capability.png)", "",
              "A standalone SVG accompanies the PNG. See [setup, reproduction and limitations](coordinated.md).", ""]
    (ROOT / "coordinated_results.md").write_text("\n".join(lines))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "PowerOptLab coordinated benchmark"})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    x = np.arange(len(selections));colors = ["#999999"]*3+["#56b4e9"]*2+["#009e73"]*2
    labels = ["Mean", "Worst phase", "Sequence", "Balanced\nno extra cap", "Balanced\n253 V", "Free I₂\nno extra cap", "Free I₂\n253 V"]
    axes[0, 0].bar(x, [r[1] for r in selections], color=colors)
    axes[0, 0].set(ylabel="Fleet P at POC (kW)", title="Different feasible sets; hardware stays fixed", xticks=x, xticklabels=labels, ylim=(0, 64))
    for i, selection in enumerate(selections):
        if selection[3] > 253+VOLTAGE_TOLERANCE:
            axes[0, 0].text(i, selection[1]+.5, "×", ha="center", color="#b2182b", fontsize=14)
    axes[0, 0].text(.02, .05, "× fails the additional 253 V criterion", transform=axes[0, 0].transAxes, color="#b2182b")
    axes[0, 1].scatter(x, [r[3] for r in selections], c=colors, s=60)
    axes[0, 1].axhline(253., color="#b2182b", ls="--", label="Illustrative upper criterion")
    axes[0, 1].set(ylabel="Maximum customer phase–neutral RMS (V)", title="Voltage axis zoomed to show criterion violations", xticks=x, xticklabels=labels)
    axes[0, 1].legend(fontsize=8)
    for allocation, label, color, style in (("balanced_converter", "Balanced", "#0072b2", "-"), ("free", "Free I₂", "#009e73", "--")):
        cid = f"s2_p3_{allocation}_253V"
        selected = sorted([r for r in devices if r["case_id"] == cid], key=lambda r: r["device_id"])
        axes[1, 0].plot(range(1, 13), [float(r["vmax_sensed_V"]) for r in selected], "o", ls=style, color=color, label=label)
        axes[1, 1].plot(range(1, 13), [float(r["converter_negative_sequence_A"]) for r in selected], "o", ls=style, color=color, label=label)
    axes[1, 0].axhline(253., color="#b2182b", ls=":")
    axes[1, 0].set(ylabel="Highest phase–neutral RMS at customer (V)", title="Both coordinated cases meet the upper criterion")
    axes[1, 1].set(ylabel="Converter negative-sequence RMS (A)", title="Free I₂ is additional actuation, not a new sensing law")
    for ax in axes[1, :]:
        ax.set(xlabel="Customer index (sorted ID; not feeder distance)", xticks=range(1, 13))
        ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(axis="y", alpha=.2)
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle("Accepted local capability references, not global upper bounds\nFull PV, 1.03 pu source, 1% upstream VUF; identical three-leg hardware", fontsize=12)
    out = ROOT / "figures"
    for extension in ("png", "svg"):
        path = out / f"coordinated_capability.{extension}"
        fig.savefig(path, dpi=180, bbox_inches="tight", metadata={"Date": None} if extension == "svg" else None)
        if extension == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    plt.close(fig)
    print("Audited 24 cases and 288 device records; generated comparison and figures")


if __name__ == "__main__":
    publish()
