"""Audit saved voltage-aware campaigns and render scientific teaching figures.

Uses CSV/JSON records, never reruns or silently discards a failed solve. Plotting
requires both complete, accepted campaigns with unchanged executed source hashes.
"""
import csv
import hashlib
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "studies/australian_pv"
MODES = ["mean_vv_vw", "worst_vv_vw", "sequence_vv_vw"]
LABELS = ["Mean", "Worst phase", "Sequence"]
COLORS = ["#0072b2", "#d55e00", "#009e73"]
MARKERS = ["o", "s", "x"]


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def number(row, key):
    return float(row[key])


def near(left, right, atol=1e-8):
    assert math.isfinite(left) and math.isfinite(right)
    assert abs(left - right) <= atol, (left, right, atol)


def audit(name):
    path = ROOT / "results" / f"voltage_aware_{name}"
    env = json.loads((path / "environment.json").read_text())
    summary = json.loads((path / "summary.json").read_text())
    for relative, digest in env["source_sha256"].items():
        assert hashlib.sha256((REPO / relative).read_bytes()).hexdigest() == digest, relative
    for file in ("Project", "Manifest"):
        assert hashlib.sha256((REPO / f"{file}.toml").read_bytes()).hexdigest() == env[f"{file}_sha256"]
    assert hashlib.sha256((ROOT / "networks/LV3_55bus.bmopf.json").read_bytes()).hexdigest() == env["network_sha256"]
    planned = {r["case_id"]: r for r in env["planned_cases"]}
    cases = read_csv(path / "cases.csv")
    devices = read_csv(path / "devices.csv")
    checks = read_csv(path / "physical_checks.csv")
    iterations = read_csv(path / "iterations.csv")
    actual = {r["case_id"]: r for r in cases}
    assert len(actual) == len(cases) and actual.keys() <= planned.keys()
    assert summary["planned"] == len(planned) and summary["completed"] == len(cases)
    failed = [r["case_id"] for r in cases if r["accepted"] != "true"]
    assert summary["unresolved_case_ids"] == failed
    assert summary["accepted"] == len(cases) - len(failed)
    complete = len(cases) == len(planned) and not failed
    assert summary["all_planned_accepted"] == complete
    assert len(devices) == 12 * len(cases)
    assert len({(r["case_id"], r["device_id"]) for r in devices}) == len(devices)
    for case in cases:
        cid = case["case_id"]
        plan = planned[cid]
        rows = [r for r in devices if r["case_id"] == cid]
        physical = [r for r in checks if r["case_id"] == cid]
        history = [r for r in iterations if r["case_id"] == cid]
        assert len(history) == int(case["iterations"])
        assert all(r["accepted"] == case["accepted"] for r in rows)
        assert all(r["operating_mode"] == plan["operating_mode"] for r in rows)
        if case["accepted"] != "true":
            continue  # retained candidates may have no finite primal quantities
        near(sum(number(r, "p_poc_W") for r in rows), number(case, "p_poc_W"))
        near(sum(number(r, "q_poc_var") for r in rows), number(case, "q_poc_var"))
        bus = json.loads((path / "network" / f"{cid}.json").read_text())
        for row in rows:
            near(number(row, "p_available_W"), 5000 * plan["pv_fraction"])
            assert row["operating_mode"] == plan["operating_mode"]
            ts = bus[row["bus"]]
            v = [complex(ts[p]["vr"], ts[p]["vi"]) - complex(ts["n"]["vr"], ts["n"]["vi"]) for p in "abc"]
            a = complex(-.5, math.sqrt(3)/2)
            u1 = (v[0]+a*v[1]+a*a*v[2])/3
            u2 = (v[0]+a*a*v[1]+a*v[2])/3
            near(number(row, "vmin_sensed_V"), min(map(abs, v)))
            near(number(row, "vmax_sensed_V"), max(map(abs, v)))
            near(number(row, "vuf_percent"), 100*abs(u2)/abs(u1))
            if row["model_present"] == "true":
                near(number(row, "p_dc_W")-number(row, "p_poc_W"),
                     number(row, "converter_loss_W")+number(row, "filter_loss_W"))
            else:
                assert row["p_dc_W"] == row["physical_checked"] == row["physical_feasible"] == ""
                near(number(row, "p_poc_W"), 0.)
                near(number(row, "q_poc_var"), 0.)
        assert case["coupling_accepted"] == case["pq_replay_accepted"] == "true"
        assert case["network_checks_passed"] == "true" and int(case["network_check_count"]) > 0
        assert number(case, "network_max_tolerance_ratio") <= 1.
        assert number(case, "pq_replay_difference_V") <= env["pq_replay_tolerance_V"]
        active = plan["operating_mode"] == "energized"
        assert int(case["physical_device_count"]) == (12 if active else 0)
        assert int(case["failed_physical_devices"]) == 0
        assert bool(physical) == active
        if active:
            assert all(r["physical_checked"] == r["physical_feasible"] == "true" for r in rows)
            assert history[-1]["reason"] == "CONVERGED"
            assert number(history[-1], "voltage_residual_inf_V") <= env["voltage_tolerance_V"]
            assert number(history[-1], "current_residual_inf_A") <= env["current_tolerance_A"]
        else:
            assert not history and case["reason"] == "DISCONNECTED_POWER_FLOW"
        for check in physical:
            value, lo, hi, tol = [number(check, k) for k in ("value", "lower", "upper", "tolerance")]
            assert math.isfinite(value) and math.isfinite(tol) and tol > 0
            assert lo-tol <= value <= hi+tol and check["passed"] == "true", check
    energy = read_csv(path / "energy.csv") if (path / "energy.csv").exists() else []
    for row in energy:
        group, mode = row["group"], row["controller"]
        plan = [r for r in planned.values() if r["group"] == group and r["controller"] == mode]
        selected_cases = [r for r in cases if r["group"] == group and r["controller"] == mode]
        accepted = len(selected_cases) == len(plan) and all(r["accepted"] == "true" for r in selected_cases)
        assert (row["accepted"] == "true") == accepted
        if not accepted:
            assert row["net_poc_kWh"] == ""
            continue
        near(sum(r["duration_h"] for r in plan), 24.)
        ids = {r["case_id"] for r in plan}
        selected = [r for r in devices if r["case_id"] in ids]
        def integrate(f):
            return sum(planned[r["case_id"]]["duration_h"]*f(r) for r in selected)/1000
        expected = {
            "available_kWh": integrate(lambda r: number(r, "p_available_W")),
            "export_kWh": integrate(lambda r: max(0., number(r, "p_poc_W"))),
            "import_kWh": integrate(lambda r: max(0., -number(r, "p_poc_W"))),
            "net_poc_kWh": integrate(lambda r: number(r, "p_poc_W")),
            "unused_dc_kWh": integrate(lambda r: number(r, "p_available_W")-(number(r, "p_dc_W") if r["model_present"] == "true" else 0.)),
            "converter_loss_kWh": integrate(lambda r: number(r, "converter_loss_W") if r["model_present"] == "true" else 0.),
            "filter_loss_kWh": integrate(lambda r: number(r, "filter_loss_W") if r["model_present"] == "true" else 0.),
            "reactive_support_kvarh": integrate(lambda r: number(r, "q_poc_var")),
        }
        for key, value in expected.items():
            near(number(row, key), value)
        near(expected["available_kWh"], 405.)
        near(expected["available_kWh"]-expected["net_poc_kWh"],
             expected["unused_dc_kWh"]+expected["converter_loss_kWh"]+expected["filter_loss_kWh"])
        assert abs(number(row, "balance_residual_Wh")) <= .001
    report = {"source_hashes_verified": len(env["source_sha256"]), "planned": len(planned),
              "record_audit_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "completed": len(cases), "accepted": len(cases)-len(failed), "all_planned_accepted": complete,
              "physical_check_rows": len(checks), "device_rows": len(devices), "energy_groups": len(energy),
              "unresolved_case_ids": failed, "missing_case_ids": sorted(planned.keys()-actual.keys()),
              "audit_scope": "saved-record consistency, source freshness, physical inequalities, terminal phasors and independent duration-weighted energy recomputation"}
    (path / "record_audit.json").write_text(json.dumps(report, indent=2)+"\n")
    print(name, report)
    return env, cases, devices, iterations, energy, complete


def figures(snapshots, chronology):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    out = ROOT / "figures"
    out.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    def save(fig, name):
        for extension in ("svg", "png"):
            fig.savefig(out / f"{name}.{extension}", dpi=180, bbox_inches="tight")
        plt.close(fig)
    _, cases, _, iterations, _, _ = snapshots
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    for mode, label, color, marker in zip(MODES, LABELS, COLORS, MARKERS):
        rows = [r for r in iterations if r["case_id"] == f"matched_s2_p3_{mode}"]
        axes[0].semilogy([int(r["iteration"]) for r in rows], [number(r, "voltage_residual_inf_V") for r in rows], marker=marker, color=color, label=label, markerfacecolor="none")
        axes[1].plot([int(r["iteration"]) for r in rows], [int(r["failed_physical_devices"]) for r in rows], marker=marker, color=color, label=label, markerfacecolor="none")
    axes[0].axhline(1e-7, color=".4", ls="--", label="Voltage stopping threshold")
    axes[0].set(ylabel="Voltage-map residual (V)", title="Numerical feedback closure")
    axes[1].set(ylabel="Devices failing the physical audit", ylim=(-.4, 12.8), yticks=[0, 4, 8, 12], title="Virtual probes can violate hardware limits")
    for ax in axes:
        ax.set_xlabel("Numerical iteration (not time)")
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    fig.suptitle("Fixed hardware, full PV, 1.03 pu source, 1% upstream VUF\nOnly the verified equilibrium is a dispatch result", fontsize=12)
    save(fig, "voltage_aware_convergence")

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
    for mode, label, color, marker in zip(MODES, LABELS, COLORS, MARKERS):
        baseline = next(r for r in cases if r["case_id"] == f"matched_s2_p3_{mode}")
        selected = [baseline]+[next(r for r in cases if r["case_id"] == f"angle{i}_{mode}") for i in (1, 2, 3)]+[baseline]
        for ax, key, scale in zip(axes.flat, ("p_poc_W", "q_poc_var", "vmax_sensed_V", "max_vuf_percent"), (1000, 1000, 1, 1)):
            ax.plot([0, 90, 180, 270, 360], [number(r, key)/scale for r in selected], marker=marker, color=color, label=label, markerfacecolor="none")
    for ax, ylabel, title in zip(axes.flat, ("Fleet P at POC (kW)", "Fleet Q at POC (kvar)", "Maximum sensed phase voltage (V)", "Maximum customer VUF (%)"), ("Active export", "Reactive exchange (absorption negative)", "Worst customer phase", "Customer unbalance")):
        ax.set(xlabel="Upstream U₂/U₁ relative angle (degrees)", ylabel=ylabel, title=title, xticks=[0, 90, 180, 270, 360])
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].text(.98, .07, "P axis zoomed; mean and sequence nearly coincide", transform=axes[0, 0].transAxes, ha="right", fontsize=8)
    fig.suptitle("1% upstream VUF has different effects as its angle changes\n|U₁|, hardware, available PV and loads stay fixed; lines guide sampled points; 360° repeats 0°", fontsize=12)
    save(fig, "voltage_aware_unbalance_angle")

    _, cases, _, _, energy, _ = chronology
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for mode, label, color in zip(MODES, LABELS, COLORS):
        for group, style in (("energized_night", "-"), ("disconnected_night", "--")):
            rows = sorted([r for r in cases if r["controller"] == mode and r["group"] == group], key=lambda r: int(r["step"]))
            edges = [number(r, "start_h") for r in rows]+[24.]
            for ax, key in zip(axes[:, 0], ("p_poc_W", "q_poc_var")):
                ax.stairs([number(r, key)/1000 for r in rows], edges, color=color, linestyle=style,
                          label=f"{label}; {'energized' if group == 'energized_night' else 'open at night'}")
    axes[0, 0].stairs([0, 9, 39, 60, 27, 0], [0, 6, 9, 12, 15, 18, 24], color=".5", linestyle=":", label="Available fleet PV")
    axes[0, 0].set(ylabel="Fleet P at POC (kW)", title="Synthetic day: constant conditions in each interval")
    axes[1, 0].set(ylabel="Fleet Q at POC (kvar)", title="Night operation retains voltage support and losses")
    axes[0, 0].legend(fontsize=7, ncol=2)
    for ax in axes[:, 0]:
        ax.set(xlabel="Hour of illustrative day", xlim=(0, 24), xticks=[0, 6, 9, 12, 15, 18, 24])
        ax.axhline(0, color=".6", linewidth=.6)
        ax.grid(alpha=.2)
    ordered = [next(r for r in energy if r["controller"] == m and r["group"] == g) for m in MODES for g in ("energized_night", "disconnected_night")]
    labels = [f"{label}\n{mode}" for label in LABELS for mode in ("Energized", "Open at night")]
    x = np.arange(6)
    bottom = np.zeros(6)
    for key, label, color in (("net_poc_kWh", "Net POC energy", "#56b4e9"), ("unused_dc_kWh", "Unused DC availability", "#bdbdbd"), ("converter_loss_kWh", "Converter loss", "#d55e00"), ("filter_loss_kWh", "Filter loss", "#e69f00")):
        y = np.array([number(r, key) for r in ordered])
        axes[0, 1].bar(x, y, bottom=bottom, label=label, color=color)
        bottom += y
    axes[0, 1].axhline(405, color=".3", ls=":")
    axes[0, 1].set(ylabel="Energy (kWh)", xticks=x, xticklabels=labels, title="405 kWh available: complete modeled energy budget")
    axes[0, 1].tick_params(axis="x", labelsize=8)
    axes[0, 1].legend(fontsize=8)
    axes[1, 1].bar(x, [number(r, "import_kWh") for r in ordered], color=[c for c in COLORS for _ in range(2)])
    axes[1, 1].set(ylabel="Inverter AC import (kWh)", xticks=x, xticklabels=labels, title="Loss imports disappear when the AC model is disconnected")
    axes[1, 1].tick_params(axis="x", labelsize=8)
    fig.suptitle("Illustrative quasi-static profile, not measured yield or a transient simulation\nDC standby outside the disconnected AC model is excluded", fontsize=12)
    save(fig, "voltage_aware_synthetic_day")


def write_findings(snapshots, chronology):
    env, cases, _, _, _, _ = snapshots
    _, day, _, _, energy, _ = chronology
    lines = ["# Generated voltage-aware study results", "",
             "Generated by `scripts/australian_pv/plot_voltage_aware.py` from the saved CSV/JSON records.", "",
             f"Executed source revision: `{env['git_revision']}`. All {len(cases)} snapshots and {len(day)} chronological intervals are accepted. "
             f"The snapshot record audit verifies {len(env['source_sha256'])} executed source/setup hashes. These are model checks, not hardware certification, regional conformity or EMT validation.", "",
             "## Matched full-PV equilibrium", "",
             "1.03 pu positive-sequence source, 1% aligned upstream VUF, fixed hardware and common reactive-priority terminal settings. "
             "Sequence denotes positive-sequence **sensing**; none of these laws requests negative-sequence cancellation.", "",
             "| Controller | Fleet P at POC (kW) | Fleet Q at POC (kvar) | Maximum sensed phase (V) | Maximum customer VUF (%) | Iterations |",
             "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for mode, label in zip(MODES, LABELS):
        row = next(r for r in cases if r["case_id"] == f"matched_s2_p3_{mode}")
        lines.append(f"| {label} | {number(row, 'p_poc_W')/1000:.4f} | {number(row, 'q_poc_var')/1000:.4f} | "
                     f"{number(row, 'vmax_sensed_V'):.4f} | {number(row, 'max_vuf_percent'):.4f} | {row['iterations']} |")
    angle = [r for r in cases if r["controller"] == "worst_vv_vw" and
             (r["group"] == "source_angle" or r["case_id"] == "matched_s2_p3_worst_vv_vw")]
    lines += ["", f"At the same 1% upstream VUF, changing its relative angle changes worst-phase fleet P from "
              f"{min(number(r, 'p_poc_W') for r in angle)/1000:.4f} to {max(number(r, 'p_poc_W') for r in angle)/1000:.4f} kW, "
              f"and Q from {min(number(r, 'q_poc_var') for r in angle)/1000:.4f} to {max(number(r, 'q_poc_var') for r in angle)/1000:.4f} kvar. "
              "Magnitude-only unbalance studies miss that dependence. The P difference is small relative to fleet output; the figure explicitly zooms that axis.", "",
              "## Synthetic 24-hour profile", "",
              "All profiles have 405 kWh of available fleet PV. Export/import are integrated per inverter before aggregation; Q is signed exchange (absorption negative). "
              "These are inverter POC quantities, excluding feeder losses and customer-load energy.", "",
              "| Controller | Night AC mode | Net POC energy (kWh) | AC import (kWh) | Converter + filter loss (kWh) | Q exchange (kvarh) |",
              "| --- | --- | ---: | ---: | ---: | ---: |"]
    for mode, label in zip(MODES, LABELS):
        for group, operation in (("energized_night", "Energized"), ("disconnected_night", "Disconnected")):
            row = next(r for r in energy if r["controller"] == mode and r["group"] == group)
            lines.append(f"| {label} | {operation} | {number(row, 'net_poc_kWh'):.4f} | {number(row, 'import_kWh'):.4f} | "
                         f"{number(row, 'converter_loss_kWh')+number(row, 'filter_loss_kWh'):.4f} | {number(row, 'reactive_support_kvarh'):.4f} |")
    lines += ["", f"The largest absolute integrated unused-DC term is {max(abs(number(r, 'unused_dc_kWh')) for r in energy):.3g} kWh "
              "(a numerical residual within the physical acceptance tolerance). This profile therefore shows loss costs rather than material DC curtailment. "
              "Energized night control retains Q exchange and imports loss power. Opening the AC contactor removes both modeled services and AC losses; "
              "the disconnected DC standby consumption is undefined and excluded, so this is not a total-product energy-saving claim.", "",
              f"Maximum fixed-PQ replay phasor difference across all {len(cases)+len(day)} cases: {max(number(r, 'pq_replay_difference_V') for r in cases+day):.3g} V. "
              f"Maximum absolute energy-identity residual: {max(abs(number(r, 'balance_residual_Wh')) for r in energy):.3g} Wh. "
              "The replay uses the same network engine with a different injection formulation.", "",
              "## Figures", "",
              "![Numerical convergence and rejected virtual probes](figures/voltage_aware_convergence.png)", "",
              "![Response to upstream negative-sequence angle](figures/voltage_aware_unbalance_angle.png)", "",
              "![Synthetic daily profile and modeled energy accounting](figures/voltage_aware_synthetic_day.png)", "",
              "SVG versions accompany each PNG. See [setup and limitations](voltage_aware.md).", ""]
    (ROOT / "voltage_aware_results.md").write_text("\n".join(lines))


if __name__ == "__main__":
    snapshots, chronology = audit("snapshots"), audit("chronology")
    if not snapshots[-1] or not chronology[-1]:
        raise SystemExit("Campaign incomplete or unresolved cases retained; matched figures withheld. See record_audit.json.")
    figures(snapshots, chronology)
    write_findings(snapshots, chronology)
