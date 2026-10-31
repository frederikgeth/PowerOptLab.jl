"""Reaudit and publish the fixed-hardware LF/carrier voltage boundary experiment."""
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from audit_emt import audit

ROOT = Path(__file__).resolve().parents[2] / "studies/australian_pv"
IDS = ("lf_only_boundary", "combined", "four_leg_combined")
LABELS = ("Three-leg: LF voltage only", "Three-leg: combined bound", "Four-leg: combined bound")


def publish(root=ROOT):
    inputs = root / "results/engineering_boundary_emt_inputs.json"
    folder = root / "results/engineering_boundary_emt"
    receipt = audit(inputs, folder)
    if not receipt["accepted"]:
        raise ValueError("Independent EMT record audit must pass before publication")
    env = json.loads(inputs.read_text())
    refs = {r["case_id"]: r for r in env["cases"]}
    result = json.loads((folder / "results.json").read_text())
    if set(refs) != set(IDS) or env["sizing_changed"]:
        raise ValueError("Unexpected boundary experiment coverage or sizing")
    for r in refs.values():
        if (r["dc_voltage_min_rating_V"], r["dc_voltage_max_rating_V"],
                r["capacitor_current_rating_A"], r["phase_peak_rating_A"]) != (747., 753., 1., 14.):
            raise ValueError("Boundary hardware/operating ceilings changed")

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11,
                         "svg.hashsalt": "engineering-boundary"})
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.7), layout="constrained",
                             gridspec_kw={"width_ratios": [1.6, 1, 1]})
    colors = {"averaged": "#7b8084", "switched": "#167d9a"}
    offsets = {"averaged": -.12, "switched": .12}
    metrics = []
    for j, ident in enumerate(IDS):
        for mode in ("averaged", "switched"):
            f = result["cases"][ident]["modes"][mode]["finest"]
            y = j + offsets[mode]
            color = "#ba3e35" if not f["all_engineering_checks_passed"] else colors[mode]
            label = mode.capitalize() if j == 1 else None
            axes[0].plot([f["v_dc_min_V"], f["v_dc_max_V"]], [y, y], "|-", color=color,
                         lw=2, markersize=11, label=label)
            axes[0].plot(f["v_dc_mean_V"], y, ".", color=color)
            axes[1].barh(y, f["capacitor_total_rms_A"], height=.19, color=color)
            axes[2].plot(max(f["current_peak_A"]), y, "o", color=color, ms=5)
            if ident == "four_leg_combined":
                axes[2].plot(f["neutral_peak_A"], y, "x", color=color, ms=7)
            metrics.append({"case_id": ident, "mode": mode,
                "dc_min_V": f["v_dc_min_V"], "dc_max_V": f["v_dc_max_V"],
                "dc_mean_V": f["v_dc_mean_V"], "capacitor_rms_A": f["capacitor_total_rms_A"],
                "largest_phase_peak_A": max(f["current_peak_A"]), "neutral_peak_A": f["neutral_peak_A"],
                "engineering_accepted": f["all_engineering_checks_passed"]})
    for ax in axes:
        ax.set_ylim(2.6, -.6)
        ax.set_yticks(range(3), LABELS if ax is axes[0] else [""]*3)
        ax.grid(axis="x", alpha=.2)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].axvspan(747, 753, color="#e3efe5", alpha=.5, zorder=0)
    for v in (747, 753):
        axes[0].axvline(v, color="#b2182b", ls="--", lw=1)
    axes[0].set(xlim=(746.7, 753.8), xlabel="DC voltage range and mean (V)",
                title="a  DC voltage extrema\nSwitching exposes the LF-only failure")
    axes[0].legend(loc="lower left", frameon=False, fontsize=9)
    axes[0].annotate("753.170 V: rejected", xy=(753.17012, .12), xytext=(750.8, .66),
                     fontsize=9, color="#a52e27", arrowprops={"arrowstyle": "->", "color": "#a52e27"})
    axes[1].axvline(1., color="#b2182b", ls="--", lw=1)
    axes[1].set(xlim=(0, 1.08), xlabel="Capacitor current (A RMS)",
                title="b  Capacitor RMS current\nFixed 1 A budget")
    axes[2].axvline(14., color="#b2182b", ls="--", lw=1)
    axes[2].set(xlim=(0, 15), xlabel="Current peak (A)",
                title="c  Phase and neutral peaks\n14 A declared operating ceiling")
    axes[2].text(.03, .03, "● largest phase    × neutral", transform=axes[2].transAxes, fontsize=8)
    fig.suptitle("Combined DC bounds leave room for carrier ripple\n"
                 "Unchanged 6 kVA / 10 A / 750 V / 300 µF hardware · 10 kHz PWM", fontsize=13)
    target = root / "figures"
    target.mkdir(exist_ok=True)
    for suffix in ("png", "svg"):
        path = target / f"engineering_boundary.{suffix}"
        fig.savefig(path, dpi=180,
                    metadata={"Date": None} if suffix == "svg" else {})
        if suffix == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    plt.close(fig)
    with (folder / "comparison.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(metrics)

    lines = ["# Combined switching and nonlinear DC-voltage boundary", "",
        "This boundary experiment compares the existing LF-only voltage guard with the combined LF/carrier guard, "
        "and exercises the four-leg implementation. It is separate from the matched feeder PF/OPF benchmark. "
        "All three retain 6 kVA, 10 A phase RMS, 750 V, 300 µF and the 1 A capacitor RMS budget. "
        "The shared 747–753 V validation band and 14 A peak ceilings are declared operating assumptions, not inferred manufacturer ratings.", "",
        "![DC extrema, capacitor RMS and phase/neutral peaks](figures/engineering_boundary.png)", "",
        "The grey rows show averaged integration; blue rows show switched integration. Red identifies the reference that fails the upper DC-voltage check. "
        "The voltage segments span measured minimum to maximum, with a dot at the measured mean. "
        "All three meshes (32, 64 and 128 steps/carrier) remain in the saved record; the figure uses the finest mesh.", "",
        "| Reference | P at POC (W) | Q at POC (var) | I₂ (A RMS) | Switched DC range (V) | Capacitor RMS (A) | Outcome |",
        "|---|---:|---:|---:|---:|---:|---|"]
    for ident, label in zip(IDS, LABELS):
        ref = refs[ident]
        v = [complex(*z) for z in ref["grid_voltage"]]
        i = [complex(*z) for z in ref["current"]]
        s = sum(x*y.conjugate() for x, y in zip(v, i))
        a = complex(-.5, math.sqrt(3)/2)
        i2 = abs((i[0]+a*a*i[1]+a*i[2])/3)
        f = result["cases"][ident]["modes"]["switched"]["finest"]
        outcome = "Pass" if f["all_engineering_checks_passed"] else "Fail: DC upper voltage"
        lines.append(f"| {label} | {s.real:.3f} | {s.imag:.3f} | {i2:.6f} | "
                     f"{f['v_dc_min_V']:.6f}–{f['v_dc_max_V']:.6f} | {f['capacitor_total_rms_A']:.6f} | {outcome} |")
    lines += ["", "The three-leg references request the same 600 W active power, −300 var reactive power and 1 A negative-sequence current. "
        "Q and P have priority; the combined guard reduces the remaining sequence service to leave voltage headroom. "
        "The four-leg reference has a different declared request (600 W total and per-phase Q of [−600, 0, 0] var); "
        "its row demonstrates that topology and is not a matched controller performance comparison.", "",
        f"The independent record audit passes **{receipt['checks']} checks** covering provenance, mechanism agreement, energy, "
        "mesh convergence, extrema witnesses and engineering flags. The LF-only switched reference fails only the upper DC-voltage check; "
        "its failure remains in the results and is not accepted as hardware-feasible.", "",
        "These are fixed-POC, feedforward L-filter circuits with a constant-power DC input and fixed mean converter loss. "
        "Finite-window DC mean drift is retained: the switched extrema include both carrier ripple and that drift. "
        "This experiment does not isolate a carrier-only error or establish periodic steady state, closed-loop regulation, "
        "feeder dynamics, LCL validity or all-time transient bounds. Analytical carrier bounds and the calibrated 10% capacitor thermal allowance remain distinct.", "",
        "Reproduce from the repository root:", "", "```sh",
        "julia --project=. scripts/australian_pv/prepare_engineering_boundary_emt.jl /tmp/boundary_inputs.json",
        "julia --project=. scripts/australian_pv/emt_circuit.jl /tmp/boundary_inputs.json /tmp/boundary_emt",
        "python3 scripts/australian_pv/audit_emt.py /tmp/boundary_inputs.json /tmp/boundary_emt",
        "python3 scripts/australian_pv/plot_engineering_boundary.py", "```", "",
        "The plotting command reaudits the committed records before regenerating the figure, table and report. "
        "New simulations retain their own input/source hashes and are not relabeled as the committed campaign.", ""]
    (root / "engineering_boundary.md").write_text("\n".join(lines))
    return receipt


if __name__ == "__main__":
    print(f"Boundary publication: {publish()['checks']} independent checks")
