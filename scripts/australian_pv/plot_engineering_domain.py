"""Reaudit the declared domain and explain conservatism of its carrier bounds."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from audit_engineering_domain import audit

ROOT = Path(__file__).resolve().parents[2]/"studies/australian_pv"


def publish():
    inputs = ROOT/"results/engineering_domain_inputs.json"
    target = ROOT/"results/engineering_domain_audit.json"
    report = audit(inputs, target)
    env = json.loads(inputs.read_text())
    refs = {r["id"]: r for r in env["cases"]}
    rows = []
    for topology in ("THREE_LEG", "FOUR_LEG"):
        for level in env["domain"]["modulation_fraction"]:
            selected = [r for r in report["records"] if refs[r["id"]]["topology"] == topology and
                        refs[r["id"]]["modulation_fraction"] == level]
            rows.append(dict(topology=topology, modulation_fraction=level, cases=len(selected),
                engineering_limit_failures=sum(not r["engineering_limits_passed"] for r in selected),
                phase_rms_ratio=max(r["maximum_phase_rms_A"]/min(refs[r["id"]]["phase_rms_bound_A"]) for r in selected),
                phase_peak_ratio=max(r["maximum_phase_peak_A"]/min(refs[r["id"]]["phase_peak_bound_A"]) for r in selected),
                dc_voltage_ratio=max(r["coupled_dc_peak_V"]/r["dc_bound_V"] for r in selected)))
    with (ROOT/"results/engineering_domain_summary.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    plt.rcParams.update({"font.size": 10, "svg.hashsalt": "engineering-domain"})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), layout="constrained")
    for topology, label, color in (("THREE_LEG", "Three-leg", "#167d9a"), ("FOUR_LEG", "Four-leg", "#b86a2d")):
        group = [r for r in rows if r["topology"] == topology]
        x = [r["modulation_fraction"] for r in group]
        axes[0].plot(x, [r["phase_rms_ratio"] for r in group], "o-", color=color, label=f"{label} RMS")
        axes[0].plot(x, [r["phase_peak_ratio"] for r in group], "s--", color=color, label=f"{label} peak")
        axes[1].plot(x, [r["dc_voltage_ratio"] for r in group], "o-", color=color, label=label)
    for ax in axes:
        ax.axhline(1., color="#b2182b", ls=":", label="Analytical envelope")
        ax.set(xlabel="Fraction of declared modulation ceiling", ylim=(0, 1.08), xticks=[.3, .7, .98])
        ax.grid(alpha=.2)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8)
    axes[0].set(title="a  Phase ripple: largest ratio in each group", ylabel="Sampled current / analytical bound")
    axes[1].set(title="b  DC charge: includes AC ripple current", ylabel="Sampled voltage peak / analytical bound")
    fig.suptitle("All-duty carrier envelopes leave conservative margin\n972 predeclared points · fixed hardware · 2 integration meshes", fontsize=13)
    for suffix in ("png", "svg"):
        path = ROOT/f"figures/engineering_domain.{suffix}"
        fig.savefig(path, dpi=180, metadata={"Date": None} if suffix == "svg" else {})
        if suffix == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    plt.close(fig)
    def maximum(key):
        return max(r[key] for r in report["records"])
    lines = ["# Declared switching-envelope validation domain", "",
        f"The independent NumPy oracle passes **{report['checks']:,} checks over {report['cases']} operating points**. "
        "The domain is declared in the preparation script and independently checked by the oracle; publication cannot shrink it to accepted points.", "",
        "| Dimension | Declared values |", "|---|---|",
        "| Topology | Three-leg; four-leg with explicit neutral |",
        "| Total converter P | 0, 600, 5000 W |", "| Total converter Q | −2000, 0, +2000 var |",
        "| Negative-sequence current | 0, 1 A RMS |", "| Negative-sequence voltage | 2% of positive sequence |",
        "| Sequence angle | 0°, 120°, 240° |", "| Fraction of 0.96 modulation ceiling | 0.30, 0.70, 0.98 |",
        "| Carrier/line angle | 0°, 60°, 120° |", "| Four-leg zero-sequence current | 0.1 A RMS |",
        "| Independent meshes | 8192, 16384 samples per carrier |", "",
        "![Conservatism of phase-current and coupled DC-charge bounds](figures/engineering_domain.png)", "",
        "The figure shows the largest ratio over P/Q, sequence and carrier-angle variations in each topology/modulation group. "
        "A ratio below one means the sampled waveform stays within its analytical carrier envelope; it does not measure total inverter utilization. "
        "The frozen bus is 750 V; AC reserves use the declared 765 V upper bound. Bounds retain their original all-duty derivation and are not fitted to this sweep.", "",
        f"The largest exact-event versus sampled phase-RMS difference is {maximum('event_rms_error_A'):.6f} A; "
        f"the largest frozen-current DC-extremum difference is {maximum('event_dc_error_V'):.6f} V. "
        f"Across meshes, the largest RMS difference is {maximum('mesh_rms_error_A'):.6f} A and the largest coupled-charge peak difference is "
        f"{maximum('mesh_dc_error_V'):.6f} V. The independent passive R/L orbit checks also enforce periodic closure and phase/neutral RMS and peak envelopes.", "",
        f"**{report['engineering_limit_failures']} points violate at least one physical engineering limit and are retained with their failure names.** "
        "These include intentionally stressed power/current combinations. Carrier-domain agreement does not turn them into acceptable dispatches. "
        "The engineering-feasibility labels are Julia SI-model diagnostics; the independent oracle validates carrier mechanisms, bounds, coverage and retention, not all those feasibility labels.", "",
        "The hardware remains 6 kVA, 10 A, 750 V, 300 µF, a 1 A capacitor budget and the same filter. "
        "This is a frozen-bus/frozen-reference, inductive carrier approximation with periodic mean-subtracted charge. "
        "It does not establish coupled nonlinear transient bounds, cumulative DC-energy drift, LCL validity or closed-loop stability. "
        "The calibrated 10% capacitor thermal allowance is separate and is not validated or retuned by this carrier sweep. "
        "Representative averaged/switched EMT is documented in the [boundary experiment](engineering_boundary.md); "
        "the [matched mixed-feeder PF/OPF benchmark](engineering_benchmark.md) separately validates the coupled thermal allocator.", "",
        "Reproduce from the repository root:", "", "```sh",
        "julia --project=. scripts/australian_pv/prepare_engineering_domain.jl /tmp/domain_inputs.json",
        "python3 scripts/australian_pv/audit_engineering_domain.py /tmp/domain_inputs.json /tmp/domain_audit.json",
        "python3 scripts/australian_pv/test_engineering_domain.py",
        "python3 scripts/australian_pv/plot_engineering_domain.py", "```", "",
        "Restore the full committed inputs using [the evidence index](evidence/README.md) before these saved-record commands; "
        "the plot reruns the complete independent audit. "
        "Eight tests cover full-domain acceptance and corrupted coverage, bounds, currents, gates, sources and failure flags.", ""]
    (ROOT/"engineering_domain.md").write_text("\n".join(lines))
    print(f"Published {report['cases']} points and {report['checks']} independent checks")


if __name__ == "__main__":
    publish()
