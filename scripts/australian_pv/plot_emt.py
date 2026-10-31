"""Publish time-domain validation and engineering failures from audited records."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from audit_emt import audit

REPO = Path(__file__).resolve().parents[2]
STUDY = REPO / "studies/australian_pv"


def generate(input_path, folder, figure_dir):
    receipt = audit(input_path, folder)
    if not receipt["accepted"]:
        raise ValueError("mechanism validation incomplete or unresolved")
    env = json.loads(input_path.read_text())
    result = json.loads((folder / "results.json").read_text())
    cross = json.loads((folder / "stress_crosscheck.json").read_text())
    if (cross["input_sha256"] != receipt["input_sha256"] or cross["results_sha256"] != receipt["results_sha256"]
            or cross["crosscheck_code_sha256"] != hashlib.sha256((REPO / "scripts/australian_pv/crosscheck_emt.jl").read_bytes()).hexdigest()):
        raise ValueError("stale stress cross-check")
    ids = [c["case_id"] for c in env["cases"]]
    labels = ["Balanced L", "DC-link boundary", "AU four-leg / 2% VUF"]
    rows = [result["cases"][i]["modes"] for i in ids]
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "svg.hashsalt": "ibr-emt-circuit"})
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.1), layout="constrained")
    x = np.arange(3)
    ax = axes[0, 0]
    ax.bar(x-.17, [r["averaged"]["finest"]["capacitor_total_rms_A"] for r in rows], .34, label="Averaged circuit", color="#4c8ab5")
    ax.bar(x+.17, [r["switched"]["finest"]["capacitor_total_rms_A"] for r in rows], .34, label="Switched circuit", color="#d67532")
    ax.axhline(1, color="#b2182b", ls="--", label="Fixed 1 A RMS budget")
    ax.set(xticks=x, xticklabels=labels, ylabel="Total capacitor current (A RMS)", title="a  Switching current exceeds the thermal budget", ylim=(0, 5.5))
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax = axes[0, 1]
    expected = [abs(complex(*r["averaged"]["finest"]["expected_dc_2omega_peak"])) for r in rows]
    ax.bar(x-.2, expected, .2, color="#777", label="Phasor reference")
    ax.bar(x, [abs(complex(*r["averaged"]["finest"]["dc_2omega_peak"])) for r in rows], .2, color="#4c8ab5", label="Averaged circuit")
    ax.bar(x+.2, [abs(complex(*r["switched"]["finest"]["dc_2omega_peak"])) for r in rows], .2, color="#d67532", label="Switched circuit")
    ax.set(xticks=x, xticklabels=labels, ylabel="DC-link 2ω voltage (V peak)", title="b  Low-frequency ripple agrees closely", ylim=(0, 10))
    ax.legend(frameon=False, fontsize=9)
    def read_csv(file):
        with (folder / file).open() as stream:
            return list(csv.DictReader(stream))
    carrier = read_csv("binding_dc_capacitor_switched_carrier.csv")
    t = np.array([float(r["time_s"]) for r in carrier])
    voltage = np.array([float(r["v_dc_V"]) for r in carrier])
    take = t >= t[-1]-.02
    d = complex(*rows[1]["switched"]["finest"]["expected_dc_2omega_peak"])
    ax = axes[1, 0]
    ax.plot((t[take]-t[take][0])*1000, 750+np.real(d*np.exp(1j*4*np.pi*50*t[take])), color="#777", ls="--", label="First-order phasor reference")
    ax.plot((t[take]-t[take][0])*1000, voltage[take], color="#d67532", label="Integrated switched voltage")
    ax.set(xlabel="Time within final line cycle (ms)", ylabel="DC-link voltage (V)", title="c  The DC capacitor is a dynamic state")
    ax.legend(frameon=False, fontsize=8)
    trace = read_csv("binding_dc_capacitor_switched_switching.csv")
    tt = np.array([float(r["time_s"]) for r in trace])
    ic = np.array([float(r["capacitor_current_A"]) for r in trace])
    ax = axes[1, 1]
    ax.plot((tt-tt[0])*1e6, ic, color="#d67532", lw=1)
    ax.set(xlabel="Time within final two carriers (µs)", ylabel="Instantaneous capacitor current (A)", title="d  Capacitor pulses carry additional RMS stress")
    ax.text(.03, .06, f"Four-cycle RMS: {rows[1]['switched']['finest']['capacitor_total_rms_A']:.3f} A\nBudget: 1 A RMS (no peak rating declared)", transform=ax.transAxes, fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})
    for ax in axes.flat:
        ax.grid(axis="y", alpha=.15)
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle("Time-domain validation exposes a switching-frequency capacitor limit\nFixed phase/DC sizing · 10 kHz centered PWM · fixed POC voltage · ideal feedforward", fontsize=13)
    figure_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        path = figure_dir / f"emt_capacitor_stress.{extension}"
        fig.savefig(path, dpi=180,
                    metadata={"Date": None} if extension == "svg" else {})
        if extension == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    plt.close(fig)
    lines = ["# Representative switched-circuit validation", "",
             "All three representative cases pass the predeclared mechanism and time-step convergence bands. All three switched cases fail the unchanged 1 A capacitor RMS budget. Agreement with a fundamental-frequency model therefore does not certify switching thermal feasibility.", "",
             "| Case | Maximum fundamental-current error (A RMS) | 2ω voltage phasor error (V peak) | Averaged capacitor (A RMS) | Switched capacitor (A RMS) |",
             "|---|---:|---:|---:|---:|"]
    for ident, r in zip(ids, rows):
        m = r["switched"]["finest"]; a = r["averaged"]["finest"]
        lines.append(f"| {ident} | {m['fundamental_current_max_error_A']:.5f} | {m['dc_2omega_error_V']:.5f} | {a['capacitor_total_rms_A']:.5f} | {m['capacitor_total_rms_A']:.5f} |")
    excess = abs(complex(*rows[1]["averaged"]["finest"]["dc_2omega_peak"])) - env["cases"][1]["dc_2omega_peak_rating_V"]
    lines += ["", f"The independent Python record audit passes {receipt['checks']} coverage, source, waveform, energy, spectral, mesh and engineering-flag checks. The Julia cross-check maps the measured non-2ω RMS stress to the existing `i_sw` bound: all three configurations are rejected because that reserve alone exceeds the 1 A budget. This point-specific measurement is not a conservative reserve over other dispatches.", "",
              f"The synthetic DC-link boundary also exceeds the 7.5 V peak 2ω limit by {excess:.7f} V in the averaged nonlinear capacitor integration with its mean DC voltage matched to 750 V. That first-order boundary approximation is retained as a failed engineering check; its tolerance is not enlarged. Practical boundary certification needs a model-error allowance or a nonlinear energy-based bound.", "",
              "The circuit integrates all phase currents, the physical neutral return and DC-capacitor voltage. Switch edges are exact integration boundaries; RK4 substeps and statistical quadrature are refined at 32/64/128 steps per carrier. The Fourier/RMS window is the final four line cycles after six cycles discarded. Energy accounting includes filter magnetic energy, capacitor energy and the declared mean converter loss.", "",
              "The input is a constant-power source at the regulated DC-link boundary, with no 2ω input power; source current may ripple. This does not validate the optional small-signal R–L source model. Converter loss is a fixed lumped mean from the accepted fundamental operating point. The centered PWM uses sampled sinusoidal feedforward and DC compensation; no PLL, current-feedback, service-controller or feeder dynamics are simulated. Finite-window DC drift is retained, so periodic steady state is not claimed. LCL/split-link dynamics, dead time, semiconductor switching losses, ESR heating, protection and stability remain outside this representative matrix.", "",
              "These results support using a dispatch-dependent PWM stress model or a certified conservative reserve during allocation before making switching-thermal claims. Increasing the capacitor rating is not used to make these points pass. See [method, source assumptions and reproduction](emt.md).", "",
              "![DC ripple and capacitor switching stress](figures/emt_capacitor_stress.png)", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--figures", type=Path, default=STUDY / "figures")
    parser.add_argument("--report", type=Path, default=STUDY / "emt_results.md")
    args = parser.parse_args()
    args.report.write_text(generate(args.input, args.folder, args.figures))
