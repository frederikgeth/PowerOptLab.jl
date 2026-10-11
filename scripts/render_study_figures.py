#!/usr/bin/env python3
"""Render scientific tutorial figures from executed Julia studies.

--check compares fresh solved data with the committed plotting inputs and checks
that each SVG was regenerated after changes to the renderer or input snapshot.
Numerical comparison uses 0.1 W, 1e-4 V/residual and 1e-5 SoC budgets.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "docs/build/plot-cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / "docs/build/plot-cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ASSETS = ROOT / "docs/src/assets/studies"
SNAPSHOT = ASSETS / "study-figure-data.json"
NAMES = ("battery-schedule", "battery-export-limit", "state-estimation")
BLUE, ORANGE, GREEN, GRAY = "#0072B2", "#D55E00", "#009E73", "#555555"
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.labelsize": 11, "axes.titlesize": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.hashsalt": "PowerOptLab-study-figures", "svg.fonttype": "path",
    "lines.linewidth": 2.0, "savefig.facecolor": "white",
})


def rounded(value):
    if isinstance(value, dict):
        return {key: rounded(item) for key, item in value.items()}
    if isinstance(value, list):
        return [rounded(item) for item in value]
    return round(value, 8) if isinstance(value, float) else value


def compare(actual, expected, path="data"):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), f"changed keys at {path}"
        for key in expected:
            compare(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert len(actual) == len(expected), f"changed length at {path}"
        for index, (a, e) in enumerate(zip(actual, expected)):
            compare(a, e, f"{path}[{index}]")
    elif isinstance(expected, (float, int)):
        tolerance = 0.1 if "_w" in path else (1e-4 if "_v" in path or "residual" in path else 1e-5)
        assert np.isfinite(actual) and abs(actual - expected) <= tolerance, (
            f"stale figure data at {path}: {actual} versus {expected}")
    else:
        assert actual == expected, f"changed value at {path}"


def source_digest():
    return hashlib.sha256(Path(__file__).read_bytes() + SNAPSHOT.read_bytes() +
        (ROOT / "scripts/requirements-study-figures.txt").read_bytes()).hexdigest()


def finish(fig, name, output, digest, previews):
    path = output / f"{name}.svg"
    fig.savefig(path, metadata={"Date": None, "Creator": "PowerOptLab / Matplotlib"})
    svg = path.read_text()
    first, rest = svg.split("\n", 1)
    normalized = first + f"\n<!-- plot-source-sha256: {digest} -->\n" + rest
    path.write_text("\n".join(line.rstrip() for line in normalized.splitlines()) + "\n")
    if previews:
        previews.mkdir(parents=True, exist_ok=True)
        fig.savefig(previews / f"{name}.png", dpi=160)
    plt.close(fig)


def battery_schedule(data):
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 5.4), sharex=True,
                             layout="constrained")
    t = data["time_h"]
    ax = axes[0]
    for key, label, color in [("ideal_power_w", "Lossless IVQ", BLUE),
                              ("resistive_power_w", "Resistive IVQ", ORANGE)]:
        power = np.asarray(data[key]) / 1000
        ax.step(t, np.r_[power, power[-1]], where="post", color=color, label=label)
    for key, label, color, marker in [
            ("ideal_pe_power_w", "Unit-efficiency PE", BLUE, "o"),
            ("calibrated_pe_power_w", "Calibrated PE", ORANGE, "s")]:
        ax.plot([0.5, 1.5], np.asarray(data[key])/1000, linestyle="none",
                marker=marker, markerfacecolor="white", color=color, label=label)
    ax.axhline(0, color=GRAY, linewidth=0.7)
    ax.axvspan(0, 1, color=GRAY, alpha=0.06)
    ax.text(0.5, 2.43, "0.30 currency/kWh", ha="center", fontsize=9)
    ax.text(1.5, 2.43, "0.10 currency/kWh", ha="center", fontsize=9)
    ax.set(ylabel="AC injection (kW)", ylim=(-2.4, 2.75),
           title="a   Discharge first, recharge second")
    ax.legend(loc="lower left", fontsize=8.5, ncols=2)
    axes[1].plot(t, data["ideal_soc"], "o-", color=BLUE, label="Lossless")
    axes[1].plot(t, data["resistive_soc"], "s-", color=ORANGE, label="Resistive")
    for bound in (0.05, 0.95):
        axes[1].axhline(bound, color=GRAY, linestyle=":", linewidth=1)
    axes[1].text(1.96, 0.965, "Usable limits", ha="right", fontsize=9, color=GRAY)
    axes[1].set(ylabel="SoC (fraction)", xlabel="Time (h)", ylim=(0, 1.03),
                title="b   Same initial and terminal charge", xlim=(0, 2),
                xticks=[0, 0.5, 1, 1.5, 2])
    axes[1].legend(loc="upper left", fontsize=9)
    for axis in axes:
        axis.grid(axis="y", alpha=0.15)
    return fig


def battery_limit(data):
    fig, ax = plt.subplots(figsize=(7.4, 3.8), layout="constrained")
    soc = np.asarray(data["soc_grid"])
    cell = np.asarray(data["cell_export_bound_w"])/1000
    reference = np.asarray(data["export_reference_w"])/1000
    rating = data["inverter_rating_w"]/1000
    crossover = data["voltage_inverter_crossover_soc"]
    ax.axvspan(soc[0], crossover, color=GREEN, alpha=0.08)
    ax.axvspan(crossover, soc[-1], color=ORANGE, alpha=0.06)
    ax.plot(soc, cell, "--", color=GRAY, linewidth=1.3, label="Cell-limited reference")
    ax.axhline(rating, color=ORANGE, linestyle=":", label="Inverter rating at Q = 0")
    ax.plot(soc, reference, color=BLUE, label="Analytic AC export limit")
    ax.plot(data["solved_soc"], np.asarray(data["linear_export_w"])/1000,
            "o", markersize=8, markerfacecolor="white", color=BLUE, label="Solved: linear OCV")
    ax.plot(data["solved_soc"], np.asarray(data["tabulated_export_w"])/1000,
            "x", markersize=7, markeredgewidth=1.5, color=GREEN, label="Solved: tabulated OCV")
    ax.text(0.28, 0.10, "Minimum cell voltage binds", ha="center", fontsize=9)
    ax.text(0.755, 0.10, "Inverter rating binds", ha="center", fontsize=9)
    ax.set(xlabel="Fixed SoC (fraction)", ylabel="Maximum AC export (kW)",
           xlim=(0.05, 0.95), ylim=(0, 3.35))
    ax.legend(loc="upper left", fontsize=8.5, ncols=2)
    ax.grid(alpha=0.15)
    return fig


def state_estimation(data):
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.5),
                             gridspec_kw={"height_ratios": [1.1, 1]}, layout="constrained")
    true, fit, weak = (np.asarray(data[key]) for key in
        ("true_neutral_v", "fitted_neutral_v", "voltage_only_neutral_v"))
    sigma = data["neutral_component_sigma_v"]
    ax = axes[0]
    ax.annotate("", xy=true, xytext=(0, 0),
                arrowprops={"arrowstyle": "->", "color": GRAY, "lw": 1.5})
    ax.plot(*true, "*", color="black", markersize=11, label="Known reference")
    ax.errorbar(*fit, xerr=sigma[0], yerr=sigma[1], fmt="o", capsize=4,
                color=BLUE, markersize=5, label="Full data, perturbed: component ±1σ")
    ax.plot(*weak, "s", markerfacecolor="white", color=ORANGE, markersize=7,
            label="Voltage only, clean: one admissible point")
    ax.axhline(0, color=GRAY, linewidth=0.6)
    ax.axvline(0, color=GRAY, linewidth=0.6)
    points = np.vstack(([0, 0], true, fit, weak))
    ax.set(xlim=(points[:, 0].min()-0.25, points[:, 0].max()+0.45),
           ylim=(points[:, 1].min()-0.25, points[:, 1].max()+0.40),
           xlabel="Real neutral-to-earth voltage (V)",
           ylabel="Imaginary component (V)", title="a   Neutral displacement at the load bus")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.15)
    ax = axes[1]
    positions = np.arange(len(data["meter_labels"]))
    ax.bar(positions, data["weighted_residuals"], color=BLUE, width=0.65,
           label="Full data, perturbed")
    ax.plot(np.asarray(data["voltage_only_positions"])-1,
            data["voltage_only_residuals"], "s", color=ORANGE,
            markerfacecolor="white", label="Voltage only, clean")
    for level in (-1, 1):
        ax.axhline(level, color=GRAY, linestyle=":", linewidth=1)
    ax.axhline(0, color=GRAY, linewidth=0.7)
    ax.set(xticks=positions, xticklabels=data["meter_labels"],
           ylabel="Residual / stated σ", ylim=(-1.2, 1.2),
           title="b   A small residual does not establish observability")
    ax.tick_params(axis="x", labelrotation=45, labelsize=9)
    ax.legend(loc="lower left", fontsize=8.5, ncols=2)
    ax.grid(axis="y", alpha=0.15)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT/"docs/build/study-figure-data")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--preview-dir", type=Path)
    args = parser.parse_args()
    data = {name: rounded(json.loads((args.data_dir/f"{name}.json").read_text()))
            for name in ("battery_storage_models", "constrained_nlls_state_estimation")}
    ASSETS.mkdir(parents=True, exist_ok=True)
    if args.check:
        compare(data, json.loads(SNAPSHOT.read_text()))
    else:
        SNAPSHOT.write_text(json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False) + "\n")
    digest = source_digest()
    if args.check:
        for name in NAMES:
            assert f"plot-source-sha256: {digest}" in (ASSETS/f"{name}.svg").read_text(), (
                f"stale figure renderer or input snapshot: {name}; regenerate the figures")
        output = ROOT / "docs/build/study-figure-check"
    else:
        output = ASSETS
    output.mkdir(parents=True, exist_ok=True)
    figures = (battery_schedule(data["battery_storage_models"]),
               battery_limit(data["battery_storage_models"]),
               state_estimation(data["constrained_nlls_state_estimation"]))
    for fig, name in zip(figures, NAMES):
        finish(fig, name, output, digest, args.preview_dir)
    print("Checked" if args.check else "Generated", "3 reproducible study figures")


if __name__ == "__main__":
    main()
