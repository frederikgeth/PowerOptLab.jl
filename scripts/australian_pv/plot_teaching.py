"""Pedagogical figures: exact local-law illustration and accepted equilibria."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2] / "studies/australian_pv"
OUT = ROOT / "figures"
PHASE_COLORS = ["#d55e00", "#0072b2", "#009e73"]
MODES = ["mean_vv_vw", "worst_vv_vw", "sequence_vv_vw"]
LABELS = ["Mean VV+VW", "Worst-phase VV+VW", "Sequence VV+VW"]
RESEARCH_PROFILE = {"pq_priority": "proportional", "current_target": "grid", "power_limit_location": "converter"}
TERMINAL_PROFILE = {"pq_priority": "var", "current_target": "converter", "power_limit_location": "poc"}


def csv_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def campaign(name=""):
    path = ROOT / "results" / name
    env = json.loads((path / "environment.json").read_text())
    cases = csv_rows(path / "cases.csv")
    expected = len(env["modes"]) * len(env["availability"]) * len(env["source_positive_sequence_pu"])
    assert len(cases) == expected and all(r["accepted"] == "true" for r in cases), name
    customers, devices = csv_rows(path / "customers.csv"), csv_rows(path / "devices.csv")
    assert len(customers) == len(devices) == 12 * expected
    return env, customers, devices


def select(rows, mode):
    return [r for r in rows if r["mode"] == mode and float(r["availability"]) == 1.
            and float(r["source_pu"]) == 1.03]


def save(fig, name):
    for extension in ("png", "svg"):
        fig.savefig(OUT / f"{name}.{extension}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def sensing_figure():
    fixture = json.loads((ROOT / "results/teaching_fixture.json").read_text())
    curves = fixture["curves"]
    magnitudes = np.array(fixture["phase_voltage_magnitude_V"])
    mean, positive = fixture["mean_magnitude_V"], fixture["positive_sequence_V"]
    q_scale = fixture["q_scale_var"] / 1000
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    ax = axes[0, 0]
    ax.axhspan(220, 240, color="#eeeeee", label="Volt-var deadband")
    ax.vlines([0, 1, 2], positive, magnitudes, color=PHASE_COLORS, lw=5)
    ax.scatter([0, 1, 2], magnitudes, color=PHASE_COLORS, s=95, zorder=3)
    ax.set(xticks=[0, 1, 2], xticklabels=["Phase A", "Phase B", "Phase C"], ylim=(205, 262),
           ylabel="Phase–neutral RMS voltage (V)", title="1. The mean hides the phase extremes")
    ax.axhline(mean, color="#7246a0", ls="--", label=f"Mean = {mean:.1f} V")
    ax.text(.98, .95, f"|U₁| = {positive:.1f} V", transform=ax.transAxes, ha="right", va="top")
    for x, value in enumerate(magnitudes):
        ax.text(x, value + .9, f"{value:.1f}", ha="center")
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[0, 1]
    ax.plot(curves["vv_v"], np.array(curves["vv_q"]) * q_scale, color="#333333")
    for phase, voltage, color in zip("ABC", magnitudes, PHASE_COLORS):
        q = np.interp(voltage, curves["vv_v"], curves["vv_q"]) * q_scale
        ax.scatter(voltage, q, edgecolor=color, facecolor="none" if phase=="B" else color,
                   s=100 if phase=="B" else 45, linewidth=1.8)
        ax.annotate(phase, (voltage, q), xytext=(5, 7 if phase != "C" else -14), textcoords="offset points")
    ax.axhline(0, color="#aaaaaa", lw=.8)
    ax.scatter([mean, positive], [0, 0], marker="x", color="#7246a0", s=65)
    ax.text(231, .3, "Mean and |U₁| lie in the deadband", fontsize=9)
    ax.set(xlabel="Voltage input (V)", ylabel="Q reference (kvar; injection positive)",
           title="2. Worst-phase VV sees opposing requests", xlim=(204, 262))

    ax = axes[1, 0]
    xx = np.linspace(210, 264, 400)
    ax.plot(xx, 100*np.interp(xx, curves["vw_v"], curves["vw_p"]), color="#333333", label="Australia A rated-power ceiling")
    availability = fixture["p_available_W"] / fixture["p_rated_W"] * 100
    ax.axhline(availability, color="#888888", ls=":", label="Available PV / rating")
    ax.axvline(mean, color="#7246a0", ls="--", label="Mean voltage")
    ax.axvline(max(magnitudes), color="#d55e00", ls="--", label="Maximum phase: worst and sequence guard")
    ax.set(xlabel="Voltage input (V)", ylabel="P ceiling (% of rated power)", ylim=(15, 107),
           title="3. The sequence policy still guards the maximum phase")
    ax.legend(loc="lower left", fontsize=8)

    ax = axes[1, 1]
    for mode, label, color in zip(MODES, LABELS, ["#7246a0", "#d55e00", "#0072b2"]):
        result = fixture["policies"][mode]
        p, q = result["p_request_W"]/1000, result["q_request_var"]/1000
        ax.scatter(p, q, color=color, s=85)
        offset = (18, 20) if mode=="mean_vv_vw" else (-10, 20 if mode=="sequence_vv_vw" else -27)
        ax.annotate(f"{label}\n({p:.2f} kW, {q:.2f} kvar)", (p, q),
                    xytext=offset, textcoords="offset points", ha="left" if mode=="mean_vv_vw" else "right", fontsize=9)
    ax.axhline(0, color="#aaaaaa", lw=.8)
    ax.set(xlabel="P request (kW)", ylabel="Q request (kvar)", xlim=(3.0, 6.4), ylim=(-3.8, 1.1),
           title="4. Same inputs and priorities, different requests")
    fig.suptitle("Why sensing policy matters: identical unbalanced voltage phasors", fontsize=15)
    fig.supxlabel("Teaching fixture: U₁ = 230 V, U₂ = 15 V, U₀ = 10 V, all real. Exact local-law requests; no feeder or plant solve.\nRequests precede plant allocation; worst-phase VV uses its declared dominant-conflict rule.", fontsize=8)
    save(fig, "sensing_and_requests")


def headroom_figure():
    env, _, devices = campaign()
    assert all(env["control_settings"].get(k)==v for k,v in RESEARCH_PROFILE.items()), "Headroom figure requires the documented research allocation profile"
    hardware = env["hardware"]
    modes = ["sequence_vv_vw", "sequence_droop", "sequence_droop_power"]
    labels = ["Sequence VV+VW\nNo I₂ service", "Sequence + I₂\nCommon scaling", "Sequence + I₂\nPower priority"]
    power, utilization = [], []
    for mode in modes:
        rows = select(devices, mode)
        total = lambda key: sum(float(r[key]) for r in rows)
        parts = [total("p_poc_W"), total("converter_loss_W"), total("filter_loss_W"), total("dc_availability_margin_W")]
        assert abs(sum(parts)-total("p_available_W")) < .01, "Power accounting must close"
        assert min(parts) >= -.01
        power.append(np.maximum(parts, 0)/1000)
        utilization.append([
            max(float(r["maximum_converter_current_A"])/hardware["i_max"] for r in rows),
            max(1-float(r["apparent_power_margin_VA"])/hardware["s_max"] for r in rows),
            max(float(r["dc_ripple_voltage_V"])/hardware["dv2_max"] for r in rows),
            max(float(r["capacitor_thermal_current_A"])/hardware["i_cap_max"] for r in rows)])
    power, utilization = np.array(power), np.array(utilization)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), layout="constrained", gridspec_kw={"width_ratios": [1, 1.2]})
    left = np.zeros(3)
    for column, label, color in zip(power.T, ["Delivered at PCC", "Converter loss", "Filter loss", "Unused PV power"],
                                    ["#009e73", "#d55e00", "#e6ab02", "#cccccc"]):
        axes[0].barh(np.arange(3), column, left=left, label=label, color=color)
        left += column
    for row, amount in enumerate(power[:, 0]):
        axes[0].text(amount/2, row, f"{amount:.2f} kW", ha="center", va="center", color="white", weight="bold")
    axes[0].set(yticks=np.arange(3), yticklabels=labels, xlabel="Fleet active-power budget (kW)", xlim=(0, 62),
                title="Available minus delivered includes losses")
    axes[0].invert_yaxis()
    axes[0].legend(loc="upper center", bbox_to_anchor=(.5, -.15), ncol=2, fontsize=8)
    image = axes[1].imshow(utilization, vmin=0, vmax=1, cmap="YlGnBu", aspect="auto")
    axes[1].set(xticks=np.arange(4), xticklabels=["Phase RMS\ncurrent", "Converter\nVA", "2ω voltage\nripple", "Capacitor\nthermal RMS"],
                yticks=np.arange(3), yticklabels=labels, title="Maximum use of each declared limit")
    for (row, column), value in np.ndenumerate(utilization):
        precision = 3 if value>.99 else 1
        axes[1].text(column, row, f"{100*value:.{precision}f}%", ha="center", va="center", color="white" if value>.65 else "black")
    fig.colorbar(image, ax=axes[1], fraction=.04, label="Utilization (1 = declared limit)")
    fig.suptitle("Sequence support consumes headroom: unchanged hardware, different allocation", fontsize=14)
    fig.supxlabel("Research profile: converter VA / grid-current command / proportional P–Q. Full PV, source 1.03 pu and 1% VUF.\nColumn maxima can occur at different customers; utilization alone does not establish which constraint determines dispatch.", fontsize=8)
    save(fig, "power_budget_and_headroom")


def terminal_figure():
    env, customers, devices = campaign("terminal_profile_validation")
    assert all(env["control_settings"].get(k)==v for k,v in TERMINAL_PROFILE.items()), "Terminal figure requires the documented reactive-priority/PCC profile"
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.7), sharey=True, layout="constrained")
    all_voltages = []
    for ax, mode, label in zip(axes, MODES, LABELS):
        rows = sorted(select(customers, mode), key=lambda r: r["bus"])
        for phase, color in zip("abc", PHASE_COLORS):
            yy = [float(r[f"v{phase}_pn_V"]) for r in rows]
            all_voltages += yy
            ax.scatter(np.arange(12), yy, color=color, s=30, label=phase.upper())
        dr = select(devices, mode)
        p = sum(float(r["p_poc_W"]) for r in dr)/1000
        q = -sum(float(r["q_poc_var"]) for r in dr)/1000
        ax.axhline(253, color="black", ls="--", lw=1)
        ax.set(title=f"{label}\n{p:.2f} kW delivered · {q:.2f} kvar absorbed",
               xticks=np.arange(12), xticklabels=[r["bus"] for r in rows], xlabel="Customer bus (identical order)")
        ax.tick_params(axis="x", rotation=90, labelsize=8)
        ax.grid(alpha=.15)
    axes[0].set_ylabel("Phase–neutral RMS voltage (V)")
    axes[0].legend(title="Phase", ncol=3, fontsize=8)
    axes[0].set_ylim(np.floor(min(all_voltages))-.5, np.ceil(max(max(all_voltages),253))+.5)
    fig.suptitle("Matched terminal controls: reactive priority, converter current command, PCC VA limit", fontsize=14)
    fig.supxlabel("Accepted equilibria: full PV, source 1.03 pu and 1% VUF; one common solver profile. Dashed line: 253 V.\nScatter points are not joined across feeder branches. These snapshots do not establish regulatory compliance.", fontsize=8)
    save(fig, "terminal_control_comparison")


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    sensing_figure()
    headroom_figure()
    terminal_figure()


if __name__ == "__main__":
    main()
