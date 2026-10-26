"""Generate the balanced, lossless P-Q tutorial figures with matplotlib.

Run: python scripts/formulations/ibr_pq_geometry.py
Requires numpy and matplotlib as development tools, not package dependencies.
Formulas are independent of JuMP; the Julia tutorial tests check solved points
against the same closed-form capability relationships.
"""
from pathlib import Path
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parents[2] / "docs/src/assets/ibr"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelcolor": "#273548", "text.color": "#273548",
    "svg.fonttype": "none", "svg.hashsalt": "PowerOptLab-IBR-PQ"})
COLORS = {"proportional": "#2463a6", "strict_q": "#168477",
    "fixed_pf": "#cb6427", "free_q": "#7d4db1"}
MARKERS = dict(zip(COLORS, ("o", "s", "^", "D")))
LABELS = {"proportional": "Proportional current backoff",
    "strict_q": "Q = 2.640 kvar, maximise P",
    "fixed_pf": "PF = −0.95, maximise P", "free_q": "Free Q, maximise P"}
imax = 1.1 * 6000 / (3 * 230)

def points(radius):
    t = radius / math.hypot(6, 2.64)
    return {"proportional": (6*t, 2.64*t),
        "strict_q": (math.sqrt(radius**2-2.64**2), 2.64),
        "fixed_pf": (.95*radius, math.sqrt(1-.95**2)*radius),
        "free_q": (radius, 0)}

def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", metadata={"Date": None})
    svg = OUT / f"{stem}.svg"
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines())+"\n")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(Path("/tmp") / f"{stem}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.5))
for ax, voltage in zip(axes, (230, 207)):
    radius = min(6, 3*voltage*imax/1000)
    theta = np.linspace(0, np.pi/2, 400)
    p, q = radius*np.cos(theta), radius*np.sin(theta)
    ax.fill_between(p[::-1], 0, q[::-1], color="#e9eef4")
    ax.plot(p, q, color="#526579", lw=2)
    ax.axhline(2.64, color=COLORS["strict_q"], lw=1.3, ls="--")
    ax.plot([0, 6], [0, 2.64], color=COLORS["proportional"], lw=1.3, ls="--")
    ax.plot([0, 6.5], [0, 6.5*math.sqrt(1-.95**2)/.95],
        color=COLORS["fixed_pf"], lw=1.3, ls="--")
    for mode, (x, y) in points(radius).items():
        ax.scatter(x, y, marker=MARKERS[mode], color=COLORS[mode], s=65,
            zorder=5, clip_on=False, edgecolor="white", linewidth=.6)
    ax.scatter(6, 2.64, color="#273548", marker="x", s=65, zorder=6)
    ax.annotate("Raw request\n6 kW + j2.64 kvar", (6, 2.64), xytext=(4.6, 4.2),
        arrowprops={"arrowstyle": "->", "color": "#526579"}, fontsize=10)
    ax.text(.55, 4.5, "Feasible P–Q region\n(reactive injection shown)", fontsize=10)
    ax.set(title=f"{voltage} V per phase · effective rating {radius:.3f} kVA",
        xlabel="Active injection P (kW)", ylabel="Reactive injection Q (kvar)",
        xlim=(0, 6.8), ylim=(0, 6.5), aspect="equal")
    ax.grid(alpha=.18)
legend = [Line2D([0], [0], color=COLORS[k], marker=MARKERS[k], lw=0,
    markersize=8, label=LABELS[k]) for k in COLORS]
fig.legend(handles=legend, loc="lower center", ncol=2, frameon=False)
fig.subplots_adjust(bottom=.22, wspace=.26)
save(fig, "pq-operating-choices")

# A close-up makes the 1% change of effective radius visible.
fig, ax = plt.subplots(figsize=(9, 5.5))
for voltage, ls in ((230, "--"), (207, "-")):
    radius = min(6, 3*voltage*imax/1000)
    q = np.linspace(1.65, 2.95, 300)
    ax.plot(np.sqrt(radius**2-q**2), q, color="#526579", ls=ls,
        label=f"{voltage} V boundary")
    for mode in ("proportional", "strict_q", "fixed_pf"):
        x, y = points(radius)[mode]
        ax.scatter(x, y, color=COLORS[mode], marker=MARKERS[mode], s=65, zorder=5)
        ax.annotate(str(voltage)+" V", (x, y), xytext=(5 if voltage == 230 else -43, 8),
            textcoords="offset points", fontsize=9, color=COLORS[mode])
ax.axhline(2.64, color=COLORS["strict_q"], lw=1, alpha=.6)
ax.set(xlabel="Active injection P (kW)", ylabel="Reactive injection Q (kvar)",
    title="Lower voltage changes the boundary; the operating rule chooses the direction",
    xlim=(5.22, 5.86), ylim=(1.7, 2.88))
ax.grid(alpha=.18)
handles, labels = ax.get_legend_handles_labels()
handles.extend(Line2D([0], [0], color=COLORS[k], marker=MARKERS[k], lw=0,
    markersize=8, label=LABELS[k]) for k in ("proportional", "strict_q", "fixed_pf"))
fig.legend(handles=handles, frameon=False, loc="lower center", ncol=2)
fig.tight_layout(rect=(0, .18, 1, 1))
save(fig, "pq-voltage-closeup")
