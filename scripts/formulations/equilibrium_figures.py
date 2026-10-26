"""Rebuild the workbook's analytical SVG/PDF figures (numpy + matplotlib).
Run: python scripts/formulations/equilibrium_figures.py
No Julia results or solver success flags are used as plotting oracles.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[2] / "docs/src/assets/learning"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "svg.hashsalt": "PowerOptLab-equilibrium-workbook"})
blue, orange, green, purple, ink = "#2463a6", "#cb6427", "#168477", "#7d4db1", "#273548"

def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    svg = OUT / f"{name}.svg"
    fig.savefig(svg, bbox_inches="tight", metadata={"Date": None})
    svg.write_text("\n".join(x.rstrip() for x in svg.read_text().splitlines())+"\n")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(Path("/tmp") / f"{name}.png", bbox_inches="tight", dpi=160)
    plt.close(fig)

def command(p):
    return np.clip(1-1.2*p, 0, 1)

def iterate(alpha):
    p = [0.0]
    for _ in range(20):
        p.append((1-alpha)*p[-1]+alpha*command(p[-1]))
    return p

root = 5/11
fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.3), layout="constrained")
p = np.linspace(0, 1, 301)
axes[0].plot(p, command(p), color=blue, lw=2.5, label="Controller output C(V(p))")
axes[0].plot(p, p, color=orange, lw=2, ls="--", label="Output equals trial power")
axes[0].scatter([root], [root], s=60, color=green, zorder=5)
axes[0].annotate("Unique equilibrium\np = 5/11", (root, root), xytext=(.51, .74),
    arrowprops={"arrowstyle": "->", "color": ink})
axes[0].set(title="An intersection defines the equilibrium", xlabel="Trial power p (pu)",
    ylabel="Power (pu)", xlim=(0, 1), ylim=(0, 1.08))
axes[0].legend(frameon=False, loc="upper center", bbox_to_anchor=(.5,-.2), fontsize=9)
axes[1].plot(range(21), iterate(1), color=orange, marker="o", markersize=3,
    label="Direct numerical iteration")
axes[1].plot(range(21), iterate(.5), color=blue, marker="s", markersize=3,
    label="Relaxed numerical iteration, α = 0.5")
axes[1].axhline(root, color=green, ls="--", label="Same equilibrium")
axes[1].set(title="The numerical method can still cycle", xlabel="Numerical iteration n (not seconds)",
    ylabel="Trial power pₙ (pu)", xlim=(0, 20), ylim=(-.04, 1.08))
axes[1].legend(frameon=False, loc="upper center", bbox_to_anchor=(.5,-.2), fontsize=9)
for ax in axes:
    ax.grid(alpha=.18)
save(fig, "equilibrium-versus-iteration")

fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.6), layout="constrained")
v = np.linspace(.985, 1.105, 400)
c = np.clip(1-12*(v-1), 0, 1)
for ax in axes:
    ax.plot(v, c, color=blue, lw=2.3, label="Controller graph p = C(v)")
    ax.plot(v, (v-1)/.1, color=orange, ls="--", label="Network v = 1 + 0.1p")
    ax.set(xlabel="Voltage v (pu)", ylabel="Delivered power p (pu)",
        xlim=(.985, 1.105), ylim=(-.05, 1.08))
    ax.grid(alpha=.18)
axes[0].scatter([1+.1*root], [root], color=green, s=65, zorder=5)
axes[0].set_title("Controller equality: one intersection")
axes[0].legend(frameon=False, loc="upper center", bbox_to_anchor=(.5,-.2), fontsize=9)
axes[1].fill_between(v, 0, c, color="#e9eef4", label="Ceiling region: 0 ≤ p ≤ C(v)")
network_p = np.linspace(0, root, 100)
axes[1].plot(1+.1*network_p, network_p, color=green, lw=4,
    label="Network-feasible dispatch choices")
for power, label, offset in ((0, "Minimise p", (8, 15)),
                              (.2, "Track p = 0.2", (-4, 25)),
                              (root, "Maximise p", (15, -4))):
    axes[1].scatter([1+.1*power], [power], color=purple, s=50, zorder=6)
    axes[1].annotate(label, (1+.1*power, power), xytext=offset,
        textcoords="offset points", fontsize=10)
axes[1].set_title("Capability inequality: the objective selects")
axes[1].legend(frameon=False, loc="upper center", bbox_to_anchor=(.5,-.2), fontsize=9)
save(fig, "controller-graph-versus-capability")

fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.3), layout="constrained")
axes[0].plot([0, 1, 2, 3], [1, 2, 1, 1], color=blue, marker="o", label="Voltage ceiling 1.030 pu")
axes[0].plot([0, 1, 2, 3], [1, 2, 1.5, 1], color=orange, marker="s", ls="--",
    label="Voltage ceiling 1.015 pu")
axes[0].set(title="Energy links consecutive physical periods", xlabel="Physical time (h)",
    ylabel="Stored energy (kWh)", xticks=[0, 1, 2, 3], ylim=(0, 2.2))
axes[0].legend(frameon=False, loc="lower left", fontsize=9)
axes[0].grid(alpha=.18)
voltage = np.linspace(.97, 1.03, 601)
axes[1].axvspan(.98, 1.02, color="#e9eef4")
axes[1].step(voltage, voltage >= 1.02, where="post", color=blue,
    label="Increasing voltage, initially off")
axes[1].step(voltage, voltage > .98, where="post", color=orange, ls="--",
    label="Decreasing voltage, initially on")
axes[1].scatter([1, 1], [0, 1], color=purple, s=45, zorder=5)
axes[1].annotate("Same voltage,\ndifferent remembered modes", (1, .5),
    xytext=(.984, .42), fontsize=10)
axes[1].set(title="Hysteresis retains physical history", xlabel="Measured voltage (pu)",
    ylabel="Curtailment mode", yticks=[0, 1], yticklabels=["Off (0)", "On (1)"],
    xlim=(.968, 1.032), ylim=(-.2, 1.2))
axes[1].legend(frameon=False, loc="upper center", bbox_to_anchor=(.5, -.19), fontsize=9)
save(fig, "physical-time-and-hysteresis")
