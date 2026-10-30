"""Plot the measured offline allocation fixture and its balanced-current slice."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "scripts/formulations/inverter_current_projection_reference.json"
OUT = ROOT / "docs/src/assets/ibr/generated-current-projection"


def phasors(rows):
    return np.array([complex(r["real"], r["imag"]) for r in rows])


def main():
    data = json.loads(DATA.read_text())
    h = data["hardware"]
    case = next(c for c in data["cases"] if c["label"] == "rail_conflict")
    assert case["accepted"] and all(c["passed"] for c in case["physical_checks"])
    # These closed-form curves are specific to this balanced, reduced-filter
    # fixture. They are not an alternate certificate for arbitrary inverters.
    assert h["topology"] == "THREE_LEG" and h["pwm_strategy"] == "NONE"
    assert all(h[k] == 0 for k in ("r_filter_grid", "x_filter_grid", "c_filter_mid",
                                  "b_filter_shunt", "a_loss", "esr_dc"))
    assert h["dc_source_r"] is None and h["r_filter_matrix"] is None
    assert h["power_limit_location"] == "converter"
    assert all(x == 0 for x in h["pwm_ac_converter_reserve"])
    v = case["voltage"]
    z = complex(h["r_filter"], h["x_filter"])
    p, q = np.meshgrid(np.linspace(-100, 5500, 700), np.linspace(-6600, 500, 800))
    current = (p - 1j*q)/(3*v)
    internal = v + z*current
    converter_power = 3*internal*np.conj(current)
    dc_power = converter_power.real + h["p_loss_fixed"] + 3*h["c_loss"]*abs(current)**2
    feasible = ((abs(current) <= h["i_max"]) & (abs(converter_power) <= h["s_max"])
                & (np.sqrt(6)*abs(internal) <= h["m_max"]*h["v_dc"])
                & (dc_power >= data["dc_source"]["p_min_W"])
                & (dc_power <= case["available"]))

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, (ax, wave) = plt.subplots(1, 2, figsize=(13, 5.5))
    ax.contourf(p/1000, q/1000, feasible.astype(float), levels=[.5, 1.5],
                colors=["#cce8dc"])
    ax.contour(p/1000, q/1000, np.sqrt(6)*abs(internal),
               levels=[h["m_max"]*h["v_dc"]], colors=["#0072b2"], linewidths=1.4)
    ax.plot([0, 4.8], [0, 0], color="#bd4b30", linestyle=":", linewidth=2)
    ax.scatter([4.8, 0], [0, 0], marker="x", s=65, color="#bd4b30", zorder=4)
    ax.annotate("Requested 4.8 kW, Q = 0", (4.8, 0), xytext=(-155, -18),
                textcoords="offset points", color="#9b3822")
    ax.annotate("Zero current also violates rails", (0, 0), xytext=(15, -35),
                textcoords="offset points", color="#9b3822")
    point = (case["p_poc_W"]/1000, case["q_poc_var"]/1000)
    ax.scatter(*point, s=65, color="#00734b", zorder=5)
    ax.annotate(f"Accepted local allocation\n{point[0]:.3f} kW, {point[1]:.3f} kvar",
                point, xytext=(25, 45), textcoords="offset points",
                arrowprops={"arrowstyle": "->", "color": "#00734b"})
    ax.text(.03, .50, "Green: feasible balanced-current slice\nBlue: continuous voltage-rail boundary",
            transform=ax.transAxes, fontsize=10)
    ax.set(xlabel="POC active injection (kW)", ylabel="POC reactive injection (kvar)",
           title="Balanced-current capability slice",
           xlim=(-.1, 5.5), ylim=(-6.6, .5))
    ax.grid(alpha=.15)

    theta = np.linspace(0, 2*np.pi, 2401)
    voltage = phasors(case["grid_voltage_V"])
    for label, key, color in (("Requested", "requested_grid_current_A", "#bd4b30"),
                              ("Allocated", "allocated_grid_current_A", "#00734b")):
        i = phasors(case[key]); u = voltage + z*i
        line = u - np.roll(u, -1)
        envelope = np.max(abs(np.sqrt(2)*np.real(line[:, None]*np.exp(1j*theta))), axis=0)
        wave.plot(np.degrees(theta), envelope, color=color, label=label)
    wave.axhline(h["m_max"]*h["v_dc"], color="black", linestyle="--", label="Available rail")
    wave.set(xlabel="Electrical angle (degrees)", ylabel="Maximum |line voltage| (V, instantaneous)",
             title="Full-cycle voltage demand", xlim=(0, 360), ylim=(490, 640))
    wave.set_xticks(np.arange(0, 361, 60)); wave.legend(loc="lower center", ncol=3, fontsize=10)
    wave.grid(alpha=.15)
    fig.suptitle("Voltage protection needs the physical filter and DC rails", fontsize=15)
    fig.text(.5, .015, "Fixed synthetic fixture: 253 V POC, 620 V DC link, 6 kVA plant; PWM disabled. "
             "Australian study sizing remains 750 V.\n"
             "The slice assumes balanced currents; the accepted allocation is a local NLP result, "
             "not a global nearest-point certificate.", ha="center", fontsize=10)
    fig.tight_layout(rect=(0, .10, 1, .99))
    for ext in ("svg", "png"):
        fig.savefig(OUT.with_suffix("."+ext), dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
