"""Independent sampled-carrier validation of the declared engineering domain.

The reference JSON comes from Julia; this oracle constructs its own inductance
map, triangular carrier and time integrals. It does not import PowerOptLab.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DECLARED_DOMAIN = dict(p_W=[0., 600., 5000.], q_var=[-2000., 0., 2000.],
    negative_sequence_A=[0., 1.], sequence_angle_rad=[0., 2*math.pi/3, 4*math.pi/3],
    modulation_fraction=[.3, .7, .98], carrier_line_angle_rad=[0., math.pi/3, 2*math.pi/3],
    voltage_negative_sequence_fraction=.02, four_leg_zero_sequence_A=.1,
    sampled_steps=[8192, 16384], event_rms_tolerance_A=.005, event_dc_tolerance_V=.001)


def independent_bounds(ref, hw, ceilings):
    four = ref["topology"] == "FOUR_LEG"
    L = hw["x_pu"] * (hw["v_nom"]**2/(hw["s_max"]/3))/(2*math.pi*50)
    B = np.column_stack((np.eye(3), -np.ones(3))) if four else np.eye(3)-np.ones((3, 3))/3
    M = np.linalg.solve(L*np.eye(3)+L*np.ones((3, 3)), B) if four else B/L
    phase_width = np.abs(M).sum(axis=1)
    neutral_width = np.abs(-M.sum(axis=0)).sum() if four else 0.
    peak = phase_width*ceilings["dc_max_V"]/(4*10000.)
    neutral_peak = neutral_width*ceilings["dc_max_V"]/(4*10000.)
    i = np.array([complex(*z) for z in ref["converter_current"]])
    leg_magnitude = np.abs(i).sum()+(abs(i.sum()) if four else 0.)
    frozen = math.sqrt(2)*leg_magnitude/(8*10000.*hw["c_dc"])
    variation = (peak.sum()+neutral_peak)/(4*10000.*hw["c_dc"])
    return dict(phase_rms_bound_A=peak/2, phase_peak_bound_A=peak,
        neutral_rms_bound_A=neutral_peak/2, neutral_peak_bound_A=neutral_peak,
        dc_frozen_bound_V=frozen, dc_ac_variation_bound_V=variation, dc_voltage_bound_V=frozen+variation)


def sampled_carrier(ref, hardware, n):
    four = ref["topology"] == "FOUR_LEG"
    period = 1e-4
    dt = period/n
    L = hardware["x_pu"] * (hardware["v_nom"]**2/(hardware["s_max"]/3))/(2*math.pi*50)
    # The fourth conductor has the same declared R/X as a phase. Inverting
    # L*I+L_n*11' here is independent of Julia's reduced-basis gate map.
    B = np.column_stack((np.eye(3), -np.ones(3))) if four else np.eye(3)-np.ones((3, 3))/3
    inverse = np.linalg.inv(L*np.eye(3)+L*np.ones((3, 3))) if four else np.eye(3)/L
    tau = (np.arange(n)+.5)/n
    carrier = np.where(tau < .5, -1+4*tau, 3-4*tau)
    duty = np.array(ref["duties"])
    gates = (2*duty[:, None]-1 >= carrier).astype(float)
    forcing = hardware["v_dc"] * inverse @ B @ (gates-duty[:, None])
    # Midpoint state for the lossless inductive primitive, then its zero-mean
    # carrier component. Sampling closure error is retained as a metric.
    wave = dt*(np.cumsum(forcing, axis=1)-forcing/2)
    wave -= wave.mean(axis=1, keepdims=True)
    neutral = -wave.sum(axis=0) if four else np.zeros(n)
    phase_rms = np.sqrt(np.mean(wave**2, axis=1))
    phase_peak = np.max(np.abs(wave), axis=1)
    frozen = np.array(ref["frozen_leg_current_A"])
    frozen_dc = (gates*frozen[:, None]).sum(axis=0)
    frozen_dc -= frozen_dc.mean()
    frozen_voltage = -dt/hardware["c_dc"]*(np.cumsum(frozen_dc)-frozen_dc/2)
    frozen_voltage -= frozen_voltage.mean()
    legs = np.vstack((wave, neutral)) if four else wave
    bridge = (gates*(frozen[:, None]+legs)).sum(axis=0)
    bridge -= bridge.mean()
    charge = -dt/hardware["c_dc"]*(np.cumsum(bridge)-bridge/2)
    charge -= charge.mean()
    # Exact periodic solution for piecewise-constant sampled forcing with the
    # declared proportional passive R/L. Remove the sampling-only mean error;
    # the analytical PWM forcing has exactly zero carrier mean.
    rate = hardware["x_pu"] * (hardware["v_nom"]**2/(hardware["s_max"]/3))/hardware["x_over_r"]/L
    decay = np.exp(-rate*dt*np.arange(1, n+1))
    f = forcing-forcing.mean(axis=1, keepdims=True)
    b = -math.expm1(-rate*dt)/rate
    initial = b*(f*decay[::-1]/decay[0]).sum(axis=1)/(-math.expm1(-rate*period))
    after = decay*(initial[:, None]+b*np.cumsum(f/decay, axis=1))
    before = np.column_stack((initial, after[:, :-1]))
    damped = before*math.exp(-rate*dt/2)+f*(-math.expm1(-rate*dt/2)/rate)
    damped_neutral = -damped.sum(axis=0) if four else np.zeros(n)
    return dict(phase_rms=phase_rms, phase_peak=phase_peak,
                neutral_rms=float(np.sqrt(np.mean(neutral**2))), neutral_peak=float(np.max(np.abs(neutral))),
                frozen_dc_min=float(frozen_voltage.min()), frozen_dc_max=float(frozen_voltage.max()),
                coupled_dc_peak=float(np.max(np.abs(charge))),
                damped_phase_rms=np.sqrt(np.mean(damped**2, axis=1)),
                damped_phase_peak=np.max(np.abs(damped), axis=1),
                damped_neutral_rms=float(np.sqrt(np.mean(damped_neutral**2))),
                damped_neutral_peak=float(np.max(np.abs(damped_neutral))),
                damped_closure_A=float(max(abs(after[:, -1]-initial))),
                sampled_closure_A=float(np.max(np.abs(dt*forcing.sum(axis=1)))))


def audit(inputs, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(dict(accepted=False, reason="AUDIT_INCOMPLETE"))+"\n")
    checks = 0
    def require(test, message):
        nonlocal checks
        checks += 1
        if not test:
            raise ValueError(message)
    try:
        env = json.loads(inputs.read_text())
        for file, digest in env["source_sha256"].items():
            require(hashlib.sha256((REPO/file).read_bytes()).hexdigest() == digest, f"stale source: {file}")
        d = env["domain"]
        require(d == DECLARED_DOMAIN, "operating domain or comparison bands changed")
        keys = ("p_W", "q_var", "negative_sequence_A", "sequence_angle_rad", "modulation_fraction", "carrier_line_angle_rad")
        plan = set(itertools.product(("THREE_LEG", "FOUR_LEG"), *(d[k] for k in keys)))
        actual = {(r["topology"], r["p_W"], r["q_var"], r["i2_A"], r["sequence_angle_rad"],
                   r["modulation_fraction"], r["carrier_line_angle_rad"]) for r in env["cases"]}
        require(actual == plan and len(actual) == len(env["cases"]), "missing or duplicated operating-domain case")
        require({r["id"] for r in env["cases"]} == set(range(1, len(plan)+1)), "invalid case IDs")
        require(not env["sizing_changed"], "sizing changed")
        hw = env["hardware"]
        require((hw["s_max"], hw["i_max"], hw["v_dc"], hw["c_dc"], hw["i_cap_max"]) ==
                (6000., 10., 750., .0003, 1.), "unexpected hardware")
        records = []
        a = complex(-.5, math.sqrt(3)/2)
        for ref in env["cases"]:
            ident = ref["id"]
            u = np.array([complex(*z) for z in ref["converter_voltage"]])
            i = np.array([complex(*z) for z in ref["converter_current"]])
            require(abs(np.dot(u, i.conj())-complex(ref["p_W"], ref["q_var"])) < 1e-8, f"{ident}: P/Q declaration")
            require(abs(abs((i[0]+a*a*i[1]+a*i[2])/3)-ref["i2_A"]) < 1e-10, f"{ident}: I2 declaration")
            four = ref["topology"] == "FOUR_LEG"
            require(abs(abs(i.sum()/3)-(d["four_leg_zero_sequence_A"] if four else 0.)) < 1e-10, f"{ident}: I0 declaration")
            rail = math.sqrt(2)*max(abs(u[j]-u[k]) for j, k in ((0, 1), (1, 2), (2, 0)))
            if four:
                rail = max(rail, math.sqrt(2)*max(abs(u)))
            require(abs(rail/(hw["v_dc"]*hw["m_max"])-ref["modulation_fraction"]) < 1e-12, f"{ident}: modulation declaration")
            instantaneous = math.sqrt(2)*np.real(u*np.exp(1j*ref["carrier_line_angle_rad"]))
            if four:
                instantaneous = np.append(instantaneous, 0.)
            duties = .5+(instantaneous-(max(instantaneous)+min(instantaneous))/2)/hw["v_dc"]
            require(np.max(np.abs(duties-ref["duties"])) < 1e-12 and np.all((duties >= 0)&(duties <= 1)), f"{ident}: clipped/inconsistent gates")
            require(bool(ref["failed_engineering_checks"]) != ref["engineering_limits_passed"], f"{ident}: engineering failure retention")
            for key, bound in independent_bounds(ref, hw, env["operating_limits"]).items():
                require(np.max(np.abs(np.asarray(bound)-ref[key])) < 1e-10, f"{ident}: incorrect analytical {key}")
            meshes = [sampled_carrier(ref, hw, n) for n in d["sampled_steps"]]
            coarse, fine = meshes[-2:]
            event_error = max(abs(fine["phase_rms"]-ref["event_phase_rms_A"]))
            dc_error = max(abs(fine["frozen_dc_min"]-ref["event_dc_min_V"]), abs(fine["frozen_dc_max"]-ref["event_dc_max_V"]))
            require(event_error <= d["event_rms_tolerance_A"], f"{ident}: event/sample AC disagreement {event_error}")
            require(dc_error <= d["event_dc_tolerance_V"], f"{ident}: event/sample DC disagreement {dc_error}")
            for sample in meshes:
                require(np.all(sample["phase_rms"] <= np.array(ref["phase_rms_bound_A"])+1e-6), f"{ident}: phase RMS envelope")
                require(np.all(sample["phase_peak"] <= np.array(ref["phase_peak_bound_A"])+1e-6), f"{ident}: phase peak envelope")
                require(sample["neutral_rms"] <= ref["neutral_rms_bound_A"]+1e-6, f"{ident}: neutral RMS envelope")
                require(sample["neutral_peak"] <= ref["neutral_peak_bound_A"]+1e-6, f"{ident}: neutral peak envelope")
                require(sample["coupled_dc_peak"] <= ref["dc_voltage_bound_V"]+1e-6, f"{ident}: coupled charge envelope")
                require(np.all(sample["damped_phase_rms"] <= np.array(ref["phase_rms_bound_A"])+1e-6), f"{ident}: passive phase RMS envelope")
                require(np.all(sample["damped_phase_peak"] <= np.array(ref["phase_peak_bound_A"])+1e-6), f"{ident}: passive phase peak envelope")
                require(sample["damped_neutral_rms"] <= ref["neutral_rms_bound_A"]+1e-6, f"{ident}: passive neutral RMS envelope")
                require(sample["damped_neutral_peak"] <= ref["neutral_peak_bound_A"]+1e-6, f"{ident}: passive neutral peak envelope")
                require(sample["damped_closure_A"] < 1e-9, f"{ident}: passive periodic closure")
            mesh_error = max(abs(fine["phase_rms"]-coarse["phase_rms"]))
            dc_mesh_error = abs(fine["coupled_dc_peak"]-coarse["coupled_dc_peak"])
            require(mesh_error <= d["event_rms_tolerance_A"], f"{ident}: AC mesh convergence")
            require(dc_mesh_error <= d["event_dc_tolerance_V"], f"{ident}: DC mesh convergence")
            records.append(dict(id=ident, engineering_limits_passed=ref["engineering_limits_passed"],
                failed_engineering_checks=ref["failed_engineering_checks"], event_rms_error_A=float(event_error),
                event_dc_error_V=dc_error, mesh_rms_error_A=float(mesh_error), mesh_dc_error_V=dc_mesh_error,
                coupled_dc_peak_V=fine["coupled_dc_peak"], dc_bound_V=ref["dc_voltage_bound_V"],
                maximum_phase_rms_A=float(max(fine["phase_rms"])), maximum_phase_peak_A=float(max(fine["phase_peak"])),
                neutral_rms_A=fine["neutral_rms"], neutral_peak_A=fine["neutral_peak"],
                sampled_closure_A=fine["sampled_closure_A"]))
        report = dict(accepted=True, checks=checks, cases=len(records), records=records,
            engineering_limit_failures=sum(not r["engineering_limits_passed"] for r in records),
            input_sha256=hashlib.sha256(inputs.read_bytes()).hexdigest(),
            audit_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), numpy_version=np.__version__,
            scope=env["scope"]+"; analytical carrier-domain validation is separate from engineering feasibility and coupled EMT")
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        return report
    except Exception as exc:
        output.write_text(json.dumps(dict(accepted=False, reason="AUDIT_FAILED", error=str(exc)))+"\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = audit(args.inputs, args.output)
    print(f"Carrier-domain audit: {result['cases']} cases, {result['checks']} checks; accepted={result['accepted']}")
