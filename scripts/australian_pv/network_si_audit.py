"""Independent SI circuit/limit reconstruction for the vendored Australian LV case.

Uses only JSON inputs and NumPy, not PowerOptLab/BMOPFTools or solver auxiliaries.
Scope: the declared delta-wye source transformer, four-conductor pi lines,
closed ideal switches, ground shunts, constant-PQ customer loads and saved PV
terminal currents. Conductor limits here concern fundamental-frequency current;
this does not simulate propagation of inverter carrier harmonics in the feeder.
"""
import math

import numpy as np


def matrix(record, prefix, size):
    return np.array([[record.get(f"{prefix}_{i+1}_{j+1}", 0.) for j in range(size)] for i in range(size)])


def audit_network(network, point, case, *, voltage_max=253.):
    checks = []

    def check(name, unit, value, lower=-math.inf, upper=math.inf, tolerance=1e-5):
        value, lower, upper, tolerance = map(float, (value, lower, upper, tolerance))
        checks.append(dict(name=name, unit=unit, value=value if math.isfinite(value) else None,
                           lower=lower if math.isfinite(lower) else None, upper=upper if math.isfinite(upper) else None,
                           tolerance=tolerance, passed=math.isfinite(value) and lower-tolerance <= value <= upper+tolerance))

    def zero(name, unit, residual, tolerance):
        check(name, unit, abs(residual), 0., 0., tolerance)

    volts = {b: {t: complex(v["vr"], v["vi"]) for t, v in ts.items()} for b, ts in point["bus"].items()}
    if set(volts) != set(network["bus"]) or any(set(volts[b]) != set(d["terminal_names"]) for b, d in network["bus"].items()):
        raise ValueError("bus/terminal coverage differs from the circuit")
    balance = {b: {t: 0j for t in ts} for b, ts in volts.items()}
    transformers = network.get("transformer", {})
    if set(transformers) != {"delta_wye"} or len(transformers["delta_wye"]) != 1 or len(network["voltage_source"]) != 1:
        raise ValueError("audit requires the declared single-source delta-wye circuit")
    tx_id, tx = next(iter(transformers["delta_wye"].items()))
    if any(tx.get(k, 0.) != 0 for k in ("g_no_load", "b_no_load", "r_series_from", "x_series_from")):
        raise ValueError("unsupported primary/no-load transformer branch")
    if tx["terminal_map_from"] != list("abc") or tx["terminal_map_to"] != list("abcn"):
        raise ValueError("unexpected transformer winding order")
    source_id, source = next(iter(network["voltage_source"].items()))
    if source["bus"] != tx["bus_from"]:
        raise ValueError("source must connect directly to the primary")
    old_source = np.array(source["v_magnitude"][:3])*np.exp(1j*np.array(source["v_angle"][:3]))
    a = np.exp(2j*np.pi/3)
    u1 = (old_source[0]+a*old_source[1]+a*a*old_source[2])/3*case["source_pu"]
    u2 = u1*case["source_vuf"]*np.exp(1j*case["source_negative_angle_rad"])
    expected_source = u1*np.array([1, a*a, a])+u2*np.array([1, a, a*a])
    for phase, value in zip("abc", expected_source):
        zero(f"source/{source_id}/{phase}", "V", volts[source["bus"]][phase]-value, 1e-6)
    zero(f"source/{source_id}/neutral", "V", volts[source["bus"]]["n"], 1e-6)

    line_loss = 0.; shunt_loss = 0.; load_power = 0j; pv_power = 0j
    conductor_ratios = []; line_currents = {}; lv_buses = {tx["bus_to"]}
    for ident, line in network["line"].items():
        if line["terminal_map_from"] != list("abcn") or line["terminal_map_to"] != list("abcn"):
            raise ValueError("audit requires four-conductor LV pi lines")
        data = network["linecode"][line["linecode"]]; length = line["length"]
        if not math.isfinite(length) or length <= 0:
            raise ValueError("invalid line length")
        impedance = length*(matrix(data, "R_series", 4)+1j*matrix(data, "X_series", 4))
        vf = np.array([volts[line["bus_from"]][t] for t in "abcn"])
        vt = np.array([volts[line["bus_to"]][t] for t in "abcn"])
        series = np.linalg.solve(impedance, vf-vt)
        yf = length*(matrix(data, "G_from", 4)+1j*matrix(data, "B_from", 4))
        yt = length*(matrix(data, "G_to", 4)+1j*matrix(data, "B_to", 4))
        end_currents = (series+yf@vf, -series+yt@vt)
        limits = line.get("i_max", data.get("i_max"))
        if limits is not None:
            limits = np.broadcast_to(np.array(limits), (4,))
            if not np.all(np.isfinite(limits)) or not np.all(limits > 0):
                raise ValueError("invalid conductor rating")
        for side, bus, current in zip(("from", "to"), (line["bus_from"], line["bus_to"]), end_currents):
            lv_buses.add(bus)
            for k, phase in enumerate("abcn"):
                balance[bus][phase] += current[k]
                if limits is not None:
                    check(f"line/{ident}/{side}/{phase}", "A", abs(current[k]), 0., limits[k])
                    conductor_ratios.append(abs(current[k])/limits[k])
        line_loss += float(np.real(vf@np.conj(end_currents[0])+vt@np.conj(end_currents[1])))
        line_currents[ident] = {side: [[z.real, z.imag] for z in flow] for side, flow in zip(("from", "to"), end_currents)}
    if source["bus"] in lv_buses:
        raise ValueError("unexpected primary-side line in LV audit")
    for ident, shunt in network["shunt"].items():
        terms = shunt["terminal_map"]; size = len(terms)
        voltage = np.array([volts[shunt["bus"]][t] for t in terms])
        current = (matrix(shunt, "G", size)+1j*matrix(shunt, "B", size))@voltage
        for terminal, value in zip(terms, current):
            balance[shunt["bus"]][terminal] += value
        shunt_loss += float(np.real(voltage@np.conj(current)))

    if set(point["devices"]) != {"pv_"+name for name in network["load"]}:
        raise ValueError("PV/customer coverage mismatch")
    shares = np.array(case["phase_shares"])
    if len(shares) != 3 or not np.all(np.isfinite(shares)) or not np.all(shares >= 0) or abs(sum(shares)-1) > 1e-12:
        raise ValueError("invalid load phase allocation")
    for ident, load in network["load"].items():
        bus = load["bus"]
        for phase, share in zip("abc", shares):
            power = 1500*case["load_scale"]*share*(1+1j*math.tan(math.acos(.95)))
            current = np.conj(power/(volts[bus][phase]-volts[bus]["n"]))
            balance[bus][phase] += current; balance[bus]["n"] -= current
            load_power += power
        current = np.array([complex(*z) for z in point["devices"]["pv_"+ident]["grid_current"]])
        if current.shape != (3,) or not np.all(np.isfinite(current)):
            raise ValueError("invalid PV current")
        for phase, value in zip("abc", current):
            balance[bus][phase] -= value
            pv_power += (volts[bus][phase]-volts[bus]["n"])*np.conj(value)
        balance[bus]["n"] += sum(current)

    parent = {b: b for b in volts}

    def find(bus):
        if parent[bus] != bus:
            parent[bus] = find(parent[bus])
        return parent[bus]

    for ident, switch in network["switch"].items():
        if switch["open_switch"] or switch["terminal_map_from"] != list("abcn") or switch["terminal_map_to"] != list("abcn"):
            raise ValueError("audit requires the declared closed four-pole switches")
        f, t = switch["bus_from"], switch["bus_to"]
        lv_buses.update((f, t))
        if find(f) == find(t):
            raise ValueError("ideal-switch loops need an additional current-sharing model")
        parent[find(f)] = find(t)
        for phase in "abcn":
            zero(f"switch/{ident}/voltage/{phase}", "V", volts[f][phase]-volts[t][phase], 1e-6)
    if lv_buses != set(volts)-{source["bus"]}:
        raise ValueError("LV planning-limit coverage differs from the circuit")
    components = {}
    for bus in volts:
        components.setdefault(find(bus), []).append(bus)
    root_members = components[find(tx["bus_to"])]
    secondary = np.array([sum(balance[b][t] for b in root_members) for t in "abc"])
    turn_ratio = math.sqrt(3)*tx["v_nom_from"]/tx["v_nom_to"]*tx["tap"]
    primary_voltage = np.array([volts[tx["bus_from"]][t] for t in "abc"])
    delta = np.array([[1., 0., -1.], [-1., 1., 0.], [0., -1., 1.]])
    winding = delta@primary_voltage/turn_ratio
    secondary_voltage = np.array([volts[tx["bus_to"]][t]-volts[tx["bus_to"]]["n"] for t in "abc"])
    ztx = complex(tx["r_series_to"], tx["x_series_to"])
    for phase, error in zip("abc", winding-secondary_voltage-ztx*secondary):
        zero(f"transformer/{tx_id}/kvl/{phase}", "V", error, 1e-6)
    for phase, value in zip("abc", secondary):
        balance[tx["bus_to"]][phase] -= value
    balance[tx["bus_to"]]["n"] += sum(secondary)
    primary = delta.T@secondary/turn_ratio
    source_power = primary_voltage@np.conj(primary)
    transformer_loss = float(ztx.real*sum(abs(secondary)**2))
    check(f"transformer/{tx_id}/apparent_power", "VA", abs(winding@np.conj(secondary)), 0., tx["s_rating"], .01)

    switch_currents = {}
    for members in components.values():
        if source["bus"] in members:
            continue # The ideal voltage source supplies the reconstructed primary current.
        edges = [(i, d) for i, d in network["switch"].items() if d["bus_from"] in members]
        incidence = np.zeros((len(members), len(edges)))
        for j, (_, switch) in enumerate(edges):
            incidence[members.index(switch["bus_from"]), j] = 1
            incidence[members.index(switch["bus_to"]), j] = -1
        for phase_index, phase in enumerate("abcn"):
            rhs = -np.array([balance[b][phase] for b in members])
            flow = np.linalg.lstsq(incidence, rhs, rcond=None)[0] if edges else np.empty(0)
            residual = incidence@flow-rhs
            for bus, error in zip(members, residual):
                zero(f"network/kcl/{bus}/{phase}", "A", error, 1e-5)
            for (ident, switch), current in zip(edges, flow):
                switch_currents.setdefault(ident, {})[phase] = [current.real, current.imag]
                limits = switch.get("i_max")
                if limits is not None:
                    rating = limits[phase_index]
                    if not math.isfinite(rating) or rating <= 0:
                        raise ValueError("invalid switch current rating")
                    check(f"switch/{ident}/current/{phase}", "A", abs(current), 0., rating)
                    conductor_ratios.append(abs(current)/rating)
    for bus in sorted(lv_buses):
        for phase in "abc":
            if voltage_max is not None:
                check(f"planning_voltage/{bus}/{phase}", "V", abs(volts[bus][phase]-volts[bus]["n"]), 0., voltage_max)
    residual = source_power.real+pv_power.real-load_power.real-line_loss-shunt_loss-transformer_loss
    zero("network/active_power_balance", "W", residual, .01)
    return dict(accepted=all(c["passed"] for c in checks), checks=checks,
                maximum_conductor_utilization=float(max(conductor_ratios)),
                source_active_power_W=float(source_power.real), pv_active_power_W=float(pv_power.real),
                load_active_power_W=float(load_power.real), line_loss_W=line_loss, grounding_loss_W=shunt_loss,
                transformer_loss_W=transformer_loss, active_power_residual_W=float(residual),
                line_currents=line_currents, switch_currents=switch_currents,
                scope="independent SI circuit reconstruction and fundamental network limits; no carrier-harmonic feeder propagation")


if __name__ == "__main__":
    import argparse
    import hashlib
    import json
    from pathlib import Path
    import sys

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("network", type=Path)
    parser.add_argument("point", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(accepted=False, reason="AUDIT_INCOMPLETE"))+"\n")
    try:
        network_bytes, point_bytes = args.network.read_bytes(), args.point.read_bytes()
        point = json.loads(point_bytes)
        result = audit_network(json.loads(network_bytes), point, point["case"],
                               voltage_max=point.get("planning_voltage_max_V", 253.))
        result.update(network_sha256=hashlib.sha256(network_bytes).hexdigest(),
                      point_sha256=hashlib.sha256(point_bytes).hexdigest(),
                      audit_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      python_version=sys.version, numpy_version=np.__version__)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        print(f"SI network audit: {len(result['checks'])} conditions; accepted={result['accepted']}")
    except Exception as exc:
        args.output.write_text(json.dumps(dict(accepted=False, reason="AUDIT_FAILED", error=str(exc)))+"\n")
        raise
    sys.exit(0 if result["accepted"] else 1)
