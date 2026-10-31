"""Independent saved-record checks for the time-stepped circuit experiment."""
import csv
import hashlib
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def audit(input_path, folder):
    target = folder / "audit.json"
    target.write_text(json.dumps({"accepted": False, "reason": "AUDIT_INCOMPLETE"}))
    checks = 0
    def require(condition, message):
        nonlocal checks
        checks += 1
        if not condition:
            raise ValueError(message)
    def close(left, right, tolerance=1e-7):
        require(math.isfinite(left) and math.isfinite(right) and abs(left-right) <= tolerance,
                f"metric mismatch: {left} vs {right}")
    try:
        env = json.loads(input_path.read_text())
        result = json.loads((folder / "results.json").read_text())
        require(hashlib.sha256(input_path.read_bytes()).hexdigest() == result["input_sha256"], "input digest mismatch")
        require(hashlib.sha256((REPO / "scripts/australian_pv/emt_circuit.jl").read_bytes()).hexdigest() == result["circuit_sha256"], "circuit digest mismatch")
        for path, digest in env["source_sha256"].items():
            require(hashlib.sha256((REPO / path).read_bytes()).hexdigest() == digest, f"stale reference source: {path}")
        for file, digest in result["waveform_sha256"].items():
            require(hashlib.sha256((folder / file).read_bytes()).hexdigest() == digest, f"waveform digest mismatch: {file}")
        refs = {c["case_id"]: c for c in env["cases"]}
        expected_waveforms = {f"{i}_{m}_{suffix}.csv" for i in refs for m in ("averaged", "switched") for suffix in ("carrier", "switching")}
        require(set(result["waveform_sha256"]) == expected_waveforms, "waveform coverage mismatch")
        require(len(refs) == len(env["cases"]) == result["planned"] == result["completed"], "incomplete/duplicate case plan")
        require(set(result["cases"]) == set(refs), "case coverage mismatch")
        timing = env["timing"]; bands = env["validation_bands"]
        expected_keys = {(i, m, n) for i in refs for m in ("averaged", "switched") for n in timing["steps_per_carrier"]}
        rows = {(r["case_id"], r["mode"], r["steps_per_carrier"]): r for r in result["metrics"]}
        require(len(rows) == len(result["metrics"]) and set(rows) == expected_keys, "mesh/mode coverage mismatch")
        evidence = {}
        all_hardware = []; all_mechanisms = []
        for ident, c in refs.items():
            vg = [complex(*p) for p in c["grid_voltage"]]
            uc = [complex(*p) for p in c["converter_voltage"]]
            ir = [complex(*p) for p in c["current"]]
            omega2 = 4*math.pi*c["frequency_Hz"]
            expected_d = 1j*sum(v*i for v, i in zip(uc, ir))/(omega2*c["dc_capacitance_F"]*c["dc_voltage_V"])
            case_modes = {}
            for mode in ("averaged", "switched"):
                for n in timing["steps_per_carrier"]:
                    r = rows[ident, mode, n]
                    fund = [complex(*p) for p in r["current_fundamental"]]
                    d = complex(*r["dc_2omega_peak"])
                    ic2 = complex(*r["capacitor_2omega_peak"])
                    close(abs(complex(*r["expected_dc_2omega_peak"])-expected_d), 0)
                    close(r["fundamental_current_max_error_A"], max(abs(x-y) for x, y in zip(fund, ir)))
                    close(r["dc_2omega_error_V"], abs(d-expected_d))
                    close(r["capacitor_2omega_error_A"], abs(ic2-1j*omega2*c["dc_capacitance_F"]*expected_d))
                    window = (timing["cycles"]-timing["discard_cycles"])/c["frequency_Hz"]
                    close(r["window_s"], window)
                    def energy(state):
                        magnetic = c["phase_L_H"]*sum(i*i for i in state[:3])/2
                        if c["topology"] == "FOUR_LEG":
                            magnetic += c["neutral_L_H"]*sum(state[:3])**2/2
                        return c["dc_capacitance_F"]*state[3]**2/2 + magnetic
                    filter_loss = c["phase_R_ohm"]*sum(i*i for i in r["current_total_rms_A"])
                    if c["topology"] == "FOUR_LEG":
                        filter_loss += c["neutral_R_ohm"]*r["neutral_total_rms_A"]**2
                    net = window*(c["p_dc_W"]-c["mean_converter_loss_W"]-r["p_grid_mean_W"]-filter_loss)
                    close(net, r["window_net_energy_J"], 1e-6)
                    residual = energy(r["window_final_state"])-energy(r["window_initial_state"])-net
                    close(residual, r["energy_residual_J"], 1e-6)
                    close(r["energy_relative_residual"], abs(residual)/max(1, abs(c["p_dc_W"])*window), 1e-8)
                    hardware = {"phase_rms": max(r["current_total_rms_A"]) <= c["phase_current_rating_A"]+1e-5,
                                "neutral_rms": c["topology"] != "FOUR_LEG" or r["neutral_total_rms_A"] <= c["neutral_current_rating_A"]+1e-5,
                                "capacitor_total_rms": r["capacitor_total_rms_A"] <= c["capacitor_current_rating_A"]+1e-5,
                                "dc_2omega_peak": abs(d) <= c["dc_2omega_peak_rating_V"]+1e-5,
                                "mean_dc_source": c["p_dc_W"] <= c["p_available_W"]+.01,
                                "modulation": r["minimum_modulation_margin"] >= -1e-8}
                    require(hardware == r["engineering_checks"], "engineering flags mismatch")
                    require(all(hardware.values()) == r["all_engineering_checks_passed"], "hardware acceptance mismatch")
                    require(r["capacitor_total_rms_A"] >= abs(ic2)/math.sqrt(2)-1e-4, "RMS spectrum inconsistency")
                fine = rows[ident, mode, timing["steps_per_carrier"][-1]]
                coarse = rows[ident, mode, timing["steps_per_carrier"][-2]]
                mesh_i = max(abs(complex(*a)-complex(*b)) for a, b in zip(fine["current_fundamental"], coarse["current_fundamental"]))
                mesh_d = abs(complex(*fine["dc_2omega_peak"])-complex(*coarse["dc_2omega_peak"]))
                mesh_cap = abs(fine["capacitor_total_rms_A"]-coarse["capacitor_total_rms_A"])
                saved = result["cases"][ident]["modes"][mode]
                require(saved["finest"] == fine, "finest result mismatch")
                for key, value in (("mesh_current_difference_A", mesh_i), ("mesh_dc_difference_V", mesh_d), ("mesh_capacitor_rms_difference_A", mesh_cap)):
                    close(value, saved[key])
                dband = bands["dc_2omega_absolute_V"]+bands["dc_2omega_relative"]*abs(expected_d)
                accepted = (fine["fundamental_current_max_error_A"] <= bands["fundamental_absolute_A"]+bands["fundamental_relative"]*max(abs(i) for i in ir)
                            and fine["dc_2omega_error_V"] <= dband
                            and fine["capacitor_2omega_error_A"] <= omega2*c["dc_capacitance_F"]*dband
                            and fine["energy_relative_residual"] <= bands["energy_relative"]
                            and mesh_i <= bands["mesh_current_A"] and mesh_d <= bands["mesh_dc_2omega_V"] and mesh_cap <= bands["mesh_current_A"])
                require(accepted == saved["mechanisms_accepted"], "mechanism acceptance mismatch")
                all_mechanisms.append(accepted)
                if mode == "switched":
                    all_hardware.append(fine["all_engineering_checks_passed"])
                case_modes[mode] = {"mechanisms_accepted": accepted, "capacitor_total_rms_A": fine["capacitor_total_rms_A"],
                                    "failed_engineering_checks": [k for k, v in fine["engineering_checks"].items() if not v]}
            require(result["cases"][ident]["mechanisms_accepted"] == all(v["mechanisms_accepted"] for v in case_modes.values()), "case acceptance mismatch")
            evidence[ident] = case_modes
        require(result["all_mechanisms_accepted"] == all(all_mechanisms), "global mechanism flag mismatch")
        require(result["all_switched_hardware_accepted"] == all(all_hardware), "global hardware flag mismatch")
        report = {"accepted": all(all_mechanisms), "checks": checks, "cases": evidence,
                  "scope": "independent metric, energy, RMS-limit, coverage and provenance reconstruction; mechanism agreement does not certify failed hardware bounds",
                  "input_sha256": result["input_sha256"], "results_sha256": hashlib.sha256((folder / "results.json").read_bytes()).hexdigest()}
        target.write_text(json.dumps(report, indent=2)+"\n")
        return report
    except Exception as exc:
        target.write_text(json.dumps({"accepted": False, "reason": "AUDIT_FAILED", "error": str(exc)})+"\n")
        raise


if __name__ == "__main__":
    import sys
    report = audit(Path(sys.argv[1]), Path(sys.argv[2]))
    print(f"EMT record audit: {report['checks']} checks; mechanism agreement accepted={report['accepted']}")
