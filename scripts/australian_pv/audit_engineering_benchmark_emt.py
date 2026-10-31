"""Check matched-case EMT selection/availability and the independent circuit records."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

from audit_emt import audit as audit_circuit

CONTROLLERS = {"mean_vv_vw", "worst_vv_vw", "sequence_vv_vw", "phase_vv_vw"}
CONDITIONS = {"nominal", "high_source", "rotated_unbalance", "night", "oversized_dc"}
MODES = ("local", "opf_equal_service", "opf_vv_released")


def validate_references(inputs, benchmark):
    env = json.loads(inputs.read_text())
    parent_env = benchmark/"environment.json"
    parent_audit = benchmark/"audit.json"
    truth = json.loads(parent_audit.read_text())
    checks = 0
    def require(test, message):
        nonlocal checks
        checks += 1
        if not test:
            raise ValueError(message)
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    require(truth["all_planned_accepted"] and truth["accepted"] == len(CONDITIONS)*len(MODES), "incomplete parent benchmark")
    require(env["benchmark_environment_sha256"] == digest(parent_env) == truth["environment_sha256"], "parent environment changed")
    require(env["benchmark_audit_sha256"] == digest(parent_audit), "parent audit changed")
    require(env["benchmark_record_sha256"] == truth["record_sha256"], "parent point plan changed")
    plan = set(itertools.product(CONDITIONS, MODES, CONTROLLERS))
    seen = {(c["parent_case_id"], c["parent_mode"], c["local_controller"]) for c in env["cases"]}
    require(seen == plan and len(env["cases"]) == len(plan), "missing or duplicated EMT selection")
    require(len({c["case_id"] for c in env["cases"]}) == len(plan), "duplicate EMT IDs")
    limits = env["operating_limits"]
    require(env["fixture_limits"] == "declared", "EMT campaign must use the declared fixture")
    for c in env["cases"]:
        key = c["parent_case_id"]+"_"+c["parent_mode"]
        parent_path = "points/"+key+".json"
        require(digest(benchmark/parent_path) == c["parent_point_sha256"] == truth["record_sha256"][parent_path], "parent phasors changed")
        point = json.loads((benchmark/parent_path).read_text())
        rows = [r for r in truth["cases"][key]["devices"] if r["local_controller"] == c["local_controller"]]
        row = next(r for r in rows if r["device_id"] == c["device_id"])
        require(row["capacitor_thermal_A"] >= max(r["capacitor_thermal_A"] for r in rows)-1e-8, "selected site is not a maximum-stress representative")
        require(abs(c["parent_capacitor_model_A"]-row["capacitor_thermal_A"]) <= 1e-8, "parent stress mismatch")
        require(c["p_available_W"] == row["p_available_W"], "site availability replaced by fixture availability")
        require(abs(c["p_dc_W"]-row["p_dc_W"]) <= 1e-8, "DC power mismatch")
        require(c["topology"] == row["topology"], "topology mismatch")
        require(c["expect_averaged_hardware_accepted"] and c["expect_switched_hardware_accepted"], "engineering acceptance expectation relaxed")
        require(c["pwm_current_envelope"] and c["pwm_dc_voltage_envelope"], "engineering guard disabled")
        require((c["dc_voltage_V"], c["dc_capacitance_F"], c["phase_current_rating_A"], c["capacitor_current_rating_A"],
                 c["phase_peak_rating_A"], c["dc_voltage_min_rating_V"], c["dc_voltage_max_rating_V"]) ==
                (750., .0003, 10., limits["capacitor_rms_A"], limits["phase_peak_A"], limits["dc_min_V"], limits["dc_max_V"]),
                "EMT hardware or operating ceiling changed")
        saved = point["devices"][c["device_id"]]
        require(c["current"] == saved["grid_current"] and c["converter_voltage"] == saved["converter_voltage"], "reference current/voltage changed")
        bus = point["bus"][row["bus"]]
        neutral = complex(bus["n"]["vr"], bus["n"]["vi"]) if c["topology"] == "FOUR_LEG" else 0j
        voltage = [complex(bus[p]["vr"], bus[p]["vi"])-neutral for p in "abc"]
        require(max(abs(complex(*a)-b) for a, b in zip(c["grid_voltage"], voltage)) < 1e-10, "fixed POC boundary changed")
        power = sum(complex(*u)*complex(*i).conjugate() for u, i in zip(c["converter_voltage"], c["current"]))
        require(abs(c["p_dc_W"]-power.real-c["mean_converter_loss_W"]) < 1e-8, "mean loss boundary mismatch")
    return dict(checks=checks, references=len(plan), input_sha256=digest(inputs),
                parent_audit_sha256=digest(parent_audit))


def audit(inputs, folder, benchmark):
    target = folder/"matched_audit.json"
    target.write_text(json.dumps(dict(accepted=False, reason="AUDIT_INCOMPLETE"))+"\n")
    try:
        refs = validate_references(inputs, benchmark)
        # Mechanism, mesh and energy agreement are mandatory. A switched hardware
        # check that fails where the model predicted a pass is a reported finding,
        # with the exact quantity and margin, never silently accepted or hidden.
        circuit = audit_circuit(inputs, folder, report_expectation_failures=True)
        results = json.loads((folder/"results.json").read_text())
        if not circuit["accepted"]:
            raise ValueError("matched EMT mechanism check failed")
        env = json.loads(inputs.read_text())
        refs_by_id = {c["case_id"]: c for c in env["cases"]}
        hardware_findings = []
        for failure in circuit["expectation_failures"]:
            fine = results["cases"][failure["case_id"]]["modes"][failure["mode"]]["finest"]
            ref = refs_by_id[failure["case_id"]]
            detail = dict(failure)
            if "dc_2omega_peak" in failure["failed_engineering_checks"]:
                detail.update(dc_2omega_peak_V=abs(complex(*fine["dc_2omega_peak"])),
                              limit_V=ref["dc_2omega_peak_rating_V"],
                              averaged_dc_2omega_peak_V=abs(complex(*results["cases"][failure["case_id"]]["modes"]["averaged"]["finest"]["dc_2omega_peak"])))
            hardware_findings.append(detail)
        predictions = {}
        for c in env["cases"]:
            fine = results["cases"][c["case_id"]]["modes"]["switched"]["finest"]
            predicted, measured = c["parent_capacitor_model_A"], fine["capacitor_total_rms_A"]
            predictions[c["case_id"]] = dict(predicted_A=predicted, switched_A=measured,
                error_A=measured-predicted,
                encloses_with_mesh_band=measured <= predicted+env["validation_bands"]["mesh_current_A"])
        missed = [key for key, row in predictions.items() if not row["encloses_with_mesh_band"]]
        report = dict(accepted=True, reference_checks=refs["checks"], references=refs["references"],
            circuit_checks=circuit["checks"], integrations=len(results["metrics"]),
            input_sha256=refs["input_sha256"], parent_audit_sha256=refs["parent_audit_sha256"],
            results_sha256=hashlib.sha256((folder/"results.json").read_bytes()).hexdigest(),
            audit_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            thermal_prediction_coverage=predictions, thermal_prediction_failures=missed,
            all_thermal_predictions_enclose=not missed,
            switched_hardware_findings=hardware_findings,
            all_switched_hardware_accepted=not hardware_findings,
            scope="declared worst-capacitor-stress selection, actual site availability, saved parent phasors, fixed mean losses and independent averaged/switched circuit records; no controller/feeder dynamics or all-time DC-drift certificate")
        target.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        return report
    except Exception as exc:
        target.write_text(json.dumps(dict(accepted=False, reason="AUDIT_FAILED", error=str(exc)))+"\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("results", type=Path)
    parser.add_argument("benchmark", type=Path)
    args = parser.parse_args()
    report = audit(args.inputs, args.results, args.benchmark)
    print(f"Matched EMT: {report['references']} references, {report['integrations']} integrations, "
          f"{report['reference_checks']} parent checks, {report['circuit_checks']} circuit checks")
