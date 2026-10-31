"""Build compact review artifacts from the complete audited campaign."""
import hashlib
import json
from pathlib import Path

from audit_engineering_benchmark_emt import validate_references
from plot_engineering_benchmark import audited_records, ROOT


def summarize():
    folder = ROOT/"results/engineering_benchmark"
    env, benchmark = audited_records(folder)
    inputs = ROOT/"results/engineering_benchmark_emt_inputs.json"
    validate_references(inputs, folder)
    refs = json.loads(inputs.read_text())
    emt = json.loads((ROOT/"results/engineering_benchmark_emt/matched_audit.json").read_text())
    domain = json.loads((ROOT/"results/engineering_domain_audit.json").read_text())
    boundary = json.loads((ROOT/"results/engineering_boundary_emt/audit.json").read_text())
    if not (emt["accepted"] and emt["all_thermal_predictions_enclose"] and domain["accepted"] and boundary["accepted"]):
        raise ValueError("complete current acceptance gates required")
    compact = [c for c in refs["cases"] if c["parent_case_id"] == "night" and c["parent_mode"] == "coordinated"
               and c["local_controller"] in ("mean_vv_vw", "phase_vv_vw")]
    if len(compact) != 2 or {c["topology"] for c in compact} != {"THREE_LEG", "FOUR_LEG"}:
        raise ValueError("compact reference coverage changed")
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    output = ROOT/"evidence"
    output.mkdir(exist_ok=True)
    (output/"compact_references.json").write_text(json.dumps(dict(
        description="Fixed-size three/four-leg night references; full records restored from current.json",
        inputs_sha256=digest(inputs), cases=compact), indent=2, sort_keys=True)+"\n")
    summary = dict(
        status="current", thermal_model=env["pwm_thermal_model"], sizing_changed=env["sizing_changed"],
        generation_revision=env["git_revision"], generation_julia=env["julia_version"],
        pf_opf=dict(points=benchmark["accepted"], inverter_rows=96,
            record_checks=benchmark["conditions"], physical_checks=benchmark["physical_checks"],
            independent_network_checks=benchmark["network_conditions"]),
        carrier_domain=dict(points=domain["cases"], checks=domain["checks"],
            retained_engineering_failures=domain["engineering_limit_failures"]),
        boundary=dict(checks=boundary["checks"], retained_lf_only_hardware_failure=True),
        emt=dict(references=emt["references"], integrations=emt["integrations"],
            reference_checks=emt["reference_checks"], circuit_checks=emt["circuit_checks"],
            all_thermal_predictions_enclose=emt["all_thermal_predictions_enclose"],
            thermal_prediction_failures=emt["thermal_prediction_failures"],
            maximum_switched_capacitor_A=max(r["switched_A"] for r in emt["thermal_prediction_coverage"].values()),
            capacitor_rating_A=1., minimum_prediction_margin_A=min(r["predicted_A"]-r["switched_A"] for r in emt["thermal_prediction_coverage"].values())),
        input_sha256={str(p.relative_to(ROOT)): digest(p) for p in (
            folder/"environment.json", folder/"audit.json", inputs,
            ROOT/"results/engineering_benchmark_emt/matched_audit.json",
            ROOT/"results/engineering_domain_inputs.json", ROOT/"results/engineering_domain_audit.json",
            ROOT/"results/engineering_boundary_emt_inputs.json", ROOT/"results/engineering_boundary_emt/audit.json")},
        scope="Service choices differ: free-Q OPF releases VV/VW. Periodic frozen-slow-quantity carrier bounds and fixed-POC EMT; no global optimality, regional conformity, feeder/controller dynamic or all-time transient certificate.")
    (output/"current_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n")
    print("Wrote compact source-linked summary and two physical reference fixtures")


if __name__ == "__main__":
    summarize()
