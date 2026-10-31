"""Build compact review artifacts from the complete audited campaign."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

from audit_engineering_benchmark_emt import audit as audit_matched_emt
from audit_engineering_domain import audit as audit_domain
from audit_emt import audit as audit_circuit
from plot_engineering_benchmark import audited_records, ROOT


def summarize(*, check=False):
    folder = ROOT/"results/engineering_benchmark"
    env, benchmark = audited_records(folder)
    inputs = ROOT/"results/engineering_benchmark_emt_inputs.json"
    refs = json.loads(inputs.read_text())
    # Replay audits separately: cross-platform floating arithmetic may alter
    # receipt bytes within the declared tolerances. Preserve immutable original
    # receipts, including in --check mode, and compare the compact semantics.
    with tempfile.TemporaryDirectory() as scratch:
        scratch = Path(scratch)
        emt_folder = shutil.copytree(ROOT/"results/engineering_benchmark_emt", scratch/"emt")
        boundary_folder = shutil.copytree(ROOT/"results/engineering_boundary_emt", scratch/"boundary")
        emt = audit_matched_emt(inputs, emt_folder, folder)
        domain = audit_domain(ROOT/"results/engineering_domain_inputs.json", scratch/"domain_audit.json")
        boundary = audit_circuit(ROOT/"results/engineering_boundary_emt_inputs.json", boundary_folder)
    boundary_results = json.loads((ROOT/"results/engineering_boundary_emt/results.json").read_text())
    retained_failure = not boundary_results["cases"]["lf_only_boundary"]["modes"]["switched"]["finest"]["all_engineering_checks_passed"]
    if not (emt["accepted"] and emt["all_thermal_predictions_enclose"] and domain["accepted"] and boundary["accepted"]):
        raise ValueError("complete current acceptance gates required")
    if not retained_failure:
        raise ValueError("retained LF-only hardware failure is missing")
    compact = [c for c in refs["cases"] if c["parent_case_id"] == "night" and c["parent_mode"] == "coordinated"
               and c["local_controller"] in ("mean_vv_vw", "phase_vv_vw")]
    if len(compact) != 2 or {c["topology"] for c in compact} != {"THREE_LEG", "FOUR_LEG"}:
        raise ValueError("compact reference coverage changed")
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    output = ROOT/"evidence"
    output.mkdir(exist_ok=True)
    artifacts = {"compact_references.json": dict(
        description="Fixed-size three/four-leg night references; full records restored from current.json",
        inputs_sha256=digest(inputs), cases=compact)}
    summary = dict(
        status="current", thermal_model=env["pwm_thermal_model"], sizing_changed=env["sizing_changed"],
        generation_revision=env["git_revision"], generation_julia=env["julia_version"],
        pf_opf=dict(points=benchmark["accepted"], inverter_rows=sum(len(c["devices"]) for c in benchmark["cases"].values()),
            record_checks=benchmark["conditions"], physical_checks=benchmark["physical_checks"],
            independent_network_checks=benchmark["network_conditions"]),
        carrier_domain=dict(points=domain["cases"], checks=domain["checks"],
            retained_engineering_failures=domain["engineering_limit_failures"]),
        boundary=dict(checks=boundary["checks"], retained_lf_only_hardware_failure=retained_failure),
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
    artifacts["current_summary.json"] = summary
    for name, value in artifacts.items():
        encoded = json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n"
        target = output/name
        if check:
            if not target.is_file() or target.read_text() != encoded:
                raise ValueError(f"stale or corrupted compact evidence: {name}")
        else:
            target.write_text(encoded)
    print("Verified" if check else "Wrote", "compact source-linked summary and two physical reference fixtures")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    summarize(check=parser.parse_args().check)
