"""Publication must match audited phasors, derived limit counts and the PF/OPF agreement."""
import csv
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from plot_engineering_benchmark import audited_records, ROOT


class BenchmarkPublicationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.folder = Path(self.tmp.name)/"benchmark"
        shutil.copytree(ROOT/"results/engineering_benchmark", self.folder)

    def change_csv(self, file, mutate):
        path = self.folder/file
        with path.open() as stream:
            reader = csv.DictReader(stream)
            fields, rows = reader.fieldnames, list(reader)
        mutate(rows)
        with path.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def test_complete_saved_point_publication(self):
        _, truth, cases = audited_records(self.folder)
        self.assertTrue(truth["all_planned_accepted"])
        self.assertEqual(len(cases), 15)

    def test_binding_summary_corrupted(self):
        self.change_csv("cases.csv", lambda rows: rows[0].update(binding_counts="capacitor_current=12"))
        with self.assertRaises(ValueError):
            audited_records(self.folder)

    def test_planning_exceedance_hidden(self):
        def hide(rows):
            row = next(r for r in rows if r["planning_253V_exceeded"] == "true")
            row["planning_253V_exceeded"] = "false"
        self.change_csv("cases.csv", hide)
        with self.assertRaises(ValueError):
            audited_records(self.folder)

    def test_agreement_rejected(self):
        path = self.folder/"audit.json"
        audit = json.loads(path.read_text())
        audit["agreement"]["nominal"]["accepted"] = False
        path.write_text(json.dumps(audit))
        with self.assertRaises(ValueError):
            audited_records(self.folder)

    def test_case_export_corrupted(self):
        self.change_csv("cases.csv", lambda rows: rows[0].update(p_poc_W="12345"))
        with self.assertRaises(ValueError):
            audited_records(self.folder)

    def test_reactive_service_corrupted(self):
        self.change_csv("devices.csv", lambda rows: rows[0].update(q_request_var="12345"))
        with self.assertRaises(ValueError):
            audited_records(self.folder)

    def test_device_row_missing(self):
        self.change_csv("devices.csv", lambda rows: rows.pop())
        with self.assertRaises(ValueError):
            audited_records(self.folder)

    def test_acceptance_changed(self):
        self.change_csv("cases.csv", lambda rows: rows[0].update(accepted="false"))
        with self.assertRaises(ValueError):
            audited_records(self.folder)

    def test_binding_constraint_corrupted(self):
        self.change_csv("devices.csv", lambda rows: rows[0].update(binding_constraints="none"))
        with self.assertRaises(ValueError):
            audited_records(self.folder)

    def test_saved_point_changed(self):
        path = self.folder/"points/nominal_local.json"
        point = json.loads(path.read_text())
        next(iter(point["devices"].values()))["grid_current"][0][0] += .01
        path.write_text(json.dumps(point))
        with self.assertRaises(ValueError):
            audited_records(self.folder)

    def test_prediction_coverage_matches_measured_errors(self):
        inputs = json.loads((ROOT/"results/engineering_benchmark_emt_inputs.json").read_text())
        folder = ROOT/"results/engineering_benchmark_emt"
        result = json.loads((folder/"results.json").read_text())
        receipt = json.loads((folder/"matched_audit.json").read_text())
        missed = []
        for ref in inputs["cases"]:
            ident = ref["case_id"]
            measured = result["cases"][ident]["modes"]["switched"]["finest"]["capacitor_total_rms_A"]
            predicted = ref["parent_capacitor_model_A"]
            self.assertAlmostEqual(receipt["thermal_prediction_coverage"][ident]["error_A"], measured-predicted)
            if measured > predicted+inputs["validation_bands"]["mesh_current_A"]:
                missed.append(ident)
        self.assertEqual(receipt["thermal_prediction_failures"], missed)
        self.assertEqual(receipt["all_thermal_predictions_enclose"], not missed)


if __name__ == "__main__":
    unittest.main()
