"""Saved-result acceptance regressions using the committed snapshot campaign."""
import contextlib
import csv
import io
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import plot_voltage_aware as study


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class SavedRecordAcceptance(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = self.root / "results/voltage_aware_snapshots"
        shutil.copytree(study.ROOT / "results/voltage_aware_snapshots", self.path)
        (self.root / "networks").mkdir()
        shutil.copy(study.ROOT / "networks/LV3_55bus.bmopf.json", self.root / "networks")

    def audit(self):
        with patch.object(study, "ROOT", self.root), contextlib.redirect_stdout(io.StringIO()):
            return study.audit("snapshots")

    def test_complete_campaign(self):
        self.assertTrue(self.audit()[-1])

    def test_missing_planned_case_withholds_comparison(self):
        path = self.path / "environment.json"
        env = json.loads(path.read_text())
        env["planned_cases"].append(dict(env["planned_cases"][0], case_id="not_completed"))
        path.write_text(json.dumps(env))
        path = self.path / "summary.json"
        summary = json.loads(path.read_text())
        summary.update(planned=42, all_planned_accepted=False)
        path.write_text(json.dumps(summary))
        self.assertFalse(self.audit()[-1])
        report = json.loads((self.path / "record_audit.json").read_text())
        self.assertEqual(report["missing_case_ids"], ["not_completed"])

    def test_failed_nonfinite_candidate_is_retained(self):
        path = self.path / "cases.csv"
        rows = study.read_csv(path)
        cid = rows[0]["case_id"]
        rows[0].update(accepted="false", reason="ALLOCATION_FAILED", p_poc_W="NaN", q_poc_var="NaN")
        write_csv(path, rows)
        path = self.path / "devices.csv"
        rows = study.read_csv(path)
        for row in rows:
            if row["case_id"] == cid:
                row.update(accepted="false", p_poc_W="NaN", q_poc_var="NaN")
        write_csv(path, rows)
        path = self.path / "summary.json"
        summary = json.loads(path.read_text())
        summary.update(accepted=40, all_planned_accepted=False, unresolved_case_ids=[cid])
        path.write_text(json.dumps(summary))
        self.assertFalse(self.audit()[-1])
        report = json.loads((self.path / "record_audit.json").read_text())
        self.assertEqual(report["unresolved_case_ids"], [cid])

    def test_false_physical_success_is_rejected(self):
        path = self.path / "physical_checks.csv"
        rows = study.read_csv(path)
        rows[0].update(value="1000000000", upper="0", passed="true")
        write_csv(path, rows)
        with self.assertRaises(AssertionError):
            self.audit()


if __name__ == "__main__":
    unittest.main()
