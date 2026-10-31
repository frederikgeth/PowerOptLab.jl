"""Publication regressions for coordinated CSV metrics and stale audit flags."""
import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import plot_coordinated as study


class CoordinatedCSVAcceptance(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = self.root / "results/coordinated"
        shutil.copytree(study.ROOT / "results/coordinated", self.path)
        (self.root / "networks").mkdir()
        shutil.copy(study.ROOT / "networks/LV3_55bus.bmopf.json", self.root / "networks")

    def audit(self):
        with patch.object(study, "ROOT", self.root):
            return study.audit_records()

    def write_rows(self, name, records):
        with (self.path / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    def assert_withheld(self):
        with self.assertRaises(AssertionError):
            self.audit()
        report = json.loads((self.path / "record_audit.json").read_text())
        self.assertFalse(report["all_planned_accepted"])

    def test_accepted_campaign(self):
        env, cases, devices = self.audit()
        self.assertEqual((len(cases), len(devices)), (24, 288))

    def test_changed_metrics_clear_previous_success(self):
        self.audit()
        records = study.rows(self.path / "cases.csv")
        records[0]["p_poc_W"] = str(float(records[0]["p_poc_W"])+1000.)
        self.write_rows("cases.csv", records)
        self.assert_withheld()

    def test_missing_device_is_rejected(self):
        records = study.rows(self.path / "devices.csv")
        self.write_rows("devices.csv", records[1:])
        self.assert_withheld()

    def test_false_inequality_success_is_rejected(self):
        records = study.rows(self.path / "network_checks.csv")
        records[0].update(value="1000", upper="0", passed="true")
        self.write_rows("network_checks.csv", records)
        self.assert_withheld()


if __name__ == "__main__":
    unittest.main()
