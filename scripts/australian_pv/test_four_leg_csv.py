import csv
from pathlib import Path
import shutil
import tempfile
import unittest
import sys

from plot_four_leg import load_evidence

SOURCE = Path(sys.argv.pop(1)) if len(sys.argv) > 1 else Path(__file__).resolve().parents[2] / "studies/australian_pv/results/four_leg"


class FourLegCSVTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.folder = Path(self.tmp.name) / "data"
        shutil.copytree(SOURCE, self.folder)

    def tearDown(self):
        self.tmp.cleanup()

    def rewrite(self, file, change):
        path = self.folder / file
        with path.open() as stream:
            rows = list(csv.DictReader(stream))
        fields = list(rows[0])
        change(rows)
        with path.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    def test_accepts_audited_records(self):
        self.assertEqual(len(load_evidence(self.folder)[0]), 11)

    def test_rejects_changed_case_voltage(self):
        self.rewrite("cases.csv", lambda rs: rs[0].update(vmax_sensed_V="999"))
        with self.assertRaises(ValueError):
            load_evidence(self.folder)

    def test_rejects_nan_metric(self):
        self.rewrite("cases.csv", lambda rs: rs[0].update(p_poc_W="nan"))
        with self.assertRaises(ValueError):
            load_evidence(self.folder)

    def test_rejects_missing_device(self):
        self.rewrite("devices.csv", lambda rs: rs.pop(0))
        with self.assertRaises(ValueError):
            load_evidence(self.folder)

    def test_rejects_false_acceptance(self):
        self.rewrite("devices.csv", lambda rs: rs[0].update(accepted="false"))
        with self.assertRaises(ValueError):
            load_evidence(self.folder)


if __name__ == "__main__":
    unittest.main()
