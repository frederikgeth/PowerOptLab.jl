import csv
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

from plot_thermal_energy import audit_records
SOURCE = Path(sys.argv.pop(1))


class ThermalEnergyCSVTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.folder = Path(self.tmp.name) / "data"
        shutil.copytree(SOURCE, self.folder)

    def tearDown(self):
        self.tmp.cleanup()

    def change(self, file, function):
        path = self.folder / file
        with path.open() as stream:
            rows = list(csv.DictReader(stream))
        fields = list(rows[0])
        function(rows)
        with path.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)

    def test_accepts_audited_mixed_records(self):
        self.assertEqual(len(audit_records(self.folder)[0]), 8)

    def test_rejects_changed_voltage(self):
        self.change("cases.csv", lambda r: r[0].update(vmax_sensed_V="999"))
        with self.assertRaises(ValueError): audit_records(self.folder)

    def test_rejects_nan(self):
        self.change("cases.csv", lambda r: r[0].update(p_poc_W="nan"))
        with self.assertRaises(ValueError): audit_records(self.folder)

    def test_rejects_missing_device(self):
        self.change("devices.csv", lambda r: r.pop())
        with self.assertRaises(ValueError): audit_records(self.folder)

    def test_rejects_wrong_controller(self):
        self.change("devices.csv", lambda r: r[0].update(controller="wrong"))
        with self.assertRaises(ValueError): audit_records(self.folder)

    def test_rejects_wrong_model(self):
        self.change("devices.csv", lambda r: r[0].update(pwm_strategy="CENTERED"))
        with self.assertRaises(ValueError): audit_records(self.folder)


if __name__ == "__main__": unittest.main()
