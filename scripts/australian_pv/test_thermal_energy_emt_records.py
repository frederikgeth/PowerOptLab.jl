import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

from audit_emt import audit
INPUT = Path(sys.argv.pop(1)); SOURCE = Path(sys.argv.pop(1))


class ThermalEnergyEMTTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.folder = Path(self.tmp.name) / "run"; shutil.copytree(SOURCE, self.folder)
        self.input = Path(self.tmp.name) / "input.json"; shutil.copyfile(INPUT, self.input)

    def tearDown(self): self.tmp.cleanup()

    def change(self, function):
        p = self.folder / "results.json"; r = json.loads(p.read_text()); function(r); p.write_text(json.dumps(r))

    def test_all_allocated_points_and_negative_control(self):
        report = audit(self.input, self.folder)
        self.assertTrue(report["accepted"])
        r = json.loads((self.folder / "results.json").read_text())
        allocated = [c for i, c in r["cases"].items() if i != "nonlinear_boundary" and not i.startswith("calibration_")]
        self.assertEqual(len(allocated), 27)
        self.assertTrue(all(c["modes"]["switched"]["finest"]["all_engineering_checks_passed"] for c in allocated))
        self.assertFalse(r["cases"]["nonlinear_boundary"]["modes"]["switched"]["finest"]["all_engineering_checks_passed"])
        self.assertTrue(r["cases"]["nonlinear_boundary"]["modes"]["averaged"]["finest"]["all_engineering_checks_passed"])

    def test_missing_mesh(self):
        self.change(lambda r: r["metrics"].pop())
        with self.assertRaises(ValueError): audit(self.input, self.folder)

    def test_corrupt_energy(self):
        self.change(lambda r: r["metrics"][0].update(window_net_energy_J=999))
        with self.assertRaises(ValueError): audit(self.input, self.folder)

    def test_false_hardware_success(self):
        self.change(lambda r: r.update(all_switched_hardware_accepted=True))
        with self.assertRaises(ValueError): audit(self.input, self.folder)

    def test_corrupt_waveform(self):
        p = self.folder / "energy_balanced_switched_switching.csv"; p.write_text(p.read_text()+"0\n")
        with self.assertRaises(ValueError): audit(self.input, self.folder)


if __name__ == "__main__": unittest.main()
