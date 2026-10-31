import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

from audit_emt import audit

INPUT = Path(sys.argv.pop(1))
SOURCE = Path(sys.argv.pop(1))


class EMTRecordTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.folder = Path(self.tmp.name) / "run"
        shutil.copytree(SOURCE, self.folder)
        self.input = Path(self.tmp.name) / "input.json"
        shutil.copyfile(INPUT, self.input)

    def tearDown(self):
        self.tmp.cleanup()

    def change(self, function):
        path = self.folder / "results.json"
        result = json.loads(path.read_text())
        function(result)
        path.write_text(json.dumps(result))

    def test_accepts_agreement_with_explicit_hardware_failures(self):
        report = audit(self.input, self.folder)
        self.assertTrue(report["accepted"])
        self.assertTrue(all("capacitor_total_rms" in c["switched"]["failed_engineering_checks"] for c in report["cases"].values()))

    def test_rejects_wrong_energy(self):
        self.change(lambda r: r["metrics"][0]["current_total_rms_A"].__setitem__(0, 99))
        with self.assertRaises(ValueError):
            audit(self.input, self.folder)

    def test_rejects_missing_mesh(self):
        self.change(lambda r: r["metrics"].pop())
        with self.assertRaises(ValueError):
            audit(self.input, self.folder)

    def test_rejects_false_hardware_acceptance(self):
        self.change(lambda r: r.update(all_switched_hardware_accepted=True))
        with self.assertRaises(ValueError):
            audit(self.input, self.folder)

    def test_rejects_changed_waveform(self):
        path = self.folder / "balanced_L_switched_switching.csv"
        path.write_text(path.read_text()+"0,0,0,0,0,0,0\n")
        with self.assertRaises(ValueError):
            audit(self.input, self.folder)

    def test_rejects_stale_circuit(self):
        self.change(lambda r: r.update(circuit_sha256="stale"))
        with self.assertRaises(ValueError):
            audit(self.input, self.folder)


if __name__ == "__main__":
    unittest.main()
