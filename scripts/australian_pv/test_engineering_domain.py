"""Corruption gates for the saved carrier-domain evidence."""
import copy
import json
import os
from pathlib import Path
import tempfile
import unittest

from audit_engineering_domain import audit, REPO

INPUT = Path(os.environ.get("ENGINEERING_DOMAIN_INPUT",
    REPO / "studies/australian_pv/results/engineering_domain_inputs.json"))


class CarrierDomainEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = json.loads(INPUT.read_text())

    def setUp(self):
        self.env = copy.deepcopy(self.fixture)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.input = Path(self.tmp.name)/"inputs.json"
        self.output = Path(self.tmp.name)/"audit.json"

    def run_audit(self):
        self.input.write_text(json.dumps(self.env, allow_nan=False))
        return audit(self.input, self.output)

    def rejected(self):
        with self.assertRaises(ValueError):
            self.run_audit()
        self.assertFalse(json.loads(self.output.read_text())["accepted"])

    def test_full_declared_domain(self):
        result = self.run_audit()
        self.assertTrue(result["accepted"])
        self.assertEqual(result["cases"], 972)
        self.assertGreater(result["checks"], 35000)
        self.assertGreater(result["engineering_limit_failures"], 0)
        self.assertLess(result["engineering_limit_failures"], 972)

    def test_missing_point(self):
        self.env["cases"].pop()
        self.rejected()

    def test_shrunk_domain(self):
        self.env["domain"]["p_W"] = [0.]
        self.env["cases"] = [r for r in self.env["cases"] if r["p_W"] == 0.]
        self.rejected()

    def test_inflated_analytical_bound(self):
        self.env["cases"][0]["phase_rms_bound_A"][0] *= 2
        self.rejected()

    def test_perturbed_reference_current(self):
        self.env["cases"][0]["converter_current"][0][0] += .01
        self.rejected()

    def test_modulation_mislabel(self):
        self.env["cases"][0]["duties"][0] += .01
        self.rejected()

    def test_masked_engineering_failure(self):
        row = next(r for r in self.env["cases"] if not r["engineering_limits_passed"])
        row["engineering_limits_passed"] = True
        self.rejected()

    def test_stale_source(self):
        self.env["source_sha256"][next(iter(self.env["source_sha256"]))] = "0"*64
        self.rejected()


if __name__ == "__main__":
    unittest.main()
