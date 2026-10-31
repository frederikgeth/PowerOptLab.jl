"""Reject altered representative selection, source availability and EMT hardware."""
import copy
import json
import os
from pathlib import Path
import tempfile
import unittest

from audit_engineering_benchmark_emt import validate_references

ROOT = Path(__file__).resolve().parents[2]/"studies/australian_pv/results"
INPUT = Path(os.environ.get("ENGINEERING_EMT_INPUT", ROOT/"engineering_benchmark_emt_inputs.json"))
BENCHMARK = Path(os.environ.get("ENGINEERING_EMT_BENCHMARK", ROOT/"engineering_benchmark"))


class MatchedEMTReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = json.loads(INPUT.read_text())

    def setUp(self):
        self.env = copy.deepcopy(self.fixture)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.input = Path(self.tmp.name)/"inputs.json"

    def check(self):
        self.input.write_text(json.dumps(self.env, allow_nan=False))
        return validate_references(self.input, BENCHMARK)

    def rejected(self):
        with self.assertRaises(ValueError):
            self.check()

    def test_complete_matched_selection(self):
        result = self.check()
        self.assertEqual(result["references"], 60)
        self.assertGreater(result["checks"], 600)

    def test_missing_representative(self):
        self.env["cases"].pop()
        self.rejected()

    def test_night_source_replaced_by_fixture(self):
        row = next(c for c in self.env["cases"] if c["parent_case_id"] == "night")
        row["p_available_W"] = 5000.
        self.rejected()

    def test_current_changed(self):
        self.env["cases"][0]["current"][0][0] += .01
        self.rejected()

    def test_hardware_changed(self):
        self.env["cases"][0]["dc_capacitance_F"] *= 2
        self.rejected()

    def test_failure_reclassified_expected(self):
        self.env["cases"][0]["expect_switched_hardware_accepted"] = False
        self.rejected()

    def test_dc_guard_disabled(self):
        self.env["cases"][0]["pwm_dc_voltage_envelope"] = False
        self.rejected()

    def test_parent_audit_changed(self):
        self.env["benchmark_audit_sha256"] = "0"*64
        self.rejected()


if __name__ == "__main__":
    unittest.main()
