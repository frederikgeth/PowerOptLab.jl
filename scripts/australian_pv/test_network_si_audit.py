"""Use a restored current-benchmark phasor point as an input fixture, not as a new run."""
import json
from pathlib import Path
import unittest

from network_si_audit import audit_network

ROOT = Path(__file__).resolve().parents[2] / "studies/australian_pv"


class NetworkSIAuditTests(unittest.TestCase):
    def setUp(self):
        self.net = json.loads((ROOT / "networks/LV3_55bus.bmopf.json").read_text())
        self.point = json.loads((ROOT / "results/engineering_benchmark/points/nominal_opf_vv_released.json").read_text())
        self.case = dict(source_pu=1., source_vuf=.01, source_negative_angle_rad=0., phase_shares=[.6, .25, .15], load_scale=1.)

    def audit(self):
        return audit_network(self.net, self.point, self.case)

    def test_circuit_and_limits(self):
        r = self.audit()
        self.assertTrue(r["accepted"])
        self.assertGreater(len(r["checks"]), 800)
        self.assertLess(abs(r["active_power_residual_W"]), .001)
        json.dumps(r, allow_nan=False)

    def test_customer_neutral_corruption(self):
        bus = next(iter(self.net["load"].values()))["bus"]
        self.point["bus"][bus]["n"]["vr"] += .05
        self.assertFalse(self.audit()["accepted"])

    def test_pv_current_corruption(self):
        next(iter(self.point["devices"].values()))["grid_current"][0][0] += .02
        self.assertFalse(self.audit()["accepted"])

    def test_conductor_limit_is_checked(self):
        baseline = self.audit()
        ident = max(baseline["line_currents"], key=lambda k: max(abs(complex(*z)) for z in baseline["line_currents"][k]["from"]))
        self.net["line"][ident]["i_max"] = [.001]*4
        r = self.audit()
        self.assertFalse(r["accepted"])
        self.assertTrue(any(c["name"].startswith("line/") and not c["passed"] for c in r["checks"]))

    def test_wrong_source_phase(self):
        self.case["source_negative_angle_rad"] = .1
        self.assertFalse(self.audit()["accepted"])

    def test_missing_pv(self):
        self.point["devices"].pop(next(iter(self.point["devices"])))
        with self.assertRaises(ValueError):
            self.audit()

    def test_transformer_tap_corruption(self):
        next(iter(self.net["transformer"]["delta_wye"].values()))["tap"] *= 1.01
        r = self.audit()
        self.assertTrue(any(c["name"].startswith("transformer/") and not c["passed"] for c in r["checks"]))

    def test_switch_voltage_corruption(self):
        switch = next(iter(self.net["switch"].values()))
        self.point["bus"][switch["bus_to"]]["a"]["vr"] += .02
        r = self.audit()
        self.assertTrue(any(c["name"].startswith("switch/") and not c["passed"] for c in r["checks"]))


if __name__ == "__main__":
    unittest.main()
