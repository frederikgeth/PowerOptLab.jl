# Shared synthetic case family for the DOE uncertainty tutorials.

using PowerOptLab
using BMOPFTools: parse_bmopf
using Dates
using Random

# Independent synthetic test draws; this distribution is a teaching assumption,
# not a calibration of the deterministic stress fixture below.
function doe_iid_test_scenarios(; count=32, seed=7319)
    rng = MersenneTwister(seed)
    return DOEScenarioSet([
        DOEScenario(id="iid-$i", role=:test,
            network=doe_uncertainty_feeder(100 + 5900rand(rng), 100 + 5900rand(rng)),
            source="independent Uniform(100,6000) W demands per bus",
            generation_method=:synthetic_iid_uniform, seed=seed)
        for i in 1:count]; dataset_id="synthetic-iid-$seed-$count")
end

function doe_uncertainty_feeder(p1_W::Real, p2_W::Real)
    return parse_bmopf("""
    {"bus":{
        "src":{"terminal_names":["1","n"],
               "perfectly_grounded_terminals":["n"]},
        "bus1":{"terminal_names":["1","n"],
                "perfectly_grounded_terminals":["n"],
                "v_min":[216.0],"v_max":[245.0]},
        "bus2":{"terminal_names":["1","n"],
                "perfectly_grounded_terminals":["n"],
                "v_min":[216.0],"v_max":[245.0]}},
     "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],
         "v_magnitude":[230.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.4}},
     "line":{
        "l1":{"bus_from":"src","bus_to":"bus1",
              "terminal_map_from":["1"],"terminal_map_to":["1"],
              "linecode":"lc","length":1.0},
        "l2":{"bus_from":"bus1","bus_to":"bus2",
              "terminal_map_from":["1"],"terminal_map_to":["1"],
              "linecode":"lc","length":1.0}},
     "load":{
        "d1":{"bus":"bus1","terminal_map":["1","n"],
              "configuration":"SINGLE_PHASE","p_nom":[$p1_W],"q_nom":[0.0]},
        "d2":{"bus":"bus2","terminal_map":["1","n"],
              "configuration":"SINGLE_PHASE","p_nom":[$p2_W],"q_nom":[0.0]}}}
    """; from_string=true)
end

function doe_uncertainty_tutorial_case()
    scenarios = DOEScenarioSet([
        DOEScenario(
            id="calibration-high-load",
            network=doe_uncertainty_feeder(5000.0, 5000.0),
            role=:calibration,
            weight=1.0,
            source="synthetic tutorial fixture",
            generation_method=:deterministic_fixture,
            seed=11,
            timestamp=DateTime(2026, 1, 1, 12),
            metadata=Dict(
                "site_id" => "synthetic-feeder",
                "aggregate_demand_kw" => 10.0,
                "uncertainty_sample_id" => "case-high")),
        DOEScenario(
            id="test-low-load",
            network=doe_uncertainty_feeder(200.0, 200.0),
            role=:stress,
            weight=0.7,
            source="synthetic tutorial fixture",
            generation_method=:held_out_fixture,
            seed=12,
            timestamp=DateTime(2026, 1, 15, 12),
            metadata=Dict(
                "site_id" => "synthetic-feeder",
                "aggregate_demand_kw" => 0.4,
                "uncertainty_sample_id" => "case-low")),
        DOEScenario(
            id="test-high-load",
            network=doe_uncertainty_feeder(5000.0, 5000.0),
            role=:test,
            weight=0.3,
            source="synthetic tutorial fixture",
            generation_method=:held_out_fixture,
            seed=13,
            timestamp=DateTime(2026, 1, 15, 12),
            metadata=Dict(
                "site_id" => "synthetic-feeder",
                "aggregate_demand_kw" => 10.0,
                "uncertainty_sample_id" => "case-high")),
    ]; dataset_id="doe-uncertainty-tutorial-v1",
       metadata=Dict("license" => "CC0-1.0", "synthetic" => true))

    connection_points = [
        ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
        ConnectionPoint(id="d2", bus="bus2", export_max=10e3),
    ]
    return (scenarios=scenarios,
            connection_points=connection_points,
            metadata=Dict(
                "case_id" => "doe-uncertainty-tutorial-v1",
                "claim_limit" =>
                    "finite synthetic scenarios do not establish probability coverage"))
end
