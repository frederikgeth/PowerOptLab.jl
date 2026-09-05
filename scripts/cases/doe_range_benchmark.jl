# Reproducible synthetic benchmark for DOE claim semantics.
#
# This case is intentionally small and uncalibrated. It demonstrates that a
# simultaneous bound-point allocation and a finite corner-tested allocation
# support different claims; it does not certify continuous AC box containment.

using PowerOptLab
using BMOPFTools: parse_bmopf

function _doe_range_benchmark_network()
    return parse_bmopf(raw"""
    {
      "bus": {
        "src": {
          "terminal_names": ["1", "2", "3", "n"],
          "perfectly_grounded_terminals": ["n"]
        },
        "b1": {
          "terminal_names": ["1", "2", "3", "n"],
          "perfectly_grounded_terminals": ["n"],
          "v_min": [200.0, 200.0, 200.0],
          "v_max": [260.0, 260.0, 260.0],
          "vneg_max": 1.0
        }
      },
      "voltage_source": {
        "vs": {
          "bus": "src",
          "terminal_map": ["1", "2", "3"],
          "v_magnitude": [230.0, 230.0, 230.0],
          "v_angle": [0.0, -2.0943951024, 2.0943951024]
        }
      },
      "linecode": {
        "lc": {
          "R_series_1_1": 0.4,
          "R_series_2_2": 0.4,
          "R_series_3_3": 0.4,
          "R_series_4_4": 0.4
        }
      },
      "line": {
        "l1": {
          "bus_from": "src",
          "bus_to": "b1",
          "terminal_map_from": ["1", "2", "3", "n"],
          "terminal_map_to": ["1", "2", "3", "n"],
          "linecode": "lc",
          "length": 1.0
        }
      }
    }
    """; from_string=true)
end

function doe_benchmark_case()
    connection_points = [
        ConnectionPoint(
            id="phase_$phase",
            bus="b1",
            phase_terminals=[string(phase)],
            neutral="n",
            export_max=20e3,
        ) for phase in 1:3
    ]
    return (
        nets=_doe_range_benchmark_network(),
        connection_points=connection_points,
        methods=[
            "bound_point_ideal_recourse" => (
                direction=:export,
                fairness=:equal,
                security=:bound_point,
                control_policy=PerfectRecourse(),
            ),
            "corners_issue_plus_local_laws" => (
                direction=:export,
                fairness=:equal,
                security=:corners,
                control_policy=IssuePlusLocalLaws(),
            ),
        ],
        seeds=Dict("case_fixture" => 20260905),
        metadata=Dict{String,Any}(
            "case_id" => "synthetic-doe-range-claim-v1",
            "dataset" => "self-contained synthetic four-wire feeder",
            "data_license" => "CC0-1.0",
            "scientific_claim" =>
                "bound-point and finite corner security support different claims",
            "claim_limit" =>
                "neither method certifies continuous nonlinear AC box containment",
        ),
    )
end
