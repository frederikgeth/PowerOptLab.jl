#!/usr/bin/env julia

# Reproducible native-static post-OPF operability study.
#
# Usage:
#   julia --project=. scripts/run_post_opf_operability_study.jl
#
# The study intentionally evaluates one solved snapshot at a time. It then
# runs two finite loading directions for constant-power and ZIP load models.
# The output is finite study evidence, not a contingency or operating-envelope
# guarantee.

using PowerOptLab
using BMOPFTools: parse_bmopf, solve_pf, SIUnitsScaling

function study_net(model::AbstractString)
    net = parse_bmopf("""
    {"bus":{
      "source":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"]},
      "loadbus":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"],
                  "v_min":[190.0,190.0,190.0],"v_max":[260.0,260.0,260.0]}},
     "voltage_source":{"vs":{"bus":"source","terminal_map":["a","b","c"],
         "v_magnitude":[230.0,230.0,230.0],"v_angle":[0.0,-2.0943951023931953,2.0943951023931953]}},
     "linecode":{"lc":{"R_series_1_1":0.08,"R_series_2_2":0.08,"R_series_3_3":0.08,
                           "R_series_4_4":0.08,"X_series_1_1":0.12,"X_series_2_2":0.12,
                           "X_series_3_3":0.12,"X_series_4_4":0.12}},
     "line":{"l1":{"bus_from":"source","bus_to":"loadbus",
         "terminal_map_from":["a","b","c","n"],"terminal_map_to":["a","b","c","n"],
         "linecode":"lc","length":1.0}},
     "load":{"ld":{"bus":"loadbus","terminal_map":["a","b","c","n"],
         "configuration":"WYE","model":"$model",
         "p_nom":[12000.0,12000.0,12000.0],"q_nom":[3000.0,3000.0,3000.0]}}
    }
    """; from_string=true)
    load = net["load"]["ld"]
    if lowercase(model) == "zip"
        load["v_nom"] = [230.0, 230.0, 230.0]
        load["alpha_z"] = [0.25, 0.25, 0.25]
        load["alpha_i"] = [0.25, 0.25, 0.25]
        load["alpha_p"] = [0.50, 0.50, 0.50]
        load["beta_z"] = [0.25, 0.25, 0.25]
        load["beta_i"] = [0.25, 0.25, 0.25]
        load["beta_p"] = [0.50, 0.50, 0.50]
    end
    net
end

const STUDY_SPEC = OperabilitySpec(
    scaling_policy=SIUnitsScaling(),
    voltage_min=190.0,
    voltage_max=260.0,
    vuf_max=0.02,
    compute_fixed_point_certificate=true)

const STUDY_DIRECTIONS = [
    OperabilityStressDirection(:uniform),
    OperabilityStressDirection(:reactive_removed; p_scale=1.0, q_scale=0.0),
]

function run_campaign(model::AbstractString)
    net = study_net(model)
    solution = solve_pf(net; per_unit=false)
    report = check_opf_operability(net, solution; spec=STUDY_SPEC)
    rows = operability_stress_rows(net, solution;
        spec=STUDY_SPEC,
        directions=STUDY_DIRECTIONS,
        lambdas=[0.0, 0.5, 1.0],
        solve=network -> solve_pf(network; per_unit=false))
    report, operability_snapshot_row(report; snapshot_id=String(model)), rows
end

campaigns = Dict{String,Vector{NamedTuple}}()
reports = Dict{String,OperabilityResult}()
snapshot_rows = NamedTuple[]
for model in ("constant_power", "zip")
    report, row, rows = run_campaign(model)
    reports[model] = report
    campaigns[model] = rows
    push!(snapshot_rows, row)
end

ensemble = operability_stress_ensemble_rows(campaigns)

println("Post-OPF operability study")
println("  scope: single_snapshot_static_ybus")
println("  models: ", join(sort!(collect(keys(campaigns))), ", "))
println("  snapshot rows:")
for row in snapshot_rows
    println("    ", row.snapshot_id,
        " status=", row.status,
        " Vmin=", round(row.minimum_terminal_voltage; digits=3),
        " VUFmax=", row.maximum_vuf,
        " certificate=", row.fixed_point_certificate_status,
        " dP/dV(high,near,low)=",
        (row.high_side_indicator_count, row.near_nose_indicator_count,
         row.low_side_indicator_count))
end
println("  finite stress summaries:")
for row in ensemble
    println("    model=", row.model,
        " direction=", row.direction,
        " status=", row.status,
        " rows=", row.row_count,
        " boundary=", row.boundary_status,
        " min_certificate_margin=", row.minimum_condition_margin)
end

@assert all(report.status == :pass for report in values(reports))
@assert all(row.status == :pass for row in ensemble)
@assert all(row.boundary_status == :not_observed for row in ensemble)
println("  result: PASS (finite declared study only)")
