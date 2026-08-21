#!/usr/bin/env julia

# Reproducible native-static post-OPF operability study.
#
# Usage:
#   julia --project=. scripts/run_post_opf_operability_study.jl
#
# The study intentionally evaluates one solved snapshot at a time. It then
# runs two finite loading directions for constant-power, ZIP, unbalanced DELTA,
# floating-neutral, and transformer-separated load cases.
# The output is finite study evidence, not a contingency or operating-envelope
# guarantee.

using PowerOptLab
using BMOPFTools: parse_bmopf, solve_pf, SIUnitsScaling

function study_net(model::AbstractString)
    if model == "single_phase_transformer"
        return parse_bmopf("""
        {"bus":{
          "hv":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
          "lv":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],
                 "v_min":[100.0],"v_max":[140.0]}},
         "voltage_source":{"vs":{"bus":"hv","terminal_map":["1"],
             "v_magnitude":[230.0],"v_angle":[0.0]}},
         "transformer":{"single_phase":{"t1":{
             "bus_from":"hv","bus_to":"lv",
             "terminal_map_from":["1","n"],"terminal_map_to":["1","n"],
             "v_nom_from":230.0,"v_nom_to":120.0,"s_rating":100000.0,
             "r_series_from":0.005,"x_series_from":0.02}}},
         "load":{"ld":{"bus":"lv","terminal_map":["1","n"],
             "configuration":"WYE","model":"constant_power",
             "p_nom":[10000.0],"q_nom":[2000.0]}}
        }
        """; from_string=true)
    end
    is_delta = model == "unbalanced_delta"
    is_floating = model == "floating_neutral"
    load_model = (is_delta || is_floating) ? "constant_power" : model
    configuration = is_delta ? "DELTA" : "WYE"
    terminal_map = is_delta ? "[\"a\",\"b\",\"c\"]" : is_floating ?
        "[\"a\",\"n\"]" :
        "[\"a\",\"b\",\"c\",\"n\"]"
    p_nom = is_delta ? "[9000.0,14000.0,11000.0]" : is_floating ? "[15000.0]" :
        "[12000.0,12000.0,12000.0]"
    q_nom = is_delta ? "[2500.0,4000.0,3000.0]" : is_floating ? "[5000.0]" :
        "[3000.0,3000.0,3000.0]"
    loadbus_grounding = is_floating ? "" :
        "\"perfectly_grounded_terminals\":[\"n\"],"
    net = parse_bmopf("""
    {"bus":{
      "source":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"]},
      "loadbus":{"terminal_names":["a","b","c","n"],$loadbus_grounding
                  "v_min":[190.0,190.0,190.0],"v_max":[260.0,260.0,260.0]}},
     "voltage_source":{"vs":{"bus":"source","terminal_map":["a","b","c"],
         "v_magnitude":[230.0,230.0,230.0],"v_angle":[0.0,-2.0943951023931953,2.0943951023931953]}},
     "linecode":{"lc":{"R_series_1_1":0.08,"R_series_2_2":0.08,"R_series_3_3":0.08,
                           "R_series_4_4":0.08,"X_series_1_1":0.12,"X_series_2_2":0.12,
                           "X_series_3_3":0.12,"X_series_4_4":0.12}},
     "line":{"l1":{"bus_from":"source","bus_to":"loadbus",
         "terminal_map_from":["a","b","c","n"],"terminal_map_to":["a","b","c","n"],
         "linecode":"lc","length":1.0}},
     "load":{"ld":{"bus":"loadbus","terminal_map":$terminal_map,
         "configuration":"$configuration","model":"$load_model",
         "p_nom":$p_nom,"q_nom":$q_nom}}
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

function study_spec(model::AbstractString)
    if model == "unbalanced_delta"
        # DELTA connection records are phase-to-phase terminal quantities;
        # retain a separate declared bound rather than comparing them with
        # phase-to-neutral WYE limits.
        return OperabilitySpec(
            scaling_policy=SIUnitsScaling(),
            voltage_min=340.0,
            voltage_max=450.0,
            vuf_max=0.10,
            compute_fixed_point_certificate=true)
    elseif model == "floating_neutral"
        return OperabilitySpec(
            scaling_policy=SIUnitsScaling(),
            voltage_min=190.0,
            voltage_max=260.0,
            vuf_max=Inf,
            compute_fixed_point_certificate=true)
    elseif model == "single_phase_transformer"
        return OperabilitySpec(
            scaling_policy=SIUnitsScaling(),
            voltage_min=100.0,
            voltage_max=140.0,
            vuf_max=Inf,
            compute_fixed_point_certificate=true)
    end
    OperabilitySpec(
        scaling_policy=SIUnitsScaling(),
        voltage_min=190.0,
        voltage_max=260.0,
        vuf_max=0.02,
        compute_fixed_point_certificate=true)
end

const STUDY_DIRECTIONS = [
    OperabilityStressDirection(:uniform),
    OperabilityStressDirection(:reactive_removed; p_scale=1.0, q_scale=0.0),
]

function run_campaign(model::AbstractString)
    net = study_net(model)
    solution = solve_pf(net; per_unit=false)
    spec = study_spec(model)
    rows = operability_stress_rows(net, solution;
        spec=spec,
        directions=STUDY_DIRECTIONS,
        lambdas=[0.0, 0.5, 1.0],
        solve=network -> solve_pf(network; per_unit=false))
    report = check_opf_operability(net, solution; spec=spec)
    report, operability_snapshot_row(report; snapshot_id=String(model)), rows
end

campaigns = Dict{String,Vector{NamedTuple}}()
reports = Dict{String,OperabilityResult}()
snapshot_rows = NamedTuple[]
for model in ("constant_power", "zip", "unbalanced_delta", "floating_neutral",
              "single_phase_transformer")
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
        " local_certificate=", row.fixed_point_local_region_status,
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
delta_report = reports["unbalanced_delta"]
@assert delta_report.provenance["operability"]["model_inventory"][
    "load_configurations"] == ["DELTA"]
@assert length(delta_report.load_connections) == 3
floating_report = reports["floating_neutral"]
neutral_voltage = abs(floating_report.node_voltages[("loadbus", "n")])
@assert isfinite(neutral_voltage) && neutral_voltage > 0.0
@assert length(floating_report.load_connections) == 1
println("  floating-neutral displacement [V]: ", round(neutral_voltage; digits=6))
transformer_report = reports["single_phase_transformer"]
@assert transformer_report.provenance["operability"]["source_buses"] == ["hv"]
@assert length(transformer_report.load_connections) == 1
transformer_voltage = transformer_report.load_connections["ld/1"]["magnitude"]
@assert isfinite(transformer_voltage) && 100.0 <= transformer_voltage <= 140.0
println("  transformer LV terminal voltage [V]: ", round(transformer_voltage; digits=6))
println("  result: PASS (finite declared study only)")
