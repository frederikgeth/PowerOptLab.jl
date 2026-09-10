# Shared LinDist3Flow test fixtures.
#
# Only fixtures used by more than one test file live here, so that every
# `lindist3flow_*_tests.jl` file can be run on its own rather than depending on
# `runtests.jl` include order to define what it uses. The guard makes repeated
# inclusion a no-op.

isdefined(@__MODULE__, :_L3F_FIXTURES_LOADED) && return
const _L3F_FIXTURES_LOADED = true

using LinearAlgebra
using JuMP
using Clarabel
using PowerOptLab

"""Single-phase source-line-load feeder, optionally with a PV generator or an
explicit neutral conductor."""
function _l3f_case(; generator=false, explicit_neutral=false)
    terminals = explicit_neutral ? ["a", "n"] : ["a"]
    source_vm = explicit_neutral ? [230.0, 0.0] : [230.0]
    source_va = zeros(length(source_vm))
    load_map = copy(terminals)
    load_p = [10_000.0]
    load_q = [2_000.0]
    lc = Dict{String,Any}(
        "R_series_1_1" => 0.2,
        "X_series_1_1" => 0.1,
    )
    if explicit_neutral
        merge!(lc, Dict{String,Any}(
            "R_series_1_2" => 0.02, "X_series_1_2" => 0.01,
            "R_series_2_2" => 0.1, "X_series_2_2" => 0.05,
        ))
    end
    net = Dict{String,Any}(
        "terminal_conventions" => Dict{String,Any}(
            "phase" => ["a"], "neutral" => explicit_neutral ? ["n"] : String[]),
        "bus" => Dict(
            "source" => Dict{String,Any}(
                "terminal_names" => copy(terminals),
                "perfectly_grounded_terminals" => explicit_neutral ? ["n"] : String[]),
            "load" => Dict{String,Any}(
                "terminal_names" => copy(terminals),
                "perfectly_grounded_terminals" => explicit_neutral ? ["n"] : String[],
                "v_min" => explicit_neutral ? [180.0, 0.0] : [180.0],
                "v_max" => explicit_neutral ? [250.0, 0.0] : [250.0])),
        "linecode" => Dict("lc" => lc),
        "line" => Dict("line" => Dict{String,Any}(
            "bus_from" => "source", "bus_to" => "load",
            "terminal_map_from" => copy(terminals),
            "terminal_map_to" => copy(terminals), "linecode" => "lc")),
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => "source", "terminal_map" => copy(terminals),
            "configuration" => "SINGLE_PHASE",
            "v_magnitude" => source_vm, "v_angle" => source_va,
            "cost" => [1.0])),
        "load" => Dict("load" => Dict{String,Any}(
            "bus" => "load", "terminal_map" => load_map,
            "configuration" => "SINGLE_PHASE", "model" => "constant_power",
            "p_nom" => load_p, "q_nom" => load_q)),
    )
    if generator
        net["generator"] = Dict("pv" => Dict{String,Any}(
            "bus" => "load", "terminal_map" => ["a"],
            "configuration" => "SINGLE_PHASE",
            "p_min" => [15_000.0], "p_max" => [15_000.0],
            "q_min" => [0.0], "q_max" => [0.0], "cost" => [0.0]))
    end
    net
end

const _L3F_DY_VPN_HV = 11_000.0 / sqrt(3)
const _L3F_DY_VPN_LV = 400.0 / sqrt(3)
const _L3F_DY_P = [30_000.0, 45_000.0, 20_000.0]
const _L3F_DY_Q = [10_000.0, 15_000.0, 5_000.0]
const _L3F_DY_D = [1.0 -1.0 0.0; 0.0 1.0 -1.0; -1.0 0.0 1.0]

_l3f_dy_opts(; kwargs...) = L3FOptions(; validate_nonlinear=false,
                                       objective=:feasibility, kwargs...)
_l3f_dy_solve(net; kwargs...) = solve_l3f_opf(net, Clarabel.Optimizer;
    options=_l3f_dy_opts(; kwargs...), solver_options=("verbose" => false,))
_l3f_dy_codes(r) = [f.code for f in r.findings]
_l3f_dy_has(r, code) = any(==(code), _l3f_dy_codes(r))

"""
Two-bus Yd/Dy bank. `source_on_delta` chooses which winding faces the source,
which is the property the formulation actually cares about.
"""
function _l3f_dy_case(subtype; source_on_delta=true, line=false, extra=Dict{String,Any}())
    delta_bus, wye_bus = "d", "y"
    bf, bt = subtype == "delta_wye" ? (delta_bus, wye_bus) : (wye_bus, delta_bus)
    vnf = subtype == "delta_wye" ? _L3F_DY_VPN_HV : _L3F_DY_VPN_LV
    vnt = subtype == "delta_wye" ? _L3F_DY_VPN_LV : _L3F_DY_VPN_HV
    source_bus = source_on_delta ? delta_bus : wye_bus
    source_v = source_on_delta ? _L3F_DY_VPN_HV : _L3F_DY_VPN_LV
    load_bus = source_on_delta ? wye_bus : delta_bus

    buses = Dict{String,Any}(
        delta_bus => Dict{String,Any}("terminal_names" => ["a","b","c"]),
        wye_bus => Dict{String,Any}("terminal_names" => ["a","b","c"]))
    lines = Dict{String,Any}()
    if line
        buses["end"] = Dict{String,Any}("terminal_names" => ["a","b","c"])
        lines["l1"] = Dict{String,Any}(
            "bus_from" => load_bus, "bus_to" => "end",
            "terminal_map_from" => ["a","b","c"], "terminal_map_to" => ["a","b","c"],
            "linecode" => "lc")
        load_bus = "end"
    end
    net = Dict{String,Any}(
        "bus" => buses,
        "linecode" => Dict("lc" => Dict{String,Any}(
            "R_series_1_1" => 0.05, "X_series_1_1" => 0.02,
            "R_series_2_2" => 0.05, "X_series_2_2" => 0.02,
            "R_series_3_3" => 0.05, "X_series_3_3" => 0.02,
            "R_series_1_2" => 0.01, "X_series_1_2" => 0.005,
            "R_series_2_3" => 0.01, "X_series_2_3" => 0.005,
            "R_series_1_3" => 0.01, "X_series_1_3" => 0.005)),
        "line" => lines,
        "transformer" => Dict(subtype => Dict("t" => merge(Dict{String,Any}(
            "bus_from" => bf, "bus_to" => bt,
            "terminal_map_from" => ["a","b","c"], "terminal_map_to" => ["a","b","c"],
            "v_nom_from" => vnf, "v_nom_to" => vnt, "s_rating" => 5.0e5), extra))),
        "voltage_source" => Dict("v" => Dict{String,Any}(
            "bus" => source_bus, "terminal_map" => ["a","b","c"],
            "configuration" => "WYE", "v_magnitude" => fill(source_v, 3),
            "v_angle" => [0.0, -2pi/3, 2pi/3])),
        "load" => Dict("l" => Dict{String,Any}(
            "bus" => load_bus, "terminal_map" => ["a","b","c"],
            "configuration" => "WYE", "model" => "constant_power",
            "p_nom" => copy(_L3F_DY_P), "q_nom" => copy(_L3F_DY_Q))))
    isempty(lines) && delete!(net, "line")
    net
end

"""Single-phase feeder used by the restricted-control and permissive suites."""
function _l3f_controls_case(; generator=Dict{String,Any}(), ibr=Dict{String,Any}())
    net = Dict{String,Any}(
        "terminal_conventions" => Dict("phase" => ["a"], "neutral" => String[]),
        "bus" => Dict(
            "source" => Dict{String,Any}("terminal_names" => ["a"]),
            "load" => Dict{String,Any}("terminal_names" => ["a"],
                "v_min" => [180.0], "v_max" => [260.0])),
        "linecode" => Dict("lc" => Dict{String,Any}(
            "R_series_1_1" => 0.2, "X_series_1_1" => 0.1)),
        "line" => Dict("line" => Dict{String,Any}(
            "bus_from" => "source", "bus_to" => "load",
            "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
            "linecode" => "lc")),
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => "source", "terminal_map" => ["a"], "configuration" => "WYE",
            "v_magnitude" => [230.0], "v_angle" => [0.0], "cost" => [1.0])),
        "load" => Dict("load" => Dict{String,Any}(
            "bus" => "load", "terminal_map" => ["a"], "configuration" => "WYE",
            "model" => "constant_power", "p_nom" => [10_000.0], "q_nom" => [2_000.0])))
    isempty(generator) || (net["generator"] = Dict("g" => generator))
    isempty(ibr) || (net["ibr"] = Dict("pv" => ibr))
    net
end
