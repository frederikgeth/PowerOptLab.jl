using Test
using Clarabel
using PowerOptLab

# OpenDSS three-winding split-phase transformer with symmetric pairwise
# leakage XHL=XHT=XLT=2%. Steinmetz conversion gives 1% on each star arm.
const _L3F_CT_S = 25_000.0
const _L3F_CT_VH = 2400.0
const _L3F_CT_VL = 120.0

function _l3f_center_tap_oracle_case()
    zbase_h = _L3F_CT_VH^2 / _L3F_CT_S
    zbase_l = _L3F_CT_VL^2 / _L3F_CT_S
    Dict{String,Any}(
        "bus" => Dict(
            "hv" => Dict{String,Any}("terminal_names" => ["h"]),
            "lv" => Dict{String,Any}("terminal_names" => ["x1", "x2"])),
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => "hv", "terminal_map" => ["h"],
            "configuration" => "SINGLE_PHASE",
            "v_magnitude" => [_L3F_CT_VH], "v_angle" => [0.0])),
        "transformer" => Dict("center_tap" => Dict("ct" => Dict{String,Any}(
            "bus_from" => "hv", "bus_to" => "lv",
            "terminal_map_from" => ["h"],
            "terminal_map_to" => ["x1", "x2"],
            "v_nom_from" => _L3F_CT_VH, "v_nom_to" => _L3F_CT_VL,
            "s_rating" => _L3F_CT_S,
            "r_series_from" => 0.005 * zbase_h,
            "x_series_from" => 0.01 * zbase_h,
            "r_series_to" => 0.005 * zbase_l,
            "x_series_to" => 0.01 * zbase_l))),
        "load" => Dict(
            "leg1" => Dict{String,Any}(
                "bus" => "lv", "terminal_map" => ["x1"],
                "configuration" => "SINGLE_PHASE", "model" => "constant_power",
                "p_nom" => [6_000.0], "q_nom" => [1_000.0]),
            "leg2" => Dict{String,Any}(
                "bus" => "lv", "terminal_map" => ["x2"],
                "configuration" => "SINGLE_PHASE", "model" => "constant_power",
                "p_nom" => [2_000.0], "q_nom" => [250.0])))
end

function _l3f_center_tap_opendss_voltages()
    dss(cmd) = OpenDSSDirect.dss(cmd)
    dss("Clear")
    dss("New Circuit.l3f_ct phases=1 bus1=hv.1.0 basekv=2.4 pu=1.0 angle=0")
    dss("Edit Vsource.source phases=1 bus1=hv.1.0 basekv=2.4 pu=1.0 angle=0")
    dss("New Transformer.ct phases=1 windings=3 xhl=2 xht=2 xlt=2")
    dss("~ wdg=1 bus=hv.1.0 conn=wye kv=2.4 kva=25 %r=0.5")
    dss("~ wdg=2 bus=lv.1.0 conn=wye kv=0.12 kva=25 %r=0.5")
    # Reversing the terminal order makes winding 3 series-aiding about ground:
    # its winding voltage is V(ground)-V(lv.2), matching BMOPF center_tap.
    dss("~ wdg=3 bus=lv.0.2 conn=wye kv=0.12 kva=25 %r=0.5")
    dss("New Load.leg1 phases=1 bus1=lv.1.0 conn=wye kv=0.12 kw=6 kvar=1 model=1")
    dss("New Load.leg2 phases=1 bus1=lv.2.0 conn=wye kv=0.12 kw=2 kvar=0.25 model=1")
    dss("Set maxiterations=200 tolerance=1e-11 controlmode=off")
    dss("Solve mode=snapshot")
    @test OpenDSSDirect.Solution.Converged()
    names = lowercase.(OpenDSSDirect.Circuit.AllNodeNames())
    values = ComplexF64.(OpenDSSDirect.Circuit.AllBusVolts())
    Dict(name => value for (name, value) in zip(names, values))
end

@testset "LinDist3Flow center-tap OpenDSS regression" begin
    net = _l3f_center_tap_oracle_case()
    report = check_l3f_applicability(net;
        options=L3FOptions(unsupported=:lower))
    @test is_l3f_applicable(report)

    result = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:feasibility,
                           unsupported=:lower),
        solver_options=("verbose" => false,))
    @test result.solve.optimal
    @test result.transformers["ct"]["p"] ≈ [6_000.0, 2_000.0] atol=1e-4
    @test result.buses["lv"]["x1"]["vm"] < result.buses["lv"]["x2"]["vm"]

    if isdefined(Main, :OpenDSSDirect)
        dss_v = _l3f_center_tap_opendss_voltages()
        @test real(dss_v["lv.1"]) > 0
        @test real(dss_v["lv.2"]) < 0
        errors = [abs(result.buses["lv"][leg]["vm"] - abs(dss_v[key])) / _L3F_CT_VL
                  for (leg, key) in (("x1", "lv.1"), ("x2", "lv.2"))]
        # L3F drops copper/leakage losses but retains the coupled first-order
        # voltage drop. At 32% aggregate loading this fixture stays well within
        # one percent of the nonlinear three-winding OpenDSS solution.
        @test maximum(errors) < 0.01
        @test maximum(errors) > 1e-5
    end
end
