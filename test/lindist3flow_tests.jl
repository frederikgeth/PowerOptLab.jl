using Test
using JuMP
using Ipopt
using Clarabel
using PowerOptLab

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

function _l3f_regulator_case(; adjustable=false)
    regulator = Dict{String,Any}(
        "bus_from" => "source", "bus_to" => "load",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
        "tap_ratio" => 1.05, "regulator_type" => "B", "s_rating" => 1e6)
    adjustable && merge!(regulator,
        Dict{String,Any}("tap_ratio_min" => 0.9, "tap_ratio_max" => 1.1))
    Dict{String,Any}(
        "bus" => Dict(
            "source" => Dict{String,Any}("terminal_names" => ["a"]),
            "load" => Dict{String,Any}("terminal_names" => ["a"],
                "v_min" => [180.0], "v_max" => [250.0])),
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => "source", "terminal_map" => ["a"],
            "configuration" => "SINGLE_PHASE",
            "v_magnitude" => [230.0], "v_angle" => [0.0])),
        "transformer" => Dict("single_phase_autotransformer" =>
            Dict("reg" => regulator)),
        "load" => Dict("load" => Dict{String,Any}(
            "bus" => "load", "terminal_map" => ["a"],
            "configuration" => "SINGLE_PHASE", "model" => "constant_power",
            "p_nom" => [10_000.0], "q_nom" => [2_000.0])))
end

@testset "LinDist3Flow coefficient oracles" begin
    vp, vq = 230cis(0.2), 218cis(-1.7)
    c = cross_voltage_coefficients(vp, vq)
    @test evaluate_cross_voltage(c, abs2(vp), abs2(vq)) ≈ vp * conj(vq)

    d = [1.0, -1.0, 0.0]
    vbar = ComplexF64[230cis(0), 230cis(-2pi / 3), 230cis(2pi / 3)]
    winding = winding_voltage_coefficients(d, vbar)
    @test evaluate_affine(winding, abs2.(vbar)) ≈ abs2(sum(d .* vbar))

    D = [1.0 -1.0 0.0; 0.0 1.0 -1.0; -1.0 0.0 1.0]
    power_map = connection_power_map(D, vbar)
    s = ComplexF64[10 + 2im, 12 - im, 8 + 3im]
    @test sum(power_map.matrix * s) ≈ sum(s)
    flipped = connection_power_map([-D[1, :]'; D[2, :]'; D[3, :]'], vbar)
    @test flipped.matrix ≈ power_map.matrix

    line = line_drop_coefficients(reshape([0.2 + 0.1im], 1, 1), [230 + 0im])
    @test line.active == reshape([0.4], 1, 1)
    @test line.reactive == reshape([0.2], 1, 1)
end

@testset "LinDist3Flow applicability and Kron boundary" begin
    net = _l3f_case(explicit_neutral=true)
    original = deepcopy(net)
    report = check_l3f_applicability(net)
    @test is_l3f_applicable(report)
    @test report.kron_reduced
    @test net == original

    report = check_l3f_applicability(net; options=L3FOptions(kron_reduce=false))
    @test !is_l3f_applicable(report)
    @test any(f -> f.code == "E.L3F.EXPLICIT_NEUTRAL_UNSUPPORTED", report.findings)

    unsupported = _l3f_case()
    unsupported["capacitor"] = Dict("s" => Dict{String,Any}())
    report = check_l3f_applicability(unsupported)
    @test any(f -> f.code == "E.L3F.COMPONENT_UNSUPPORTED", report.findings)

    meshed = _l3f_case()
    meshed["bus"]["third"] = Dict{String,Any}("terminal_names" => ["a"])
    meshed["line"]["l2"] = Dict{String,Any}(
        "bus_from" => "load", "bus_to" => "third",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"], "linecode" => "lc")
    meshed["line"]["l3"] = Dict{String,Any}(
        "bus_from" => "third", "bus_to" => "source",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"], "linecode" => "lc")
    report = check_l3f_applicability(meshed)
    @test any(f -> f.code == "E.L3F.TOPOLOGY_NOT_RADIAL", report.findings)

    unbounded = _l3f_case(generator=true)
    delete!(unbounded["generator"]["pv"], "q_max")
    report = check_l3f_applicability(unbounded)
    @test any(f -> f.code == "E.L3F.DEVICE_ARITY", report.findings)
end

@testset "LinDist3Flow radial LP" begin
    options = L3FOptions(validate_nonlinear=false, objective=:feasibility)
    build = build_l3f_opf(_l3f_case(), Ipopt.Optimizer; options)
    @test l3f_model_class(build) == :LP
    @test Set(keys(build.constraints)) == Set((:source_voltage, :line_voltage_drop,
                                               :nodal_active_balance, :nodal_reactive_balance))
    result = solve_l3f_opf(_l3f_case(), Ipopt.Optimizer; options,
        solver_options=("print_level" => 0,))
    @test result.solve.optimal
    @test result.lines["line"]["p"] ≈ [10_000.0] atol=1e-4
    @test result.lines["line"]["q"] ≈ [2_000.0] atol=1e-4
    @test result.buses["load"]["a"]["w"] ≈ 48_500.0 atol=1e-3

    reverse = solve_l3f_opf(_l3f_case(generator=true), Ipopt.Optimizer; options,
        solver_options=("print_level" => 0,))
    @test reverse.solve.optimal
    @test reverse.lines["line"]["p"] ≈ [-5_000.0] atol=1e-4
    @test reverse.buses["load"]["a"]["w"] ≈ 54_500.0 atol=1e-3

    shunted = _l3f_case()
    shunted["shunt"] = Dict("g" => Dict{String,Any}(
        "bus" => "load", "terminal_map" => ["a"], "G_1_1" => 0.01))
    shunt_result = solve_l3f_opf(shunted, Ipopt.Optimizer; options,
        solver_options=("print_level" => 0,))
    @test shunt_result.solve.optimal
    @test shunt_result.buses["load"]["a"]["w"] ≈ 48_500 / 1.004 atol=1e-3
end

@testset "LinDist3Flow fixed regulator and exact SOC bounds" begin
    options = L3FOptions(validate_nonlinear=false, objective=:feasibility)
    regulator = solve_l3f_opf(_l3f_regulator_case(), Clarabel.Optimizer; options,
        solver_options=("verbose" => false,))
    @test regulator.solve.optimal
    @test regulator.buses["load"]["a"]["vm"] ≈ 230 / 1.05 atol=1e-5
    @test regulator.transformers["reg"]["p"] ≈ [10_000.0] atol=1e-4

    report = check_l3f_applicability(_l3f_regulator_case(adjustable=true))
    @test any(f -> f.code == "E.L3F.ADJUSTABLE_TAP_UNSUPPORTED", report.findings)

    bounded = _l3f_case(generator=true)
    gen = bounded["generator"]["pv"]
    gen["p_min"] = [0.0]; gen["p_max"] = [15_000.0]
    gen["s_max"] = [8_000.0]
    result = solve_l3f_opf(bounded, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:cost),
        solver_options=("verbose" => false,))
    @test result.solve.optimal
    @test result.generators["pv"]["pg"] ≈ [8_000.0] atol=1e-3
    build = build_l3f_opf(bounded, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false))
    @test l3f_model_class(build) == :SOCP
end

@testset "LinDist3Flow delta constant-power allocation" begin
    v = 230.0
    net = Dict{String,Any}(
        "bus" => Dict("b" => Dict{String,Any}("terminal_names" => ["a", "b", "c"])),
        "voltage_source" => Dict("s" => Dict{String,Any}(
            "bus" => "b", "terminal_map" => ["a", "b", "c"],
            "configuration" => "WYE", "v_magnitude" => fill(v, 3),
            "v_angle" => [0.0, -2pi / 3, 2pi / 3])),
        "load" => Dict("d" => Dict{String,Any}(
            "bus" => "b", "terminal_map" => ["a", "b", "c"],
            "configuration" => "DELTA", "model" => "constant_power",
            "p_nom" => [10_000.0, 12_000.0, 8_000.0],
            "q_nom" => [2_000.0, -1_000.0, 3_000.0])))
    result = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:feasibility),
        solver_options=("verbose" => false,))
    @test result.solve.optimal
    @test sum(result.sources["s"]["pg"]) ≈ 30_000.0 atol=1e-5
    @test sum(result.sources["s"]["qg"]) ≈ 4_000.0 atol=1e-5
end
