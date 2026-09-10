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

function _l3f_controls_solve(net; per_unit=false, objective=:source_import)
    solve_l3f_opf(net, Clarabel.Optimizer; options=L3FOptions(
        unsupported=:lower, objective=objective, per_unit=per_unit, s_base=25_000.0),
        solver_options=("verbose" => false,))
end

function _l3f_controls_three_phase(; phase_cap=[1_500.0, 1_500.0, 1_500.0],
                                   aggregate_cap=3_000.0)
    tm = ["a", "b", "c"]
    Dict{String,Any}(
        "terminal_conventions" => Dict("phase" => tm, "neutral" => String[]),
        "bus" => Dict("b" => Dict{String,Any}("terminal_names" => tm)),
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => "b", "terminal_map" => tm, "configuration" => "WYE",
            "v_magnitude" => fill(230.0, 3),
            "v_angle" => [0.0, -2pi/3, 2pi/3], "cost" => ones(3))),
        "generator" => Dict("g" => Dict{String,Any}(
            "bus" => "b", "terminal_map" => tm, "configuration" => "WYE",
            "p_min" => fill(1_000.0, 3), "p_max" => fill(1_000.0, 3),
            "q_min" => [1_000.0, -1_000.0, 0.0],
            "q_max" => [1_000.0, -1_000.0, 0.0], "s_max" => phase_cap,
            "aggregate_s_max" => aggregate_cap, "cost" => zeros(3))))
end

@testset "LinDist3Flow restricted controls" begin
    @testset "native IBR lowering preserves input and signed PF" begin
        original = Dict{String,Any}(
            "generator" => Dict{String,Any}(),
            "control_profile" => Dict("lead" => Dict(
                "power_factor" => Dict("pf" => -0.8))),
            "ibr" => Dict("pv" => Dict{String,Any}(
                "bus" => "b", "terminal_map" => ["a", "b", "c"],
                "_l3f_neutral_reduced" => true,
                "topology" => "FOUR_LEG", "p_min" => fill(0.0, 3),
                "p_max" => fill(10.0, 3), "q_min" => fill(-10.0, 3),
                "q_max" => fill(10.0, 3), "s_max" => fill(12.0, 3),
                "control_profile" => "lead")))
        net = deepcopy(original); findings = PowerOptLab.L3FFinding[]
        PowerOptLab._l3f_lower_restricted_controls!(findings, net)
        @test original["ibr"]["pv"]["terminal_map"] == ["a", "b", "c"]
        @test isempty(net["ibr"])
        @test net["generator"]["pv"]["configuration"] == "WYE"
        @test net["generator"]["pv"]["terminal_map"] == ["a", "b", "c"]
        @test net["generator"]["pv"]["fixed_pf"] == -0.8
        @test any(f -> f.code == "L.L3F.IBR_TO_GENERATOR", findings)
    end

    @testset "unsupported control laws are rejected rather than ignored" begin
        net = Dict{String,Any}(
            "control_profile" => Dict("vv" => Dict("volt_var" => Dict())),
            "ibr" => Dict("pv" => Dict{String,Any}(
                "bus" => "b", "terminal_map" => ["a", "n"],
                "topology" => "SINGLE_PHASE", "p_min" => [0.0], "p_max" => [1.0],
                "q_min" => [-1.0], "q_max" => [1.0], "control_profile" => "vv")))
        findings = PowerOptLab.L3FFinding[]
        PowerOptLab._l3f_lower_restricted_controls!(findings, net)
        @test any(f -> f.code == "E.L3F.IBR_UNSUPPORTED", findings)
        @test !haskey(get(net, "generator", Dict()), "pv")
    end

    @testset "aggregate and hard-target validation" begin
        valid = Dict{String,Any}("generator" => Dict("g" => Dict{String,Any}(
            "p_min" => [5.0], "p_max" => [5.0], "q_min" => [-2.0],
            "q_max" => [2.0], "v_target" => 230.0,
            "aggregate_s_max" => 8.0)))
        findings = PowerOptLab.L3FFinding[]
        PowerOptLab._l3f_validate_restricted_controls!(findings, valid)
        @test isempty(findings)
        invalid = deepcopy(valid)
        invalid["generator"]["g"]["p_max"] = [6.0]
        PowerOptLab._l3f_validate_restricted_controls!(findings, invalid)
        @test any(f -> f.code == "E.L3F.GENERATOR_CONTROL_INVALID", findings)
    end

    @testset "aggregate apparent power permits phase-Q cancellation" begin
        for per_unit in (false, true)
            result = _l3f_controls_solve(_l3f_controls_three_phase(); per_unit)
            @test result.solve.optimal
            @test result.generators["g"]["qg"] ≈ [1_000.0, -1_000.0, 0.0] atol=1e-3
            @test !_l3f_controls_solve(_l3f_controls_three_phase(
                phase_cap=[1_200.0, 1_500.0, 1_500.0]); per_unit).solve.optimal
            @test !_l3f_controls_solve(_l3f_controls_three_phase(
                aggregate_cap=2_999.0); per_unit).solve.optimal
        end
    end

    @testset "public generator controls in SI and PU" begin
        base = Dict{String,Any}(
            "bus" => "load", "terminal_map" => ["a"], "configuration" => "WYE",
            "p_min" => [5_000.0], "p_max" => [5_000.0],
            "q_min" => [-5_000.0], "q_max" => [5_000.0], "cost" => [0.0])
        pq = deepcopy(base); pq["q_min"] = [1_000.0]; pq["q_max"] = [1_000.0]
        for per_unit in (false, true)
            result = _l3f_controls_solve(_l3f_controls_case(generator=pq); per_unit)
            @test result.solve.optimal
            @test result.generators["g"]["pg"] ≈ [5_000.0] atol=1e-4
            @test result.generators["g"]["qg"] ≈ [1_000.0] atol=1e-4

            pf = deepcopy(base); pf["fixed_pf"] = -0.8
            result = _l3f_controls_solve(_l3f_controls_case(generator=pf); per_unit)
            @test result.solve.optimal
            @test result.generators["g"]["qg"][1] ≈ 3_750.0 atol=1e-3

            aggregate = deepcopy(base)
            aggregate["aggregate_p_min"] = 5_000.0
            aggregate["aggregate_p_max"] = 5_000.0
            aggregate["aggregate_q_min"] = 1_500.0
            aggregate["aggregate_q_max"] = 1_500.0
            aggregate["aggregate_s_max"] = 5_500.0
            result = _l3f_controls_solve(_l3f_controls_case(generator=aggregate); per_unit)
            @test result.solve.optimal
            @test result.generators["g"]["qg"] ≈ [1_500.0] atol=1e-3
            phase_limited = deepcopy(aggregate); phase_limited["s_max"] = [5_100.0]
            @test !_l3f_controls_solve(
                _l3f_controls_case(generator=phase_limited); per_unit).solve.optimal
        end
    end

    @testset "hard voltage target and reactive cap" begin
        generator = Dict{String,Any}(
            "bus" => "load", "terminal_map" => ["a"], "configuration" => "WYE",
            "p_min" => [5_000.0], "p_max" => [5_000.0],
            "q_min" => [-5_000.0], "q_max" => [5_000.0], "cost" => [0.0],
            "v_target" => [sqrt(51_000.0)])
        for per_unit in (false, true)
            result = _l3f_controls_solve(_l3f_controls_case(generator=generator); per_unit)
            @test result.solve.optimal
            @test result.generators["g"]["qg"] ≈ [2_500.0] atol=2e-3
            @test result.validation["status"] == "not_requested"
        end
        tight = deepcopy(generator); tight["q_min"] = [-2_000.0]; tight["q_max"] = [2_000.0]
        @test !_l3f_controls_solve(_l3f_controls_case(generator=tight)).solve.optimal
    end

    @testset "native IBR economic dispatch, trace, and preservation" begin
        for pf in (-0.8, 0.8), objective in (:cost, :source_import), per_unit in (false, true)
            ibr = Dict{String,Any}(
                "bus" => "load", "terminal_map" => ["a"],
                "topology" => "SINGLE_PHASE", "s_max" => [4_000.0],
                "fixed_pf" => pf)
            net = _l3f_controls_case(ibr=ibr); original = deepcopy(net)
            result = _l3f_controls_solve(net; per_unit, objective)
            @test result.solve.optimal
            @test net == original
            @test result.generators["pv"]["pg"][1] ≈ 3_200.0 atol=2e-3
            @test result.generators["pv"]["qg"][1] ≈ -sign(pf) * 2_400.0 atol=2e-3
            @test result.generators["pv"]["original_component"] == "ibr"
            @test result.generators["pv"]["original_id"] == "pv"
            @test any(f -> f.code == "L.L3F.IBR_TO_GENERATOR",
                      result.applicability.findings)
        end
    end

    @testset "unsupported native controls fail publicly" begin
        base = Dict{String,Any}("bus" => "load", "terminal_map" => ["a"],
            "topology" => "SINGLE_PHASE", "s_max" => [4_000.0])
        for (field, value) in (("grid_forming", true), ("r_filter", [0.1]),
                               ("dc_control", "droop"), ("time_series", Dict()))
            bad = deepcopy(base); bad[field] = value
            report = check_l3f_applicability(_l3f_controls_case(ibr=bad);
                options=L3FOptions(unsupported=:lower))
            @test !is_l3f_applicable(report)
            @test any(f -> f.code == "E.L3F.IBR_UNSUPPORTED", report.findings)
        end
        unknown = deepcopy(base); unknown["mystery_physics"] = 1.0
        conflict = deepcopy(base); conflict["fixed_pf"] = 0.9
        conflict["control_profile"] = "pf"
        conflict_net = _l3f_controls_case(ibr=conflict)
        conflict_net["control_profile"] = Dict("pf" =>
            Dict("power_factor" => Dict("pf" => 0.9)))
        intersected = deepcopy(base)
        intersected["dc_link_coupled"] = true
        intersected["p_dc_min"] = 0.0; intersected["p_dc_max"] = 0.0
        intersected["aggregate_p_min"] = 1.0
        for badnet in (_l3f_controls_case(ibr=unknown), conflict_net,
                       _l3f_controls_case(ibr=intersected))
            report = check_l3f_applicability(badnet;
                options=L3FOptions(unsupported=:lower))
            @test !is_l3f_applicable(report)
        end
    end

    @testset "permissive mode drops controls but retains dispatch" begin
        ibr = Dict{String,Any}(
            "bus" => "load", "terminal_map" => ["a"],
            "topology" => "SINGLE_PHASE", "s_max" => [4_000.0],
            "p_min" => [2_000.0], "p_max" => [2_000.0],
            "q_min" => [500.0], "q_max" => [500.0],
            "r_filter" => [0.1], "x_filter" => [0.2],
            "grid_forming" => true, "v_ref_internal" => 230.0,
            "p_avail" => 3_000.0)
        net = _l3f_controls_case(ibr=ibr); original = deepcopy(net)
        for per_unit in (false, true)
            result = solve_l3f_opf(net, Clarabel.Optimizer; options=L3FOptions(
                unsupported=:permissive, objective=:source_import, per_unit=per_unit,
                s_base=25_000.0), solver_options=("verbose" => false,))
            @test result.solve.optimal
            @test result.generators["pv"]["pg"] ≈ [2_000.0] atol=1e-3
            @test result.generators["pv"]["qg"] ≈ [500.0] atol=1e-3
            @test net == original
            warning = only(filter(f -> f.code == "A.L3F.IBR_FIELDS_DROPPED",
                                  result.applicability.findings))
            fields = warning.evidence["original_fields"]
            @test all(haskey(fields, key) for key in
                      ("r_filter", "x_filter", "grid_forming", "v_ref_internal"))
            @test !haskey(fields, "p_avail")
        end

        droop = deepcopy(ibr)
        delete!(droop, "r_filter"); delete!(droop, "x_filter")
        delete!(droop, "grid_forming"); delete!(droop, "v_ref_internal")
        droop["control_profile"] = "vv"
        droopnet = _l3f_controls_case(ibr=droop)
        droopnet["control_profile"] = Dict("vv" => Dict("volt_var" => Dict(
            "breakpoints" => [0.9, 0.95, 1.05, 1.1], "q_limits" => [-1.0, 1.0])))
        result = solve_l3f_opf(droopnet, Clarabel.Optimizer; options=L3FOptions(
            unsupported=:permissive, objective=:source_import),
            solver_options=("verbose" => false,))
        @test result.solve.optimal
        @test result.generators["pv"]["pg"] ≈ [2_000.0] atol=1e-3

        malformed = deepcopy(droopnet)
        malformed["control_profile"]["vv"]["volt_var"]["breakpoints"][2] = NaN
        report = check_l3f_applicability(malformed;
            options=L3FOptions(unsupported=:permissive))
        @test !is_l3f_applicable(report)

        dcgraph = _l3f_controls_case(ibr=Dict{String,Any}(
            "bus" => "load", "terminal_map" => ["a"],
            "topology" => "SINGLE_PHASE", "s_max" => [4_000.0]))
        dcgraph["dc_bus"] = Dict("dc" => Dict{String,Any}())
        report = check_l3f_applicability(dcgraph;
            options=L3FOptions(unsupported=:permissive))
        @test !is_l3f_applicable(report)
        @test any(f -> f.code == "E.L3F.DC_SUBSYSTEM_UNSUPPORTED", report.findings)
    end

    @testset "permissive neutral current handling uses declared role" begin
        net = Dict{String,Any}(
            "terminal_conventions" => Dict("phase" => ["a", "b"], "neutral" => ["n"]),
            "ibr" => Dict("x" => Dict{String,Any}(
                "bus" => "b", "terminal_map" => ["a", "n", "b"],
                "topology" => "FOUR_LEG", "s_max" => [1.0, 1.0],
                "i_max" => [10.0, 5.0, 20.0])))
        findings = PowerOptLab.L3FFinding[]
        PowerOptLab._l3f_permissive_restricted_controls!(findings, net)
        @test net["ibr"]["x"]["i_max"] == [10.0, 20.0]
        @test only(findings).evidence["original_fields"]["i_max_neutral"] == 5.0

        malformed = deepcopy(net)
        malformed["ibr"]["x"]["r_filter"] = [NaN]
        findings = PowerOptLab.L3FFinding[]
        PowerOptLab._l3f_permissive_restricted_controls!(findings, malformed)
        @test haskey(malformed["ibr"]["x"], "r_filter")
    end
end

@testset "LinDist3Flow projection result contract" begin
    for per_unit in (false, true)
        net = _l3f_controls_case()
        net["bus"]["load"]["vuf_max"] = 0.02
        original = deepcopy(net)
        options = L3FOptions(unsupported=:permissive, per_unit=per_unit,
                             objective=:feasibility)
        result = solve_l3f_opf(net; options, solver_options=("verbose" => false,))
        @test result.solve.optimal
        @test net == original
        @test result.formulation["network_semantics"] == "projected"
        @test result.formulation["physical_feasibility_certified"] === false
        @test result.validation["status"] == "not_requested"
        projection = only(filter(p -> get(p["evidence"], "field", nothing) == "vuf_max",
                                 result.formulation["projections"]))
        @test projection["evidence"]["original_value"] == 0.02
        @test projection["component"] == "bus" && projection["id"] == "load"
        net["bus"]["load"]["v_min"] = [240.0]
        @test !solve_l3f_opf(net; options,
                            solver_options=("verbose" => false,)).solve.optimal
    end
end
