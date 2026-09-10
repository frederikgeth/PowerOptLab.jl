using Test
using JuMP
using Clarabel
using PowerOptLab

include("lindist3flow_fixtures.jl")

@testset "LinDist3Flow permissive preprocessing" begin
    solve_perm(net; pu=false) = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(unsupported=:permissive, objective=:feasibility,
                           per_unit=pu), solver_options=("verbose" => false,))

    for pu in (false, true)
        net = _l3f_controls_case()
        net["bus"]["load"]["vuf_max"] = 0.02
        net["line"]["line"]["va_diff_max"] = 0.1
        original = deepcopy(net)
        result = solve_perm(net; pu)
        @test result.solve.optimal
        @test net == original
        findings = filter(f -> startswith(f.code, "A.L3F.PERMISSIVE_"),
                          result.applicability.findings)
        @test any(f -> f.code == "A.L3F.PERMISSIVE_VOLTAGE_LIMIT_DROPPED" &&
                       f.evidence["original_value"] == 0.02, findings)
        @test any(f -> f.code == "A.L3F.PERMISSIVE_LINE_LIMIT_DROPPED" &&
                       f.evidence["original_value"] == 0.1, findings)
        @test result.applicability.lowered

        bad = deepcopy(net); bad["bus"]["load"]["vuf_max"] = NaN
        @test !is_l3f_applicable(check_l3f_applicability(bad;
            options=L3FOptions(unsupported=:permissive)))
        strict = check_l3f_applicability(net;
            options=L3FOptions(unsupported=:approximate))
        @test !is_l3f_applicable(strict)
    end

    for pu in (false, true), kind in ("delta_delta", "closed_delta_regulator")
        net = _l3f_dy_case("delta_wye")
        tx = deepcopy(net["transformer"]["delta_wye"]["t"])
        net["transformer"] = Dict(kind => Dict("t" => tx))
        if kind == "closed_delta_regulator"
            tx["tap_ratio"] = ones(3)
            delete!(tx, "v_nom_from"); delete!(tx, "v_nom_to")
        end
        tx["r_series_from"] = 0.2
        tx["g_no_load"] = 0.01
        tx["s_rating"] = 1e5
        original = deepcopy(net)
        result = solve_perm(net; pu)
        @test result.solve.optimal
        @test net == original
        @test any(f -> f.code == "A.L3F.PERMISSIVE_TRANSFORMER_IDEALIZED" &&
                       haskey(f.evidence, "original_value"), result.applicability.findings)
        @test any(f -> f.code == "A.L3F.PERMISSIVE_TRANSFORMER_LIMIT_DROPPED" &&
                       f.evidence["original_value"] == 1e5,
                  result.applicability.findings)
    end

    malformed = _l3f_dy_case("delta_wye")
    tx = pop!(malformed["transformer"]["delta_wye"], "t")
    malformed["transformer"] = Dict("delta_delta" => Dict("t" => tx))
    tx["r_series_from"] = "bad"
    report = check_l3f_applicability(malformed;
        options=L3FOptions(unsupported=:permissive))
    @test !is_l3f_applicable(report)

    conflicting = _l3f_dy_case("delta_wye")
    ctx = pop!(conflicting["transformer"]["delta_wye"], "t")
    conflicting["transformer"] = Dict("delta_delta" => Dict("t" => ctx))
    ctx["g_no_load"] = 0.01
    ctx["no_load_shunt"] = Dict("winding" => 2, "g" => 0.01, "b" => -0.02)
    conflict_report = check_l3f_applicability(conflicting;
        options=L3FOptions(unsupported=:permissive))
    @test !is_l3f_applicable(conflict_report)
    @test any(f -> f.code == "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED" &&
                   occursin("conflicts", f.message), conflict_report.findings)

    disconnected = _l3f_controls_case()
    disconnected["load"]["load"]["bus"] = "missing"
    @test !is_l3f_applicable(check_l3f_applicability(disconnected;
        options=L3FOptions(unsupported=:permissive)))
end
