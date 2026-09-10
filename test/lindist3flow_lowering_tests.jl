using Test
using LinearAlgebra
using JuMP
using Ipopt
using Clarabel
using PowerOptLab

# `L3FOptions(unsupported=...)` widens the admissible input without widening the
# formulation. The tests below separate the two claims it makes:
#
#   :lower        the rewrite preserves the canonical L3F component data. Each
#                 case may be compared with a hand-written supported equivalent,
#                 but neither side is an exact AC model.
#   :approximate  the rewrite is a projection. The tests pin what is preserved
#                 (tangency at nominal voltage, total ZIP fractions) and what is
#                 not, rather than asserting an accuracy the construction cannot
#                 deliver.

_l3f_low_opts(; kwargs...) = L3FOptions(; validate_nonlinear=false,
                                        objective=:feasibility, kwargs...)
_l3f_low_solve(net; kwargs...) = solve_l3f_opf(net, Clarabel.Optimizer;
    options=_l3f_low_opts(; kwargs...), solver_options=("verbose" => false,))
_l3f_low_codes(report) = [f.code for f in report.findings]
_l3f_low_has(report, code) = any(==(code), _l3f_low_codes(report))
_l3f_low_finding(report, code) =
    only(f for f in report.findings if f.code == code)

"""
Three-bus single-phase feeder `s -[l1]- m -[l2]- l`.

Pass `bridge=false` to leave the `m`-to-`l` gap open for a test that inserts the
device under study there.
"""
function _l3f_low_case(; bridge::Bool=true)
    Dict{String,Any}(
        "bus" => Dict(
            "s" => Dict{String,Any}("terminal_names" => ["a"]),
            "m" => Dict{String,Any}("terminal_names" => ["a"]),
            "l" => Dict{String,Any}("terminal_names" => ["a"])),
        "linecode" => Dict("lc" => Dict{String,Any}(
            "R_series_1_1" => 0.2, "X_series_1_1" => 0.1)),
        "line" => Dict{String,Any}(
            "l1" => Dict{String,Any}(
                "bus_from" => "s", "bus_to" => "m",
                "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
                "linecode" => "lc"),
            (bridge ? ("l2" => Dict{String,Any}(
                "bus_from" => "m", "bus_to" => "l",
                "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
                "linecode" => "lc"),) : ())...),
        "voltage_source" => Dict("v" => Dict{String,Any}(
            "bus" => "s", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
            "v_magnitude" => [230.0], "v_angle" => [0.0], "cost" => [1.0])),
        "load" => Dict("d" => Dict{String,Any}(
            "bus" => "l", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
            "model" => "constant_power",
            "p_nom" => [10_000.0], "q_nom" => [2_000.0])))
end

"""Assert two canonical L3F networks agree at every terminal."""
function _l3f_low_agree(a, b; rtol=1e-8)
    for bus in intersect(keys(a.buses), keys(b.buses))
        for terminal in keys(a.buses[bus])
            @test a.buses[bus][terminal]["w"] ≈ b.buses[bus][terminal]["w"] rtol=rtol
        end
    end
    for (sid, source) in a.sources
        @test source["pg"] ≈ b.sources[sid]["pg"] rtol=rtol
        @test source["qg"] ≈ b.sources[sid]["qg"] rtol=rtol
    end
end

@testset "LinDist3Flow lowering policy ladder" begin
    net = _l3f_low_case(bridge=false)
    net["switch"] = Dict("sw" => Dict{String,Any}(
        "bus_from" => "m", "bus_to" => "l",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
        "open_switch" => false))

    # :reject is the default and is unchanged by any of this.
    strict = check_l3f_applicability(net)
    @test !is_l3f_applicable(strict)
    @test !strict.lowered
    @test _l3f_low_has(strict, "E.L3F.SWITCH_UNSUPPORTED")
    @test L3FOptions().unsupported == :reject

    lowered = check_l3f_applicability(net; options=L3FOptions(unsupported=:lower))
    @test is_l3f_applicable(lowered)
    @test lowered.lowered
    @test _l3f_low_has(lowered, "L.L3F.SWITCH_LOWERED")

    # A canonical lowering is :info, never :warning — it makes no AC accuracy claim.
    @test _l3f_low_finding(lowered, "L.L3F.SWITCH_LOWERED").severity == :info
    @test all(f -> f.severity != :error, lowered.findings)

    # The policy and its effect reach the published result contract.
    result = _l3f_low_solve(net; unsupported=:lower)
    @test result.formulation["unsupported_policy"] == "lower"
    @test result.formulation["lowered"] == true
    @test _l3f_low_solve(_l3f_low_case()).formulation["lowered"] == false

    @test_throws ArgumentError L3FOptions(unsupported=:whatever)
end

@testset "LinDist3Flow canonical lowering: switches" begin
    switched = _l3f_low_case(bridge=false)
    switched["switch"] = Dict("sw" => Dict{String,Any}(
        "bus_from" => "m", "bus_to" => "l",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
        "open_switch" => false, "s_max" => [100_000.0], "i_max" => [500.0]))
    # The hand-written equivalent: a branch with no series impedance.
    manual = _l3f_low_case(bridge=false)
    manual["line"]["sw"] = Dict{String,Any}(
        "bus_from" => "m", "bus_to" => "l",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
        "R_series_1_1" => 0.0, "X_series_1_1" => 0.0)
    _l3f_low_agree(_l3f_low_solve(switched; unsupported=:lower), _l3f_low_solve(manual))

    # A closed ideal switch imposes equal squared voltage on both sides.
    result = _l3f_low_solve(switched; unsupported=:lower)
    @test result.buses["m"]["a"]["w"] ≈ result.buses["l"]["a"]["w"] rtol=1e-10
    prepared = PowerOptLab._l3f_prepare(switched;
        options=L3FOptions(validate_nonlinear=false, unsupported=:lower))
    @test prepared.network["line"]["_l3f_switch_sw"]["s_max"] == [100_000.0]
    @test prepared.network["line"]["_l3f_switch_sw"]["i_max"] == [500.0]
    switch_build = build_l3f_opf(switched, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, unsupported=:lower))
    @test length(switch_build.constraints[:line_apparent_power]) == 2
    @test length(switch_build.constraints[:line_current]) == 2

    # An open switch carries no current, so its far side is a separate island
    # and is reported as such rather than being quietly energized.
    opened = deepcopy(switched)
    opened["switch"]["sw"]["open_switch"] = true
    report = check_l3f_applicability(opened; options=L3FOptions(unsupported=:lower))
    @test _l3f_low_has(report, "L.L3F.SWITCH_OPEN_REMOVED")
    @test _l3f_low_has(report, "E.L3F.SOURCE_MISSING")
    @test !is_l3f_applicable(report)
end

@testset "LinDist3Flow canonical lowering: capacitors" begin
    for (configuration, terminals, q_rated, v_nom) in (
            ("SINGLE_PHASE", ["a"], [5_000.0], 230.0),
            ("WYE", ["a","b","c"], [5_000.0, 4_000.0, 6_000.0], 230.0),
            ("DELTA", ["a","b","c"], [5_000.0, 4_000.0, 6_000.0], 400.0))
        n = length(terminals)
        base = Dict{String,Any}(
            "bus" => Dict("s" => Dict{String,Any}("terminal_names" => copy(terminals)),
                          "l" => Dict{String,Any}("terminal_names" => copy(terminals))),
            "linecode" => Dict("lc" => Dict{String,Any}(
                Dict("R_series_$(k)_$(k)" => 0.2 for k in 1:n)...,
                Dict("X_series_$(k)_$(k)" => 0.1 for k in 1:n)...)),
            "line" => Dict("l1" => Dict{String,Any}(
                "bus_from" => "s", "bus_to" => "l",
                "terminal_map_from" => copy(terminals),
                "terminal_map_to" => copy(terminals), "linecode" => "lc")),
            "voltage_source" => Dict("v" => Dict{String,Any}(
                "bus" => "s", "terminal_map" => copy(terminals),
                "configuration" => n == 1 ? "SINGLE_PHASE" : "WYE",
                "v_magnitude" => fill(230.0, n),
                "v_angle" => n == 1 ? [0.0] : [0.0, -2pi/3, 2pi/3])),
            "load" => Dict("d" => Dict{String,Any}(
                "bus" => "l", "terminal_map" => copy(terminals),
                "configuration" => n == 1 ? "SINGLE_PHASE" : "WYE",
                "model" => "constant_power",
                "p_nom" => fill(10_000.0, n), "q_nom" => fill(2_000.0, n))))

        capacitive = deepcopy(base)
        capacitive["capacitor"] = Dict("c" => Dict{String,Any}(
            "bus" => "l", "terminal_map" => copy(terminals),
            "configuration" => configuration, "q_rated" => copy(q_rated),
            "v_nom" => v_nom))

        # Hand-written equivalent: the terminal admittance matrix of the same
        # coils, D' * diag(q / v_nom^2) * D.
        D = configuration == "DELTA" ?
            [1.0 -1.0 0.0; 0.0 1.0 -1.0; -1.0 0.0 1.0] : Matrix{Float64}(I, n, n)
        Y = transpose(D) * Diagonal(q_rated ./ v_nom^2) * D
        manual = deepcopy(base)
        shunt = Dict{String,Any}("bus" => "l", "terminal_map" => copy(terminals))
        for i in 1:n, j in i:n
            iszero(Y[i, j]) || (shunt["B_$(i)_$(j)"] = Y[i, j])
        end
        manual["shunt"] = Dict("c" => shunt)

        _l3f_low_agree(_l3f_low_solve(capacitive; unsupported=:lower),
                       _l3f_low_solve(manual))
    end

    # Invalid capacitor data is an error, not a silent skip.
    broken = _l3f_low_case()
    broken["capacitor"] = Dict("c" => Dict{String,Any}(
        "bus" => "l", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
        "q_rated" => [5_000.0], "v_nom" => 0.0))
    @test _l3f_low_has(check_l3f_applicability(broken;
        options=L3FOptions(unsupported=:lower)), "E.L3F.CAPACITOR_INVALID")
end

@testset "LinDist3Flow canonical lowering: line shunts" begin
    # BMOPF requires exactly one impedance source. Retain deterministic inline
    # precedence internally, but reject the ambiguous public declaration.
    shunted = _l3f_low_case()
    merge!(shunted["line"]["l2"], Dict{String,Any}(
        "R_series_1_1" => 0.2, "X_series_1_1" => 0.1,
        "G_from_1_1" => 2e-4, "B_from_1_1" => 3e-4,
        "G_to_1_1" => 1e-4, "B_to_1_1" => 5e-4))
    shunted["linecode"]["lc2"] = Dict{String,Any}(
        "R_series_1_1" => 0.2, "X_series_1_1" => 0.1,
        "G_from_1_1" => 9.0)
    shunted["line"]["l2"]["linecode"] = "lc2"
    dual_source = check_l3f_applicability(shunted;
        options=L3FOptions(unsupported=:lower))
    @test _l3f_low_has(dual_source, "E.L3F.LINE_IMPEDANCE_SOURCE")
    delete!(shunted["line"]["l2"], "linecode")

    # Hand-written pi equivalent: each declared half on its own bus.
    manual = _l3f_low_case()
    delete!(manual["line"]["l2"], "linecode")
    merge!(manual["line"]["l2"], Dict{String,Any}(
        "R_series_1_1" => 0.2, "X_series_1_1" => 0.1))
    manual["shunt"] = Dict(
        "from" => Dict{String,Any}("bus" => "m", "terminal_map" => ["a"],
            "G_1_1" => 2e-4, "B_1_1" => 3e-4),
        "to" => Dict{String,Any}("bus" => "l", "terminal_map" => ["a"],
            "G_1_1" => 1e-4, "B_1_1" => 5e-4))

    _l3f_low_agree(_l3f_low_solve(shunted; unsupported=:lower), _l3f_low_solve(manual))
    @test _l3f_low_has(check_l3f_applicability(shunted;
        options=L3FOptions(unsupported=:lower)), "L.L3F.LINE_SHUNT_LOWERED")
    # Still rejected under the default policy.
    @test _l3f_low_has(check_l3f_applicability(shunted), "E.L3F.LINE_SHUNT_UNSUPPORTED")

    # The unused linecode is never merged into the selected inline source.
    prepared = PowerOptLab._l3f_prepare(shunted;
        options=L3FOptions(validate_nonlinear=false, unsupported=:lower))
    @test prepared.network["shunt"]["_l3f_lineshunt_l2_from"]["G_1_1"] == 2e-4
    rated_inline = deepcopy(shunted)
    rated_inline["line"]["l2"]["i_max"] = [100.0]
    rated_report = check_l3f_applicability(rated_inline;
        options=L3FOptions(unsupported=:lower))
    @test is_l3f_applicable(rated_report)
    rated_build = build_l3f_opf(rated_inline, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, unsupported=:lower))
    @test length(rated_build.constraints[:line_current]) == 2
    @test all(JuMP.constraint_object(c).set isa JuMP.MOI.RotatedSecondOrderCone
              for c in values(rated_build.constraints[:line_current]))

    # The line-owned pi shunt belongs inside the line's endpoint ampacity cone;
    # an electrically identical standalone bank belongs only in nodal balance.
    manual_rated = deepcopy(manual)
    manual_rated["line"]["l2"]["i_max"] = [100.0]
    manual_build = build_l3f_opf(manual_rated, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false))
    lowered_cone = JuMP.constraint_object(
        rated_build.constraints[:line_current][("l2", "from", 1)]).func
    standalone_cone = JuMP.constraint_object(
        manual_build.constraints[:line_current][("l2", "from", 1)]).func
    lowered_w = rated_build.variables[:w][("m", "a")]
    standalone_w = manual_build.variables[:w][("m", "a")]
    @test any(!iszero(JuMP.coefficient(lowered_cone[i], lowered_w)) for i in 3:4)
    @test all(iszero(JuMP.coefficient(standalone_cone[i], standalone_w)) for i in 3:4)

    # A linecode carries per-unit-length shunt, so it scales with `length`.
    coded = _l3f_low_case()
    coded["linecode"]["shunted"] = Dict{String,Any}(
        "R_series_1_1" => 0.2, "X_series_1_1" => 0.1, "B_from_1_1" => 1e-4)
    merge!(coded["line"]["l2"],
        Dict{String,Any}("linecode" => "shunted", "length" => 3.0))
    scaled = _l3f_low_case()
    delete!(scaled["line"]["l2"], "linecode")
    merge!(scaled["line"]["l2"],
        Dict{String,Any}("R_series_1_1" => 0.6, "X_series_1_1" => 0.3))
    scaled["shunt"] = Dict("s" => Dict{String,Any}(
        "bus" => "m", "terminal_map" => ["a"], "B_1_1" => 3e-4))
    _l3f_low_agree(_l3f_low_solve(coded; unsupported=:lower), _l3f_low_solve(scaled))
    inherited_rated = deepcopy(coded)
    inherited_rated["linecode"]["shunted"]["i_max"] = [100.0]
    inherited_report = check_l3f_applicability(inherited_rated;
        options=L3FOptions(unsupported=:lower))
    @test is_l3f_applicable(inherited_report)

    # A from-end pi shunt must be inside the line rating even when the series
    # path carries zero power. This distinguishes endpoint total power from a
    # cone stamped on the series-flow variable alone.
    endpoint_rated = _l3f_low_case()
    empty!(endpoint_rated["load"])
    delete!(endpoint_rated["line"]["l2"], "linecode")
    merge!(endpoint_rated["line"]["l2"], Dict{String,Any}(
        "R_series_1_1" => 0.2, "X_series_1_1" => 0.1,
        "B_from_1_1" => 1.0e-2, "s_max" => [100.0]))
    @test !_l3f_low_solve(endpoint_rated; unsupported=:lower).solve.optimal
    endpoint_rated["line"]["l2"]["s_max"] = [1_000.0]
    @test _l3f_low_solve(endpoint_rated; unsupported=:lower).solve.optimal
    delete!(endpoint_rated["line"]["l2"], "s_max")
    endpoint_rated["line"]["l2"]["i_max"] = [0.1]
    @test !_l3f_low_solve(endpoint_rated; unsupported=:lower).solve.optimal
    @test !_l3f_low_solve(endpoint_rated; unsupported=:lower,
                          per_unit=false).solve.optimal
    endpoint_rated["line"]["l2"]["i_max"] = [5.0]
    @test _l3f_low_solve(endpoint_rated; unsupported=:lower).solve.optimal
    @test _l3f_low_solve(endpoint_rated; unsupported=:lower,
                         per_unit=false).solve.optimal

    # A fully stored matrix preserves both orientations, while a triangular
    # source is mirrored during reconstruction. Direct entries win when both
    # orientations are deliberately present.
    full = Dict{String,Any}(
        "bus" => Dict("a" => Dict("terminal_names" => ["a", "b"]),
                      "b" => Dict("terminal_names" => ["a", "b"])),
        "line" => Dict("l" => Dict{String,Any}(
            "bus_from" => "a", "bus_to" => "b",
            "terminal_map_from" => ["a", "b"],
            "terminal_map_to" => ["a", "b"],
            "R_series_1_1" => 0.2, "X_series_1_1" => 0.1,
            "R_series_2_2" => 0.2, "X_series_2_2" => 0.1,
            "G_from_1_2" => 1.0, "G_from_2_1" => 2.0,
            "G_from_1_1" => 3.0)),
        "linecode" => Dict{String,Any}(), "shunt" => Dict{String,Any}())
    triangular = deepcopy(full)
    delete!(triangular["line"]["l"], "G_from_2_1")
    full_findings = PowerOptLab.L3FFinding[]
    PowerOptLab._l3f_lower_line_shunts!(full_findings, full)
    lowered_full = full["shunt"]["_l3f_lineshunt_l_from"]
    @test lowered_full["G_1_2"] == 1.0
    @test lowered_full["G_2_1"] == 2.0
    @test lowered_full["G_1_1"] == 3.0
    @test any(f -> f.code == "W.L3F.LINE_SHUNT_ASYMMETRIC", full_findings)

    triangular_findings = PowerOptLab.L3FFinding[]
    PowerOptLab._l3f_lower_line_shunts!(triangular_findings, triangular)
    lowered_triangular = triangular["shunt"]["_l3f_lineshunt_l_from"]
    @test lowered_triangular["G_1_2"] == 1.0
    @test lowered_triangular["G_2_1"] == 1.0
    @test !any(f -> f.code == "W.L3F.LINE_SHUNT_ASYMMETRIC", triangular_findings)

    invalid = deepcopy(coded)
    invalid["linecode"]["shunted"]["G_from_2_2"] = 1e-4
    report = check_l3f_applicability(invalid; options=L3FOptions(unsupported=:lower))
    @test _l3f_low_has(report, "E.L3F.LINE_SHUNT_INVALID")
end

@testset "LinDist3Flow canonical lowering: transformer impedance" begin
    function transformer_case(; extra=Dict{String,Any}())
        net = _l3f_low_case(bridge=false)
        net["transformer"] = Dict("single_phase" => Dict("t" => merge(Dict{String,Any}(
            "bus_from" => "m", "bus_to" => "l",
            "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
            "v_nom_from" => 230.0, "v_nom_to" => 115.0, "s_rating" => 1.0e5), extra)))
        net["load"]["d"]["p_nom"] = [4_000.0]; net["load"]["d"]["q_nom"] = [800.0]
        net
    end

    # From-side leakage: an ideal unit behind a series impedance.
    leaky = transformer_case(extra=Dict{String,Any}(
        "r_series_from" => 0.05, "x_series_from" => 0.15))
    manual = transformer_case()
    manual["bus"]["int"] = Dict{String,Any}("terminal_names" => ["a"])
    manual["transformer"]["single_phase"]["t"]["bus_from"] = "int"
    manual["line"]["leak"] = Dict{String,Any}(
        "bus_from" => "m", "bus_to" => "int",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
        "R_series_1_1" => 0.05, "X_series_1_1" => 0.15)
    _l3f_low_agree(_l3f_low_solve(leaky; unsupported=:lower), _l3f_low_solve(manual))

    # No-load admittance: a shunt across the to-side coil.
    magnetising = transformer_case(extra=Dict{String,Any}(
        "g_no_load" => 1.0e-4, "b_no_load" => -3.0e-4))
    manual_shunt = transformer_case()
    manual_shunt["shunt"] = Dict("nl" => Dict{String,Any}(
        "bus" => "l", "terminal_map" => ["a"],
        "G_1_1" => 1.0e-4, "B_1_1" => -3.0e-4))
    _l3f_low_agree(_l3f_low_solve(magnetising; unsupported=:lower),
                   _l3f_low_solve(manual_shunt))

    # Both windings and the shunt together.
    full = transformer_case(extra=Dict{String,Any}(
        "r_series_from" => 0.05, "x_series_from" => 0.15,
        "r_series_to" => 0.01, "x_series_to" => 0.03,
        "g_no_load" => 1.0e-4, "b_no_load" => -3.0e-4))
    report = check_l3f_applicability(full; options=L3FOptions(unsupported=:lower))
    @test is_l3f_applicable(report)
    @test count(==("L.L3F.TRANSFORMER_LEAKAGE_LOWERED"), _l3f_low_codes(report)) == 2
    @test _l3f_low_has(report, "L.L3F.TRANSFORMER_NO_LOAD_LOWERED")

    result = _l3f_low_solve(full; unsupported=:lower)
    @test result.solve.optimal
    # The internal buses and leakage branches are visible in the result, which is
    # what makes the rewrite auditable rather than hidden.
    @test haskey(result.buses, "_l3f_xfmr_t_from")
    @test haskey(result.buses, "_l3f_xfmr_t_to")
    @test haskey(result.lines, "_l3f_leakage_t_from")
    @test haskey(result.lines, "_l3f_leakage_t_to")
    @test result.buses["m"]["a"]["w"] > result.buses["_l3f_xfmr_t_from"]["a"]["w"]

    # An ideal transformer is untouched: nothing to lower, nothing reported.
    ideal = check_l3f_applicability(transformer_case();
        options=L3FOptions(unsupported=:lower))
    @test !ideal.lowered
    @test !any(f -> startswith(f.code, "L.L3F."), ideal.findings)

    neutral = transformer_case(extra=Dict{String,Any}("r_neutral_from" => 0.1))
    @test _l3f_low_has(check_l3f_applicability(neutral;
        options=L3FOptions(unsupported=:lower)),
        "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED")
    rated_no_load = transformer_case(extra=Dict{String,Any}(
        "g_no_load" => 1.0e-4, "i_max_to" => [100.0]))
    rated_no_load_report = check_l3f_applicability(rated_no_load;
        options=L3FOptions(unsupported=:lower))
    @test is_l3f_applicable(rated_no_load_report)
    rated_no_load_build = build_l3f_opf(rated_no_load, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, unsupported=:lower))
    @test length(rated_no_load_build.constraints[:transformer_current]) == 1
    to_current = only(values(rated_no_load_build.constraints[:transformer_current]))
    to_current_function = JuMP.constraint_object(to_current).func
    to_voltage = rated_no_load_build.variables[:w][("l", "a")]
    @test !iszero(JuMP.coefficient(to_current_function[3], to_voltage))

    # Leakage must not move a winding current limit to the synthetic internal
    # endpoint. Pin the cone's live-voltage coordinate directly.
    endpoint = transformer_case(extra=Dict{String,Any}(
        "v_nom_from" => 230.0, "v_nom_to" => 230.0,
        "r_series_to" => 0.2, "x_series_to" => 0.1,
        "i_max_to" => [45.0]))
    for per_unit in (false, true)
        build = build_l3f_opf(endpoint, Clarabel.Optimizer;
            options=L3FOptions(unsupported=:lower, per_unit=per_unit))
        cone = JuMP.constraint_object(build.constraints[:transformer_current][
            ("single_phase", "t", "to", 1)]).func
        external_w = build.variables[:w][("l", "a")]
        internal_w = build.variables[:w][("_l3f_xfmr_t_to", "a")]
        @test JuMP.coefficient(cone[1], external_w) > 0.0
        @test iszero(JuMP.coefficient(cone[1], internal_w))
    end

    # The exciting shunt remains at the original winding-2 bus. Its endpoint
    # current expression follows BMOPFTools' terminal-current convention even
    # when a to-side leakage line lies between that bus and the ideal core.
    shunt_endpoint = transformer_case(extra=Dict{String,Any}(
        "v_nom_from" => 230.0, "v_nom_to" => 230.0,
        "r_series_to" => 0.2, "x_series_to" => 0.1,
        "g_no_load" => 0.02, "i_max_to" => [47.5]))
    for per_unit in (false, true)
        build = build_l3f_opf(shunt_endpoint, Clarabel.Optimizer;
            options=L3FOptions(unsupported=:lower, per_unit=per_unit))
        JuMP.set_silent(build.model); JuMP.optimize!(build.model)
        @test JuMP.is_solved_and_feasible(build.model)
        cone = JuMP.constraint_object(build.constraints[:transformer_current][
            ("single_phase", "t", "to", 1)]).func
        scale = per_unit ? build.options.s_base : 1.0
        @test JuMP.value(cone[3]) * scale ≈ -4_000.0 atol=1e-3
        @test JuMP.value(cone[4]) * scale ≈ -800.0 atol=1e-3
    end

    # And under :reject the leakage is still refused rather than dropped.
    @test _l3f_low_has(check_l3f_applicability(leaky),
                       "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED")
end

@testset "LinDist3Flow canonical lowering: coupled center tap" begin
    function center_tap_case(; extra=Dict{String,Any}())
        Dict{String,Any}(
            "bus" => Dict(
                "hv" => Dict{String,Any}("terminal_names" => ["h"]),
                "lv" => Dict{String,Any}("terminal_names" => ["x1", "x2"])),
            "voltage_source" => Dict("source" => Dict{String,Any}(
                "bus" => "hv", "terminal_map" => ["h"],
                "configuration" => "SINGLE_PHASE",
                "v_magnitude" => [2400.0], "v_angle" => [0.0])),
            "transformer" => Dict("center_tap" => Dict("ct" => merge(Dict{String,Any}(
                "bus_from" => "hv", "bus_to" => "lv",
                "terminal_map_from" => ["h"],
                "terminal_map_to" => ["x1", "x2"],
                "v_nom_from" => 2400.0, "v_nom_to" => 120.0,
                "s_rating" => 25_000.0), extra))),
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

    leaky = center_tap_case(extra=Dict{String,Any}(
        "r_series_from" => 0.4, "x_series_from" => 0.8,
        "r_series_to" => 0.002, "x_series_to" => 0.004,
        "g_no_load" => 2.0e-5, "b_no_load" => -5.0e-5))
    strict = check_l3f_applicability(leaky)
    @test _l3f_low_has(strict, "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED")

    report = check_l3f_applicability(leaky;
        options=L3FOptions(unsupported=:lower))
    @test is_l3f_applicable(report)
    @test count(==("L.L3F.TRANSFORMER_LEAKAGE_LOWERED"),
                _l3f_low_codes(report)) == 2
    @test _l3f_low_has(report, "L.L3F.TRANSFORMER_NO_LOAD_LOWERED")

    build = build_l3f_opf(leaky, Clarabel.Optimizer;
        options=_l3f_low_opts(unsupported=:lower))
    shunt = build.network["shunt"]["_l3f_noload_ct"]
    @test shunt["terminal_map"] == ["x1", "x2"]
    @test shunt["G_1_1"] == 2.0e-5
    @test shunt["B_1_1"] == -5.0e-5
    @test !haskey(shunt, "G_2_2") && !haskey(shunt, "B_2_2")

    result = _l3f_low_solve(leaky; unsupported=:lower)
    @test result.solve.optimal
    @test haskey(result.buses, "_l3f_xfmr_ct_from")
    @test haskey(result.buses, "_l3f_xfmr_ct_to")
    @test haskey(result.lines, "_l3f_leakage_ct_from")
    @test haskey(result.lines, "_l3f_leakage_ct_to")

    # The primary star arm carries the aggregate of both legs. The two
    # secondary arms then carry their individual powers. These equalities pin
    # the coupled three-winding LinDistFlow construction directly.
    primary = result.lines["_l3f_leakage_ct_from"]
    secondary = result.lines["_l3f_leakage_ct_to"]
    wh = result.buses["hv"]["h"]["w"]
    wc = result.buses["_l3f_xfmr_ct_from"]["h"]["w"]
    @test wc ≈ wh - 2 * (0.4 * primary["p"][1] + 0.8 * primary["q"][1]) atol=1e-4
    for k in 1:2
        wint = result.buses["_l3f_xfmr_ct_to"]["x$k"]["w"]
        wout = result.buses["lv"]["x$k"]["w"]
        @test wint ≈ wc / 20.0^2 atol=1e-6
        @test wout ≈ wint - 2 * (0.002 * secondary["p"][k] +
                                0.004 * secondary["q"][k]) atol=1e-5
    end
    @test result.buses["lv"]["x1"]["vm"] < result.buses["lv"]["x2"]["vm"]

    rated = deepcopy(leaky)
    tx = rated["transformer"]["center_tap"]["ct"]
    tx["i_max_from"], tx["i_max_to"] = [10.0], [100.0, 100.0]
    for per_unit in (false, true)
        b = build_l3f_opf(rated, Clarabel.Optimizer;
            options=L3FOptions(unsupported=:lower, per_unit=per_unit))
        for (side, phi, physical, internal) in (
                ("from", 1, ("hv", "h"), ("_l3f_xfmr_ct_from", "h")),
                ("to", 1, ("lv", "x1"), ("_l3f_xfmr_ct_to", "x1")),
                ("to", 2, ("lv", "x2"), ("_l3f_xfmr_ct_to", "x2")))
            cone = JuMP.constraint_object(b.constraints[:transformer_current][
                ("center_tap", "ct", side, phi)]).func
            @test JuMP.coefficient(cone[1], b.variables[:w][physical]) > 0.0
            @test iszero(JuMP.coefficient(cone[1], b.variables[:w][internal]))
        end
    end
end

@testset "LinDist3Flow projection: load laws" begin
    # The ZP tangent matches a voltage-exponent law in value and slope at v_nom.
    # Exponent 0 and 2 coincide with native classes; 1 is the half-to-Z, half-to-P
    # split of a constant-current term.
    v_nom = 230.0
    for (gamma, expected_z, expected_p) in ((0.0, 0.0, 1.0), (1.0, 0.5, 0.5),
                                            (2.0, 1.0, 0.0), (1.4, 0.7, 0.3))
        net = _l3f_low_case()
        merge!(net["load"]["d"], Dict{String,Any}("model" => "exponential",
            "v_nom" => [v_nom], "gamma_p" => [gamma], "gamma_q" => [gamma]))
        report = check_l3f_applicability(net; options=L3FOptions(unsupported=:approximate))
        @test is_l3f_applicable(report)
        finding = _l3f_low_finding(report, "A.L3F.LOAD_LAW_PROJECTED")
        @test finding.severity == :warning

        # Tangency: at nominal voltage the two laws agree, and their
        # first derivatives in (V/v_nom) agree.
        f_true(u) = u^gamma
        f_zp(u) = expected_p + expected_z * u^2
        @test f_zp(1.0) ≈ f_true(1.0)
        @test 2 * expected_z ≈ gamma
        for u in (0.9, 1.1)
            @test abs(f_zp(u) - f_true(u)) < 0.02      # second order in (u-1)
        end
    end

    # At the fixed source voltage, arbitrary finite exponential tangents retain
    # nominal P and Q in the actual solved balance, including gamma outside the
    # common [0,2] range and in both coordinate systems.
    for gamma in (-1.0, 0.0, 1.0, 1.4, 2.0, 4.0), per_unit in (false, true)
        at_source = _l3f_low_case()
        empty!(at_source["line"]); empty!(at_source["linecode"])
        delete!(at_source["bus"], "m"); delete!(at_source["bus"], "l")
        load = at_source["load"]["d"]
        merge!(load, Dict{String,Any}("bus" => "s", "model" => "exponential",
            "v_nom" => [230.0], "gamma_p" => gamma, "gamma_q" => 2 - gamma))
        result = _l3f_low_solve(at_source; unsupported=:approximate, per_unit)
        @test result.solve.optimal
        @test result.sources["v"]["pg"][1] ≈ 10_000.0 atol=1e-3
        @test result.sources["v"]["qg"][1] ≈ 2_000.0 atol=1e-3
    end

    # A ZIP current fraction is split evenly and the total fractions survive.
    net = _l3f_low_case()
    merge!(net["load"]["d"], Dict{String,Any}("model" => "zip", "v_nom" => [v_nom],
        "alpha_z" => [0.2], "alpha_i" => [0.5], "alpha_p" => [0.3],
        "beta_z" => [0.1], "beta_i" => [0.4], "beta_p" => [0.5]))
    @test _l3f_low_has(check_l3f_applicability(net), "E.L3F.ZIP_CURRENT_UNSUPPORTED")
    @test _l3f_low_has(check_l3f_applicability(net;
        options=L3FOptions(unsupported=:lower)), "E.L3F.ZIP_CURRENT_UNSUPPORTED")
    report = check_l3f_applicability(net; options=L3FOptions(unsupported=:approximate))
    @test is_l3f_applicable(report)
    @test _l3f_low_has(report, "A.L3F.LOAD_LAW_PROJECTED")

    # Solved, the projected load equals the hand-written ZP load it becomes.
    projected = deepcopy(net)
    merge!(projected["load"]["d"], Dict{String,Any}("model" => "zip",
        "alpha_z" => [0.45], "alpha_i" => [0.0], "alpha_p" => [0.55],
        "beta_z" => [0.3], "beta_i" => [0.0], "beta_p" => [0.7]))
    _l3f_low_agree(_l3f_low_solve(net; unsupported=:approximate),
                   _l3f_low_solve(projected))

    # Omitted P and Q coefficient families each retain their independent
    # constant-power default. Projecting alpha_i must not erase q_nom, and the
    # symmetric beta_i case must not erase p_nom.
    for active_current in (true, false), per_unit in (false, true)
        missing = _l3f_low_case()
        merge!(missing["load"]["d"], Dict{String,Any}(
            "model" => "zip", "v_nom" => [v_nom],
            (active_current ? "alpha_i" : "beta_i") => [1.0]))
        result = _l3f_low_solve(missing; unsupported=:approximate, per_unit)
        @test result.solve.optimal
        @test result.sources["v"][active_current ? "qg" : "pg"][1] ≈
              (active_current ? 2_000.0 : 10_000.0) atol=1e-3
    end

    # A constant-current load becomes the 0.5/0.5 tangent.
    current = _l3f_low_case()
    merge!(current["load"]["d"],
        Dict{String,Any}("model" => "constant_current", "v_nom" => [v_nom]))
    @test _l3f_low_has(check_l3f_applicability(current), "E.L3F.LOAD_MODEL_UNSUPPORTED")
    half = _l3f_low_case()
    merge!(half["load"]["d"], Dict{String,Any}("model" => "zip", "v_nom" => [v_nom],
        "alpha_z" => [0.5], "alpha_p" => [0.5], "beta_z" => [0.5], "beta_p" => [0.5]))
    _l3f_low_agree(_l3f_low_solve(current; unsupported=:approximate),
                   _l3f_low_solve(half))
end

@testset "LinDist3Flow projection: taps and bound rejection" begin
    net = _l3f_low_case(bridge=false)
    net["transformer"] = Dict("single_phase_autotransformer" => Dict(
        "r" => Dict{String,Any}(
            "bus_from" => "m", "bus_to" => "l",
            "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
            "tap_ratio" => 1.04, "regulator_type" => "B",
            "tap_ratio_min" => 0.95, "tap_ratio_max" => 1.10)))

    @test _l3f_low_has(check_l3f_applicability(net), "E.L3F.ADJUSTABLE_TAP_UNSUPPORTED")
    @test _l3f_low_has(check_l3f_applicability(net;
        options=L3FOptions(unsupported=:lower)), "E.L3F.ADJUSTABLE_TAP_UNSUPPORTED")

    report = check_l3f_applicability(net; options=L3FOptions(unsupported=:approximate))
    @test is_l3f_applicable(report)
    finding = _l3f_low_finding(report, "A.L3F.ADJUSTABLE_TAP_PROJECTED")
    @test finding.severity == :warning
    @test finding.evidence["fixed"] == [1.04]     # the declared tap lies inside

    # Solved, it equals the same regulator declared at that fixed tap.
    fixed = deepcopy(net)
    delete!(fixed["transformer"]["single_phase_autotransformer"]["r"], "tap_ratio_min")
    delete!(fixed["transformer"]["single_phase_autotransformer"]["r"], "tap_ratio_max")
    _l3f_low_agree(_l3f_low_solve(net; unsupported=:approximate), _l3f_low_solve(fixed))
    @test _l3f_low_solve(net; unsupported=:approximate).buses["l"]["a"]["vm"] ≈
        _l3f_low_solve(fixed).buses["l"]["a"]["vm"] rtol=1e-9

    # With no declared tap inside the interval, the midpoint is used and said so.
    midpoint = deepcopy(net)
    midpoint["transformer"]["single_phase_autotransformer"]["r"]["tap_ratio"] = 1.5
    mid_report = check_l3f_applicability(midpoint;
        options=L3FOptions(unsupported=:approximate))
    @test _l3f_low_finding(mid_report,
        "A.L3F.ADJUSTABLE_TAP_PROJECTED").evidence["fixed"] == [1.025]

    for (lo, hi) in ((1.1, 0.9), ("bad", 1.1), (-0.9, 1.1))
        invalid = deepcopy(net)
        tx = invalid["transformer"]["single_phase_autotransformer"]["r"]
        tx["tap_ratio_min"], tx["tap_ratio_max"] = lo, hi
        report = check_l3f_applicability(invalid;
            options=L3FOptions(unsupported=:approximate))
        @test !is_l3f_applicable(report)
        @test _l3f_low_has(report, "E.L3F.TAP_INTERVAL_INVALID")
    end
    fixed_conflict = deepcopy(net)
    tx = fixed_conflict["transformer"]["single_phase_autotransformer"]["r"]
    tx["tap_ratio"], tx["tap_ratio_min"], tx["tap_ratio_max"] = 1.05, 1.0, 1.0
    for policy in (:reject, :lower, :approximate)
        @test _l3f_low_has(check_l3f_applicability(fixed_conflict;
            options=L3FOptions(unsupported=policy)), "E.L3F.TAP_INTERVAL_INVALID")
    end
    fixed_only = deepcopy(net)
    tx = fixed_only["transformer"]["single_phase_autotransformer"]["r"]
    delete!(tx, "tap_ratio")
    tx["tap_ratio_min"] = tx["tap_ratio_max"] = 1.05
    fixed_only_report = check_l3f_applicability(fixed_only)
    @test is_l3f_applicable(fixed_only_report)
    fixed_declared = deepcopy(fixed_only)
    tx2 = fixed_declared["transformer"]["single_phase_autotransformer"]["r"]
    tx2["tap_ratio"] = 1.05
    _l3f_low_agree(_l3f_low_solve(fixed_only), _l3f_low_solve(fixed_declared))

    # Bus limits the formulation does not assess remain errors under every
    # policy; approximate mode never silently relaxes an engineering bound.
    bus_limits = Dict{String,Any}(
        "vpn_min" => [200.0], "vpn_max" => [250.0],
        "vn_max" => 250.0,
        "vpos_min" => 0.9, "vpos_max" => 1.1,
        "vuf_max" => 0.02, "vneg_max" => 0.02,
        "vzero_max" => 0.02, "va_diff_min" => -0.1,
        "va_diff_max" => 0.1)
    for (field, value) in bus_limits, policy in (:reject, :lower, :approximate)
        limited = _l3f_low_case()
        limited["bus"]["l"][field] = value
        report = check_l3f_applicability(limited;
            options=L3FOptions(unsupported=policy))
        @test !is_l3f_applicable(report)
        @test any(f -> f.code == "E.L3F.LIMIT_UNSUPPORTED" &&
                       occursin(field, f.message), report.findings)
        @test !any(f -> startswith(f.code, "A.L3F.BUS_"), report.findings)
    end
    for field in ("va_diff_min", "va_diff_max"), policy in (:reject, :lower, :approximate)
        limited = _l3f_low_case()
        limited["line"]["l2"][field] = field == "va_diff_min" ? -0.1 : 0.1
        report = check_l3f_applicability(limited;
            options=L3FOptions(unsupported=policy))
        @test !is_l3f_applicable(report)
        @test any(f -> f.code == "E.L3F.LIMIT_UNSUPPORTED" &&
                       occursin(field, f.message), report.findings)
    end
    nominal = _l3f_low_case()
    nominal["bus"]["l"]["va_nom"] = [0.0]
    @test !_l3f_low_has(check_l3f_applicability(nominal), "E.L3F.LIMIT_UNSUPPORTED")
end

@testset "LinDist3Flow projection is flagged in the replay" begin
    exact = _l3f_low_case()
    exact["capacitor"] = Dict("c" => Dict{String,Any}(
        "bus" => "l", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
        "q_rated" => [3_000.0], "v_nom" => 230.0))
    result = solve_l3f_opf(exact, Ipopt.Optimizer;
        options=L3FOptions(unsupported=:lower, objective=:feasibility,
                           validate_nonlinear=true),
        solver_options=("print_level" => 0,))
    # A canonical lowering leaves the replay meaningful for the supported model.
    @test result.validation["replayed_network"] == "as_supplied"

    projected = _l3f_low_case()
    merge!(projected["load"]["d"], Dict{String,Any}(
        "model" => "constant_current", "v_nom" => [230.0]))
    approximate = solve_l3f_opf(projected, Ipopt.Optimizer;
        options=L3FOptions(unsupported=:approximate, objective=:feasibility,
                           validate_nonlinear=true),
        solver_options=("print_level" => 0,))
    # Here both the model and the replay describe the substituted load, so the
    # reported error does not include the projection error. The flag says so.
    @test approximate.validation["replayed_network"] == "projected"
    @test approximate.validation["status"] == "replayed"
end

@testset "LinDist3Flow lowering leaves genuine obstacles alone" begin
    # Projection must never invent physics. A meshed island, DC subsystem, and
    # active IBR stay errors at every policy level because there is no
    # defensible substitution. Center taps left this list once their coupled
    # one-to-two map and canonical three-winding-star lowering were added.
    for (name, mutate!) in (
        "meshed" => net -> (net["line"]["l3"] = Dict{String,Any}(
                "bus_from" => "l", "bus_to" => "s", "terminal_map_from" => ["a"],
                "terminal_map_to" => ["a"], "linecode" => "lc")),
        "dc" => net -> (net["dc_bus"] = Dict("d" => Dict{String,Any}())),
        "ibr" => net -> (net["ibr"] = Dict("i" => Dict{String,Any}())),
    )
        net = _l3f_low_case(); mutate!(net)
        for policy in (:reject, :lower, :approximate)
            report = check_l3f_applicability(net;
                options=L3FOptions(unsupported=policy))
            @test !is_l3f_applicable(report)
        end
    end
end
