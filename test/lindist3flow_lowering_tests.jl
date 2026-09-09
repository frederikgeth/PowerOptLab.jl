using Test
using LinearAlgebra
using JuMP
using Ipopt
using Clarabel
using PowerOptLab

# `L3FOptions(unsupported=...)` widens the admissible input without widening the
# formulation. The tests below separate the two claims it makes:
#
#   :lower        the rewrite is EXACT. Each case is solved twice — once from the
#                 richer component, once from a hand-written equivalent built
#                 only from supported components — and the two must agree to
#                 solver tolerance. That is the whole content of "exact".
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

"""Assert two solved networks describe the same physics at every terminal."""
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

    # An exact rewrite is :info, never :warning — it makes no accuracy claim.
    @test _l3f_low_finding(lowered, "L.L3F.SWITCH_LOWERED").severity == :info
    @test all(f -> f.severity != :error, lowered.findings)

    # The policy and its effect reach the published result contract.
    result = _l3f_low_solve(net; unsupported=:lower)
    @test result.formulation["unsupported_policy"] == "lower"
    @test result.formulation["lowered"] == true
    @test _l3f_low_solve(_l3f_low_case()).formulation["lowered"] == false

    @test_throws ArgumentError L3FOptions(unsupported=:whatever)
end

@testset "LinDist3Flow exact lowering: switches" begin
    switched = _l3f_low_case(bridge=false)
    switched["switch"] = Dict("sw" => Dict{String,Any}(
        "bus_from" => "m", "bus_to" => "l",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
        "open_switch" => false))
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

    # An open switch carries no current, so its far side is a separate island
    # and is reported as such rather than being quietly energized.
    opened = deepcopy(switched)
    opened["switch"]["sw"]["open_switch"] = true
    report = check_l3f_applicability(opened; options=L3FOptions(unsupported=:lower))
    @test _l3f_low_has(report, "L.L3F.SWITCH_OPEN_REMOVED")
    @test _l3f_low_has(report, "E.L3F.SOURCE_MISSING")
    @test !is_l3f_applicable(report)
end

@testset "LinDist3Flow exact lowering: capacitors" begin
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

@testset "LinDist3Flow exact lowering: line shunts" begin
    shunted = _l3f_low_case()
    merge!(shunted["line"]["l2"], Dict{String,Any}(
        "G_from_1_1" => 2e-4, "B_from_1_1" => 3e-4,
        "G_to_1_1" => 1e-4, "B_to_1_1" => 5e-4))

    # Hand-written pi equivalent: each declared half on its own bus.
    manual = _l3f_low_case()
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
end

@testset "LinDist3Flow exact lowering: transformer impedance" begin
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

    # And under :reject the leakage is still refused rather than dropped.
    @test _l3f_low_has(check_l3f_applicability(leaky),
                       "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED")
end

@testset "LinDist3Flow projection: load laws" begin
    # The ZP tangent matches a voltage-exponent law in value and slope at v_nom.
    # Exponent 0 and 2 are reproduced exactly; 1 is the half-to-Z, half-to-P
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

        # Tangency: at exactly nominal voltage the two laws agree, and their
        # first derivatives in (V/v_nom) agree.
        f_true(u) = u^gamma
        f_zp(u) = expected_p + expected_z * u^2
        @test f_zp(1.0) ≈ f_true(1.0)
        @test 2 * expected_z ≈ gamma
        for u in (0.9, 1.1)
            @test abs(f_zp(u) - f_true(u)) < 0.02      # second order in (u-1)
        end
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

@testset "LinDist3Flow projection: taps and dropped limits" begin
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

    # Bus limits the formulation does not assess are dropped, and the finding
    # says the solved problem is a relaxation.
    limited = _l3f_low_case()
    limited["bus"]["l"]["vpp_min"] = 200.0
    limited["bus"]["l"]["vneg_max"] = 0.02
    @test _l3f_low_has(check_l3f_applicability(limited), "E.L3F.LIMIT_UNSUPPORTED")
    dropped = check_l3f_applicability(limited; options=L3FOptions(unsupported=:approximate))
    @test is_l3f_applicable(dropped)
    evidence = _l3f_low_finding(dropped, "A.L3F.BUS_LIMIT_DROPPED").evidence
    @test sort(evidence["limits"]) == ["vneg_max", "vpp_min"]
end

@testset "LinDist3Flow projection is flagged in the replay" begin
    exact = _l3f_low_case()
    exact["capacitor"] = Dict("c" => Dict{String,Any}(
        "bus" => "l", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
        "q_rated" => [3_000.0], "v_nom" => 230.0))
    result = solve_l3f_opf(exact, Ipopt.Optimizer;
        options=L3FOptions(unsupported=:lower, objective=:feasibility),
        solver_options=("print_level" => 0,))
    # An exact lowering leaves the replay meaningful: same physics both sides.
    @test result.validation["replayed_network"] == "as_supplied"

    projected = _l3f_low_case()
    merge!(projected["load"]["d"], Dict{String,Any}(
        "model" => "constant_current", "v_nom" => [230.0]))
    approximate = solve_l3f_opf(projected, Ipopt.Optimizer;
        options=L3FOptions(unsupported=:approximate, objective=:feasibility),
        solver_options=("print_level" => 0,))
    # Here both the model and the replay describe the substituted load, so the
    # reported error does not include the projection error. The flag says so.
    @test approximate.validation["replayed_network"] == "projected"
    @test approximate.validation["status"] == "replayed"
end

@testset "LinDist3Flow lowering leaves genuine obstacles alone" begin
    # Projection must never invent physics. A transformer subtype with no
    # supported voltage map, a meshed island, and a DC subsystem stay errors at
    # every policy level, because there is no defensible substitution. (Yd/Dy
    # banks left this list once their maps were written; center_tap has not.)
    for (name, mutate!) in (
        "center_tap" => net -> (delete!(net["line"], "l2"); net["transformer"] = Dict("center_tap" => Dict(
            "t" => Dict{String,Any}("bus_from" => "m", "bus_to" => "l",
                "terminal_map_from" => ["a"], "terminal_map_to" => ["a"])))),
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
