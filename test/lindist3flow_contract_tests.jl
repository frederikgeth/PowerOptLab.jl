using Test
using LinearAlgebra
using JuMP
using Ipopt
using Clarabel
using PowerOptLab

# Contract-level coverage for the LinDist3Flow builder: every documented option,
# every published result field, and every diagnostic code. The physics oracles
# live in `lindist3flow_tests.jl`; the literature evidence in
# `lindist3flow_literature_tests.jl`.

const _L3F_FEASIBLE = L3FOptions(validate_nonlinear=false, objective=:feasibility)
_l3f_clarabel() = ("verbose" => false,)
_l3f_ipopt() = ("print_level" => 0,)

"""Minimal single-phase two-bus feeder: source -- line -- load."""
function _l3f_two_bus(; line_extra=Dict{String,Any}(), load_extra=Dict{String,Any}(),
                      source_extra=Dict{String,Any}())
    net = Dict{String,Any}(
        "bus" => Dict(
            "source" => Dict{String,Any}("terminal_names" => ["a"]),
            "load" => Dict{String,Any}("terminal_names" => ["a"],
                "v_min" => [180.0], "v_max" => [260.0])),
        "linecode" => Dict("lc" => Dict{String,Any}(
            "R_series_1_1" => 0.2, "X_series_1_1" => 0.1)),
        "line" => Dict("line" => merge(Dict{String,Any}(
            "bus_from" => "source", "bus_to" => "load",
            "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
            "linecode" => "lc"), line_extra)),
        "voltage_source" => Dict("source" => merge(Dict{String,Any}(
            "bus" => "source", "terminal_map" => ["a"],
            "configuration" => "SINGLE_PHASE",
            "v_magnitude" => [230.0], "v_angle" => [0.0], "cost" => [1.0]), source_extra)),
        "load" => Dict("load" => merge(Dict{String,Any}(
            "bus" => "load", "terminal_map" => ["a"],
            "configuration" => "SINGLE_PHASE", "model" => "constant_power",
            "p_nom" => [10_000.0], "q_nom" => [2_000.0]), load_extra)))
    net
end

"""Three-phase bus pair carrying a bank of three single-phase regulator units."""
function _l3f_wye_bank(; taps=[1.05, 1.00, 0.95], v=2400.0)
    units = Dict{String,Any}()
    for (k, phase) in enumerate(("a", "b", "c"))
        units["reg_$phase"] = Dict{String,Any}(
            "bus_from" => "source", "bus_to" => "regulated",
            "terminal_map_from" => [phase], "terminal_map_to" => [phase],
            "tap_ratio" => taps[k], "regulator_type" => "B")
    end
    Dict{String,Any}(
        "bus" => Dict(
            "source" => Dict{String,Any}("terminal_names" => ["a", "b", "c"]),
            "regulated" => Dict{String,Any}("terminal_names" => ["a", "b", "c"])),
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => "source", "terminal_map" => ["a", "b", "c"],
            "configuration" => "WYE", "v_magnitude" => fill(v, 3),
            "v_angle" => [0.0, -2pi / 3, 2pi / 3])),
        "transformer" => Dict("single_phase_autotransformer" => units),
        "load" => Dict("load" => Dict{String,Any}(
            "bus" => "regulated", "terminal_map" => ["a", "b", "c"],
            "configuration" => "WYE", "model" => "constant_power",
            "p_nom" => fill(100_000.0, 3), "q_nom" => fill(20_000.0, 3))))
end

_l3f_codes(report) = [f.code for f in report.findings]
_l3f_has(report, code) = any(==(code), _l3f_codes(report))

@testset "LinDist3Flow per-conductor radiality" begin
    # Three single-phase regulator units on one bus pair occupy disjoint
    # conductors, so the bank is radial even though the bus graph shows three
    # parallel edges. This is the IEEE 13/34/123 wye-bank topology.
    net = _l3f_wye_bank()
    report = check_l3f_applicability(net)
    @test is_l3f_applicable(report)
    result = solve_l3f_opf(net, Clarabel.Optimizer; options=_L3F_FEASIBLE,
        solver_options=_l3f_clarabel())
    @test result.solve.optimal
    for (k, phase) in enumerate(("a", "b", "c"))
        @test result.buses["regulated"][phase]["vm"] ≈
            2400.0 * [1.05, 1.00, 0.95][k] atol=1e-6
    end
    @test length(result.transformers) == 3
    @test sum(result.sources["source"]["pg"]) ≈ 300_000.0 atol=1e-4

    # Genuinely parallel conductors on the same phase remain non-radial.
    parallel = _l3f_wye_bank()
    parallel["transformer"]["single_phase_autotransformer"]["dup"] = Dict{String,Any}(
        "bus_from" => "source", "bus_to" => "regulated",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
        "tap_ratio" => 1.0, "regulator_type" => "B")
    @test _l3f_has(check_l3f_applicability(parallel), "E.L3F.TOPOLOGY_NOT_RADIAL")

    # A conductor-level cycle through a third bus is still non-radial.
    looped = _l3f_two_bus()
    looped["bus"]["mid"] = Dict{String,Any}("terminal_names" => ["a"])
    looped["line"]["l2"] = Dict{String,Any}(
        "bus_from" => "load", "bus_to" => "mid",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"], "linecode" => "lc")
    looped["line"]["l3"] = Dict{String,Any}(
        "bus_from" => "mid", "bus_to" => "source",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"], "linecode" => "lc")
    @test _l3f_has(check_l3f_applicability(looped), "E.L3F.TOPOLOGY_NOT_RADIAL")
end

@testset "LinDist3Flow ideal single-phase transformer" begin
    net = Dict{String,Any}(
        "bus" => Dict(
            "hv" => Dict{String,Any}("terminal_names" => ["a"]),
            "lv" => Dict{String,Any}("terminal_names" => ["a"])),
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => "hv", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
            "v_magnitude" => [11_000.0], "v_angle" => [0.0])),
        "transformer" => Dict("single_phase" => Dict("t" => Dict{String,Any}(
            "bus_from" => "hv", "bus_to" => "lv",
            "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
            "v_nom_from" => 11_000.0, "v_nom_to" => 230.0, "s_rating" => 1e5))),
        "load" => Dict("load" => Dict{String,Any}(
            "bus" => "lv", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
            "model" => "constant_power", "p_nom" => [10_000.0], "q_nom" => [2_000.0])))

    si = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:feasibility, per_unit=false),
        solver_options=_l3f_clarabel())
    pu = solve_l3f_opf(net, Clarabel.Optimizer; options=_L3F_FEASIBLE,
        solver_options=_l3f_clarabel())
    @test si.solve.optimal && pu.solve.optimal
    # Ideal fixed ratio: the secondary sits exactly at the nominal turns ratio.
    @test si.buses["lv"]["a"]["vm"] ≈ 230.0 atol=1e-6
    @test pu.buses["lv"]["a"]["vm"] ≈ si.buses["lv"]["a"]["vm"] rtol=1e-9
    @test si.transformers["t"]["effective_ratio_from_to"] ≈ 11_000.0 / 230.0
    @test si.transformers["t"]["p"] ≈ [10_000.0] atol=1e-4
    @test si.transformers["t"]["reversed_from_input"] == false

    # An off-nominal `tap` multiplies the nominal ratio (BMOPF `N = N0 * tap`).
    tapped = deepcopy(net)
    tapped["transformer"]["single_phase"]["t"]["tap"] = 1.025
    tapped_result = solve_l3f_opf(tapped, Clarabel.Optimizer; options=_L3F_FEASIBLE,
        solver_options=_l3f_clarabel())
    @test tapped_result.buses["lv"]["a"]["vm"] ≈ 230.0 / 1.025 atol=1e-6

    # Leakage and no-load admittance are outside the ideal contract.
    nonideal = deepcopy(net)
    nonideal["transformer"]["single_phase"]["t"]["r_series_from"] = 0.5
    @test _l3f_has(check_l3f_applicability(nonideal),
                   "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED")

    # An adjustable `tap` interval is rejected as a non-affine decision.
    adjustable = deepcopy(net)
    merge!(adjustable["transformer"]["single_phase"]["t"],
        Dict{String,Any}("tap_min" => 0.95, "tap_max" => 1.05))
    @test _l3f_has(check_l3f_applicability(adjustable),
                   "E.L3F.ADJUSTABLE_TAP_UNSUPPORTED")
end

@testset "LinDist3Flow lateral transformers below a three-phase trunk" begin
    # A phase-b transformer feeds a one-phase lateral after an unbalanced
    # three-phase section.  This is deliberately narrower than the three-unit
    # bank in the IEEE 13-bus regression: it pins the 3φ -> 1φ transition,
    # conductor-specific power balance, fixed turns ratio, and SI/pu agreement.
    net = Dict{String,Any}(
        "bus" => Dict(
            "source" => Dict{String,Any}("terminal_names" => ["a", "b", "c"]),
            "trunk" => Dict{String,Any}("terminal_names" => ["a", "b", "c"]),
            "lv" => Dict{String,Any}("terminal_names" => ["x"])),
        "linecode" => Dict("trunk" => Dict{String,Any}(
            "R_series_1_1" => 0.02, "X_series_1_1" => 0.01,
            "R_series_2_2" => 0.02, "X_series_2_2" => 0.01,
            "R_series_3_3" => 0.02, "X_series_3_3" => 0.01)),
        "line" => Dict("trunk" => Dict{String,Any}(
            "bus_from" => "source", "bus_to" => "trunk",
            "terminal_map_from" => ["a", "b", "c"],
            "terminal_map_to" => ["a", "b", "c"], "linecode" => "trunk")),
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => "source", "terminal_map" => ["a", "b", "c"],
            "configuration" => "WYE", "v_magnitude" => fill(2300.0, 3),
            "v_angle" => [0.0, -2pi / 3, 2pi / 3])),
        "transformer" => Dict("single_phase" => Dict("phase_b" => Dict{String,Any}(
            "bus_from" => "trunk", "bus_to" => "lv",
            "terminal_map_from" => ["b"], "terminal_map_to" => ["x"],
            "v_nom_from" => 2300.0, "v_nom_to" => 230.0, "s_rating" => 20_000.0))),
        "load" => Dict(
            "unbalanced_trunk" => Dict{String,Any}(
                "bus" => "trunk", "terminal_map" => ["a", "b", "c"],
                "configuration" => "WYE", "model" => "constant_power",
                "p_nom" => [1_000.0, 2_000.0, 3_000.0],
                "q_nom" => [100.0, 200.0, 300.0]),
            "lateral" => Dict{String,Any}(
                "bus" => "lv", "terminal_map" => ["x"],
                "configuration" => "SINGLE_PHASE", "model" => "constant_power",
                "p_nom" => [4_000.0], "q_nom" => [500.0])))

    si = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:feasibility,
                           per_unit=false),
        solver_options=_l3f_clarabel())
    pu = solve_l3f_opf(net, Clarabel.Optimizer;
        options=_L3F_FEASIBLE, solver_options=_l3f_clarabel())
    @test si.solve.optimal && pu.solve.optimal
    @test si.sources["source"]["pg"] ≈ [1_000.0, 6_000.0, 3_000.0] atol=1e-4
    @test si.sources["source"]["qg"] ≈ [100.0, 700.0, 300.0] atol=1e-4
    @test si.transformers["phase_b"]["p"] ≈ [4_000.0] atol=1e-4
    @test si.transformers["phase_b"]["q"] ≈ [500.0] atol=1e-4
    @test si.buses["lv"]["x"]["vm"] ≈
        si.buses["trunk"]["b"]["vm"] / 10.0 rtol=1e-9
    @test pu.sources["source"]["pg"] ≈ si.sources["source"]["pg"] rtol=1e-9
    @test pu.buses["lv"]["x"]["vm"] ≈ si.buses["lv"]["x"]["vm"] rtol=1e-9

    # A BMOPF center-tap transformer is not two independent single-phase
    # transformers: its one primary drives two oppositely oriented secondary
    # half-windings. Exercise that rectangular map with unequal leg loads below
    # the same 3φ trunk.
    split = deepcopy(net)
    split["bus"]["lv"]["terminal_names"] = ["x1", "x2"]
    split["transformer"] = Dict("center_tap" => Dict("split_phase" => Dict{String,Any}(
        "bus_from" => "trunk", "bus_to" => "lv",
        "terminal_map_from" => ["b"], "terminal_map_to" => ["x1", "x2"],
        "v_nom_from" => 2300.0, "v_nom_to" => 115.0, "s_rating" => 20_000.0)))
    split["load"] = Dict(
        "unbalanced_trunk" => deepcopy(net["load"]["unbalanced_trunk"]),
        "leg_1" => Dict{String,Any}(
            "bus" => "lv", "terminal_map" => ["x1"],
            "configuration" => "SINGLE_PHASE", "model" => "constant_power",
            "p_nom" => [4_000.0], "q_nom" => [500.0]),
        "leg_2" => Dict{String,Any}(
            "bus" => "lv", "terminal_map" => ["x2"],
            "configuration" => "SINGLE_PHASE", "model" => "constant_power",
            "p_nom" => [1_500.0], "q_nom" => [100.0]))
    split_si = solve_l3f_opf(split, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:feasibility,
                           per_unit=false), solver_options=_l3f_clarabel())
    split_pu = solve_l3f_opf(split, Clarabel.Optimizer;
        options=_L3F_FEASIBLE, solver_options=_l3f_clarabel())
    @test split_si.solve.optimal && split_pu.solve.optimal
    @test split_si.transformers["split_phase"]["p"] ≈ [4_000.0, 1_500.0] atol=1e-4
    @test split_si.transformers["split_phase"]["q"] ≈ [500.0, 100.0] atol=1e-4
    @test split_si.sources["source"]["pg"] ≈ [1_000.0, 7_500.0, 3_000.0] atol=1e-4
    @test split_si.buses["lv"]["x1"]["vm"] ≈
        split_si.buses["trunk"]["b"]["vm"] / 20.0 rtol=1e-9
    @test split_si.buses["lv"]["x2"]["vm"] ≈
        split_si.buses["lv"]["x1"]["vm"] rtol=1e-9
    @test mod(split_si.buses["lv"]["x1"]["reference_angle"] -
              split_si.buses["lv"]["x2"]["reference_angle"], 2pi) ≈ pi atol=1e-12
    @test split_pu.transformers["split_phase"]["p"] ≈
        split_si.transformers["split_phase"]["p"] rtol=1e-9
    @test split_pu.buses["lv"]["x1"]["vm"] ≈
        split_si.buses["lv"]["x1"]["vm"] rtol=1e-9

    # A 240 V leg-to-leg channel uses the same fixed 180-degree reference and
    # splits its terminal power equally between the two hot legs.
    across = deepcopy(split)
    across["load"]["across_legs"] = Dict{String,Any}(
        "bus" => "lv", "terminal_map" => ["x1", "x2"],
        "configuration" => "SINGLE_PHASE", "model" => "constant_power",
        "p_nom" => [2_300.0], "q_nom" => [200.0])
    across_result = solve_l3f_opf(across, Clarabel.Optimizer;
        options=_L3F_FEASIBLE, solver_options=_l3f_clarabel())
    @test across_result.solve.optimal
    @test across_result.transformers["split_phase"]["p"] ≈
        [5_150.0, 2_650.0] atol=1e-4
    @test across_result.transformers["split_phase"]["q"] ≈
        [600.0, 200.0] atol=1e-4

    # The scalar nameplate constrains aggregate primary VA. Current ratings
    # bind the aggregate primary and each retained secondary leg separately.
    rated = deepcopy(split)
    merge!(rated["transformer"]["center_tap"]["split_phase"], Dict{String,Any}(
        "i_max_from" => [10.0], "i_max_to" => [100.0, 100.0]))
    rated_build = build_l3f_opf(rated, Clarabel.Optimizer;
        options=_L3F_FEASIBLE)
    @test length(rated_build.constraints[:transformer_apparent_power]) == 1
    @test length(rated_build.constraints[:transformer_current]) == 3
    tight_nameplate = deepcopy(rated)
    tight_nameplate["transformer"]["center_tap"]["split_phase"]["s_rating"] = 5_000.0
    @test !solve_l3f_opf(tight_nameplate, Clarabel.Optimizer;
        options=_L3F_FEASIBLE, solver_options=_l3f_clarabel()).solve.optimal
    tight_primary_current = deepcopy(rated)
    tight_primary_current["transformer"]["center_tap"]["split_phase"]["i_max_from"] = [2.0]
    @test !solve_l3f_opf(tight_primary_current, Clarabel.Optimizer;
        options=_L3F_FEASIBLE, solver_options=_l3f_clarabel()).solve.optimal
    tight_leg_current = deepcopy(rated)
    tight_leg_current["transformer"]["center_tap"]["split_phase"]["i_max_to"] = [20.0, 100.0]
    @test !solve_l3f_opf(tight_leg_current, Clarabel.Optimizer;
        options=_L3F_FEASIBLE, solver_options=_l3f_clarabel()).solve.optimal

    # Reverse traversal is underdetermined as a one-to-two voltage map and is
    # therefore a typed applicability error, independent of power-flow sign.
    reversed = deepcopy(split)
    reversed["voltage_source"]["source"]["bus"] = "lv"
    reversed["voltage_source"]["source"]["terminal_map"] = ["x1", "x2"]
    reversed["voltage_source"]["source"]["configuration"] = "WYE"
    reversed["voltage_source"]["source"]["v_magnitude"] = [115.0, 115.0]
    reversed["voltage_source"]["source"]["v_angle"] = [0.0, pi]
    @test _l3f_has(check_l3f_applicability(reversed),
                   "E.L3F.CENTER_TAP_ORIENTATION_UNSUPPORTED")
end

@testset "LinDist3Flow reversed branch orientation" begin
    # BMOPF `bus_from`/`bus_to` need not point away from the source. The model
    # reorients, and the published flow is parent->child, flagged as reversed.
    net = _l3f_two_bus()
    net["line"]["line"]["bus_from"] = "load"
    net["line"]["line"]["bus_to"] = "source"
    result = solve_l3f_opf(net, Clarabel.Optimizer; options=_L3F_FEASIBLE,
        solver_options=_l3f_clarabel())
    @test result.solve.optimal
    @test result.lines["line"]["reversed_from_input"] == true
    @test result.lines["line"]["parent"] == "source"
    @test result.lines["line"]["child"] == "load"
    @test result.lines["line"]["p"] ≈ [10_000.0] atol=1e-4
    @test result.buses["load"]["a"]["w"] ≈ 48_500.0 atol=1e-3
end

@testset "LinDist3Flow reference policies" begin
    net = _l3f_two_bus()
    flat = Dict(("source", "a") => 230.0 + 0im, ("load", "a") => 230.0 + 0im)

    # Internal option derivation preserves every field except explicit
    # overrides, so adding a future field has one centralized copy path.
    original_options = L3FOptions(validate_nonlinear=false,
        reference_policy=:explicit, kron_reduce=false,
        require_neutral_provenance=true, unsupported=:approximate,
        objective=:source_import, per_unit=false, s_base=12_345.0)
    copied_options = PowerOptLab._l3f_with_options(original_options;
        reference_policy=:source_propagated)
    @test copied_options.reference_policy == :source_propagated
    for name in setdiff(collect(fieldnames(L3FOptions)), [:reference_policy])
        @test getfield(copied_options, name) == getfield(original_options, name)
    end

    auto = build_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, reference_policy=:auto),
        reference=flat)
    @test auto.reference.provenance == :explicit

    propagated = build_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, reference_policy=:auto))
    @test propagated.reference.provenance == :source_propagated
    # The default reference is the flat no-load profile: no line drop is applied.
    @test propagated.reference.voltage[("load", "a")] ≈ 230.0 + 0im

    # `:source_propagated` names the propagated profile, so a supplied reference
    # must not silently displace it.
    forced = build_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, reference_policy=:source_propagated),
        reference=Dict(("source", "a") => 230.0 + 0im, ("load", "a") => 100.0 + 0im))
    @test forced.reference.provenance == :source_propagated
    @test forced.reference.voltage[("load", "a")] ≈ 230.0 + 0im

    # `:explicit` refuses to invent one.
    @test _l3f_has(check_l3f_applicability(net;
            options=L3FOptions(reference_policy=:explicit)), "E.L3F.REFERENCE_MISSING")
    @test is_l3f_applicable(check_l3f_applicability(net;
            options=L3FOptions(reference_policy=:explicit), reference=flat))

    # A reference that contradicts the declared source is rejected, and the
    # nested bus => terminal => phasor spelling is accepted.
    @test_throws L3FInapplicableError build_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false),
        reference=Dict(("source", "a") => 100.0 + 0im, ("load", "a") => 230.0 + 0im))
    nested = build_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false),
        reference=Dict("source" => Dict("a" => Dict("vr" => 230.0, "vi" => 0.0)),
                       "load" => Dict("a" => Dict("vr" => 225.0, "vi" => 0.0))))
    @test nested.reference.voltage[("load", "a")] ≈ 225.0 + 0im

    # The provenance hash is part of the published contract and tracks content.
    result = solve_l3f_opf(net, Clarabel.Optimizer; options=_L3F_FEASIBLE,
        solver_options=_l3f_clarabel())
    @test result.formulation["reference_provenance"] == "source_propagated"
    @test length(result.formulation["reference_hash"]) == 64
    @test result.formulation["reference_hash"] == propagated.reference.source_hash
    @test nested.reference.source_hash != propagated.reference.source_hash
end

@testset "LinDist3Flow power-flow reference" begin
    net = _l3f_two_bus()
    reference = l3f_reference_from_powerflow(net; solver_options=_l3f_ipopt())
    @test reference isa L3FReferenceState
    @test reference.provenance == :power_flow
    @test reference.nonlinear_status in (:OPTIMAL, :LOCALLY_SOLVED)
    @test length(reference.source_hash) == 64

    # The helper can create an explicit reference independently of whether a
    # caller elects to require neutral-reduction provenance on a later build.
    explicit_required = l3f_reference_from_powerflow(net;
        options=L3FOptions(reference_policy=:explicit),
        solver_options=_l3f_ipopt())
    @test explicit_required.provenance == :power_flow
    explicit_build = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, reference_policy=:explicit),
        reference=explicit_required, solver_options=_l3f_clarabel())
    @test explicit_build.solve.optimal
    @test explicit_build.formulation["reference_provenance"] == "power_flow"

    # The source is fixed data and is restored exactly; the load bus carries the
    # real drop, which the flat profile by construction does not.
    @test reference.voltage[("source", "a")] ≈ 230.0 + 0im
    @test abs(reference.voltage[("load", "a")]) < 230.0
    @test abs(reference.voltage[("load", "a")]) ≈ 220.2 atol=0.5

    flat = build_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false))
    @test abs(flat.reference.voltage[("load", "a")]) ≈ 230.0

    # It feeds straight back in and is recorded in the published provenance.
    refined = solve_l3f_opf(net, Clarabel.Optimizer; options=_L3F_FEASIBLE,
        reference=reference, solver_options=_l3f_clarabel())
    @test refined.solve.optimal
    @test refined.formulation["reference_provenance"] == "power_flow"
    @test refined.formulation["reference_hash"] == reference.source_hash
    # Linearizing elsewhere changes the answer; both remain close to the truth.
    @test refined.buses["load"]["a"]["w"] != flat.reference.voltage[("load", "a")]
    @test abs(refined.buses["load"]["a"]["vm"] - 220.2) < 1.0

    # `:source_propagated` still refuses to be displaced, even by a real
    # power-flow reference — the policy names the profile, not the quality.
    forced = build_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, reference_policy=:source_propagated),
        reference=reference)
    @test forced.reference.provenance == :source_propagated

    # A power flow is determined, so an unpinned generator range is refused with
    # an actionable message rather than a solver failure.
    ranged = _l3f_two_bus()
    ranged["generator"] = Dict("pv" => Dict{String,Any}(
        "bus" => "load", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
        "p_min" => [0.0], "p_max" => [4_000.0], "q_min" => [0.0], "q_max" => [0.0],
        "cost" => [0.0]))
    err = try
        l3f_reference_from_powerflow(ranged; solver_options=_l3f_ipopt()); nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("dispatch", sprint(showerror, err))

    # Pinning it at a previous solution is the successive-linearization loop.
    first_pass = solve_l3f_opf(ranged, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:cost),
        solver_options=_l3f_clarabel())
    @test first_pass.generators["pv"]["pg"] ≈ [4_000.0] atol=1e-3
    looped = l3f_reference_from_powerflow(ranged; dispatch=first_pass,
        solver_options=_l3f_ipopt())
    @test looped.provenance == :power_flow
    # The generator lifts the load-bus voltage, so this reference sits above the
    # one taken without any dispatch.
    @test abs(looped.voltage[("load", "a")]) > abs(reference.voltage[("load", "a")])
    second_pass = solve_l3f_opf(ranged, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:cost),
        reference=looped, solver_options=_l3f_clarabel())
    @test second_pass.solve.optimal
    @test second_pass.formulation["reference_provenance"] == "power_flow"

    # An inapplicable network is refused before any solver runs.
    @test_throws L3FInapplicableError l3f_reference_from_powerflow(
        _l3f_two_bus(load_extra=Dict{String,Any}("configuration" => "ZIGZAG")))
end

@testset "LinDist3Flow objectives and cost hygiene" begin
    net = _l3f_two_bus()
    net["generator"] = Dict("pv" => Dict{String,Any}(
        "bus" => "load", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
        "p_min" => [0.0], "p_max" => [4_000.0], "q_min" => [0.0], "q_max" => [0.0],
        "cost" => [0.0]))

    # :source_import minimizes total source active injection and is reported in W.
    imports = [solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:source_import,
                           per_unit=pu, s_base=17_000.0),
        solver_options=_l3f_clarabel()) for pu in (false, true)]
    @test imports[1].objective ≈ imports[2].objective rtol=1e-8
    @test imports[1].objective ≈ 6_000.0 atol=1e-3
    @test all(r -> r.objective ≈ sum(r.sources["source"]["pg"]), imports)

    # :feasibility is a zero objective.
    feasible = solve_l3f_opf(net, Clarabel.Optimizer; options=_L3F_FEASIBLE,
        solver_options=_l3f_clarabel())
    @test feasible.objective == 0.0

    # :cost prices the zero-cost generator down to the source's marginal price.
    cost = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:cost),
        solver_options=_l3f_clarabel())
    @test cost.generators["pv"]["pg"] ≈ [4_000.0] atol=1e-3
    @test cost.objective ≈ 6.0 atol=1e-4   # 6 kW imported at 1 currency/kWh

    # BMOPF requires `cost`; defaulting it to zero silently would flatten the
    # objective, so a dispatchable unit without one is reported.
    missing_cost = deepcopy(net)
    delete!(missing_cost["generator"]["pv"], "cost")
    report = check_l3f_applicability(missing_cost; options=L3FOptions(objective=:cost))
    @test _l3f_has(report, "W.L3F.COST_MISSING")
    @test is_l3f_applicable(report)          # a warning never blocks the build
    # It is objective-specific: :feasibility never prices anything.
    @test !_l3f_has(check_l3f_applicability(missing_cost;
            options=L3FOptions(objective=:feasibility)), "W.L3F.COST_MISSING")
    # A non-dispatchable unit (p_min == p_max) has no dispatch direction to flatten.
    fixed = deepcopy(missing_cost)
    fixed["generator"]["pv"]["p_min"] = [4_000.0]
    @test !_l3f_has(check_l3f_applicability(fixed;
            options=L3FOptions(objective=:cost)), "W.L3F.COST_MISSING")
end

@testset "LinDist3Flow rating constraints bind" begin
    # Apparent-power and live-voltage current ratings are native second-order cones
    # on every rated element. Each is checked by making it the binding limit.
    demand = hypot(10_000.0, 2_000.0)

    for (field, radius) in ("s_max" => demand / 2, "i_max" => demand / (2 * 230.0))
        net = _l3f_two_bus(line_extra=Dict{String,Any}(field => [radius]))
        @test is_l3f_applicable(check_l3f_applicability(net))
        result = solve_l3f_opf(net, Clarabel.Optimizer; options=_L3F_FEASIBLE,
            solver_options=_l3f_clarabel())
        @test !result.solve.optimal      # the feeder cannot serve the load
        # The same rating, doubled, admits the load.
        relaxed = _l3f_two_bus(line_extra=Dict{String,Any}(field => [4 * radius]))
        @test solve_l3f_opf(relaxed, Clarabel.Optimizer; options=_L3F_FEASIBLE,
            solver_options=_l3f_clarabel()).solve.optimal
    end

    # A linecode rating is inherited when the line does not override it.
    inherited = _l3f_two_bus()
    inherited["linecode"]["lc"]["i_max"] = [demand / (2 * 230.0)]
    @test !solve_l3f_opf(inherited, Clarabel.Optimizer; options=_L3F_FEASIBLE,
        solver_options=_l3f_clarabel()).solve.optimal

    # Source ratings must bind identically in both coordinate systems. Upstream
    # per-unit preparation leaves a source's `s_max`/`i_max` in SI, so stamping
    # the working copy verbatim would silently widen the nameplate by `s_base`.
    for (field, radius) in ("s_max" => demand / 2, "i_max" => demand / (2 * 230.0))
        limited = _l3f_two_bus(source_extra=Dict{String,Any}(field => [radius]))
        for s_base in (1e6, 25_000.0)
            @test !solve_l3f_opf(limited, Clarabel.Optimizer;
                options=L3FOptions(validate_nonlinear=false, objective=:feasibility,
                                   per_unit=true, s_base=s_base),
                solver_options=_l3f_clarabel()).solve.optimal
        end
        @test !solve_l3f_opf(limited, Clarabel.Optimizer;
            options=L3FOptions(validate_nonlinear=false, objective=:feasibility,
                               per_unit=false),
            solver_options=_l3f_clarabel()).solve.optimal
        ample = _l3f_two_bus(source_extra=Dict{String,Any}(field => [4 * radius]))
        @test solve_l3f_opf(ample, Clarabel.Optimizer; options=_L3F_FEASIBLE,
            solver_options=_l3f_clarabel()).solve.optimal
    end

    # Generator apparent-power and live-voltage current cones.
    net = _l3f_two_bus()
    net["generator"] = Dict("pv" => Dict{String,Any}(
        "bus" => "load", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
        "p_min" => [0.0], "p_max" => [1e5], "q_min" => [0.0], "q_max" => [0.0],
        "s_max" => [3_000.0], "cost" => [0.0]))
    capped = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:cost),
        solver_options=_l3f_clarabel())
    @test capped.generators["pv"]["pg"] ≈ [3_000.0] atol=1e-3
    current_capped = deepcopy(net)
    delete!(current_capped["generator"]["pv"], "s_max")
    current_capped["generator"]["pv"]["i_max"] = [10.0]      # 10 A at 230 V
    result = solve_l3f_opf(current_capped, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:cost),
        solver_options=_l3f_clarabel())
    result_si = solve_l3f_opf(current_capped, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:cost,
                           per_unit=false),
        solver_options=_l3f_clarabel())
    @test result.solve.optimal
    @test result_si.solve.optimal
    # The current cone uses the solved terminal voltage, not the 230 V
    # reference; voltage drop therefore lowers the available real power.
    live_radius = 10.0 * sqrt(result.buses["load"]["a"]["w"])
    @test result.generators["pv"]["pg"] ≈ [live_radius] atol=1e-2
    @test result_si.generators["pv"]["pg"] ≈ result.generators["pv"]["pg"] atol=1e-2
    # The SI cone is the exact same inequality, reciprocally scaled about the
    # 230 V reference so Clarabel does not see axes spanning several orders of
    # magnitude: [I*w/Vref, I*Vref/2, p, q] ∈ Qr.
    current_build_si = build_l3f_opf(current_capped, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:cost,
                           per_unit=false))
    current_cone_si = JuMP.constraint_object(
        current_build_si.constraints[:generator_current][("pv", 1)]).func
    @test JuMP.coefficient(current_cone_si[1],
        current_build_si.variables[:w][("load", "a")]) ≈ 10.0 / 230.0
    @test JuMP.constant(current_cone_si[2]) ≈ 10.0 * 230.0 / 2

    # BMOPF types transformer current limits as `number[]`; both shapes are
    # accepted on a single-conductor device and mean the same thing.
    bank = _l3f_wye_bank()
    units = bank["transformer"]["single_phase_autotransformer"]
    units["reg_a"]["i_max_from"] = [20.0]
    units["reg_b"]["i_max_from"] = 20.0
    @test is_l3f_applicable(check_l3f_applicability(bank))
    build = build_l3f_opf(bank, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false))
    @test length(build.constraints[:transformer_current]) == 2
    @test all(JuMP.constraint_object(c).set isa JuMP.MOI.RotatedSecondOrderCone
              for c in values(build.constraints[:transformer_current]))
    @test !_l3f_has(check_l3f_applicability(bank), "E.L3F.LIMIT_INVALID")
    # A per-conductor vector that does not match the device arity is still invalid.
    units["reg_c"]["i_max_from"] = [20.0, 20.0]
    @test _l3f_has(check_l3f_applicability(bank), "E.L3F.LIMIT_INVALID")
end

@testset "LinDist3Flow nonlinear replay contract" begin
    net = _l3f_two_bus()
    result = solve_l3f_opf(net, Ipopt.Optimizer;
        options=L3FOptions(objective=:feasibility, validate_nonlinear=true),
        solver_options=_l3f_ipopt())
    validation = result.validation
    # `status` describes the replay, not the accuracy: the linearization omits
    # series losses, so a converged replay still differs from the linear answer.
    @test validation["status"] == "replayed"
    @test validation["nonlinear_solve_status"] in ("OPTIMAL", "LOCALLY_SOLVED")
    @test validation["compared_terminals"] == 2
    @test validation["physical_limits"] == "unassessed"
    @test isfinite(validation["maximum_voltage_magnitude_error"])
    @test validation["maximum_voltage_magnitude_error"] > 0.0
    @test !haskey(validation, "within_tolerance")
    @test result.reference.nonlinear_status == :replayed
    @test solve_diagnostics(result).validation == "replayed"

    # An explicit tolerance is what turns the replay into a judgement.
    loose = solve_l3f_opf(net, Ipopt.Optimizer;
        options=L3FOptions(objective=:feasibility, validate_nonlinear=true), solver_options=_l3f_ipopt(),
        voltage_tolerance=5.0)
    @test loose.validation["within_tolerance"] == true
    tight = solve_l3f_opf(net, Ipopt.Optimizer;
        options=L3FOptions(objective=:feasibility, validate_nonlinear=true), solver_options=_l3f_ipopt(),
        voltage_tolerance=1e-6)
    @test tight.validation["within_tolerance"] == false
    @test tight.validation["voltage_tolerance"] == 1e-6

    # The replay is skipped, not faked, when the linear solve is not optimal.
    # An unreachable voltage floor keeps the model an LP so Ipopt can take it.
    infeasible = _l3f_two_bus()
    infeasible["bus"]["load"]["v_min"] = [259.0]
    skipped = solve_l3f_opf(infeasible, Ipopt.Optimizer;
        options=L3FOptions(objective=:feasibility, validate_nonlinear=true),
        solver_options=_l3f_ipopt())
    @test !skipped.solve.optimal
    @test skipped.validation["status"] == "not_run"
    @test isnan(skipped.objective)
    @test all(isnan, skipped.lines["line"]["p"])

    # `validate_l3f_solution` operates on the result's own snapshot.
    standalone = validate_l3f_solution(result; voltage_tolerance=5.0)
    @test standalone["status"] == "replayed"
    @test standalone["within_tolerance"] == true
end

@testset "LinDist3Flow default is a one-shot affine/SOC solve" begin
    build = build_l3f_opf(_l3f_two_bus(), Clarabel.Optimizer)
    @test all(JuMP.start_value(variable) === nothing
              for variable in JuMP.all_variables(build.model))
    @test JuMP.objective_function_type(build.model) <: JuMP.GenericAffExpr
    @test l3f_model_class(build) in (:LP, :SOCP)
    result = solve_l3f_opf(_l3f_two_bus(), Clarabel.Optimizer;
        solver_options=_l3f_clarabel())
    @test result.solve.optimal
    @test result.validation["status"] == "not_requested"
    @test result.reference.provenance == :source_propagated
end

@testset "LinDist3Flow islands and source assignment" begin
    two = Dict{String,Any}(
        "bus" => Dict("s1" => Dict{String,Any}("terminal_names" => ["a"]),
                      "s2" => Dict{String,Any}("terminal_names" => ["a"])),
        "voltage_source" => Dict(
            "v1" => Dict{String,Any}("bus" => "s1", "terminal_map" => ["a"],
                "configuration" => "SINGLE_PHASE",
                "v_magnitude" => [230.0], "v_angle" => [0.0]),
            "v2" => Dict{String,Any}("bus" => "s2", "terminal_map" => ["a"],
                "configuration" => "SINGLE_PHASE",
                "v_magnitude" => [400.0], "v_angle" => [0.5])),
        "load" => Dict("l" => Dict{String,Any}("bus" => "s2", "terminal_map" => ["a"],
            "configuration" => "SINGLE_PHASE", "model" => "constant_power",
            "p_nom" => [1_000.0], "q_nom" => [0.0])))
    report = check_l3f_applicability(two)
    @test is_l3f_applicable(report)
    @test report.islands == [["s1"], ["s2"]]
    @test report.roots == ["s1", "s2"]
    result = solve_l3f_opf(two, Clarabel.Optimizer; options=_L3F_FEASIBLE,
        solver_options=_l3f_clarabel())
    @test result.buses["s1"]["a"]["vm"] ≈ 230.0 atol=1e-6
    @test result.buses["s2"]["a"]["vm"] ≈ 400.0 atol=1e-6

    unsourced = _l3f_two_bus()
    delete!(unsourced, "voltage_source")
    @test _l3f_has(check_l3f_applicability(unsourced), "E.L3F.SOURCE_MISSING")

    doubled = _l3f_two_bus()
    doubled["voltage_source"]["extra"] = Dict{String,Any}(
        "bus" => "load", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
        "v_magnitude" => [230.0], "v_angle" => [0.0])
    @test _l3f_has(check_l3f_applicability(doubled), "E.L3F.MULTIPLE_SOURCES")

    # Bus-level connectivity is insufficient: every retained conductor must be
    # reachable from a source terminal, even when an explicit phasor is supplied.
    isolated = _l3f_two_bus()
    isolated["bus"]["load"]["terminal_names"] = ["a", "b"]
    explicit = Dict(("source", "a") => 230.0 + 0im,
                    ("load", "a") => 225.0 + 0im,
                    ("load", "b") => 225.0cis(-2pi / 3))
    isolated_report = check_l3f_applicability(isolated; reference=explicit)
    @test !is_l3f_applicable(isolated_report)
    @test _l3f_has(isolated_report, "E.L3F.CONDUCTOR_UNREACHABLE")
end

@testset "LinDist3Flow neutral grounding survives reduction" begin
    # A winding grounded through a nonzero impedance is not a solidly grounded
    # one. Kron reduction deletes `r_neutral_*`/`x_neutral_*` before the
    # transformer ideality check can see them, so without an explicit
    # provenance check the same device would be accepted when supplied with an
    # explicit neutral and rejected when supplied already reduced. The verdict
    # must depend on the device, not on which form it arrived in.
    explicit(rn) = Dict{String,Any}(
        "terminal_conventions" => Dict("phase" => ["a"], "neutral" => ["n"]),
        "bus" => Dict(
            "hv" => Dict{String,Any}("terminal_names" => ["a", "n"],
                "perfectly_grounded_terminals" => ["n"]),
            "lv" => Dict{String,Any}("terminal_names" => ["a", "n"],
                "perfectly_grounded_terminals" => ["n"])),
        "transformer" => Dict("single_phase" => Dict("t" => merge(Dict{String,Any}(
            "bus_from" => "hv", "bus_to" => "lv",
            "terminal_map_from" => ["a", "n"], "terminal_map_to" => ["a", "n"],
            "v_nom_from" => 11_000.0, "v_nom_to" => 230.0, "s_rating" => 1.0e5), rn))),
        "voltage_source" => Dict("v" => Dict{String,Any}(
            "bus" => "hv", "terminal_map" => ["a", "n"], "configuration" => "WYE",
            "v_magnitude" => [11_000.0, 0.0], "v_angle" => [0.0, 0.0])),
        "load" => Dict("l" => Dict{String,Any}(
            "bus" => "lv", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
            "model" => "constant_power", "p_nom" => [1.0e4], "q_nom" => [2.0e3])))
    reduced(rn) = Dict{String,Any}(
        "bus" => Dict("hv" => Dict{String,Any}("terminal_names" => ["a"]),
                      "lv" => Dict{String,Any}("terminal_names" => ["a"])),
        "transformer" => Dict("single_phase" => Dict("t" => merge(Dict{String,Any}(
            "bus_from" => "hv", "bus_to" => "lv",
            "terminal_map_from" => ["a"], "terminal_map_to" => ["a"],
            "v_nom_from" => 11_000.0, "v_nom_to" => 230.0, "s_rating" => 1.0e5), rn))),
        "voltage_source" => Dict("v" => Dict{String,Any}(
            "bus" => "hv", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
            "v_magnitude" => [11_000.0], "v_angle" => [0.0])),
        "load" => Dict("l" => Dict{String,Any}(
            "bus" => "lv", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
            "model" => "constant_power", "p_nom" => [1.0e4], "q_nom" => [2.0e3])))

    impedance = Dict{String,Any}("r_neutral_to" => 5.0, "x_neutral_to" => 2.0)
    solid = Dict{String,Any}("r_neutral_to" => 0.0, "x_neutral_to" => 0.0)

    for policy in (:reject, :lower, :approximate)
        options = L3FOptions(unsupported=policy)
        for build in (explicit, reduced)
            report = check_l3f_applicability(build(impedance); options=options)
            @test !is_l3f_applicable(report)
            @test _l3f_has(report, "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED")
        end
        # A solidly grounded neutral really is ideal; it must not be flagged.
        for build in (explicit, reduced)
            @test is_l3f_applicable(check_l3f_applicability(build(solid); options=options))
            @test is_l3f_applicable(check_l3f_applicability(
                build(Dict{String,Any}()); options=options))
        end
    end

    # `:permissive` is the mode that accepts the idealization, and it records
    # the discarded value rather than dropping it silently.
    permissive = L3FOptions(unsupported=:permissive)
    for build in (explicit, reduced)
        report = check_l3f_applicability(build(impedance); options=permissive)
        @test is_l3f_applicable(report)
        idealized = [f for f in report.findings
                     if f.code == "A.L3F.PERMISSIVE_TRANSFORMER_IDEALIZED"]
        # One per discarded field, each carrying the value it discarded.
        @test length(idealized) == 2
        @test all(f -> f.severity == :warning, idealized)
        @test Set(f.evidence["field"] for f in idealized) ==
              Set(["r_neutral_to", "x_neutral_to"])
        @test Set(Float64(f.evidence["original_value"]) for f in idealized) ==
              Set([5.0, 2.0])
    end
end

@testset "LinDist3Flow neutral provenance and Kron boundary" begin
    # A perfectly grounded terminal is an explicit return conductor whatever it
    # is named, so it must route to Kron reduction rather than to a confusing
    # zero-reference-phasor error.
    grounded = Dict{String,Any}(
        "bus" => Dict("s" => Dict{String,Any}(
            "terminal_names" => ["a", "g"], "perfectly_grounded_terminals" => ["g"])),
        "voltage_source" => Dict("v" => Dict{String,Any}(
            "bus" => "s", "terminal_map" => ["a", "g"], "configuration" => "WYE",
            "v_magnitude" => [230.0, 0.0], "v_angle" => [0.0, 0.0])),
        "load" => Dict("l" => Dict{String,Any}(
            "bus" => "s", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
            "model" => "constant_power", "p_nom" => [1_000.0], "q_nom" => [0.0])))
    report = check_l3f_applicability(grounded)
    @test report.kron_reduced
    # `kron_reduce_bmopf` only eliminates conductors declared as neutrals, so a
    # grounded terminal under another label survives. The diagnostic must name
    # that obstacle rather than the zero reference phasor it later causes.
    @test _l3f_has(report, "E.L3F.GROUNDED_TERMINAL_RETAINED")
    @test occursin("perfectly grounded",
        only(f.message for f in report.findings
             if f.code == "E.L3F.GROUNDED_TERMINAL_RETAINED"))
    # Declaring the same conductor a neutral lets the reducer eliminate it.
    named = deepcopy(grounded)
    named["terminal_conventions"] = Dict{String,Any}("phase" => ["a"], "neutral" => ["g"])
    named_report = check_l3f_applicability(named)
    @test named_report.kron_reduced
    @test is_l3f_applicable(named_report)
    @test _l3f_has(check_l3f_applicability(grounded;
            options=L3FOptions(kron_reduce=false)), "E.L3F.EXPLICIT_NEUTRAL_UNSUPPORTED")

    # Provenance can be demanded when a case claims to be already reduced.
    reduced = _l3f_two_bus()
    strict = L3FOptions(require_neutral_provenance=true)
    @test _l3f_has(check_l3f_applicability(reduced; options=strict),
                   "E.L3F.NEUTRAL_REDUCTION_UNDECLARED")
    declared = deepcopy(reduced)
    declared["_meta"] = Dict{String,Any}("kron_reduction" => Dict{String,Any}("applied" => true))
    @test is_l3f_applicable(check_l3f_applicability(declared; options=strict))
    # A phasor reference says nothing about how an explicit neutral was removed.
    @test !is_l3f_applicable(check_l3f_applicability(reduced; options=strict,
        reference=Dict(("source", "a") => 230.0 + 0im, ("load", "a") => 230.0 + 0im)))

    forced = deepcopy(named)
    delete!(forced["bus"]["s"], "perfectly_grounded_terminals")
    strict_forced = check_l3f_applicability(forced;
        options=L3FOptions(unsupported=:lower))
    @test _l3f_has(strict_forced, "E.L3F.NEUTRAL_GROUNDING_PROJECTION")
    approximate_forced = check_l3f_applicability(forced;
        options=L3FOptions(unsupported=:approximate))
    @test is_l3f_applicable(approximate_forced)
    @test _l3f_has(approximate_forced, "A.L3F.NEUTRAL_GROUNDING_PROJECTED")
end

@testset "LinDist3Flow diagnostic code inventory" begin
    # Every code the compiler can emit must have a reachable case. A new code
    # without one fails this test rather than shipping untested.
    emitted = Set{String}()
    record!(report) = union!(emitted, _l3f_codes(report))

    base = _l3f_two_bus()
    record!(check_l3f_applicability(base; options=L3FOptions(objective=:cost)))

    for (family, payload) in ("switch" => Dict{String,Any}(),
                              "capacitor" => Dict{String,Any}(),
                              "ibr" => Dict{String,Any}(),
                              "dc_bus" => Dict{String,Any}())
        net = _l3f_two_bus(); net[family] = Dict("x" => payload)
        record!(check_l3f_applicability(net))
    end

    let net = _l3f_two_bus()
        net["control_profile"] = Dict("p" => Dict{String,Any}())
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["load"]["load"]["time_series"] = "profile"
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["bus"]["load"]["vpos_min"] = 0.9
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["bus"]["load"]["v_min"] = [300.0]      # exceeds v_max
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus(); net["bus"] = Dict{String,Any}()
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["load"]["load"]["bus"] = "nowhere"
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["load"]["load"]["terminal_map"] = ["z"]
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["load"]["load"]["configuration"] = "ZIGZAG"
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["load"]["load"]["model"] = "exponential"
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        merge!(net["load"]["load"], Dict{String,Any}("model" => "zip",
            "v_nom" => [230.0], "alpha_z" => [0.5], "alpha_i" => [0.5]))
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        merge!(net["load"]["load"], Dict{String,Any}("model" => "zip"))  # no v_nom
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["load"]["load"]["p_nom"] = [10_000.0, 1.0]
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["linecode"]["lc"]["G_from_1_1"] = 1e-6
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        merge!(net["line"]["line"], Dict{String,Any}(
            "R_series_1_1" => 0.2, "X_series_1_1" => 0.1))
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        delete!(net["line"]["line"], "linecode")
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["line"]["line"]["i_max"] = [-1.0]
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["shunt"] = Dict("s" => Dict{String,Any}(
            "bus" => "load", "terminal_map" => ["a"], "G_1_1" => 0.01, "G_2_2" => 9.0))
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["transformer"] = Dict("n_winding" => Dict("t" => Dict{String,Any}(
            "bus_from" => "source", "bus_to" => "load",
            "terminal_map_from" => ["a"], "terminal_map_to" => ["a"])))
        record!(check_l3f_applicability(net))
    end
    let net = _l3f_two_bus()
        net["voltage_source"]["source"]["configuration"] = "DELTA"
        net["voltage_source"]["source"]["v_magnitude"] = [0.0]
        record!(check_l3f_applicability(net))
    end
    record!(check_l3f_applicability(_l3f_wye_bank();
        options=L3FOptions(reference_policy=:explicit)))
    let bank = _l3f_wye_bank()
        bank["transformer"]["single_phase_autotransformer"]["reg_a"]["tap_ratio"] = 0.0
        record!(check_l3f_applicability(bank))
    end
    let bank = _l3f_wye_bank()
        bank["transformer"]["single_phase_autotransformer"]["reg_a"]["terminal_map_from"] = ["a", "b"]
        bank["transformer"]["single_phase_autotransformer"]["reg_a"]["terminal_map_to"] = ["a", "b"]
        record!(check_l3f_applicability(bank))
    end
    let net = _l3f_two_bus()
        net["generator"] = Dict("g" => Dict{String,Any}(
            "bus" => "load", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
            "p_min" => [1.0], "p_max" => [0.0], "q_min" => [0.0], "q_max" => [0.0],
            "cost" => [0.0]))
        record!(check_l3f_applicability(net))
    end

    let net = _l3f_two_bus()
        net["bus"]["source"]["terminal_names"] = ["a", "n"]
        net["bus"]["source"]["perfectly_grounded_terminals"] = ["n"]
        record!(check_l3f_applicability(net; options=L3FOptions(kron_reduce=false)))
    end
    let net = _l3f_two_bus()
        net["bus"]["source"]["perfectly_grounded_terminals"] = ["a"]
        record!(check_l3f_applicability(net))
    end
    record!(check_l3f_applicability(_l3f_two_bus();
        options=L3FOptions(require_neutral_provenance=true)))
    let bank = _l3f_wye_bank()
        merge!(bank["transformer"]["single_phase_autotransformer"]["reg_a"],
               Dict{String,Any}("tap_ratio_min" => 0.9, "tap_ratio_max" => 1.1))
        record!(check_l3f_applicability(bank))
    end
    let bank = _l3f_wye_bank()
        tx = bank["transformer"]["single_phase_autotransformer"]["reg_a"]
        tx["tap_ratio_min"], tx["tap_ratio_max"] = 1.1, 0.9
        record!(check_l3f_applicability(bank;
            options=L3FOptions(unsupported=:approximate)))
    end
    let net = _l3f_two_bus()
        net["_meta"] = Dict{String,Any}("kron_reduction" => Dict{String,Any}(
            "forced_ground_buses" => ["load"], "changes" => Any[]))
        record!(check_l3f_applicability(net))
    end

    expected = Set([
        "E.L3F.ADJUSTABLE_TAP_UNSUPPORTED", "E.L3F.BUS_UNKNOWN",
        "E.L3F.CENTER_TAP_ORIENTATION_UNSUPPORTED", "E.L3F.CONDUCTOR_UNREACHABLE",
        "E.L3F.COMPONENT_UNSUPPORTED", "E.L3F.CONNECTION_UNSUPPORTED",
        "E.L3F.CONTROL_PROFILE_UNSUPPORTED", "E.L3F.DC_SUBSYSTEM_UNSUPPORTED",
        "E.L3F.DEVICE_ARITY", "E.L3F.DEVICE_DATA_INVALID",
        "E.L3F.EMPTY_NETWORK", "E.L3F.EXPLICIT_NEUTRAL_UNSUPPORTED",
        "E.L3F.GROUNDED_TERMINAL_RETAINED",
        "E.L3F.LIMIT_INVALID", "E.L3F.LIMIT_UNSUPPORTED",
        "E.L3F.LINE_IMPEDANCE_SOURCE", "E.L3F.LINE_MATRIX_INVALID",
        "E.L3F.LINE_SHUNT_UNSUPPORTED",
        "E.L3F.LOAD_MODEL_UNSUPPORTED", "E.L3F.MULTIPLE_SOURCES",
        "E.L3F.NEUTRAL_REDUCTION_UNDECLARED", "E.L3F.NEUTRAL_GROUNDING_PROJECTION",
        "E.L3F.NEUTRAL_LIMIT_DISCARDED", "E.L3F.REFERENCE_MISSING",
        "E.L3F.REFERENCE_ZERO_WINDING", "E.L3F.SHUNT_INVALID",
        "E.L3F.SOURCE_MISSING", "E.L3F.SWITCH_UNSUPPORTED",
        "E.L3F.TERMINAL_MAP_INVALID", "E.L3F.TIME_SERIES_UNSUPPORTED",
        "E.L3F.TOPOLOGY_NOT_RADIAL", "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED",
        "E.L3F.TRANSFORMER_RATIO_INVALID", "E.L3F.TRANSFORMER_UNSUPPORTED",
        "E.L3F.TAP_INTERVAL_INVALID", "E.L3F.VOLTAGE_BOUND_INVALID",
        "E.L3F.ZIP_CURRENT_UNSUPPORTED",
        "W.L3F.COST_MISSING",
    ])
    # Codes that only the lowering pass can emit; reachable cases live in
    # lindist3flow_lowering_tests.jl, which owns that policy surface.
    lowering = Set([
        "E.L3F.CAPACITOR_INVALID", "L.L3F.SWITCH_LOWERED",
        "E.L3F.LINE_SHUNT_INVALID", "W.L3F.LINE_SHUNT_ASYMMETRIC",
        "L.L3F.SWITCH_OPEN_REMOVED", "L.L3F.CAPACITOR_LOWERED",
        "L.L3F.LINE_SHUNT_LOWERED", "L.L3F.TRANSFORMER_LEAKAGE_LOWERED",
        "L.L3F.TRANSFORMER_NO_LOAD_LOWERED", "A.L3F.LOAD_LAW_PROJECTED",
        "A.L3F.ADJUSTABLE_TAP_PROJECTED", "A.L3F.NEUTRAL_GROUNDING_PROJECTED",
    ])
    # Yd/Dy orientation codes; reachable cases live in
    # lindist3flow_delta_transformer_tests.jl.
    delta = Set(["E.L3F.DELTA_ORIENTATION_UNSUPPORTED",
                 "A.L3F.DELTA_ZERO_SEQUENCE_GAUGE",
                 "E.L3F.DELTA_DELTA_REQUIRES_APPROXIMATION",
                 "A.L3F.DELTA_DELTA_ZERO_SEQUENCE_PROJECTED",
                 "E.L3F.LOCAL_BANK_RATING_UNSUPPORTED"])
    controls = Set(["E.L3F.IBR_INVALID", "E.L3F.IBR_UNSUPPORTED",
                    "E.L3F.GENERATOR_CONTROL_INVALID",
                    "L.L3F.IBR_TO_GENERATOR"])
    # Codes whose reachable case lives in a dedicated testset above.
    covered = union(emitted, Set([
        "E.L3F.MULTIPLE_SOURCES", "E.L3F.SOURCE_MISSING",
        "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED", "E.L3F.KRON_REDUCTION_FAILED",
        "E.L3F.CENTER_TAP_ORIENTATION_UNSUPPORTED",
        "E.L3F.NEUTRAL_LIMIT_DISCARDED",
    ]))
    @test setdiff(expected, covered) == Set{String}()
    # Anything emitted here that the inventory does not name is a new code.
    @test setdiff(emitted,
        union(expected, lowering, delta, controls,
              Set(["E.L3F.KRON_REDUCTION_FAILED"]))) == Set{String}()
end

@testset "LinDist3Flow build refuses inapplicable input" begin
    meshed = _l3f_two_bus()
    meshed["bus"]["third"] = Dict{String,Any}("terminal_names" => ["a"])
    meshed["line"]["l2"] = Dict{String,Any}(
        "bus_from" => "load", "bus_to" => "third",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"], "linecode" => "lc")
    meshed["line"]["l3"] = Dict{String,Any}(
        "bus_from" => "third", "bus_to" => "source",
        "terminal_map_from" => ["a"], "terminal_map_to" => ["a"], "linecode" => "lc")

    err = try
        build_l3f_opf(meshed, Clarabel.Optimizer;
            options=L3FOptions(validate_nonlinear=false))
        nothing
    catch e
        e
    end
    @test err isa L3FInapplicableError
    @test !is_l3f_applicable(err.report)
    message = sprint(showerror, err)
    @test occursin("L3F-BMOPF is inapplicable", message)
    @test occursin("E.L3F.TOPOLOGY_NOT_RADIAL", message)
    # A warning alone never blocks a build.
    warned = _l3f_two_bus()
    warned["generator"] = Dict("g" => Dict{String,Any}(
        "bus" => "load", "terminal_map" => ["a"], "configuration" => "SINGLE_PHASE",
        "p_min" => [0.0], "p_max" => [1.0], "q_min" => [0.0], "q_max" => [0.0]))
    @test build_l3f_opf(warned, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:cost)) isa L3FBuild
end

@testset "LinDist3Flow delta generator channel allocation" begin
    v = 230.0
    net = Dict{String,Any}(
        "bus" => Dict("b" => Dict{String,Any}("terminal_names" => ["a", "b", "c"])),
        "voltage_source" => Dict("s" => Dict{String,Any}(
            "bus" => "b", "terminal_map" => ["a", "b", "c"], "configuration" => "WYE",
            "v_magnitude" => fill(v, 3), "v_angle" => [0.0, -2pi / 3, 2pi / 3],
            "cost" => fill(1.0, 3))),
        "generator" => Dict("g" => Dict{String,Any}(
            "bus" => "b", "terminal_map" => ["a", "b", "c"], "configuration" => "DELTA",
            "p_min" => zeros(3), "p_max" => fill(5_000.0, 3),
            "q_min" => zeros(3), "q_max" => zeros(3), "cost" => fill(0.0, 3))),
        "load" => Dict("l" => Dict{String,Any}(
            "bus" => "b", "terminal_map" => ["a", "b", "c"], "configuration" => "WYE",
            "model" => "constant_power", "p_nom" => fill(10_000.0, 3), "q_nom" => zeros(3))))
    result = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:cost),
        solver_options=_l3f_clarabel())
    @test result.solve.optimal
    # A delta generator displaces source import one-for-one in aggregate.
    @test sum(result.generators["g"]["pg"]) ≈ 15_000.0 atol=1e-2
    @test sum(result.sources["s"]["pg"]) ≈ 15_000.0 atol=1e-2

    # The channel-to-terminal map of a closed delta has rank 2, so the published
    # per-channel split is determined only up to a circulating component. Total
    # complex power, by contrast, is conserved exactly for any reference.
    D = [1.0 -1.0 0.0; 0.0 1.0 -1.0; -1.0 0.0 1.0]
    vbar = ComplexF64[v, v * cis(-2pi / 3), v * cis(2pi / 3)]
    H = connection_power_map(D, vbar)
    @test rank(H.matrix) == 2
    @test H.matrix * H.reference_winding_voltage ≈ zeros(ComplexF64, 3) atol=1e-9
    @test sum(H.reference_winding_voltage) ≈ 0 atol=1e-9
    for s in (ComplexF64[3 + 1im, -2 + 4im, 5 - 2im], ComplexF64[1, 1, 1])
        @test sum(H.matrix * s) ≈ sum(s)
    end
end
