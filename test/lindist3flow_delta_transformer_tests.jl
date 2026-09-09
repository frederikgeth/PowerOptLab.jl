using Test
using LinearAlgebra
using JuMP
using Clarabel
using PowerOptLab

# Yd / Dy transformer banks.
#
# The coil relation of an ideal three-phase delta-wye bank is
#
#     D v_delta = g v_wye,      g = sqrt(3) * v_nom_delta / v_nom_wye
#
# with `D` the delta incidence and `v_wye` measured against the (grounded,
# reduced) star point. This is BMOPFTools' executable convention: its
# `wye_delta` uses n_eff = sqrt(3)/N and its `delta_wye` uses n_eff = N*sqrt(3),
# both with N = v_nom_from/v_nom_to and both v_nom quoted phase-to-neutral. The
# two spellings collapse to the same statement.
#
# `D` has rank 2, and that is the whole story of this component. The relation
# determines the wye voltages from the delta ones but not the reverse: a delta
# winding neither imposes nor carries a zero-sequence terminal voltage. So the
# bank is exact with the delta winding upstream, and needs a gauge — an
# assumption, offered only under `unsupported=:approximate` — with the wye
# winding upstream.

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

@testset "LinDist3Flow Yd/Dy coil relation" begin
    # Both spellings state the same coil relation, so their gains agree whenever
    # the same physical bank is described.
    for (subtype, source_on_delta) in (("delta_wye", true), ("wye_delta", false))
        net = _l3f_dy_case(subtype; source_on_delta=source_on_delta)
        build = build_l3f_opf(net, Clarabel.Optimizer;
            options=L3FOptions(validate_nonlinear=false,
                               unsupported=source_on_delta ? :reject : :approximate))
        reference = build.reference
        v_delta = ComplexF64[reference.voltage[("d", t)] for t in ("a","b","c")]
        v_wye = ComplexF64[reference.voltage[("y", t)] for t in ("a","b","c")]
        g = sqrt(3.0) * _L3F_DY_VPN_HV / _L3F_DY_VPN_LV
        @test _L3F_DY_D * v_delta ≈ g .* v_wye rtol=1e-10
    end

    # The wye side of an ideal bank carries no zero-sequence voltage: D's rows
    # sum to zero, so anything it produces does too. That is the physical
    # content of an ideal delta blocking zero sequence.
    build = build_l3f_opf(_l3f_dy_case("delta_wye"), Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false))
    v_wye = ComplexF64[build.reference.voltage[("y", t)] for t in ("a","b","c")]
    @test abs(sum(v_wye)) < 1e-9 * sum(abs, v_wye)

    # The bank's map is singular, which is why orientation is part of the model
    # and not a bookkeeping detail.
    @test rank(_L3F_DY_D) == 2

    # Power is conserved exactly through a non-diagonal channel map.
    vbar = ComplexF64[_L3F_DY_VPN_HV * cis(a) for a in (0.0, -2pi/3, 2pi/3)]
    g = sqrt(3.0) * _L3F_DY_VPN_HV / _L3F_DY_VPN_LV
    H = connection_power_map(_L3F_DY_D ./ g, vbar)
    @test rank(H.matrix) == 2
    for s in (ComplexF64[3+1im, -2+4im, 5-2im], ComplexF64[1, 1, 1])
        @test sum(H.matrix * s) ≈ sum(s)
    end
end

@testset "LinDist3Flow Yd/Dy solves with the delta upstream" begin
    for subtype in ("delta_wye", "wye_delta")
        net = _l3f_dy_case(subtype; source_on_delta=true, line=true)
        report = check_l3f_applicability(net)
        @test is_l3f_applicable(report)            # exact at the default policy
        @test !_l3f_dy_has(report, "E.L3F.DELTA_ORIENTATION_UNSUPPORTED")

        result = _l3f_dy_solve(net)
        @test result.solve.optimal
        # The turns ratio lands the secondary at its nominal magnitude.
        for t in ("a","b","c")
            @test result.buses["y"][t]["vm"] ≈ _L3F_DY_VPN_LV atol=1e-6
        end
        # Lossless balance: the source carries exactly the load.
        @test sum(result.sources["v"]["pg"]) ≈ sum(_L3F_DY_P) atol=1e-3
        @test sum(result.sources["v"]["qg"]) ≈ sum(_L3F_DY_Q) atol=1e-3
        @test result.transformers["t"]["effective_ratio_from_to"] ≈
            (subtype == "delta_wye" ? _L3F_DY_VPN_HV / _L3F_DY_VPN_LV :
                                      _L3F_DY_VPN_LV / _L3F_DY_VPN_HV)

        # SI and per-unit coordinates agree through a non-diagonal map.
        si = _l3f_dy_solve(net; per_unit=false)
        @test si.solve.optimal
        for bus in keys(result.buses), terminal in keys(result.buses[bus])
            @test result.buses[bus][terminal]["w"] ≈ si.buses[bus][terminal]["w"] rtol=1e-7
        end
        @test result.transformers["t"]["p"] ≈ si.transformers["t"]["p"] rtol=1e-7
    end
end

@testset "LinDist3Flow Yd/Dy orientation gate" begin
    for subtype in ("delta_wye", "wye_delta")
        net = _l3f_dy_case(subtype; source_on_delta=false)

        # With the wye winding upstream the delta terminals are undetermined in
        # their common component, so :reject and :lower both refuse.
        for policy in (:reject, :lower)
            report = check_l3f_applicability(net; options=L3FOptions(unsupported=policy))
            @test !is_l3f_applicable(report)
            @test _l3f_dy_has(report, "E.L3F.DELTA_ORIENTATION_UNSUPPORTED")
        end
        @test_throws L3FInapplicableError build_l3f_opf(net, Clarabel.Optimizer;
            options=L3FOptions(validate_nonlinear=false))

        # :approximate accepts the zero-zero-sequence gauge, as a warning.
        report = check_l3f_applicability(net; options=L3FOptions(unsupported=:approximate))
        @test is_l3f_applicable(report)
        gauge = only(f for f in report.findings
                     if f.code == "A.L3F.DELTA_ZERO_SEQUENCE_GAUGE")
        @test gauge.severity == :warning
        @test gauge.evidence["delta_bus"] == "d"
        @test gauge.evidence["wye_bus"] == "y"

        result = _l3f_dy_solve(net; unsupported=:approximate)
        @test result.solve.optimal
        for t in ("a","b","c")
            @test result.buses["d"][t]["vm"] ≈ _L3F_DY_VPN_HV atol=1e-6
        end
        @test sum(result.sources["v"]["pg"]) ≈ sum(_L3F_DY_P) atol=1e-3

        # The gauge is exactly "no zero-sequence voltage at the delta bus".
        build = build_l3f_opf(net, Clarabel.Optimizer;
            options=L3FOptions(validate_nonlinear=false, unsupported=:approximate))
        v_delta = ComplexF64[build.reference.voltage[("d", t)] for t in ("a","b","c")]
        @test abs(sum(v_delta)) < 1e-9 * sum(abs, v_delta)
    end
end

@testset "LinDist3Flow Yd/Dy data and rating contract" begin
    # A Kron-reduced bank carries three phase conductors per side; the wye star
    # point is the eliminated conductor.
    short = _l3f_dy_case("delta_wye")
    short["transformer"]["delta_wye"]["t"]["terminal_map_from"] = ["a","b"]
    short["transformer"]["delta_wye"]["t"]["terminal_map_to"] = ["a","b"]
    @test _l3f_dy_has(check_l3f_applicability(short), "E.L3F.DEVICE_ARITY")

    invalid = _l3f_dy_case("delta_wye")
    invalid["transformer"]["delta_wye"]["t"]["v_nom_to"] = 0.0
    @test _l3f_dy_has(check_l3f_applicability(invalid),
                      "E.L3F.TRANSFORMER_RATIO_INVALID")

    # Leakage is still refused by default and lowered exactly under :lower,
    # exactly as for every other subtype.
    leaky = _l3f_dy_case("delta_wye"; extra=Dict{String,Any}(
        "r_series_from" => 0.5, "x_series_from" => 1.5))
    @test _l3f_dy_has(check_l3f_applicability(leaky),
                      "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED")
    lowered = check_l3f_applicability(leaky; options=L3FOptions(unsupported=:lower))
    @test is_l3f_applicable(lowered)
    @test _l3f_dy_has(lowered, "L.L3F.TRANSFORMER_LEAKAGE_LOWERED")
    @test _l3f_dy_solve(leaky; unsupported=:lower).solve.optimal

    # A rating binds the winding it is declared on. The two sides of a Yd bank
    # carry different terminal power, so a from-side limit must reach the
    # from-side expression — a distinction that is vacuous for a diagonal map.
    base = _l3f_dy_case("delta_wye")
    @test _l3f_dy_solve(base).solve.optimal
    # 11 kV side: ~95 kVA total, so ~33 kVA per terminal at 6351 V is ~5.2 A.
    tight_from = _l3f_dy_case("delta_wye";
        extra=Dict{String,Any}("i_max_from" => fill(1.0, 3)))
    @test !_l3f_dy_solve(tight_from).solve.optimal
    ample_from = _l3f_dy_case("delta_wye";
        extra=Dict{String,Any}("i_max_from" => fill(100.0, 3)))
    @test _l3f_dy_solve(ample_from).solve.optimal
    # 400 V side carries the same power at a much lower voltage, so the same
    # ampere figure is generous there and restrictive here.
    tight_to = _l3f_dy_case("delta_wye";
        extra=Dict{String,Any}("i_max_to" => fill(20.0, 3)))
    @test !_l3f_dy_solve(tight_to).solve.optimal
    ample_to = _l3f_dy_case("delta_wye";
        extra=Dict{String,Any}("i_max_to" => fill(1_000.0, 3)))
    @test _l3f_dy_solve(ample_to).solve.optimal
end

@testset "LinDist3Flow Yd/Dy against OpenDSS" begin
    # The oracle must implement the SAME vector group, or it is measuring the
    # convention rather than the formulation. BMOPFTools pairs delta coil k with
    # wye phase k; OpenDSS's default delta node order `d.1.2.3` pairs them
    # differently, so the delta winding is written `d.2.3.1` here. See the
    # vector-group testset below for what the default order does instead.
    function opendss_wye_voltages(vm, va, delta_nodes; line=false)
        dss(cmd) = OpenDSSDirect.dss(cmd)
        dss("Clear")
        dss("New Circuit.l3fdy bus1=dummy basekv=11 phases=3")
        dss("Edit Vsource.source enabled=no")
        for phase in 1:3
            dss("New Vsource.v$phase phases=1 bus1=d.$phase.0 " *
                "basekv=$(vm[phase] / 1000) pu=1 angle=$(rad2deg(va[phase])) " *
                "mvasc3=1e9 mvasc1=1e9")
        end
        dss("New Transformer.t phases=3 windings=2 XHL=0.000001 %loadloss=0 " *
            "%noloadloss=0 buses=[$delta_nodes y.1.2.3.0] conns=[delta wye] " *
            "kvs=[11 0.4] kvas=[500 500]")
        load_bus = "y"
        if line
            R = [0.05 0.01 0.01; 0.01 0.05 0.01; 0.01 0.01 0.05]
            X = [0.02 0.005 0.005; 0.005 0.02 0.005; 0.005 0.005 0.02]
            rows(M) = join((join(M[i, 1:i], " ") for i in 1:3), " | ")
            dss("New LineCode.lc nphases=3 units=km rmatrix=[$(rows(R))] " *
                "xmatrix=[$(rows(X))] cmatrix=[0 | 0 0 | 0 0 0]")
            dss("New Line.l1 phases=3 bus1=y.1.2.3 bus2=end.1.2.3 linecode=lc " *
                "length=1 units=km")
            load_bus = "end"
        end
        for phase in 1:3
            dss("New Load.d$phase bus1=$load_bus.$phase.0 phases=1 conn=wye model=1 " *
                "kv=$(_L3F_DY_VPN_LV / 1000) kw=$(_L3F_DY_P[phase] / 1000) " *
                "kvar=$(_L3F_DY_Q[phase] / 1000) vminpu=0")
        end
        dss("Set maxiterations=200 tolerance=1e-10")
        dss("Set controlmode=off")
        dss("Solve mode=snapshot")
        @test OpenDSSDirect.Solution.Converged()
        names = lowercase.(OpenDSSDirect.Circuit.AllNodeNames())
        volts = ComplexF64.(OpenDSSDirect.Circuit.AllBusVolts())
        Dict(n => v for (n, v) in zip(names, volts))
    end

    if isdefined(Main, :OpenDSSDirect)
        balanced_vm, balanced_va = fill(_L3F_DY_VPN_HV, 3), [0.0, -2pi/3, 2pi/3]
        # Deliberately unbalanced in both magnitude and angle: this is the case
        # that can tell two vector groups apart, and the one a balanced test
        # would silently pass under either.
        unbalanced_vm = [1.06, 0.94, 1.00] .* _L3F_DY_VPN_HV
        unbalanced_va = [0.0, -2.05, 2.2]

        for (tag, vm, va) in (("balanced", balanced_vm, balanced_va),
                              ("unbalanced", unbalanced_vm, unbalanced_va))
            net = _l3f_dy_case("delta_wye"; source_on_delta=true)
            net["voltage_source"]["v"]["v_magnitude"] = collect(vm)
            net["voltage_source"]["v"]["v_angle"] = collect(va)
            result = _l3f_dy_solve(net)
            @test result.solve.optimal
            dss_v = opendss_wye_voltages(vm, va, "d.2.3.1")
            # Magnitudes agreeing under an UNBALANCED reference is the real
            # content: it says the two decks pair the same delta coil with the
            # same wye phase. An ideal bank feeding its own bus carries no
            # linearization error, so this is an identity up to solver tolerance.
            for (phase, terminal) in enumerate(("a","b","c"))
                @test result.buses["y"][terminal]["vm"] ≈ abs(dss_v["y.$phase"]) rtol=1e-7
            end
            # What remains is a uniform rotation (180 degrees here: the two decks
            # traverse the same coil in opposite senses). Uniform is the property
            # that matters -- a per-phase difference would mean a different
            # pairing, which the magnitudes above would already have caught.
            offsets = [rad2deg(result.buses["y"][t]["reference_angle"]) -
                       rad2deg(angle(dss_v["y.$k"])) for (k, t) in enumerate(("a","b","c"))]
            @test all(o -> abs(rem(o - offsets[1], 360, RoundNearest)) < 1e-6, offsets)
        end

        # OpenDSS's DEFAULT delta node order is a different vector group, and a
        # balanced comparison cannot see it while an unbalanced one can.
        let net = _l3f_dy_case("delta_wye"; source_on_delta=true)
            net["voltage_source"]["v"]["v_magnitude"] = collect(unbalanced_vm)
            net["voltage_source"]["v"]["v_angle"] = collect(unbalanced_va)
            result = _l3f_dy_solve(net)
            matched = opendss_wye_voltages(unbalanced_vm, unbalanced_va, "d.2.3.1")
            default = opendss_wye_voltages(unbalanced_vm, unbalanced_va, "d.1.2.3")
            rel(d) = maximum(abs(result.buses["y"][t]["vm"] - abs(d["y.$k"]))
                             for (k, t) in enumerate(("a","b","c"))) / _L3F_DY_VPN_LV
            @test rel(matched) < 1e-7
            @test rel(default) > 0.05        # about 11%, and not linearization error
        end

        # With a line downstream the linearization error appears, and only then.
        net = _l3f_dy_case("delta_wye"; source_on_delta=true, line=true)
        result = _l3f_dy_solve(net)
        dss_v = opendss_wye_voltages(balanced_vm, balanced_va, "d.2.3.1"; line=true)
        errors = [abs(result.buses[bus][t]["vm"] - abs(dss_v["$bus.$k"])) / _L3F_DY_VPN_LV
                  for bus in ("y", "end") for (k, t) in enumerate(("a","b","c"))]
        @test length(errors) == 6
        @test maximum(errors[1:3]) < 1e-7      # the bank itself is exact
        @test 1e-5 < maximum(errors) < 0.004   # the line carries the omitted losses
    else
        @test_skip "Requires OpenDSSDirect"
    end
end

@testset "LinDist3Flow Yd/Dy vector group is a real choice" begin
    # BMOPFTools pairs delta coil k with wye phase k. Writing the delta winding
    # in OpenDSS's default node order pairs them differently, and the difference
    # is NOT a modelling detail that washes out:
    #
    #   balanced reference    the two groups differ by a uniform 60 degree
    #                         rotation, which this formulation cannot observe
    #   unbalanced reference  they differ by a cyclic relabelling of which delta
    #                         pair drives which wye phase, so the magnitudes
    #                         themselves move
    #
    # A balanced-only comparison therefore passes under either convention. This
    # pins the distinction so that a future "fix" to match OpenDSS's default
    # cannot quietly change the physics.
    V = _L3F_DY_VPN_HV
    D_next = _L3F_DY_D                                     # coil k spans (k, k+1)
    D_prev = [1.0 0.0 -1.0; -1.0 1.0 0.0; 0.0 -1.0 1.0]    # coil k spans (k, k-1)

    balanced = ComplexF64[V, V * cis(-2pi/3), V * cis(2pi/3)]
    a, b = D_next * balanced, D_prev * balanced
    @test abs.(a) ≈ abs.(b)                                # indistinguishable
    offsets = rad2deg.(angle.(a) .- angle.(b))
    @test all(o -> abs(rem(o - offsets[1], 360, RoundNearest)) < 1e-9, offsets)
    @test abs(rem(offsets[1], 360, RoundNearest)) ≈ 60.0 atol=1e-9

    unbalanced = ComplexF64[1.06V, 0.94V * cis(-2.05), V * cis(2.2)]
    a, b = D_next * unbalanced, D_prev * unbalanced
    @test !isapprox(abs.(a), abs.(b); rtol=1e-3)           # now distinguishable
    # The other group is the same magnitudes on different phases: coil k of one
    # is the negated coil k-1 of the other.
    @test abs.(b) ≈ circshift(abs.(a), 1)
    @test maximum(abs.(abs.(a) .- abs.(b))) / V > 0.05     # ~11% here, not noise
end

@testset "LinDist3Flow Yd/Dy reference rotation invariance" begin
    # The claim relied on above, stated directly: rotating every reference
    # phasor downstream of the bank leaves the solved squared magnitudes
    # unchanged, because the formulation only ever forms ratios within a bus.
    net = _l3f_dy_case("delta_wye"; source_on_delta=true, line=true)
    base = _l3f_dy_solve(net)
    build = build_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false))
    rotated = Dict{Tuple{String,String},ComplexF64}()
    for (key, value) in build.reference.voltage
        rotated[key] = key[1] == "d" ? value : value * cis(deg2rad(60.0))
    end
    turned = _l3f_dy_solve(net; reference_policy=:auto)
    turned = solve_l3f_opf(net, Clarabel.Optimizer; options=_l3f_dy_opts(),
        reference=rotated, solver_options=("verbose" => false,))
    @test turned.solve.optimal
    for bus in keys(base.buses), terminal in keys(base.buses[bus])
        @test base.buses[bus][terminal]["w"] ≈ turned.buses[bus][terminal]["w"] rtol=1e-9
    end
end
