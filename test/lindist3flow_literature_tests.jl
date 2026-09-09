using Test
using LinearAlgebra
using JuMP
using Ipopt
using PowerOptLab

# Independent transcription of the IEEE 37 topology and impedance data used in
# Bazrafshan, Gatsis & Zhu (PSCC 2018), repository commit
# 79cb2c13980880ce2176238c7e85560a2ff79566.  The maintained OpenDSS IEEE 37
# deck at electricdss-tst commit 3b208397160213cae4a9e2d0a7d1aa3528ce26e1
# was used to cross-check names and lengths; no source code from either project
# is vendored here.  The paper replaces delta loads by constant-power grounded-wye
# equivalents, removes line shunts, and reports the fixed open-delta ratios
# 0.9062 (AB) and 0.9062 (CB).
const _L3F_IEEE37_LINES = [
    ("L1", "701", "702", 722, 0.96), ("L2", "702", "705", 724, 0.40),
    ("L3", "702", "713", 723, 0.36), ("L4", "702", "703", 722, 1.32),
    ("L5", "703", "727", 724, 0.24), ("L6", "703", "730", 723, 0.60),
    ("L7", "704", "714", 724, 0.08), ("L8", "704", "720", 723, 0.80),
    ("L9", "705", "742", 724, 0.32), ("L10", "705", "712", 724, 0.24),
    ("L11", "706", "725", 724, 0.28), ("L12", "707", "724", 724, 0.76),
    ("L13", "707", "722", 724, 0.12), ("L14", "708", "733", 723, 0.32),
    ("L15", "708", "732", 724, 0.32), ("L16", "709", "731", 723, 0.60),
    ("L17", "709", "708", 723, 0.32), ("L18", "710", "735", 724, 0.20),
    ("L19", "710", "736", 724, 1.28), ("L20", "711", "741", 723, 0.40),
    ("L21", "711", "740", 724, 0.20), ("L22", "713", "704", 723, 0.52),
    ("L23", "714", "718", 724, 0.52), ("L24", "720", "707", 724, 0.92),
    ("L25", "720", "706", 723, 0.60), ("L26", "727", "744", 723, 0.28),
    ("L27", "730", "709", 723, 0.20), ("L28", "733", "734", 723, 0.56),
    ("L29", "734", "737", 723, 0.64), ("L30", "734", "710", 724, 0.52),
    ("L31", "737", "738", 723, 0.40), ("L32", "738", "711", 723, 0.40),
    ("L33", "744", "728", 724, 0.20), ("L34", "744", "729", 724, 0.28),
    ("L35", "799r", "701", 721, 1.85),
]

const _L3F_IEEE37_Z_PER_KFT = Dict(
    721 => ComplexF64[
        0.2926+0.1973im 0.0673-0.0368im 0.0337-0.0417im
        0.0673-0.0368im 0.2646+0.1900im 0.0673-0.0368im
        0.0337-0.0417im 0.0673-0.0368im 0.2926+0.1973im] ./ 5.28,
    722 => ComplexF64[
        0.4751+0.2973im 0.1629-0.0326im 0.1234-0.0607im
        0.1629-0.0326im 0.4488+0.2678im 0.1629-0.0326im
        0.1234-0.0607im 0.1629-0.0326im 0.4751+0.2973im] ./ 5.28,
    723 => ComplexF64[
        1.2936+0.6713im 0.4871+0.2111im 0.4585+0.1521im
        0.4871+0.2111im 1.3022+0.6326im 0.4871+0.2111im
        0.4585+0.1521im 0.4871+0.2111im 1.2936+0.6713im] ./ 5.28,
    724 => ComplexF64[
        2.0952+0.7758im 0.5204+0.2738im 0.4926+0.2123im
        0.5204+0.2738im 2.1068+0.7398im 0.5204+0.2738im
        0.4926+0.2123im 0.5204+0.2738im 2.0952+0.7758im] ./ 5.28,
)

# (bus, delta leg, kW, kvar).  "abc" is a balanced three-phase delta total.
const _L3F_IEEE37_DELTA_LOADS = [
    ("701", "ab", 140.0, 70.0), ("701", "bc", 140.0, 70.0),
    ("701", "ca", 350.0, 175.0), ("712", "ca", 85.0, 40.0),
    ("713", "ca", 85.0, 40.0), ("714", "ab", 17.0, 8.0),
    ("714", "bc", 21.0, 10.0), ("718", "ab", 85.0, 40.0),
    ("720", "ca", 85.0, 40.0), ("722", "bc", 140.0, 70.0),
    ("722", "ca", 21.0, 10.0), ("724", "bc", 42.0, 21.0),
    ("725", "bc", 42.0, 21.0), ("727", "ca", 42.0, 21.0),
    ("728", "abc", 126.0, 63.0), ("729", "ab", 42.0, 21.0),
    ("730", "ca", 85.0, 40.0), ("731", "bc", 85.0, 40.0),
    ("732", "ca", 42.0, 21.0), ("733", "ab", 85.0, 40.0),
    ("734", "ca", 42.0, 21.0), ("735", "ca", 85.0, 40.0),
    ("736", "bc", 42.0, 21.0), ("737", "ab", 140.0, 70.0),
    ("738", "ab", 126.0, 62.0), ("740", "ca", 85.0, 40.0),
    ("741", "ca", 42.0, 21.0), ("742", "ab", 8.0, 4.0),
    ("742", "bc", 85.0, 40.0), ("744", "ab", 42.0, 21.0),
]

function _l3f_ieee37_wye_loads()
    loads = Dict{String,Vector{ComplexF64}}()
    for (bus, leg, p, q) in _L3F_IEEE37_DELTA_LOADS
        s = (p + q * im) * 1_000
        phase = get!(loads, bus, zeros(ComplexF64, 3))
        if leg == "abc"
            phase .+= s / 3
        elseif leg == "ab"
            phase[1] += s / 2; phase[2] += s / 2
        elseif leg == "bc"
            phase[2] += s / 2; phase[3] += s / 2
        else
            phase[3] += s / 2; phase[1] += s / 2
        end
    end
    loads
end

function _l3f_ieee37_linecode(code, Z)
    out = Dict{String,Any}()
    for i in axes(Z, 1), j in i:size(Z, 2)
        out["R_series_$(i)_$(j)"] = real(Z[i, j])
        out["X_series_$(i)_$(j)"] = imag(Z[i, j])
    end
    out
end

function _l3f_ieee37_case(; include_regulator=true)
    phases = ["a", "b", "c"]
    line_buses = unique(vcat([x[2] for x in _L3F_IEEE37_LINES],
                             [x[3] for x in _L3F_IEEE37_LINES]))
    source_bus = include_regulator ? "799" : "799r"
    buses = Dict{String,Any}(bus => Dict{String,Any}(
        "terminal_names" => copy(phases)) for bus in line_buses)
    buses[source_bus] = Dict{String,Any}("terminal_names" => copy(phases))
    lines = Dict{String,Any}(id => Dict{String,Any}(
        "bus_from" => from, "bus_to" => to,
        "terminal_map_from" => copy(phases), "terminal_map_to" => copy(phases),
        "linecode" => string(code), "length" => len)
        for (id, from, to, code, len) in _L3F_IEEE37_LINES)
    vbase = 4800 / sqrt(3)
    parent_v = vbase .* ComplexF64[1, cis(-2pi / 3), cis(2pi / 3)]
    effective_ratios = [0.9062, 0.9062]
    # The paper reports effective r.  BMOPF Type B stores its reciprocal tap.
    taps = inv.(effective_ratios)
    source_v = include_regulator ? parent_v :
        regulator_gain_matrix("open-delta", taps; connection="ABBC",
                              regulator_type="B") \ parent_v
    net = Dict{String,Any}(
        "terminal_conventions" => Dict("phase" => phases, "neutral" => String[]),
        "bus" => buses,
        "linecode" => Dict(string(k) => _l3f_ieee37_linecode(k, Z)
                           for (k, Z) in _L3F_IEEE37_Z_PER_KFT),
        "line" => lines,
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => source_bus, "terminal_map" => copy(phases),
            "configuration" => "WYE", "v_magnitude" => abs.(source_v),
            "v_angle" => angle.(source_v))),
        "load" => Dict("load_$bus" => Dict{String,Any}(
            "bus" => bus, "terminal_map" => copy(phases),
            "configuration" => "WYE", "model" => "constant_power",
            "p_nom" => real.(s), "q_nom" => imag.(s))
            for (bus, s) in _l3f_ieee37_wye_loads()),
    )
    if include_regulator
        net["transformer"] = Dict("open_delta_regulator" => Dict("reg" => Dict{String,Any}(
            "bus_from" => "799", "bus_to" => "799r",
            "terminal_map_from" => copy(phases), "terminal_map_to" => copy(phases),
            "connection" => "ABBC", "tap_ratio" => taps, "regulator_type" => "B")))
    end
    net
end

function _l3f_ieee37_opendss_voltages(source_v)
    OpenDSSDirect.dss("Clear")
    OpenDSSDirect.dss("New Circuit.l3f37 bus1=dummy basekv=4.8 phases=3")
    OpenDSSDirect.dss("Edit Vsource.source enabled=no")
    for (k, phase) in enumerate((1, 2, 3))
        OpenDSSDirect.dss("New Vsource.v$phase phases=1 bus1=799r.$phase.0 " *
            "basekv=$(abs(source_v[k]) / 1000) pu=1 angle=$(rad2deg(angle(source_v[k]))) " *
            "mvasc3=1e9 mvasc1=1e9")
    end
    for (code, Z) in sort!(collect(_L3F_IEEE37_Z_PER_KFT); by=first)
        rows(M) = join((join(M[i, 1:i], " ") for i in axes(M, 1)), " | ")
        OpenDSSDirect.dss("New LineCode.$code nphases=3 units=kft " *
            "rmatrix=[$(rows(real.(Z)))] xmatrix=[$(rows(imag.(Z)))] " *
            "cmatrix=[0 | 0 0 | 0 0 0]")
    end
    for (id, from, to, code, len) in _L3F_IEEE37_LINES
        OpenDSSDirect.dss("New Line.$id phases=3 bus1=$from.1.2.3 bus2=$to.1.2.3 " *
            "linecode=$code length=$len units=kft")
    end
    for (bus, s) in _l3f_ieee37_wye_loads(), phase in 1:3
        abs(s[phase]) == 0 && continue
        OpenDSSDirect.dss("New Load.load_$(bus)_$phase phases=1 bus1=$bus.$phase.0 " *
            "conn=wye model=1 kv=$(4800 / sqrt(3) / 1000) " *
            "kw=$(real(s[phase]) / 1000) kvar=$(imag(s[phase]) / 1000) vminpu=0")
    end
    OpenDSSDirect.dss("Set maxiterations=100 tolerance=1e-10")
    OpenDSSDirect.dss("Solve mode=snapshot controlmode=off")
    @test OpenDSSDirect.Solution.Converged()
    names = lowercase.(OpenDSSDirect.Circuit.AllNodeNames())
    values = ComplexF64.(OpenDSSDirect.Circuit.AllBusVolts())
    Dict(name => value for (name, value) in zip(names, values))
end

@testset "LinDist3Flow Bazrafshan IEEE 37 literature regression" begin
    loads = _l3f_ieee37_wye_loads()
    phase_load = reduce(+, values(loads))
    # The repository workbook plus its published half-adjacent conversion code
    # yields this split.  It does not reproduce the per-phase split printed in
    # the paper, although both agree on the feeder total below.
    @test real.(phase_load) ./ 2.5e6 ≈ [0.3636, 0.2732, 0.3460] atol=5e-5
    @test imag.(phase_load) ./ 2.5e6 ≈ [0.1774, 0.1342, 0.1688] atol=5e-5
    @test sum(phase_load) / 2.5e6 ≈ 0.9828 + 0.4804im atol=5e-5

    net = _l3f_ieee37_case()
    report = check_l3f_applicability(net; options=L3FOptions(kron_reduce=false))
    @test is_l3f_applicable(report)
    @test length(net["line"]) == 35
    @test length(net["bus"]) == 37
    result = solve_l3f_opf(net, Ipopt.Optimizer;
        options=L3FOptions(kron_reduce=false, validate_nonlinear=false,
                           objective=:feasibility),
        solver_options=("print_level" => 0,))
    @test result.solve.optimal
    @test sum(result.sources["source"]["pg"]) ≈ real(sum(phase_load)) atol=1e-3
    @test sum(result.sources["source"]["qg"]) ≈ imag(sum(phase_load)) atol=1e-3

    # OpenDSS is used as a nonlinear feeder oracle downstream of the regulator.
    # Supplying the exact ideal-regulator secondary phasors avoids introducing
    # OpenDSS's discrete RegControl and finite transformer impedance into this
    # fixed, continuous component contract.
    source = net["voltage_source"]["source"]
    parent_v = source["v_magnitude"] .* cis.(source["v_angle"])
    secondary_v = regulator_gain_matrix("open-delta", inv.([0.9062, 0.9062]);
        connection="ABBC", regulator_type="B") \ parent_v
    if isdefined(Main, :OpenDSSDirect)
        dss_v = _l3f_ieee37_opendss_voltages(secondary_v)
        errors = Float64[]
        for (bus, terminal_results) in result.buses,
            (phase, terminal) in enumerate(("a", "b", "c"))
            bus == "799" && continue
            key = "$(lowercase(bus)).$phase"
            haskey(dss_v, key) || continue
            push!(errors, abs(sqrt(terminal_results[terminal]["w"]) - abs(dss_v[key])) /
                          (4800 / sqrt(3)))
        end
        @test length(errors) == 108
        @test maximum(errors) < 0.012
        @test sqrt(sum(abs2, errors) / length(errors)) < 0.006
    else
        @test_skip "Requires OpenDSSDirect"
    end
end
