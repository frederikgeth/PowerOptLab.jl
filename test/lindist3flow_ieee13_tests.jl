using Test
using LinearAlgebra
using JuMP
using Clarabel
using PowerOptLab

# Modified IEEE 13 node test feeder.
#
# Data transcribed from the IEEE PES Distribution Test Feeder Working Group
# archive (feeder13.zip, "IEEE 13 Node Test Feeder", 2004) and cross-checked
# against the maintained OpenDSS deck `IEEE13Nodeckt.dss` / `IEEELineCodes.DSS`.
# The impedance matrices below are ohms per 1000 ft; multiplying a diagonal by
# 5.28 reproduces the published ohms-per-mile values exactly (0.065625*5.28 =
# 0.3465 for configuration 601).
#
# IEEE 13 is the feeder whose regulator is a bank of three single-phase units on
# one bus pair. That topology is radial per conductor but shows three parallel
# edges on the bus graph, so it exercises the conductor-level radiality rule
# directly. It also brings coverage IEEE 37 does not: buses with only one or two
# retained phases, a single delta leg between two terminals, constant-impedance
# loads in both wye and delta connection, an in-line transformer, and fixed
# shunt capacitors.
#
# Modifications, each forced by the LinDist3Flow supported slice and mirrored in
# the OpenDSS oracle so the two decks describe identical physics:
#
#   1. The 115 kV source and substation transformer are removed; the source sits
#      at bus 650 at exactly nominal 4.16 kV (2401.78 V line-to-neutral).
#   2. The regulator bank is ideal and at fixed taps 1.0625 / 1.0500 / 1.06875,
#      with leakage and load loss removed. RegControl is not used.
#   3. XFM-1 (633-634) becomes three ideal single-phase units with leakage and
#      load loss removed. Each carries the bank nameplate rather than a third of
#      it, so splitting the bank introduces no per-phase limit the original does
#      not have; OpenDSS ignores transformer ratings in a power flow either way.
#   4. Loads 692 and 611 are constant-current in the original (`Model=5`). The
#      formulation rejects a constant-current law as non-affine in squared
#      voltage, so both become constant power, and the OpenDSS deck uses
#      `Model=1` for them to match.
#   5. Line shunt capacitance is dropped (the original deck also omits it).
#
# Everything else is the published feeder: all line lengths and configurations,
# the phasings, both capacitor banks, and every other load model and value.

const _L3F_IEEE13_Z_PER_KFT = Dict(
    "601" => ComplexF64[
        0.065625+0.192784091im  0.029545455+0.095018939im 0.029924242+0.080227273im
        0.029545455+0.095018939im 0.063920455+0.19844697im 0.02907197+0.072897727im
        0.029924242+0.080227273im 0.02907197+0.072897727im 0.064659091+0.195984848im],
    "602" => ComplexF64[
        0.142537879+0.22375im    0.029924242+0.080227273im 0.029545455+0.095018939im
        0.029924242+0.080227273im 0.14157197+0.226950758im 0.02907197+0.072897727im
        0.029545455+0.095018939im 0.02907197+0.072897727im 0.140833333+0.229393939im],
    "603" => ComplexF64[
        0.251780303+0.255132576im 0.039128788+0.086950758im
        0.039128788+0.086950758im 0.250719697+0.256988636im],
    "604" => ComplexF64[
        0.250719697+0.256988636im 0.039128788+0.086950758im
        0.039128788+0.086950758im 0.251780303+0.255132576im],
    "605" => ComplexF64[0.251742424+0.255208333im;;],
    "606" => ComplexF64[
        0.151174242+0.084526515im 0.060454545+0.006212121im 0.053958333-0.002708333im
        0.060454545+0.006212121im 0.149450758+0.076534091im 0.060454545+0.006212121im
        0.053958333-0.002708333im 0.060454545+0.006212121im 0.151174242+0.084526515im],
    "607" => ComplexF64[0.254261364+0.097045455im;;],
    # The 671-692 switch, kept as the deck's near-zero series impedance.
    "sw"  => ComplexF64[1e-4 0 0; 0 1e-4 0; 0 0 1e-4],
)

# (id, from bus, from terminals, to bus, to terminals, configuration, kft)
const _L3F_IEEE13_LINES = [
    ("650632", "RG60", ["a","b","c"], "632", ["a","b","c"], "601", 2.000),
    ("632670", "632",  ["a","b","c"], "670", ["a","b","c"], "601", 0.667),
    ("670671", "670",  ["a","b","c"], "671", ["a","b","c"], "601", 1.333),
    ("671680", "671",  ["a","b","c"], "680", ["a","b","c"], "601", 1.000),
    ("632633", "632",  ["a","b","c"], "633", ["a","b","c"], "602", 0.500),
    ("632645", "632",  ["c","b"],     "645", ["c","b"],     "603", 0.500),
    ("645646", "645",  ["c","b"],     "646", ["c","b"],     "603", 0.300),
    ("692675", "692",  ["a","b","c"], "675", ["a","b","c"], "606", 0.500),
    ("671684", "671",  ["a","c"],     "684", ["a","c"],     "604", 0.300),
    ("684611", "684",  ["c"],         "611", ["c"],         "605", 0.300),
    ("684652", "684",  ["a"],         "652", ["a"],         "607", 0.800),
    ("671692", "671",  ["a","b","c"], "692", ["a","b","c"], "sw",  1.000),
]

const _L3F_IEEE13_BUS_TERMINALS = Dict(
    "650" => ["a","b","c"], "RG60" => ["a","b","c"], "632" => ["a","b","c"],
    "633" => ["a","b","c"], "634" => ["a","b","c"], "645" => ["b","c"],
    "646" => ["b","c"], "670" => ["a","b","c"], "671" => ["a","b","c"],
    "680" => ["a","b","c"], "684" => ["a","c"], "611" => ["c"],
    "652" => ["a"], "692" => ["a","b","c"], "675" => ["a","b","c"],
)

const _L3F_IEEE13_VLN = 4160.0 / sqrt(3)          # 2401.777 V
const _L3F_IEEE13_VLN_634 = 480.0 / sqrt(3)       # 277.128 V
const _L3F_IEEE13_TAPS = [1.0625, 1.0500, 1.06875]

function _l3f_ieee13_linecode(Z)
    out = Dict{String,Any}()
    for i in axes(Z, 1), j in i:size(Z, 2)
        out["R_series_$(i)_$(j)"] = real(Z[i, j])
        out["X_series_$(i)_$(j)"] = imag(Z[i, j])
    end
    out
end

_l3f_ieee13_wye(bus, terminals, p, q; extra...) = Dict{String,Any}(
    "bus" => bus, "terminal_map" => collect(terminals), "configuration" => "WYE",
    "model" => "constant_power", "p_nom" => collect(float.(p)),
    "q_nom" => collect(float.(q)), extra...)

function _l3f_ieee13_case()
    buses = Dict{String,Any}(bus => Dict{String,Any}("terminal_names" => copy(t))
                             for (bus, t) in _L3F_IEEE13_BUS_TERMINALS)
    lines = Dict{String,Any}(id => Dict{String,Any}(
        "bus_from" => bf, "terminal_map_from" => copy(mf),
        "bus_to" => bt, "terminal_map_to" => copy(mt),
        "linecode" => code, "length" => len)
        for (id, bf, mf, bt, mt, code, len) in _L3F_IEEE13_LINES)

    # Regulator bank: three single-phase units sharing the 650 -> RG60 bus pair.
    regulators = Dict{String,Any}()
    for (k, phase) in enumerate(("a", "b", "c"))
        regulators["reg_$phase"] = Dict{String,Any}(
            "bus_from" => "650", "bus_to" => "RG60",
            "terminal_map_from" => [phase], "terminal_map_to" => [phase],
            "tap_ratio" => _L3F_IEEE13_TAPS[k], "regulator_type" => "B",
            "s_rating" => 1.666e6)
    end
    # XFM-1 idealized as three single-phase units, 4.16 kV -> 480 V.
    transformers = Dict{String,Any}()
    for phase in ("a", "b", "c")
        transformers["xfm_$phase"] = Dict{String,Any}(
            "bus_from" => "633", "bus_to" => "634",
            "terminal_map_from" => [phase], "terminal_map_to" => [phase],
            "v_nom_from" => _L3F_IEEE13_VLN, "v_nom_to" => _L3F_IEEE13_VLN_634,
            "s_rating" => 5.0e5)
    end

    loads = Dict{String,Any}(
        # Delta spot load at 671: 1155 kW / 660 kvar split across the three legs.
        "671" => Dict{String,Any}(
            "bus" => "671", "terminal_map" => ["a","b","c"], "configuration" => "DELTA",
            "model" => "constant_power",
            "p_nom" => fill(385_000.0, 3), "q_nom" => fill(220_000.0, 3)),
        "634" => _l3f_ieee13_wye("634", ["a","b","c"],
            (160_000, 120_000, 120_000), (110_000, 90_000, 90_000)),
        "645" => _l3f_ieee13_wye("645", ["b"], (170_000,), (125_000,)),
        # Single delta leg B-C, constant impedance; v_nom is the line-to-line
        # nominal because the coil spans two phases.
        "646" => Dict{String,Any}(
            "bus" => "646", "terminal_map" => ["b","c"], "configuration" => "DELTA",
            "model" => "constant_impedance", "v_nom" => 4160.0,
            "p_nom" => [230_000.0], "q_nom" => [132_000.0]),
        # Single delta leg C-A; constant current in the original, see header.
        "692" => Dict{String,Any}(
            "bus" => "692", "terminal_map" => ["c","a"], "configuration" => "DELTA",
            "model" => "constant_power",
            "p_nom" => [170_000.0], "q_nom" => [151_000.0]),
        "675" => _l3f_ieee13_wye("675", ["a","b","c"],
            (485_000, 68_000, 290_000), (190_000, 60_000, 212_000)),
        # Constant current in the original, see header.
        "611" => _l3f_ieee13_wye("611", ["c"], (170_000,), (80_000,)),
        "652" => Dict{String,Any}(
            "bus" => "652", "terminal_map" => ["a"], "configuration" => "WYE",
            "model" => "constant_impedance", "v_nom" => 2400.0,
            "p_nom" => [128_000.0], "q_nom" => [86_000.0]),
        # The 632-671 distributed load, lumped at 670 as in the OpenDSS deck.
        "670" => _l3f_ieee13_wye("670", ["a","b","c"],
            (17_000, 66_000, 117_000), (10_000, 38_000, 68_000)),
    )

    # Capacitor banks as fixed shunt susceptance, B = Q_rated / V_rated^2.
    shunts = Dict{String,Any}(
        "cap1" => Dict{String,Any}("bus" => "675", "terminal_map" => ["a","b","c"]),
        "cap2" => Dict{String,Any}("bus" => "611", "terminal_map" => ["c"],
            "B_1_1" => 100_000.0 / 2400.0^2),
    )
    for k in 1:3
        shunts["cap1"]["B_$(k)_$(k)"] = 200_000.0 / _L3F_IEEE13_VLN^2
    end

    Dict{String,Any}(
        "terminal_conventions" => Dict("phase" => ["a","b","c"], "neutral" => String[]),
        "bus" => buses,
        "linecode" => Dict{String,Any}(code => _l3f_ieee13_linecode(Z)
                                       for (code, Z) in _L3F_IEEE13_Z_PER_KFT),
        "line" => lines,
        "transformer" => Dict{String,Any}(
            "single_phase_autotransformer" => regulators,
            "single_phase" => transformers),
        "shunt" => shunts,
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => "650", "terminal_map" => ["a","b","c"], "configuration" => "WYE",
            "v_magnitude" => fill(_L3F_IEEE13_VLN, 3),
            "v_angle" => [0.0, -2pi / 3, 2pi / 3])),
        "load" => loads,
    )
end

const _L3F_IEEE13_PHASE_INDEX = Dict("a" => 1, "b" => 2, "c" => 3)

"""Build the matching OpenDSS deck and return `"bus.phase" => complex volts`."""
function _l3f_ieee13_opendss_voltages()
    dss(cmd) = OpenDSSDirect.dss(cmd)
    dss("Clear")
    dss("New Circuit.l3f13 bus1=sourcedummy basekv=4.16 phases=3")
    dss("Edit Vsource.source enabled=no")
    for phase in 1:3
        angle = rad2deg([0.0, -2pi / 3, 2pi / 3][phase])
        dss("New Vsource.v$phase phases=1 bus1=650.$phase.0 " *
            "basekv=$(_L3F_IEEE13_VLN / 1000) pu=1 angle=$angle mvasc3=1e9 mvasc1=1e9")
    end
    for (code, Z) in sort!(collect(_L3F_IEEE13_Z_PER_KFT); by=first)
        n = size(Z, 1)
        rows(M) = join((join(M[i, 1:i], " ") for i in 1:n), " | ")
        zeros_row = join((join(zeros(i), " ") for i in 1:n), " | ")
        dss("New LineCode.lc$code nphases=$n units=kft " *
            "rmatrix=[$(rows(real.(Z)))] xmatrix=[$(rows(imag.(Z)))] " *
            "cmatrix=[$zeros_row]")
    end
    # Ideal regulator bank and ideal XFM-1: negligible leakage, no loss.
    for (k, phase) in enumerate(("a", "b", "c"))
        p = _L3F_IEEE13_PHASE_INDEX[phase]
        dss("New Transformer.reg$p phases=1 windings=2 XHL=0.000001 " *
            "%loadloss=0 %noloadloss=0 " *
            "buses=[650.$p RG60.$p] conns=[wye wye] " *
            "kvs=[$(_L3F_IEEE13_VLN / 1000) $(_L3F_IEEE13_VLN / 1000)] kvas=[1666 1666]")
        dss("Transformer.reg$p.wdg=2 tap=$(_L3F_IEEE13_TAPS[k])")
        dss("New Transformer.xfm$p phases=1 windings=2 XHL=0.000001 " *
            "%loadloss=0 %noloadloss=0 " *
            "buses=[633.$p 634.$p] conns=[wye wye] " *
            "kvs=[$(_L3F_IEEE13_VLN / 1000) $(_L3F_IEEE13_VLN_634 / 1000)] kvas=[500 500]")
    end
    for (id, bf, mf, bt, mt, code, len) in _L3F_IEEE13_LINES
        nodes(m) = join(("." * string(_L3F_IEEE13_PHASE_INDEX[t]) for t in m))
        dss("New Line.l$id phases=$(length(mf)) bus1=$bf$(nodes(mf)) " *
            "bus2=$bt$(nodes(mt)) linecode=lc$code length=$len units=kft")
    end
    dss("New Load.l671 bus1=671.1.2.3 phases=3 conn=delta model=1 kv=4.16 " *
        "kw=1155 kvar=660 vminpu=0")
    for (k, phase) in enumerate((1, 2, 3))
        kw = (160.0, 120.0, 120.0)[k]; kvar = (110.0, 90.0, 90.0)[k]
        dss("New Load.l634$phase bus1=634.$phase.0 phases=1 conn=wye model=1 " *
            "kv=$(_L3F_IEEE13_VLN_634 / 1000) kw=$kw kvar=$kvar vminpu=0")
        kw675 = (485.0, 68.0, 290.0)[k]; kvar675 = (190.0, 60.0, 212.0)[k]
        dss("New Load.l675$phase bus1=675.$phase.0 phases=1 conn=wye model=1 " *
            "kv=$(_L3F_IEEE13_VLN / 1000) kw=$kw675 kvar=$kvar675 vminpu=0")
        kw670 = (17.0, 66.0, 117.0)[k]; kvar670 = (10.0, 38.0, 68.0)[k]
        dss("New Load.l670$phase bus1=670.$phase.0 phases=1 conn=wye model=1 " *
            "kv=$(_L3F_IEEE13_VLN / 1000) kw=$kw670 kvar=$kvar670 vminpu=0")
    end
    dss("New Load.l645 bus1=645.2.0 phases=1 conn=wye model=1 " *
        "kv=$(_L3F_IEEE13_VLN / 1000) kw=170 kvar=125 vminpu=0")
    dss("New Load.l646 bus1=646.2.3 phases=1 conn=delta model=2 kv=4.16 " *
        "kw=230 kvar=132 vminpu=0")
    dss("New Load.l692 bus1=692.3.1 phases=1 conn=delta model=1 kv=4.16 " *
        "kw=170 kvar=151 vminpu=0")
    dss("New Load.l611 bus1=611.3.0 phases=1 conn=wye model=1 kv=2.4 " *
        "kw=170 kvar=80 vminpu=0")
    dss("New Load.l652 bus1=652.1.0 phases=1 conn=wye model=2 kv=2.4 " *
        "kw=128 kvar=86 vminpu=0")
    dss("New Capacitor.cap1 bus1=675 phases=3 kvar=600 kv=4.16")
    dss("New Capacitor.cap2 bus1=611.3 phases=1 kvar=100 kv=2.4")
    dss("Set maxiterations=200 tolerance=1e-10")
    dss("Set controlmode=off")
    dss("Solve mode=snapshot")
    @test OpenDSSDirect.Solution.Converged()
    names = lowercase.(OpenDSSDirect.Circuit.AllNodeNames())
    values = ComplexF64.(OpenDSSDirect.Circuit.AllBusVolts())
    Dict(name => value for (name, value) in zip(names, values))
end

@testset "LinDist3Flow modified IEEE 13 feeder" begin
    net = _l3f_ieee13_case()
    report = check_l3f_applicability(net)
    @test is_l3f_applicable(report)
    @test length(net["bus"]) == 15
    @test length(net["line"]) == 12
    @test report.islands == [sort(collect(keys(_L3F_IEEE13_BUS_TERMINALS)))]
    @test report.roots == ["650"]

    # The regulator bank is three units on one bus pair, and XFM-1 is another
    # three. Both are radial per conductor and rejected by a bus-level test.
    @test length(net["transformer"]["single_phase_autotransformer"]) == 3
    @test length(net["transformer"]["single_phase"]) == 3

    result = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:feasibility),
        solver_options=("verbose" => false,))
    @test result.solve.optimal
    @test result.formulation["problem_class"] == "SOCP"

    # Fixed ideal bank: RG60 is the source magnitude times each unit's tap.
    for (k, phase) in enumerate(("a", "b", "c"))
        @test result.buses["RG60"][phase]["vm"] ≈
            _L3F_IEEE13_VLN * _L3F_IEEE13_TAPS[k] atol=1e-6
    end
    # Buses keep only their retained phases all the way through the result.
    @test sort(collect(keys(result.buses["611"]))) == ["c"]
    @test sort(collect(keys(result.buses["684"]))) == ["a", "c"]
    @test sort(collect(keys(result.buses["645"]))) == ["b", "c"]

    # SI/per-unit equivalence on a feeder with mixed phasing, delta ZP loads,
    # shunts, two transformer subtypes and a regulator bank.
    si = solve_l3f_opf(net, Clarabel.Optimizer;
        options=L3FOptions(validate_nonlinear=false, objective=:feasibility,
                           per_unit=false),
        solver_options=("verbose" => false,))
    @test si.solve.optimal
    for bus in keys(result.buses), terminal in keys(result.buses[bus])
        @test result.buses[bus][terminal]["w"] ≈ si.buses[bus][terminal]["w"] rtol=1e-7
    end

    if isdefined(Main, :OpenDSSDirect)
        dss_v = _l3f_ieee13_opendss_voltages()
        errors = Float64[]
        signed_errors = Float64[]
        per_bus = Dict{String,Float64}()
        for (bus, terminals) in result.buses
            bus == "650" && continue
            base = bus == "634" ? _L3F_IEEE13_VLN_634 : _L3F_IEEE13_VLN
            for (terminal, values) in terminals
                key = "$(lowercase(bus)).$(_L3F_IEEE13_PHASE_INDEX[terminal])"
                haskey(dss_v, key) || continue
                signed = (values["vm"] - abs(dss_v[key])) / base
                push!(signed_errors, signed)
                push!(errors, abs(signed))
                per_bus[bus] = max(get(per_bus, bus, 0.0), abs(signed))
            end
        end
        # 14 buses downstream of the source: 611 and 652 carry one phase, 645,
        # 646 and 684 carry two, the other nine carry three. 9*3 + 3*2 + 2 = 35.
        @test length(errors) == 35
        # Bracketed both ways. The upper bound is the published LinDist3Flow
        # envelope for this feeder: Table I of arXiv:2210.08550 reports 0.007 to
        # 0.01 pu maximum deviation from a Z-bus power flow on IEEE 13, and
        # Sankur et al. (arXiv:1606.04492) report under 1% on a modified IEEE 13.
        # Measured here: 0.51% maximum, 0.37% RMS. The lower bound fails loudly
        # if the comparison ever becomes vacuous rather than accurate.
        @test 0.003 < maximum(errors) < 0.008
        @test 0.002 < sqrt(sum(abs2, errors) / length(errors)) < 0.005

        # Every deviation is an overestimate: omitting series losses can only
        # make the linear model optimistic about voltage. A sign flip would mean
        # a different error source had appeared.
        @test all(>=(0.0), signed_errors)
        @test sum(signed_errors) / length(signed_errors) > 0.002

        # The loss gap itself, from the oracle: OpenDSS delivers what the feeder
        # actually costs, the lossless balance delivers only the load.
        # TotalPower is complex kW + j kvar, negative into the circuit.
        opendss_import = -real(sum(OpenDSSDirect.Circuit.TotalPower())) * 1_000
        l3f_import = sum(result.sources["source"]["pg"])
        @test l3f_import < opendss_import
        @test (opendss_import - l3f_import) / l3f_import ≈ 0.030 atol=5e-3

        # The regulated bus is matched by the ideal OpenDSS bank, which is what
        # makes this an independent check of the regulator rather than of the
        # secondary phasors being fed in.
        @test per_bus["RG60"] < 1e-4
        # The error grows with electrical distance, as the omitted losses imply.
        @test per_bus["632"] < per_bus["675"]
    else
        @test_skip "Requires OpenDSSDirect"
    end
end
