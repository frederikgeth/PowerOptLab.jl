using Test
using LinearAlgebra
using JuMP
using PowerOptLab

# Component builders below derive from BMOPFTools powerflow_comparison_tests.jl
# at revision 72e6cec22a66cf376c37ec4d64aef350b9f1100d; see data/state_estimation/license.md.
# Modified IEEE 13 node test feeder.
#
# Data transcribed from the IEEE PES Distribution Test Feeder Working Group
# archive (feeder13.zip, "IEEE 13 Node Test Feeder", 2004) and cross-checked
# against the maintained OpenDSS deck `IEEE13Nodeckt.dss` / `IEEELineCodes.DSS`.
# The impedance matrices below are ohms per 1000 ft; multiplying a diagonal by
# 5.28 reproduces the published ohms-per-mile values exactly (0.065625*5.28 =
# 0.3465 for configuration 601).
#
# SE-IEEE13-finite: derived from this repository's modified IEEE 13 fixture.
# The transformer/regulator losses below are deliberately specified finite
# parameters, not an approximation to the former ideal units. Both models use
# the same 0.5% winding resistances; regulator XHL=1%, transformer XHL=2%.
# Fixed taps, constant-power loads, no line shunts, source at 650. This is a
# named modified benchmark, not the original IEEE operating case.
const _SE_IEEE13_Z_PER_KFT = Dict(
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
const _SE_IEEE13_LINES = [
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

const _SE_IEEE13_BUS_TERMINALS = Dict(
    "650" => ["a","b","c"], "RG60" => ["a","b","c"], "632" => ["a","b","c"],
    "633" => ["a","b","c"], "634" => ["a","b","c"], "645" => ["b","c"],
    "646" => ["b","c"], "670" => ["a","b","c"], "671" => ["a","b","c"],
    "680" => ["a","b","c"], "684" => ["a","c"], "611" => ["c"],
    "652" => ["a"], "692" => ["a","b","c"], "675" => ["a","b","c"],
)

const _SE_IEEE13_VLN = 4160.0 / sqrt(3)          # 2401.777 V
const _SE_IEEE13_VLN_634 = 480.0 / sqrt(3)       # 277.128 V
const _SE_IEEE13_TAPS = [1.0625, 1.0500, 1.06875]

function _se_ieee13_linecode(Z)
    out = Dict{String,Any}()
    for i in axes(Z, 1), j in i:size(Z, 2)
        out["R_series_$(i)_$(j)"] = real(Z[i, j])
        out["X_series_$(i)_$(j)"] = imag(Z[i, j])
    end
    out
end

_se_ieee13_wye(bus, terminals, p, q; extra...) = Dict{String,Any}(
    "bus" => bus, "terminal_map" => collect(terminals), "configuration" => "WYE",
    "model" => "constant_power", "p_nom" => collect(float.(p)),
    "q_nom" => collect(float.(q)), extra...)

function _se_ieee13_case()
    buses = Dict{String,Any}(bus => Dict{String,Any}("terminal_names" => copy(t))
                             for (bus, t) in _SE_IEEE13_BUS_TERMINALS)
    lines = Dict{String,Any}(id => Dict{String,Any}(
        "bus_from" => bf, "terminal_map_from" => copy(mf),
        "bus_to" => bt, "terminal_map_to" => copy(mt),
        "linecode" => code, "length" => len)
        for (id, bf, mf, bt, mt, code, len) in _SE_IEEE13_LINES)

    # Regulator bank: three single-phase units sharing the 650 -> RG60 bus pair.
    regulators = Dict{String,Any}()
    for (k, phase) in enumerate(("a", "b", "c"))
        regulators["reg_$phase"] = Dict{String,Any}(
            "bus_from" => "650", "bus_to" => "RG60",
            "terminal_map_from" => [phase], "terminal_map_to" => [phase],
            "tap_ratio" => _SE_IEEE13_TAPS[k], "regulator_type" => "B",
            "s_rating" => 1.666e6,
            "r_series_from" => .005 * _SE_IEEE13_VLN^2 / 1.666e6,
            "r_series_to" => .005 * _SE_IEEE13_VLN^2 / 1.666e6,
            "x_series_from" => .01 * _SE_IEEE13_VLN^2 / 1.666e6)
    end
    # XFM-1 represented as three finite single-phase units, 4.16 kV -> 480 V.
    transformers = Dict{String,Any}()
    for phase in ("a", "b", "c")
        transformers["xfm_$phase"] = Dict{String,Any}(
            "bus_from" => "633", "bus_to" => "634",
            "terminal_map_from" => [phase], "terminal_map_to" => [phase],
            "v_nom_from" => _SE_IEEE13_VLN, "v_nom_to" => _SE_IEEE13_VLN_634,
            "s_rating" => 5.0e5,
            "r_series_from" => .005 * _SE_IEEE13_VLN^2 / 5e5,
            "r_series_to" => .005 * _SE_IEEE13_VLN_634^2 / 5e5,
            "x_series_from" => .02 * _SE_IEEE13_VLN^2 / 5e5)
    end

    loads = Dict{String,Any}(
        # Delta spot load at 671: 1155 kW / 660 kvar split across the three legs.
        "671" => Dict{String,Any}(
            "bus" => "671", "terminal_map" => ["a","b","c"], "configuration" => "DELTA",
            "model" => "constant_power",
            "p_nom" => fill(385_000.0, 3), "q_nom" => fill(220_000.0, 3)),
        "634" => _se_ieee13_wye("634", ["a","b","c"],
            (160_000, 120_000, 120_000), (110_000, 90_000, 90_000)),
        "645" => _se_ieee13_wye("645", ["b"], (170_000,), (125_000,)),
        # Single delta leg B-C, constant power; v_nom is the line-to-line
        # nominal because the coil spans two phases.
        "646" => Dict{String,Any}(
            "bus" => "646", "terminal_map" => ["b","c"], "configuration" => "DELTA",
            "model" => "constant_power", "v_nom" => 4160.0,
            "p_nom" => [230_000.0], "q_nom" => [132_000.0]),
        # Single delta leg C-A; constant current in the original, see header.
        "692" => Dict{String,Any}(
            "bus" => "692", "terminal_map" => ["c","a"], "configuration" => "DELTA",
            "model" => "constant_power",
            "p_nom" => [170_000.0], "q_nom" => [151_000.0]),
        "675" => _se_ieee13_wye("675", ["a","b","c"],
            (485_000, 68_000, 290_000), (190_000, 60_000, 212_000)),
        # Constant current in the original, see header.
        "611" => _se_ieee13_wye("611", ["c"], (170_000,), (80_000,)),
        "652" => Dict{String,Any}(
            "bus" => "652", "terminal_map" => ["a"], "configuration" => "WYE",
            "model" => "constant_power", "v_nom" => 2400.0,
            "p_nom" => [128_000.0], "q_nom" => [86_000.0]),
        # The 632-671 distributed load, lumped at 670 as in the OpenDSS deck.
        "670" => _se_ieee13_wye("670", ["a","b","c"],
            (17_000, 66_000, 117_000), (10_000, 38_000, 68_000)),
    )

    # Capacitor banks as fixed shunt susceptance, B = Q_rated / V_rated^2.
    shunts = Dict{String,Any}(
        "cap1" => Dict{String,Any}("bus" => "675", "terminal_map" => ["a","b","c"]),
        "cap2" => Dict{String,Any}("bus" => "611", "terminal_map" => ["c"],
            "B_1_1" => 100_000.0 / 2400.0^2),
    )
    for k in 1:3
        shunts["cap1"]["B_$(k)_$(k)"] = 200_000.0 / _SE_IEEE13_VLN^2
    end

    Dict{String,Any}(
        "terminal_conventions" => Dict("phase" => ["a","b","c"], "neutral" => String[]),
        "bus" => buses,
        "linecode" => Dict{String,Any}(code => _se_ieee13_linecode(Z)
                                       for (code, Z) in _SE_IEEE13_Z_PER_KFT),
        "line" => lines,
        "transformer" => Dict{String,Any}(
            "single_phase_autotransformer" => regulators,
            "single_phase" => transformers),
        "shunt" => shunts,
        "voltage_source" => Dict("source" => Dict{String,Any}(
            "bus" => "650", "terminal_map" => ["a","b","c"], "configuration" => "WYE",
            "v_magnitude" => fill(_SE_IEEE13_VLN, 3),
            "v_angle" => [0.0, -2pi / 3, 2pi / 3])),
        "load" => loads,
    )
end

const _SE_IEEE13_PHASE_INDEX = Dict("a" => 1, "b" => 2, "c" => 3)

"""Build the matching OpenDSS deck and return `"bus.phase" => complex volts`."""
function _se_ieee13_opendss_voltages()
    dss(cmd) = OpenDSSDirect.dss(cmd)
    dss("Clear")
    dss("New Circuit.se13finite bus1=sourcedummy basekv=4.16 phases=3")
    dss("Edit Vsource.source enabled=no")
    for phase in 1:3
        angle = rad2deg([0.0, -2pi / 3, 2pi / 3][phase])
        dss("New Vsource.v$phase phases=1 bus1=650.$phase.0 " *
            "basekv=$(_SE_IEEE13_VLN / 1000) pu=1 angle=$angle mvasc3=1e9 mvasc1=1e9")
    end
    for (code, Z) in sort!(collect(_SE_IEEE13_Z_PER_KFT); by=first)
        n = size(Z, 1)
        rows(M) = join((join(M[i, 1:i], " ") for i in 1:n), " | ")
        zeros_row = join((join(zeros(i), " ") for i in 1:n), " | ")
        dss("New LineCode.lc$code nphases=$n units=kft " *
            "rmatrix=[$(rows(real.(Z)))] xmatrix=[$(rows(imag.(Z)))] " *
            "cmatrix=[$zeros_row]")
    end
    # Finite regulator bank and XFM-1 with the specified winding losses.
    for (k, phase) in enumerate(("a", "b", "c"))
        p = _SE_IEEE13_PHASE_INDEX[phase]
        dss("New Transformer.reg$p phases=1 windings=2 XHL=1 " *
            "%loadloss=1 %noloadloss=0 " *
            "buses=[650.$p RG60.$p] conns=[wye wye] " *
            "kvs=[$(_SE_IEEE13_VLN / 1000) $(_SE_IEEE13_VLN / 1000)] kvas=[1666 1666]")
        dss("Transformer.reg$p.wdg=2 tap=$(_SE_IEEE13_TAPS[k])")
        dss("New Transformer.xfm$p phases=1 windings=2 XHL=2 " *
            "%loadloss=1 %noloadloss=0 " *
            "buses=[633.$p 634.$p] conns=[wye wye] " *
            "kvs=[$(_SE_IEEE13_VLN / 1000) $(_SE_IEEE13_VLN_634 / 1000)] kvas=[500 500]")
    end
    for (id, bf, mf, bt, mt, code, len) in _SE_IEEE13_LINES
        nodes(m) = join(("." * string(_SE_IEEE13_PHASE_INDEX[t]) for t in m))
        dss("New Line.l$id phases=$(length(mf)) bus1=$bf$(nodes(mf)) " *
            "bus2=$bt$(nodes(mt)) linecode=lc$code length=$len units=kft")
    end
    dss("New Load.l671 bus1=671.1.2.3 phases=3 conn=delta model=1 kv=4.16 " *
        "kw=1155 kvar=660 vminpu=0")
    for (k, phase) in enumerate((1, 2, 3))
        kw = (160.0, 120.0, 120.0)[k]; kvar = (110.0, 90.0, 90.0)[k]
        dss("New Load.l634$phase bus1=634.$phase.0 phases=1 conn=wye model=1 " *
            "kv=$(_SE_IEEE13_VLN_634 / 1000) kw=$kw kvar=$kvar vminpu=0")
        kw675 = (485.0, 68.0, 290.0)[k]; kvar675 = (190.0, 60.0, 212.0)[k]
        dss("New Load.l675$phase bus1=675.$phase.0 phases=1 conn=wye model=1 " *
            "kv=$(_SE_IEEE13_VLN / 1000) kw=$kw675 kvar=$kvar675 vminpu=0")
        kw670 = (17.0, 66.0, 117.0)[k]; kvar670 = (10.0, 38.0, 68.0)[k]
        dss("New Load.l670$phase bus1=670.$phase.0 phases=1 conn=wye model=1 " *
            "kv=$(_SE_IEEE13_VLN / 1000) kw=$kw670 kvar=$kvar670 vminpu=0")
    end
    dss("New Load.l645 bus1=645.2.0 phases=1 conn=wye model=1 " *
        "kv=$(_SE_IEEE13_VLN / 1000) kw=170 kvar=125 vminpu=0")
    dss("New Load.l646 bus1=646.2.3 phases=1 conn=delta model=1 kv=4.16 " *
        "kw=230 kvar=132 vminpu=0")
    dss("New Load.l692 bus1=692.3.1 phases=1 conn=delta model=1 kv=4.16 " *
        "kw=170 kvar=151 vminpu=0")
    dss("New Load.l611 bus1=611.3.0 phases=1 conn=wye model=1 kv=2.4 " *
        "kw=170 kvar=80 vminpu=0")
    dss("New Load.l652 bus1=652.1.0 phases=1 conn=wye model=1 kv=2.4 " *
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


function _se_net_1ph_xfmr()
    # pf_1ph_xfmr.dss: hv ──[1-ph YY xfmr, 11 kV / 0.24 kV, 50 kVA]── lv
    # %r=1.0 per winding, xhl=4.0%, %noloadloss=0.3, %imag=1.5
    # load: 30 kW + 10 kVAr on lv
    #
    # Impedance conversion (single_phase subtype, WYE-WYE):
    #   v_nom_from = 11 000 V, v_nom_to = 240 V, s_rating = 50 000 VA
    #   z_base_from = 11000² / 50000 = 2420 Ω
    #   z_base_to   =   240² / 50000 = 1.152 Ω
    #   r_series_from = 0.01 × 2420 = 24.2 Ω
    #   r_series_to   = 0.01 × 1.152 = 0.01152 Ω
    #   x_series_from = 0.04 × 2420 = 96.8 Ω  (all xsc on winding 1)
    #   x_series_to   = 0.0
    #
    # No-load branch. This hand-built fixture keeps the no-load admittance tiny
    # (its purpose is to validate the leakage/voltage behaviour; the no-load
    # branch on the correct winding-2 base is validated end-to-end by the
    # `from_dss`-based fidelity test, whose fixture carries the real core loss).
    s   = 50_000.0
    vf  = 11_000.0
    vt  =    240.0
    zbf = vf^2 / s
    zbt = vt^2 / s
    yb  = s / vf^2
    nl  = 0.003
    cm  = 0.015
    G_nl = nl * yb
    Y_nl = cm * yb
    B_nl = sqrt(max(Y_nl^2 - G_nl^2, 0.0))

    Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "hv" => Dict{String,Any}(
                "terminal_names"  => ["1", "n"],
                "neutral_terminal"=> "n"),
            "lv" => Dict{String,Any}(
                "terminal_names"  => ["1", "n"],
                "neutral_terminal"=> "n")),
        "voltage_source" => Dict{String,Any}(
            "source" => Dict{String,Any}(
                "bus"          => "hv",
                "terminal_map" => ["1"],
                "v_magnitude"  => [11_000.0],
                "v_angle"      => [0.0])),
        "shunt" => Dict{String,Any}(
            "grnd_hv" => Dict{String,Any}(
                "bus"          => "hv",
                "terminal_map" => ["n"],
                "G_1_1"        => 1000.0,
                "B_1_1"        => 0.0),
            "grnd_lv" => Dict{String,Any}(
                "bus"          => "lv",
                "terminal_map" => ["n"],
                "G_1_1"        => 1000.0,
                "B_1_1"        => 0.0)),
        "transformer" => Dict{String,Any}(
            "single_phase" => Dict{String,Any}(
                "t1" => Dict{String,Any}(
                    "bus_from"         => "hv",
                    "bus_to"           => "lv",
                    "terminal_map_from"=> ["1", "n"],
                    "terminal_map_to"  => ["1", "n"],
                    "v_nom_from"       => vf,
                    "v_nom_to"         => vt,
                    "s_rating"         => s,
                    "r_series_from"    => 0.01 * zbf,
                    "r_series_to"      => 0.01 * zbt,
                    "x_series_from"    => 0.04 * zbf,
                    "x_series_to"      => 0.0,
                    "g_no_load"        => G_nl,
                    "b_no_load"        => B_nl))),
        "load" => Dict{String,Any}(
            "ld1" => Dict{String,Any}(
                "bus"           => "lv",
                "terminal_map"  => ["1", "n"],
                "configuration" => "WYE",
                "p_nom"         => [30_000.0],
                "q_nom"         => [10_000.0])))
end


function _se_net_yd_xfmr()
    # pf_yd_xfmr.dss: hv ──[Yd xfmr, 11 kV wye / 0.415 kV delta, 500 kVA]── lv
    # %r=1.0 per winding, xhl=4.0%, %noloadloss=0.3, %imag=1.5
    # HV neutral earthed via R=0.001 Ω reactor → shunt G=1000 S at hv.n
    # LV delta: near-solid ground (R=0.001 Ω, G≈1000 S) at lv.1 to anchor
    #   common mode deterministically — mirrors ODS Reactor.grnd_lv.
    # loads: unbalanced (3:1:2) delta loads to expose per-phase errors
    #
    # Impedance conversion (wye_delta subtype):
    #   v_nom_from = 11 000 V (wye, line voltage)
    #   v_nom_to   =    415 V (delta, line voltage)
    #   s_rating   = 500 000 VA
    #   zbf = 11000² / 500000 = 242 Ω
    #   zbt =   415² / 500000 = 0.34445 Ω
    #
    #   _add_yd_transformer! variables: Iw = wye line current, Id = delta arm current.
    #   Current constraint: n_eff·Id = Iw_k - Iw_{k-1}, so at rated load |Id|=|Iw|/√3.
    #   Power loss matching to OpenDSS (%r per winding on winding kVA base):
    #     Rw × |Iw|² = %r × (S/3)  →  Rw = %r × zbf       (no √3)
    #     Rd × |Id|² = %r × (S/3)  →  Rd = %r × zbt × 3   (|Id|=rated_LV_arm_current)
    #   Wait — for the delta winding (to-side here), |Id| at rated = (S/3)/vt (arm current).
    #   So Rd × ((S/3)/vt)² = %r × S/3  →  Rd = %r × vt²/S = %r × zbt.
    #   Result: r_series = %r × z_base, no √3 on either side.
    #   xhl split 50/50 between windings (a valid T-model choice summing to xhl):
    #   r_series_from = 0.01 × 242 = 2.42 Ω    (wye winding)
    #   r_series_to   = 0.01 × 0.34445 = 3.4445e-3 Ω  (delta winding)
    #   x_series_from = 0.04/2 × 242 = 4.84 Ω  (half xhl on wye side)
    #   x_series_to   = 0.04/2 × 0.34445 = 6.889e-3 Ω  (half xhl on delta side)
    s   = 500_000.0
    vf  = 11_000.0
    vt  =    415.0
    zbf = vf^2 / s
    zbt = vt^2 / s
    nl  = 0.003
    cm  = 0.015
    # No-load branch kept tiny here (this hand-built fixture validates the
    # leakage/voltage behaviour and the near-solid delta-common-mode anchor; the
    # no-load branch on the correct winding-2 base is validated end-to-end by the
    # `from_dss` fidelity test). Referred to V_LN = vf/√3 as a small placeholder.
    yb_nl = s / (vf / sqrt(3))^2   # = 3·s/vf²
    G_nl = nl * yb_nl
    Y_nl = cm * yb_nl
    B_nl = sqrt(max(Y_nl^2 - G_nl^2, 0.0))

    Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "hv" => Dict{String,Any}(
                "terminal_names"  => ["1", "2", "3", "n"],
                "neutral_terminal"=> "n"),
            "lv" => Dict{String,Any}(
                "terminal_names"  => ["1", "2", "3"])),
        "voltage_source" => Dict{String,Any}(
            "source" => Dict{String,Any}(
                "bus"          => "hv",
                "terminal_map" => ["1", "2", "3"],
                "v_magnitude"  => [vf/sqrt(3), vf/sqrt(3), vf/sqrt(3)],
                "v_angle"      => [0.0, -2π/3, 2π/3])),
        "shunt" => Dict{String,Any}(
            "grnd_hv" => Dict{String,Any}(
                "bus"          => "hv",
                "terminal_map" => ["n"],
                "G_1_1"        => 1000.0,
                "B_1_1"        => 0.0),
            "grnd_lv" => Dict{String,Any}(
                "bus"          => "lv",
                "terminal_map" => ["1"],
                "G_1_1"        => 1000.0,
                "B_1_1"        => 0.0)),
        "transformer" => Dict{String,Any}(
            "wye_delta" => Dict{String,Any}(
                "t1" => Dict{String,Any}(
                    "bus_from"         => "hv",
                    "bus_to"           => "lv",
                    "terminal_map_from"=> ["1", "2", "3", "n"],
                    "terminal_map_to"  => ["1", "2", "3"],
                    "v_nom_from"       => vf,
                    "v_nom_to"         => vt,
                    "s_rating"         => s,
                    "r_series_from"    => 0.01 * zbf,
                    "r_series_to"      => 0.01 * zbt,
                    "x_series_from"    => 0.04 / 2 * zbf,
                    "x_series_to"      => 0.04 / 2 * zbt,
                    "g_no_load"        => G_nl,
                    "b_no_load"        => B_nl))),
        "load" => Dict{String,Any}(
            "ld1" => Dict{String,Any}(
                "bus"           => "lv",
                "terminal_map"  => ["1", "2"],
                "configuration" => "DELTA",
                "p_nom"         => [180_000.0],
                "q_nom"         => [60_000.0]),
            "ld2" => Dict{String,Any}(
                "bus"           => "lv",
                "terminal_map"  => ["2", "3"],
                "configuration" => "DELTA",
                "p_nom"         => [60_000.0],
                "q_nom"         => [20_000.0]),
            "ld3" => Dict{String,Any}(
                "bus"           => "lv",
                "terminal_map"  => ["3", "1"],
                "configuration" => "DELTA",
                "p_nom"         => [120_000.0],
                "q_nom"         => [40_000.0])))
end


function _se_net_dy_xfmr()
    # pf_dy_xfmr.dss: hv ──[Dy xfmr, 11 kV delta / 0.415 kV wye, 500 kVA]── lv
    # %r=1.0 per winding, xhl=4.0%, %noloadloss=0.3, %imag=1.5
    # LV neutral earthed via R=0.3 Ω reactor → shunt G=1/0.3 S
    # loads: unbalanced wye-to-neutral 120 kW + 40 kVAr each phase (3:3:3 here, balanced)
    #
    # Impedance conversion (delta_wye subtype):
    #   v_nom_from = 11 000 V (delta, line voltage)
    #   v_nom_to   =    415 V (wye, line voltage)
    #   s_rating   = 500 000 VA
    #   zbf = 11000² / 500000 = 242 Ω
    #   zbt =   415² / 500000 = 0.34445 Ω
    #
    #   _add_yd_transformer! variables: Id = delta arm current (from-side), Iw = wye line current.
    #   Current constraint: n_eff·Id = Iw_k - Iw_{k+1}, so |Id| = |Iw|/√3 balanced.
    #   Power loss matching to OpenDSS (%r per winding on winding kVA base):
    #     Rd × |Id|² = %r × (S/3)  →  Rd = %r × zbf       (no √3: |Id|_rated = (S/3)/(vf/√3·√3))
    #     Rw × |Iw|² = %r × (S/3)  →  Rw = %r × zbt       (no √3)
    #   Result: r_series = %r × z_base, no √3 on either side.
    #   xhl split 50/50 between windings:
    #   r_series_from = 0.01 × 242 = 2.42 Ω         (delta winding, from-side)
    #   r_series_to   = 0.01 × 0.34445 = 3.4445e-3 Ω (wye winding, to-side)
    #   x_series_from = 0.04/2 × 242 = 4.84 Ω       (half xhl on delta side)
    #   x_series_to   = 0.04/2 × 0.34445 = 6.889e-3 Ω (half xhl on wye side)
    s   = 500_000.0
    vf  = 11_000.0
    vt  =    415.0
    zbf = vf^2 / s
    zbt = vt^2 / s
    nl  = 0.003
    cm  = 0.015
    # No-load branch kept tiny here (this hand-built fixture validates the
    # leakage/voltage behaviour; the no-load branch on the correct winding-2 base
    # is validated end-to-end by the `from_dss` fidelity test). Small placeholder.
    yb_nl = s / (vf / sqrt(3))^2   # = 3·s/vf²
    G_nl = nl * yb_nl
    Y_nl = cm * yb_nl
    B_nl = sqrt(max(Y_nl^2 - G_nl^2, 0.0))

    Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "hv" => Dict{String,Any}(
                "terminal_names" => ["1", "2", "3"]),
            "lv" => Dict{String,Any}(
                "terminal_names"  => ["1", "2", "3", "n"],
                "neutral_terminal"=> "n")),
        "voltage_source" => Dict{String,Any}(
            "source" => Dict{String,Any}(
                "bus"          => "hv",
                "terminal_map" => ["1", "2", "3"],
                "v_magnitude"  => [vf/sqrt(3), vf/sqrt(3), vf/sqrt(3)],
                "v_angle"      => [0.0, -2π/3, 2π/3])),
        "shunt" => Dict{String,Any}(
            "grnd_lv" => Dict{String,Any}(
                "bus"          => "lv",
                "terminal_map" => ["n"],
                "G_1_1"        => 1.0 / 0.3,
                "B_1_1"        => 0.0)),
        "transformer" => Dict{String,Any}(
            "delta_wye" => Dict{String,Any}(
                "t1" => Dict{String,Any}(
                    "bus_from"         => "hv",
                    "bus_to"           => "lv",
                    "terminal_map_from"=> ["1", "2", "3"],
                    "terminal_map_to"  => ["1", "2", "3", "n"],
                    "v_nom_from"       => vf,
                    "v_nom_to"         => vt,
                    "s_rating"         => s,
                    "r_series_from"    => 0.01 * zbf,
                    "r_series_to"      => 0.01 * zbt,
                    "x_series_from"    => 0.04 / 2 * zbf,
                    "x_series_to"      => 0.04 / 2 * zbt,
                    "g_no_load"        => G_nl,
                    "b_no_load"        => B_nl))),
        "load" => Dict{String,Any}(
            "ld1" => Dict{String,Any}(
                "bus"           => "lv",
                "terminal_map"  => ["1", "n"],
                "configuration" => "WYE",
                "p_nom"         => [120_000.0],
                "q_nom"         => [40_000.0]),
            "ld2" => Dict{String,Any}(
                "bus"           => "lv",
                "terminal_map"  => ["2", "n"],
                "configuration" => "WYE",
                "p_nom"         => [120_000.0],
                "q_nom"         => [40_000.0]),
            "ld3" => Dict{String,Any}(
                "bus"           => "lv",
                "terminal_map"  => ["3", "n"],
                "configuration" => "WYE",
                "p_nom"         => [120_000.0],
                "q_nom"         => [40_000.0])))
end


function _se_net_autotransformer()
    # pf_autotransformer.dss: src ──[1-ph regulator, 2.4 kV, 500 kVA, tap 1.05]── reg
    # ANSI Type B at tap a=1.05 → n_eff = 1/a; lossless boost V_reg = a·V_src.
    # %r=0.5 per winding, xhl=1.0 (both kv equal → same z_base both sides).
    #   z_base = 2400² / 500000 = 11.52 Ω
    #   r_series_from = r_series_to = 0.005 × 11.52   (winding %r on own base)
    #   x_series_from = 0.01 × 11.52  (all xhl on winding 1, Γ-model)
    s   = 500_000.0
    v   = 2400.0
    zb  = v^2 / s
    Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "src" => Dict{String,Any}(
                "terminal_names"  => ["1", "n"],
                "neutral_terminal"=> "n"),
            "reg" => Dict{String,Any}(
                "terminal_names"  => ["1", "n"],
                "neutral_terminal"=> "n")),
        "voltage_source" => Dict{String,Any}(
            "vs" => Dict{String,Any}(
                "bus"          => "src",
                "terminal_map" => ["1"],
                "v_magnitude"  => [2400.0],
                "v_angle"      => [0.0])),
        "shunt" => Dict{String,Any}(
            "g_src" => Dict{String,Any}(
                "bus" => "src", "terminal_map" => ["n"],
                "G_1_1" => 1000.0, "B_1_1" => 0.0),
            "g_reg" => Dict{String,Any}(
                "bus" => "reg", "terminal_map" => ["n"],
                "G_1_1" => 1000.0, "B_1_1" => 0.0)),
        "transformer" => Dict{String,Any}(
            "single_phase_autotransformer" => Dict{String,Any}(
                "reg" => Dict{String,Any}(
                    "bus_from"         => "src",
                    "bus_to"           => "reg",
                    "terminal_map_from"=> ["1", "n"],
                    "terminal_map_to"  => ["1", "n"],
                    "tap_ratio"        => 1.05,
                    "regulator_type"   => "B",
                    "s_rating"         => s,
                    "r_series_from"    => 0.005 * zb,
                    "r_series_to"      => 0.005 * zb,
                    "x_series_from"    => 0.01  * zb))),
        "load" => Dict{String,Any}(
            "ld" => Dict{String,Any}(
                "bus"           => "reg",
                "terminal_map"  => ["1", "n"],
                "configuration" => "SINGLE_PHASE",
                "p_nom"         => [400_000.0],
                "q_nom"         => [100_000.0])))
end


function _se_net_open_delta_reg_original()
    # pf_open_delta_reg.dss: src ──[open-delta regulator, ABBC]── reg
    # Two single-phase line-to-line regulators (reg1 across A-B tap 1.05,
    # reg2 across B-C tap 1.025), 500 kVA each, 4.157 kV L-L, %r=0.5, xhl=1.0.
    # The monolithic BMOPF open_delta_regulator reproduces this with two cores.
    #   v_LL = 2400·√3 = 4156.9 V; z_base = v_LL² / S referred to the L-L coil.
    s    = 500_000.0
    vll  = 2400.0 * sqrt(3)
    zb   = vll^2 / s
    Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "src" => Dict{String,Any}(
                "terminal_names"  => ["1", "2", "3", "n"],
                "neutral_terminal"=> "n"),
            "reg" => Dict{String,Any}(
                "terminal_names"  => ["1", "2", "3", "n"],
                "neutral_terminal"=> "n")),
        "voltage_source" => Dict{String,Any}(
            "vs" => Dict{String,Any}(
                "bus"          => "src",
                "terminal_map" => ["1", "2", "3"],
                "v_magnitude"  => [2400.0, 2400.0, 2400.0],
                "v_angle"      => [0.0, -2π/3, 2π/3])),
        "shunt" => Dict{String,Any}(
            "g_src" => Dict{String,Any}(
                "bus" => "src", "terminal_map" => ["n"],
                "G_1_1" => 1000.0, "B_1_1" => 0.0),
            "g_reg" => Dict{String,Any}(
                "bus" => "reg", "terminal_map" => ["n"],
                "G_1_1" => 1000.0, "B_1_1" => 0.0)),
        "transformer" => Dict{String,Any}(
            "open_delta_regulator" => Dict{String,Any}(
                "od" => Dict{String,Any}(
                    "bus_from"         => "src",
                    "bus_to"           => "reg",
                    "terminal_map_from"=> ["1", "2", "3", "n"],
                    "terminal_map_to"  => ["1", "2", "3", "n"],
                    "connection"       => "ABBC",
                    "tap_ratio"        => [1.05, 1.025],
                    "regulator_type"   => "B",
                    "s_rating"         => s,
                    "r_series_from"    => 0.005 * zb,
                    "r_series_to"      => 0.005 * zb,
                    "x_series_from"    => 0.01  * zb))),
        "load" => Dict{String,Any}(
            "lda" => Dict{String,Any}(
                "bus" => "reg", "terminal_map" => ["1", "n"],
                "configuration" => "SINGLE_PHASE", "p_nom" => [20_000.0], "q_nom" => [0.0]),
            "ldb" => Dict{String,Any}(
                "bus" => "reg", "terminal_map" => ["2", "n"],
                "configuration" => "SINGLE_PHASE", "p_nom" => [20_000.0], "q_nom" => [0.0]),
            "ldc" => Dict{String,Any}(
                "bus" => "reg", "terminal_map" => ["3", "n"],
                "configuration" => "SINGLE_PHASE", "p_nom" => [20_000.0], "q_nom" => [0.0])))
end


function _se_net_center_tap()
    zh = 2400.0^2 / 25000.; zl = 120.0^2 / 25000.
    Dict{String,Any}(
        "bus" => Dict("hv"=>Dict("terminal_names"=>["1","n"],"perfectly_grounded_terminals"=>["n"]),
                      "lv"=>Dict("terminal_names"=>["1","n","2"],"perfectly_grounded_terminals"=>["n"])),
        "voltage_source" => Dict("s"=>Dict("bus"=>"hv","terminal_map"=>["1"],"v_magnitude"=>[2400.],"v_angle"=>[0.])),
        "transformer" => Dict("center_tap"=>Dict("ct"=>Dict{String,Any}(
            "bus_from"=>"hv","bus_to"=>"lv","terminal_map_from"=>["1","n"],"terminal_map_to"=>["1","n","2"],
            "v_nom_from"=>2400.,"v_nom_to"=>120.,"s_rating"=>25000.,
            "r_series_from"=>.005zh,"x_series_from"=>.01zh,"r_series_to"=>.005zl,"x_series_to"=>.01zl))),
        "load"=>Dict("a"=>Dict("bus"=>"lv","terminal_map"=>["1","n"],"configuration"=>"SINGLE_PHASE","p_nom"=>[6000.],"q_nom"=>[1000.]),
                     "b"=>Dict("bus"=>"lv","terminal_map"=>["2","n"],"configuration"=>"SINGLE_PHASE","p_nom"=>[2000.],"q_nom"=>[250.])))
end

# A finite grounded shunt provides a physical common-mode path for the
# open-delta oracle; constant-power wye loads alone are an unstable reference.
function _se_net_open_delta_reg()
    net = _se_net_open_delta_reg_original()
    net["shunt"]["reference"] = Dict{String,Any}(
        "bus"=>"reg", "terminal_map"=>["1","2","3"],
        "G_1_1"=>.001,"G_2_2"=>.001,"G_3_3"=>.001)
    net
end
