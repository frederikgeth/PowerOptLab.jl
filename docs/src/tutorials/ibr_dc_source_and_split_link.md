# [DC-source impedance and split-link carrier stress](@id ibr-dc-source-and-split-link)

This tutorial shows how bridge switching current divides between the DC-link
capacitor and a finite upstream source impedance. It then resolves the resulting
voltage ripple across unequal upper and lower split-link banks.

The default model assumes the source is effectively open at carrier frequency.
Supplying `pwm_dc_source_r` and/or `pwm_dc_source_l` replaces that assumption with
a series R–L source branch in parallel with the link capacitor.

## 1. Reuse a complete balanced network

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf

network = parse_bmopf("""
{"bus":{
  "grid":{"terminal_names":["a","b","c","n"],
          "perfectly_grounded_terminals":["n"]},
  "poc":{"terminal_names":["a","b","c","n"],
         "perfectly_grounded_terminals":["n"],
         "v_min":[180.0,180.0,180.0],"v_max":[270.0,270.0,270.0]}},
 "voltage_source":{"vs":{"bus":"grid","terminal_map":["a","b","c"],
   "v_magnitude":[230.0,230.0,230.0],"v_angle":[0.0,-2.0944,2.0944]}},
 "linecode":{"lc":{"R_series_1_1":0.05,"R_series_2_2":0.05,
   "R_series_3_3":0.05,"R_series_4_4":0.05}},
 "line":{"l1":{"bus_from":"grid","bus_to":"poc",
   "terminal_map_from":["a","b","c","n"],
   "terminal_map_to":["a","b","c","n"],
   "linecode":"lc","length":1.0}}}
"""; from_string=true)
```

## 2. Compare open, high-, and finite-impedance sources

```julia
function four_leg(; source_r=nothing, source_l=0.0, harmonics=64)
    AdvancedInverter(
        id="four-leg", bus="poc",
        phase_terminals=["a","b","c"], neutral="n",
        topology=:FOUR_LEG,
        s_max=20e3, i_max=40.0, In_max=40.0,
        v_dc=700.0, c_dc=1.1e-3, i_cap_max=15.0,
        r_filter=0.05, x_filter=0.15, m_max=0.96,
        f_sw=10e3, pwm_strategy=:CENTERED,
        pwm_fundamental_samples=72, pwm_carrier_samples=128,
        pwm_dc_harmonics=harmonics,
        pwm_dc_source_r=source_r, pwm_dc_source_l=source_l,
    )
end

cases = [
    "open" => four_leg(),
    "high Z" => four_leg(source_r=1e6),
    "finite R" => four_leg(source_r=0.01),
    "finite RL" => four_leg(source_r=0.02, source_l=0.2e-6),
]

function dc_summary(label, device)
    r = solve_advanced_inverter(network, device)
    @assert r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    return (
        case=label,
        bridge_a=r.i_dc_bridge_switching_rms,
        capacitor_a=r.i_cap_switching,
        source_a=r.i_dc_source_switching_rms,
        ripple_rms_v=r.dv_switching_rms,
        source_loss_w=r.p_dc_source_switching_loss,
        network_margin=r.pwm_dc_network_margin,
        reserve_margin_a=r.pwm_reserve_margin,
    )
end

dc_results = [dc_summary(label, device) for (label, device) in cases]
foreach(println, dc_results)
```

The high-impedance case should converge to the open-source capacitor current and
voltage ripple. A low *resistive* source impedance diverts part of each retained
harmonic away from the capacitor. An inductive source can instead create
circulating-current magnification near parallel antiresonance, so capacitor RMS
may exceed bridge RMS even though harmonic KCL is exact. RMS magnitudes do not
generally add arithmetically: KCL holds for the complex current at every harmonic
before RMS aggregation.

For a resistive source branch, verify the independent loss identity:

```julia
r_finite = solve_advanced_inverter(network, four_leg(source_r=0.01))
@assert isapprox(r_finite.p_dc_source_switching_loss,
    0.01*r_finite.i_dc_source_switching_rms^2; rtol=1e-8)
```

This loss is upstream of the converter and is deliberately excluded from
`p_dc`. Add it to a system boundary only when that boundary includes the source
impedance.

## 3. Check spectral convergence

The source-current and voltage diagnostics use the retained Fourier series.
Unretained bridge-current energy is assigned to the capacitor, conservatively,
for thermal reserve closure.

```julia
for h in (16, 32, 64)
    r = solve_advanced_inverter(network,
        four_leg(source_r=0.02, source_l=0.2e-6, harmonics=h))
    @assert r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    println((harmonics=h,
        cap=r.i_cap_switching,
        source=r.i_dc_source_switching_rms,
        ripple=r.dv_switching_pp,
        margin=r.pwm_dc_network_margin))
end
```

Increase `pwm_carrier_samples` when increasing the harmonic ceiling; always
maintain `pwm_dc_harmonics ≤ pwm_carrier_samples/2`. Peak-to-peak voltage is
usually more truncation-sensitive than RMS current.

## 4. Resolve unequal split-bank voltages

```julia
cu, cl = 2.5e-3, 3.5e-3
split = AdvancedInverter(
    id="split", bus="poc",
    phase_terminals=["a","b","c"], neutral="n",
    topology=:SPLIT_DC,
    s_max=20e3, i_max=40.0, In_max=40.0,
    v_dc=900.0, c_dc=3e-3,
    c_dc_upper=cu, c_dc_lower=cl, i_cap_max=25.0,
    r_filter=0.05, x_filter=0.15, m_max=0.96,
    f_sw=12e3, pwm_strategy=:SPWM,
    pwm_fundamental_samples=72, pwm_carrier_samples=128,
    pwm_dc_harmonics=64, pwm_dc_source_r=0.05,
)

r_split = solve_advanced_inverter(network, split)
@assert r_split.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
ceq = cu*cl/(cu+cl)

@assert isapprox(r_split.dv_switching_upper_rms,
    ceq/cu*r_split.dv_switching_rms; rtol=1e-10)
@assert isapprox(r_split.dv_switching_lower_rms,
    ceq/cl*r_split.dv_switching_rms; rtol=1e-10)

@show r_split.dv_switching_upper_rms r_split.dv_switching_lower_rms
@show r_split.dv_switching_upper_pp r_split.dv_switching_lower_pp
```

The smaller upper bank sees the larger voltage ripple. This voltage allocation
is distinct from the neutral-current sharing used in the low-frequency thermal
budget; both should be inspected when the half-banks differ.

## 5. Screen source-capacitor antiresonance

For a lossless source inductor, parallel cancellation occurs near
`2πf = 1/sqrt(Ls*Ceq)`. The public diagnostic is
`pwm_dc_network_margin`; zero is singular and values close to zero are suspect.

```julia
f_sw = 10e3
ceq = 1.1e-3
l_resonant = 1 / ((2pi*f_sw)^2 * ceq)
@show l_resonant
```

Do not solve exactly at this ideal lossless point. Add physically justified
damping and use a measured impedance spectrum when the source is a battery,
DC/DC converter, cable, or controlled rectifier. The scalar margin screens one
passive R–L/C network; it is not a Nyquist stability margin.

## 6. Publication checks

Require all of the following:

- `pwm_reserve_margin ≥ 0` within the selected closure tolerance;
- `pwm_modulation_margin ≥ 0`;
- finite `pwm_dc_network_margin` comfortably above zero;
- stable RMS and peak-to-peak diagnostics under a sampling/harmonic refinement;
- documented source-impedance frequency range and system boundary for its loss;
- individual split-bank voltage and thermal checks when `Cu != Cl`.
