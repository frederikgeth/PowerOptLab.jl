# [Carrier harmonics through L and LCL filters](@id ibr-ac-harmonics-lcl)

<!-- study-executable -->

The network and device data are synthetic. After [installing the pinned
dependencies](../getting_started.md), paste the Julia blocks into a fresh
session in order, or run `julia --project=. scripts/run_study_tutorials.jl`
from the repository root. CI executes these published blocks on Julia 1.10
and current stable.

This tutorial follows the ideal PWM pole-voltage error through the physical AC
filter. It demonstrates why converter-side current, grid-side current, midpoint-
capacitor current, and neutral current must be reported separately for an LCL
converter.

The calculation is a local carrier-harmonic audit around a solved fundamental
operating point. It is more detailed than an RMS allowance, but it is not EMT or
a controller impedance scan.

## 1. Create an unbalanced four-wire operating point

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
   "v_magnitude":[245.0,215.0,230.0],"v_angle":[0.05,-2.15,2.0]}},
 "linecode":{"lc":{"R_series_1_1":0.05,"R_series_2_2":0.05,
   "R_series_3_3":0.05,"R_series_4_4":0.05}},
 "line":{"l1":{"bus_from":"grid","bus_to":"poc",
   "terminal_map_from":["a","b","c","n"],
   "terminal_map_to":["a","b","c","n"],
   "linecode":"lc","length":1.0}}}
"""; from_string=true)
```

## 2. Compare a reduced inductor with an explicit LCL circuit

The two series reactances sum to the reduced model's reactance. Therefore, with
no midpoint capacitor the two representations are identical.

Hold POC delivery at 15 kW and Q at zero for every filter and harmonic-grid
case. Use per-unit coordinates on a 20 kVA base, with bound relaxation disabled.
The synthetic `c_loss=0.02` W/A² coefficient selects a low-current solution.
This compares filters at a declared operating point rather than allowing their
maximum-export solves to select different powers. Raw-SI convergence is less
reliable for the coupled PWM-reserve equations; per-unit scaling changes the
numerical coordinates, while the physical checks remain in W, var and A.

```julia
common = (
    id="inverter", bus="poc",
    phase_terminals=["a","b","c"], neutral="n",
    topology=:FOUR_LEG,
    s_max=20e3, i_max=45.0, i_grid_max=45.0, In_max=100.0,
    v_dc=700.0, c_dc=1.1e-3, i_cap_max=50.0, c_loss=0.02,
    f_sw=10e3, pwm_strategy=:CENTERED, pwm_ac_ripple=true,
    pwm_fundamental_samples=36, pwm_carrier_samples=128,
    pwm_ac_harmonics=64,
)

operating_point = (
    objective=:min_loss, p_set=15e3, q_set=0.0,
    per_unit=true, s_base=20e3,
    solver_options=(tol=1e-7, acceptable_iter=0, max_iter=500,
                    bound_relax_factor=0.0),
)

function solve_filter(device)
    result = solve_advanced_inverter(network, device; operating_point...)
    @assert solve_status(result).publishable "$(device.id): $(result.termination_status)"
    @assert abs(result.p_poc - 15000.0) < 0.1
    @assert abs(result.q_poc) < 0.1
    @assert maximum(result.i_ac_total_rms) <= 45.0 + 0.01
    @assert maximum(result.i_grid_total_rms) <= 45.0 + 0.01
    @assert result.i_neutral_total_rms <= 100.0 + 0.01
    @assert all(result.i_ac_switching_reserved .>= result.i_ac_switching_rms .- 1e-5)
    @assert all(result.i_grid_switching_reserved .>= result.i_grid_switching_rms .- 1e-5)
    return result
end

reduced = AdvancedInverter(; r_filter=0.05, x_filter=0.15, common...)
series_only = AdvancedInverter(; r_filter=0.02, x_filter=0.05,
    r_filter_grid=0.03, x_filter_grid=0.10, common...)
lcl = AdvancedInverter(; r_filter=0.02, x_filter=0.05,
    r_filter_grid=0.03, x_filter_grid=0.10,
    c_filter_mid=20e-6, r_filter_damping=1.0, common...)

r_reduced = solve_filter(reduced)
r_series = solve_filter(series_only)
r_lcl = solve_filter(lcl)
@assert all(r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
            for r in (r_reduced, r_series, r_lcl))

@show r_reduced.i_ac_switching_rms r_reduced.i_grid_switching_rms
@show r_series.i_ac_switching_rms r_series.i_grid_switching_rms
@show r_lcl.i_ac_switching_rms r_lcl.i_grid_switching_rms
@show r_lcl.i_filter_shunt_switching_rms
@assert maximum(r_lcl.i_grid_switching_rms) < maximum(r_lcl.i_ac_switching_rms)
@assert all(r_lcl.i_ac_switching_reserved .>= r_lcl.i_ac_switching_rms .- 1e-5)
@assert all(r_lcl.i_grid_switching_reserved .>= r_lcl.i_grid_switching_rms .- 1e-5)
@assert maximum(r_lcl.i_ac_total_rms) <= 45.0 + 0.01
@assert maximum(r_lcl.i_grid_total_rms) <= 45.0 + 0.01
```

Check the limiting case numerically:

```julia
@assert all(isapprox.(r_reduced.i_ac_switching_rms,
                      r_series.i_ac_switching_rms; rtol=1e-6))
@assert all(isapprox.(r_reduced.i_grid_switching_rms,
                      r_series.i_grid_switching_rms; rtol=1e-6))
```

With the midpoint capacitor present, converter current splits between the shunt
branch and grid arm. A low grid-side ripple does not imply low semiconductor or
filter-capacitor stress.

## 3. Compare fundamental, switching, and total RMS currents

```julia
for phase in eachindex(r_lcl.i_mag)
    @assert isapprox(r_lcl.i_ac_total_rms[phase],
        hypot(r_lcl.i_mag[phase], r_lcl.i_ac_switching_rms[phase]);
        rtol=1e-10)
    @assert isapprox(r_lcl.i_grid_total_rms[phase],
        hypot(r_lcl.i_grid_mag[phase], r_lcl.i_grid_switching_rms[phase]);
        rtol=1e-10)
end

@show r_lcl.i_mag r_lcl.i_ac_switching_rms r_lcl.i_ac_total_rms
@show r_lcl.i_grid_mag r_lcl.i_grid_switching_rms r_lcl.i_grid_total_rms
@show r_lcl.i_ac_switching_reserved r_lcl.i_grid_switching_reserved
```

The sequential closure allocates the predicted carrier RMS inside `i_max`,
`i_grid_max`, and `In_max`. Confirm that each `*_reserved` diagnostic covers its
prediction before publishing the point.

## 4. Show why the neutral inductor matters

```julia
no_neutral_l = AdvancedInverter(; r_filter=0.05, x_filter=0.15,
    x_filter_neutral=0.0, common...)
with_neutral_l = AdvancedInverter(; r_filter=0.05, x_filter=0.15,
    x_filter_neutral=0.30, common...)

r_no_ln = solve_filter(no_neutral_l)
r_with_ln = solve_filter(with_neutral_l)

@show r_no_ln.i_neutral_switching_rms r_with_ln.i_neutral_switching_rms
@show r_no_ln.i_neutral_total_rms r_with_ln.i_neutral_total_rms
@assert solve_status(r_no_ln).publishable && solve_status(r_with_ln).publishable
@assert r_with_ln.i_neutral_switching_rms < r_no_ln.i_neutral_switching_rms
```

The fourth-wire impedance participates in a shared return path, so it changes
all phase-to-neutral harmonic equations. It must not be added as an independent
scalar correction after solving three uncoupled phase circuits.

For coupled or unequal conductors, provide the complete primitive matrices in
`a,b,c,n` order:

```julia
Xprimitive = [0.15 0.01 0.01 0.02;
              0.01 0.15 0.01 0.02;
              0.01 0.01 0.15 0.02;
              0.02 0.02 0.02 0.30]
```

Use `x_filter_matrix=Xprimitive` instead of the scalar phase/neutral reactances;
the two parameterisations cannot be mixed.

## 5. Refine the retained harmonic series

```julia
for (samples, harmonics) in ((128,32), (128,64), (256,128))
    parameters = merge(common, (r_filter=0.02, x_filter=0.05,
        r_filter_grid=0.03, x_filter_grid=0.10,
        c_filter_mid=20e-6, r_filter_damping=1.0,
        pwm_carrier_samples=samples, pwm_ac_harmonics=harmonics,
    ))
    device = AdvancedInverter(; parameters...)
    r = solve_filter(device)
    @assert r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    println((samples=samples, harmonics=harmonics,
        converter=maximum(r.i_ac_switching_rms),
        grid=maximum(r.i_grid_switching_rms),
        shunt=maximum(r.i_filter_shunt_switching_rms),
        peak_to_peak=maximum(r.i_ac_switching_pp)))
end
```

RMS quantities normally converge sooner than peak-to-peak reconstructions.
Refine both fundamental-angle and carrier grids when an operating point lies
close to a physical rating.

## 6. Interpret the resonance diagnostic correctly

`filter_resonance_hz` is the undamped scalar LCL natural-frequency estimate. It
helps catch implausible component selections, but it is not a stability margin.
The carrier audit reuses constant R/L/C primitives at harmonic frequencies and
assumes a stiff high-frequency POC. It omits:

- controller and digital-delay impedance;
- frequency-dependent winding resistance, core loss, and saturation;
- capacitor ESR/ESL unless represented separately in a future model;
- common-mode capacitance and protective-earth paths; and
- switching dead time and device capacitances.

Use an impedance scan or EMT model when those mechanisms determine the result.

## 7. Publication checks

Report converter, grid, shunt, and neutral paths separately; demonstrate
harmonic-grid convergence; verify reserve coverage and total-RMS ratings; state
the assumed high-frequency grid boundary; and avoid describing the scalar
resonance estimate as a closed-loop stability result.
