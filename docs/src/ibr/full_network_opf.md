# Run a full network OPF with advanced inverters

**Yes: the constrained inverters can participate in a shared network OPF.**
The one-device `solve_constrained_inverter` function is a convenience entry
point. The underlying `stamp_device!` interface lets several devices share
network equations, operational limits and an objective selected by the caller.
The network voltages and all inverter operating points are solved together.

This worked case connects the [controller-to-equilibrium workbook](../tutorials/controller_program_to_equilibrium.md)
to a circuit-aware model. It answers three separate validation questions:

| Question | Comparison |
|:--|:--|
| Do the two code paths encode the same equations when no limit binds? | Forward controls versus maximum export with matching constraints |
| Does the solved Q follow the declared volt-var law? | Exact and smoothed laws evaluated independently at the replayed voltages |
| Does coordinated dispatch meet the additional network requirements? | Free-Q OPF versus local controls at the same high availability and a feasible benchmark |
| Is the electrical solution consistent with an independent implementation? | Exported P/Q replay through a four-wire radial power-flow solver |

## The seven-bus case

```text
source ── b1 ── b2 [PV1] ── b4 [PV2]
           │
           └── b3 ── b5
                │
                └── b6 [PV3]
```

Each line carries three phase conductors and a neutral, with a full complex
primitive impedance matrix including mutual coupling. Only the source neutral
is grounded. The six downstream buses have unequal constant-PQ phase-to-neutral
loads, so their neutral voltages and phase magnitudes differ. There are no line
shunts or transformers in this particular validation fixture.

Three three-leg inverters have ratings of 6, 6 and 5 kVA, nonzero series
filters, fixed and current-quadratic losses, DC links and physical capability
constraints. Their converter currents are explicitly balanced; the topology
alone would require only zero current sum. Average-voltage, positive-sequence
and worst-phase voltage policies share a sloping volt-var curve. Every model
retains the nonlinear constant-PQ load and network equations.

The fixture is in `scripts/formulations/inverter_network_opf.jl`. Run:

```sh
julia --project=. scripts/formulations/inverter_network_opf.jl
```

## Compose several inverters in one OPF

This is the essential public-API pattern. The example utilities supply the
network and device data; the composition uses the ordinary staged builder.

```@example full_compose
using PowerOptLab, JuMP, Ipopt
import BMOPFTools
include(joinpath(pkgdir(PowerOptLab), "scripts", "formulations", "inverter_network_opf.jl"))
case = inverter_network_case()

devices = [ConstrainedInverter(plant, VoltVarWattConstraints(policy);
               dc_source=DCSourceBounds(p_min=0.))
           for (plant, policy) in zip(case.plants, case.policies)]
requests = [InverterControlRequest(p_available=available,
                p_rated=plant.s_max, q_scale=plant.s_max)
            for (plant, available) in zip(case.plants, case.availability)]
handles = Dict{Int,Any}()

hook! = ctx -> begin
    for i in eachindex(devices)
        handles[i] = stamp_device!(ctx, devices[i]; request=requests[i])
    end
    physical = [inverter_handles(devices[i], handles[i]) for i in eachindex(devices)]
    model = BMOPFTools.opf_model(ctx)
    set_optimizer_attribute(model, "tol", 1e-9)
    @objective(model, Max, sum(h.p_poc*h.sb/10000 for h in physical))
end
ctx = BMOPFTools.build_opf_model(case.net; per_unit=true, s_base=1e6,
    add_objective=false, model_hook! = hook!, optimizer=Ipopt.Optimizer)
BMOPFTools.optimize_opf!(ctx)
status = solve_status(BMOPFTools.opf_model(ctx))
results = [extract_device(devices[i], handles[i], status) for i in eachindex(devices)]
@assert status.publishable && length(results)==3
[(P_W=r.plant.p_poc,Q_var=r.plant.q_poc) for r in results]
```

The handles use model units; multiplying P by `h.sb` converts it to watts.
Dividing this entire objective by 10000 only changes numerical scaling. The
stamp adds no objective. A caller can instead select a cost, tracking or other
dispatch criterion and combine it with other network devices and constraints.

`solve_status(model)` applies the package's normal publication contract to a
custom JuMP model. Passing it to `extract_device` masks numerical outputs when
the solve is unsuccessful. `inverter_handles` works with the returned handles
of both controlled and constrained inverter wrappers.

This fixture has **no native IBR ports**. Stamping a custom plant does not
automatically replace an imported native IBR with the same identifier. When
adapting an imported network, establish component ownership explicitly or
remove the intended native components from a copy before adding their custom
replacements. Preserve unrelated native devices; count each physical port once.

## Check 1: matching controls and operating constraints

First give the inverters DC availability of 3500, 4500 and 3000 W. The
forward model uses `SequenceController` with converter-current targets. The
optimisation model uses the same smooth voltage laws, sensing reference,
balanced-current allocation and DC availability, then maximises total POC P.
The plant has no LCL midpoint or POC shunt here, so both Q definitions use the
same grid current at the POC. Inverter S/current saturation is inactive.

The broad native terminal-to-ground bounds are 200–280 V in this comparison.
It does **not** impose the tighter network requirements introduced in Check 2.

```@example full_inverter_opf
using PowerOptLab
include(joinpath(pkgdir(PowerOptLab), "scripts", "formulations", "inverter_network_opf.jl"))
validation = inverter_full_opf_validation()
@assert validation.forward.status.publishable && validation.matched.status.publishable
@assert validation.agreement.p_error_W < .2
@assert validation.agreement.q_error_var < .1
@assert validation.agreement.voltage_error_V < .001
validation.agreement
```

Maximum per-inverter discrepancies are approximately **0.12 W and 0.04 var**;
the largest complex bus-voltage difference is **0.00018 V**. Branch and converter
currents also agree within the declared tolerances.

Read this as a **consistency check between two code paths**, not as physical
evidence. With saturation inactive, maximum export under `p_dc ≤ p_available`
binds the source budget, which is exactly what the forward controller
commands; both models then reduce to the same equations by construction. The
residual is dominated by the forward limiter's smooth common scale
(``1-s \approx 2.6\times10^{-5}`` here), which multiplies the whole command. The
regime in which the formulations are meant to differ, saturation, is covered by
the single-device oracles and the [reproduction cases](formulation_reproductions.md).

### Does Q follow the declared law?

Satisfying the encoded Q equality is automatic, so the audit evaluates the
volt-var law **without the package's controller code**, at the independently
replayed voltages (local-neutral reference). It reports two comparisons:

- against a re-derivation of the smooth encoding: agreement to about
  ``10^{-11}`` var for the optimisation model, and exactly the limiter scale
  ``s`` for the forward model;
- against the **exact** piecewise-linear law with hard extrema and hard
  positive/negative parts.

```@example full_inverter_opf
[(inverter=i, Q_var=round(validation.matched.devices[i].plant.q_poc,digits=2),
  exact_law_var=round(x.exact_q_var,digits=2),
  deviation_var=round(x.exact_error_var,sigdigits=3),
  analytic_bound_var=x.smoothing_bound_var)
 for (i,x) in enumerate(validation.law_audit)]
```

The average-voltage and positive-sequence inverters follow the exact law,
because their sensed voltages lie several volts from any curve breakpoint. The
worst-phase `:net` inverter deviates by about **11 var (0.23% of rating)**. At
this operating point its low-voltage branch requests only about 0.006 p.u.,
inside the `conflict_epsilon = 0.01` width of the smooth positive part. The
algebraic smooth max/min errs by at most ``\varepsilon_c/2`` per branch, and
the two-stage smooth phase extrema by at most ``\varepsilon_v``, giving the
bound ``S_{rated}(\varepsilon_c/2 + |\mathrm{d}q/\mathrm{d}V|\,\varepsilon_v)
= 32.5`` var. The forward controller shares this deviation; reduce the widths
when the exact firmware law matters more than solver robustness.

## Check 2: coordinate dispatch under binding network limits

Now increase availability to each inverter's rating. The autonomous comparison
retains its local voltage laws without adding the new service requirements.
The coordinated OPF instead uses free Q within ±40% of each rating and imposes:

- phase-to-neutral service voltage between 220 and **244 V** at every downstream bus;
- a **12 A** total fundamental conductor-current rating on all four feeder-head conductors;
- the same inverter filter, loss, source and hardware constraints.

It maximises an illustrative weighted-export objective with weights
``[1,1.25,0.9]`` and a quadratic reactive-use penalty. These are declared
planning preferences, not firmware or market rules. The objective is

```math
\sum_i w_i\frac{P_i}{10000\,\mathrm W}
-0.05\sum_i\left(\frac{Q_i}{S_{rated,i}}\right)^2.
```

The code normally uses native `vpn_min` / `vpn_max` bus fields. **`v_min` /
`v_max` bound terminal-to-ground voltage**, which differs when the neutral
moves. An alternative hook formulation writes the phase-to-neutral squared
magnitudes directly. Tests compare their dispatch, objective and complex
voltages; both formulations select the same solution within tolerance.

```@example full_inverter_opf
opf = validation.dispatch
audit = validation.audits.dispatch
@assert opf.status.publishable
@assert abs(audit.phase_neutral_max_V-244.) < 1e-5
@assert abs(audit.feeder_max_A-12.) < 1e-6
@assert validation.audits.unmanaged.phase_neutral_max_V > 247.
[(inverter=i, P_W=d.plant.p_poc, Q_var=d.plant.q_poc)
 for (i,d) in enumerate(opf.devices)]
```

The resulting P/Q values are approximately:

| Inverter | P (W) | Q (var) |
|:--|--:|--:|
| PV1 | 5502.17 | −2400.00 |
| PV2 | 2400.54 | −2400.00 |
| PV3 | 4959.28 | −600.13 |

Negative Q means reactive absorption. Both the 244 V and 12 A network limits
bind; converter apparent-power bounds bind too. Downstream neutral voltage
reaches about 1.6 V, making the voltage reference material to the result.

At the **same high availability**, the local-controller comparison reaches
about 247.94 V and 13.29 A. It violates the newly declared network requirements.
The coordinated solution is deliberately a different operating problem, so
agreement with that autonomous dispatch would be the wrong validation target.

The example also solves a feasible benchmark fixing each inverter at 2000 W
and zero Q under the same coordinated constraints. The optimised objective is
better than that benchmark. This is a useful feasibility/improvement check;
it does not establish a globally optimal solution of the nonconvex AC OPF.

## Check 3: independently replay the exported dispatch

The replay has no JuMP, BMOPFTools network equations or controller evaluator.
It takes the physical impedances, loads and solved P/Q as inputs. At each trial
voltage it calculates constant-PQ load currents and balanced inverter current:

```math
I_{load,k}=\overline{\frac{S_{load,k}}{V_k-V_n}},\qquad
I_{1,inv}=\overline{\frac{S_{inv}}{3U_1}}.
```

It reconstructs phase and neutral demand, sums branch currents backward through
the radial tree, then updates four-conductor voltages forward using
``V_{child}=V_{parent}-ZI_{branch}``. Iteration continues until the voltage
change is below the declared tolerance. Replaying fixed P/Q checks electrical
consistency, rather than reproducing an inverter's autonomous law.

```@example full_inverter_opf
for name in (:forward,:matched,:unmanaged,:dispatch,:benchmark)
    a = getproperty(validation.audits,name)
    @assert a.replay.converged
    @assert a.voltage_error_V < 1e-5 && a.current_error_A < 1e-5
    @assert a.kcl_error_A < 1e-5 && a.kvl_error_V < 1e-6
    @assert a.power_balance_error_VA < .01
end
[(case=name, delta_V_V=getproperty(validation.audits,name).voltage_error_V,
  delta_I_A=getproperty(validation.audits,name).current_error_A,
  power_residual_VA=getproperty(validation.audits,name).power_balance_error_VA)
 for name in (:forward,:matched,:unmanaged,:dispatch,:benchmark)]
```

The residuals are near floating-point precision for this fixture. Tests use
the wider physical-unit tolerances above. They additionally check filter and
semiconductor losses, converter-terminal power, POC power, source budgets and
capability margins. Changing the per-unit power base from 1 MVA to 10 kVA
preserves the solution. A deliberately interrupted solve cannot publish a
dispatch or a successful audit.

## What is ready, and what remains explicit?

| Capability | Current route |
|:--|:--|
| Several inverters, full network equations, custom objective | Shared staged model; executable example above |
| THREE_LEG, FOUR_LEG, SPLIT_DC capability constraints | Supported by `ConstrainedInverter`; this larger agreement fixture validates THREE_LEG |
| Native network limits and additional operating constraints | Native fields or explicit hook constraints, with the intended reference/units |
| Multiple periods | `build_multi_context` plus per-period requests and explicit linking constraints; an inverter alone has no stored-energy state |
| Dedicated constrained-fleet/import/scheduling convenience API | Not supplied; the current `ControlledInverterFleetSpec` describes forward controls |
| Automatic PWM reserve closure, single-phase wrapper, GFM operating laws, complete exact-PWL controller graphs | Remain outside this wrapper's supported scope |

This example validates classic per-unit coordinates and a radial four-wire
constant-PQ network. Its independent replay does not cover transformers, shunts,
meshes, LCL current allocation or unbalanced inverter sequence dispatch. Those
need their corresponding oracle assumptions. The existing
[LCL/sequence reproduction cases](formulation_reproductions.md) and
[feature comparison table](../components/inverter_operating_constraints.md#Simulation-and-optimisation-formulations)
document the separate coverage. Successful local solves and agreement checks
provide no dynamic-stability, uniqueness or global-optimality guarantee.
