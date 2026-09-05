# Select smoothing in a physical inverter model

Start with [control intent, compilation and solver pathways](compilation.md) for
the separation between physical semantics, encodings and backend configuration.

`PiecewiseLinearLaw` now accepts `formulation=...`. The legacy
`smoothing_epsilon=...` constructor still selects BMOPFTools softplus. Each
volt-watt, volt-var or negative-sequence gain curve may select its own smoothing
family. This changes the curve representation; the controller's extrema, limiter,
voltage floor and plant-capability equations retain their separately configured
behavior. A C2 curve therefore does not make every other operation a compact
spline or eliminate every source of numerical difficulty.

## Define the physical case

This executable example uses a deliberately small balanced feeder, an existing
three-leg inverter plant and an average-voltage volt-watt law. It illustrates
integration, not a claim about representative distribution-network performance.

```@example controller_formulation
using PowerOptLab, JuMP, Ipopt
using BMOPFTools: parse_bmopf
network = parse_bmopf("""
{"bus":{
 "grid":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"]},
 "poc":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"],
        "v_min":[180,180,180],"v_max":[270,270,270]}},
 "voltage_source":{"vs":{"bus":"grid","terminal_map":["a","b","c"],
   "v_magnitude":[230,230,230],"v_angle":[0,-2.0943951023931953,2.0943951023931953]}},
 "linecode":{"lc":{"R_series_1_1":0.05,"R_series_2_2":0.05,
                     "R_series_3_3":0.05,"R_series_4_4":0.05}},
 "line":{"l":{"bus_from":"grid","bus_to":"poc",
   "terminal_map_from":["a","b","c","n"],"terminal_map_to":["a","b","c","n"],
   "linecode":"lc","length":1}}}
""";from_string=true)
request = InverterControlRequest(p_available=2e3,q_scale=4e3)
intent = VoltVarWattIntent(
    volt_watt=PWLFunction([220.,240.,250.,270.],[1.,1.,.2,.2]),
    sensing=:average_voltage)
function make_device(family,config)
    rep = smoothing_for_error(intent.volt_watt,family,config.fraction_error)
    encoding = VoltVarWattEncoding(volt_watt=rep)
    controller = SequenceController(lower_positive_policy(intent,encoding))
    plant = AdvancedInverter(id="research",bus="poc",phase_terminals=["a","b","c"],
        neutral="n",topology=:THREE_LEG,s_max=20e3,i_max=40.,v_dc=700.,
        c_dc=1.1e-3,r_filter=.05,x_filter=.15)
    return ControlledDevice(plant,controller)
end
case = controlled_inverter_case(network,make_device,request)
methods = [FormulationMethod(string(family),family,Ipopt.Optimizer;
    options=(tol=1e-8,max_iter=400),configure! = set_silent,
    metadata=(optimizer="Ipopt",family=string(family)))
    for family in (SoftplusFormulation,LocalC2Formulation)]
configs = [(fraction_error=1e-4,s_base=1e6,selection_objective=:loss)]
rows = run_formulation_experiment([case],methods;configurations=configs,on_error=:throw)
@assert all(r -> r["strict_solver_success"],rows)
[(r["method"],r["metrics"].pre_capability_target_gap_A) for r in rows]
```

The assertion checks this tutorial's representative solve; it is not the
runner's acceptance policy. The adapter builds a fresh staged model and copies
the network for each run. `network` and `request` can also be configuration
callbacks, so voltage, impedance, available power or hardware can vary without
editing the runner. `device_builder` can read any user-defined configuration keys.
Choose a custom `FormulationCase` when you need fleet/DOE construction, arbitrary
starts, a different objective, or extra model constraints.

## Interpret what was measured

The recorded phasors are candidate POC quantities in volts and amperes. They can
be present even for an unsuccessful solve. `pre_capability_target_gap_A` compares
numeric exact and smooth controller targets at the **same candidate voltage**,
before the plant-specific capability backoff. It is neither the distance between
two equilibria nor a residual of the complete plant/controller equations.
Requested active/reactive powers are also reported. A custom metric callback
`metrics(ctx, handles, device, request)` can add quantities relevant to your study;
its record is kept under `metrics.custom`.

For the established full solved-result interface, pass the same device to
`solve_controlled_inverter`. Its `exact_smooth_current_residual` includes the
plant-capability comparison at the solved phasors and is published only after a
strictly successful outcome. See [the control-study methodology](../ibr/control_study_methodology.md)
for the distinction between same-state replay, hardware feasibility and network
fixed-point evidence. The existing network fixed-point/sensitivity interfaces
remain available; this branch introduces no new continuation or re-solve method.

## Tail semantics and faithful evaluation

A controller curve keeps its endpoint value beyond its breakpoint range.
`PWLFunction` defaults to bounded numeric evaluation and `formulate_pwl!` imposes
that bounded domain. Controllers instead use `smooth_pwl_expression`, which
builds the smoothing of the **flat extension** without artificial voltage bounds.
The numeric `evaluate_smooth` path uses the same selected family as stamping.

```@example controller_formulation
law = PiecewiseLinearLaw([220.,240.,250.],[1.,1.,.2];
    formulation=LocalC2Formulation(.05))
canonical = PWLFunction(law;input_unit=:V,output_unit=:pu)
@assert primitive_value(canonical,280.;domain_policy=:flat_extension)==.2
@assert isapprox(primitive_value(canonical,280.,law.formulation;
    domain_policy=:flat_extension),.2)
nothing # hide
```

Compact smoothing equals the canonical tail sufficiently far from the endpoint;
softplus approaches it asymptotically. Close to a breakpoint, mixed slope changes
can make a complete surrogate either exceed or undershoot the canonical curve.
Preserve physical hardware constraints and examine controller composition rather
than interpreting the curve's error budget as a safety certificate.

Exact, convex-hull and complementarity graph objects are intentionally not
accepted as `PiecewiseLinearLaw` smooth families. Their use changes the problem
class. The scalar case already supports these objects; a researcher extending
an AC controller should build and document the graph coupling explicitly and
select a backend capable of the resulting nonlinear/integer/MPCC model.

## What the advanced IBR refactor shares

The production controller equations now call the same primitive layer used by
standalone experiments. The following table identifies the physical role of each
construction; these roles explain why one smoothing family should not simply be
substituted for every square root in a model.

| Operation in the existing IBR model | Shared construction | Preserved interpretation |
|:--|:--|:--|
| Volt-watt, volt-var and unbalance gain curves | `smooth_pwl_expression` / `primitive_value` | Selected family and flat tails; legacy softplus default |
| Phase extrema, available-power minimum and positive/negative conflict branches | `selector_expression` / `selector_value` with `AlgebraicFormulation` | Existing physical widths and binary operation order |
| Symmetric P/Q priority clipping | `symmetric_clip_expression` / `symmetric_clip_value` | Existing clipping formula |
| Composition of successive capability scale factors | `selector_expression(...; kind=:nonnegative_min)` | Zero-preserving, nonnegative under-approximation of the minimum |
| Current, apparent-power and ripple capability denominators | `MagnitudeApproximation(...; direction=:upper)` | Avoid creating headroom through magnitude underestimation |
| Negative-sequence voltage signal | `MagnitudeApproximation(...; direction=:lower)` | Zero at balance, existing eta error budget |
| Advanced inverter linear-current loss | `MagnitudeApproximation(...; direction=:lower)` | Existing per-leg ampere budget and loss deficit bound |

The physical squared capability constraints, sequence transformations, voltage
floors, plant-aware allocation and controller tie/conflict policies remain
controller/plant logic. The shared layer does not infer these semantics. Named
BMOPFTools differentiability annotations are retained for the shifted norms, and
normalized auxiliary rows still prevent repeated expression expansion.

Per-curve smoothing families and existing width/rating settings remain
configurable. Selectors in the established controller retain their algebraic
family. The public selector builders support softplus, local C2 and custom hinge
families for experiments; changing the entire controller's selector family would
require an explicit policy change and assessment of its composed behavior.

### Why some roots remain exact and implicit

Phase/positive-sequence voltage magnitudes and P/Q priority headroom retain their
nonnegative implicit roots. Replacing `u²=radicand, u≥0` by a shifted norm changes
the quantity; replacing it by `positive_root_expression` imposes a positive
radicand floor. Neither is an algebraically neutral refactor. At zero the implicit
root may be degenerate, so this preservation is not a claim that all numerical
conditioning issues are solved. Use a domain restriction only when it is part of
the intended physical/control model, and document it as such.
