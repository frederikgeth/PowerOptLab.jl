# From control intent to solver equations

For illustrated bound-dependent equality/inequality specialization, continue to
[bounds and relation-aware lowering](bounds_and_relations.md). The companion
[tolerance and solver guide](tolerances_and_solvers.md) distinguishes OpenDSS
control-loop criteria, approximation budgets and whole-model solver capabilities.

A volt-var/watt characteristic specifies what the controller is meant to do.
A softplus, local spline, exact graph or relaxation specifies how we represent
that characteristic in an optimization model. Keep those decisions separate:
changing the representation should not silently change voltage sensing, power
bases, sign conventions or firmware conflict handling.

This page follows one intent through numeric evaluation, full smooth IBR
construction, and bounded curve-port experiments with different solvers. The
examples execute during the documentation build. Optional solver examples also
run in CI on the supported Julia versions.

## Declare the intent once

```@example control_compilation
using PowerOptLab, JuMP, Ipopt
intent = VoltVarWattIntent(
    volt_watt=PWLFunction([220.,240.,250.,270.],[1.,1.,0.,0.]),
    volt_var=PWLFunction([210.,220.,240.,250.],[.5,0.,0.,-.5]),
    sensing=:average_voltage,
    volt_watt_basis=:available)
encoding = VoltVarWattEncoding(
    volt_watt=smoothing_for_error(intent.volt_watt,SoftplusFormulation,1e-3),
    volt_var=smoothing_for_error(intent.volt_var,LocalC2Formulation,1e-3),
    extrema_epsilon=.05)
policy = lower_positive_policy(intent,encoding)
controller = SequenceController(policy)
(typeof(policy),encoding.volt_watt.width,encoding.volt_var.width)
```

`VoltVarWattIntent` stores immutable physical curves, sensing, the active-power
basis and the established conflict/guard semantics. Inputs are volts; outputs
are dimensionless fractions. Volt-watt scales available or rated power;
volt-var scales `InverterControlRequest.q_scale`, with positive reactive power
meaning injection. Curves must be non-increasing; watt fractions lie in `[0,1]`.
Unit symbols are checked metadata: convert kV to V before construction.
The intent extends each curve flat beyond its outer knots. It contains no solver,
smoothing family, smoothing width, integer encoding or per-unit base.

`VoltVarWattEncoding` selects a representation independently for each curve.
Widths are in **physical volts**, including when the optimization variables are
per unit. Equal widths across families do not mean equal error: the example
chooses conservative **output-fraction** error budgets instead. This is a bound
on the curve approximation, not the composed controller or network solution.

`lower_positive_policy` builds the corresponding existing average-, worst-phase
or positive-sequence policy. It adds no model variables or constraints and picks
no optimizer. Construct it once per intent/encoding combination and reuse it
across snapshots. Existing `PiecewiseLinearLaw` constructors remain available;
the intent interface is an additional semantic frontend.

## What “compilation” means here

There are three distinct operations:

1. **Semantic lowering:** intent plus encoding becomes a control policy or a
   bounded curve graph. The full policy selects voltage measurements, evaluates
   curves, resolves conflicts, forms P/Q requests and applies sequence/current
   allocation and limiting rules.
2. **Model lowering:** stamping adds physical device equations and JuMP
   variables, expressions and constraints. JuMP and MathOptInterface (MOI)
   provide the nonlinear derivative interface or constraint bridges required by
   the attached optimizer. A bridge translates a supported mathematical
   representation; it is not permission to replace a graph by its convex hull.
3. **Numerical solution:** an external optimizer solves those equations. Julia's
   first-use machine-code compilation is an additional runtime cost, separate
   from both model construction and solver iterations.

The implemented routes are:

| Route | What is constructed | Backend and scope |
|:--|:--|:--|
| `evaluate_exact(controller, measurement, request, ratings)` | Numeric replay of the canonical curves and established controller semantics | No optimization; a same-state reference, not a network equilibrium solve |
| `evaluate_smooth(...)` | Numeric replay using the policy's chosen approximations | No optimization; intended to match the stamped smooth controller |
| Full smooth IBR policy | `SequenceController` → `ControlledDevice` → staged device/controller equations → KCL → nonlinear program | Ipopt, or another compatible NLP optimizer such as MadNLP; convergence is case-dependent |
| Smooth curve port | An output variable and a smooth equality, with explicit voltage bounds | Ipopt or MadNLP in this tutorial |
| `ExactPWLGraph` curve port | PiecewiseLinearOpt segment-selection graph with discrete structure | HiGHS for the otherwise affine example; nonlinear AC coupling would require a suitable MINLP approach |
| `PWLConvexHull` curve port | Convex combinations of all bounded graph vertices, without adjacency | LP in this example; HiGHS solves the relaxation |
| `ComplementarityGraph` curve port | Nonnegative hinge pairs with a difference equation and `MOI.Complements` | MathOptComplements bridges → CCOpt's MOI/NLPModelsJuMP interface → MPCC relaxation/homotopy with MadNLP |

**Full IBR policy lowering currently accepts smooth curve representations only.**
Passing a MIP, hull or complementarity encoding to `lower_positive_policy` raises
`UnsupportedFormulation`. Those representations are available through
`formulate_control_curve!` at an already sensed voltage. That operation does not
encode phase extrema, conflict resolution, P/Q conversion, hardware capability
or the network. Extending their use to a full AC controller requires explicit
formulations for that surrounding logic and a compatible solver.

### How the smooth equations and derivatives reach the solver

A continuous PWL curve with flat tails can be written

```math
f(V)=y_1+\sum_j c_j\max(V-k_j,0),\qquad
f_\epsilon(V)=y_1+\sum_j c_j h_\epsilon(V-k_j).
```

The slope changes `c_j` and knot locations `k_j` are prepared in the immutable
`PWLFunction`. Zero slope changes need no hinge. The representation determines
`h`: BMOPFTools telescoping softplus, a compact C2 quartic patch to the hinge,
or the algebraic square-root hinge. This is not interpolation by an unconstrained
cubic spline, and “sigmoid” alone does not specify the curve construction. The
logistic sigmoid is the derivative of the softplus hinge.

For **plain JuMP** models, PowerOptLab registers a scalar operator for the whole
smooth curve and supplies its analytic first and second derivatives. It caches
operator registration by curve, family, width and coordinate scales. The softplus
value callback uses BMOPFTools' numeric helper; its input vectors are captured
once. The explicit derivatives avoid differentiating through that helper's
Float64 conversion. For the other families, values and derivatives use the
shared hinge interface.

For **staged BMOPFTools** construction, softplus expressions instead use
`BMOPFTools.opf_piecewise_linear_expression`. Its default
`softplus=:user_defined` registers stable `log1pexp`/logistic hinge operators and
reuses them by normalized width. The telescope remains visible as a sum of
hinges. Local C2, algebraic and custom PWL smoothings use PowerOptLab's scalar
whole-curve operator also in the staged path. These routes represent the same
selected smooth law, although their expression graphs and floating-point
operation orders differ.

JuMP combines the supplied scalar derivatives with those of the surrounding
expressions to provide the NLP's derivative information. We do not hide an
entire network inside one large user-defined operator. That would obscure its
structure and generally make derivative work less efficient. See JuMP's
[nonlinear modeling and operator documentation](https://jump.dev/JuMP.jl/stable/manual/nonlinear/).
The curve family's smoothness also does not change every other controller
operation: see the [IBR primitive integration map](controllers.md#What-the-advanced-IBR-refactor-shares).

A specialized BMOPFTools option, `build_opf_model(...; softplus=:builtin)`, emits
native `log1p(exp(...))` expressions for wrappers that reject user-defined
operators, notably the pinned DiffOpt nonlinear route. It has a narrower
numerically safe range than stable `log1pexp`. It is **not** a general speed or
stability setting, and it does not convert PowerOptLab's other custom operators
into native expressions. Use a custom staged `FormulationCase` when selecting
this upstream option; `controlled_inverter_case` does not expose that keyword.
Neither this option nor a successful solve establishes sensitivity validity at
a degenerate point.

## Worked example: expose all the encodings

Consider the intent above at an already measured voltage of 245 V. Its canonical
volt-watt fraction is 0.5. Use the explicit physical domain `[210,280]` V. The
following builds a normalized smooth curve port:

```@example control_compilation
model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model,v,start=245/230)
@constraint(model,v == 245/230)
h = formulate_control_curve!(model,intent,:volt_watt,v,encoding;
    domain=(210.,280.),input_scale=230.,output_scale=1.)
@objective(model,Max,h.output)
optimize!(model)
@assert is_solved_and_feasible(model)
audit_pwl(h)
```

`input_scale*v` is in volts and `output_scale*h.output` is the fraction. Maximizing
at a fixed input is deliberate: the exact graph fixes the output, whereas the
hull permits an interval. This isolates the meaning of the encoding from network
feedback and solver-dependent selection among multiple equilibria.

The reusable script runs all three smoothings and the hull with Ipopt. It uses
the same semantic intent and a fresh model per method:

```@example control_compilation
include(joinpath(dirname(pathof(PowerOptLab)),"..","scripts","formulations","control_lowering.jl"))
rows = ControlLoweringExample.run_core()
@assert all(r -> r["strict_solver_success"],rows)
[(r["method"],r["metrics"].fraction,only(r["observations"]).exact_graph_error) for r in rows]
```

The smooth runs should agree with their selected surrogate; they need not agree
exactly with the PWL graph. Here symmetry at 245 V makes some surrogate errors
particularly small, so this single point is not an error-bound experiment.
The independently imposed fraction budget is `1e-3`.

For the hull, maximizing returns **0.875**, not 0.5: the convex combination of
`(240 V,1)` and `(280 V,0)` gives 0.875 at 245 V. It is a valid relaxed point with
a canonical graph error of 0.375. At 10 kW available power, these fractions would
request 8.75 kW and 5 kW respectively, before capability limiting. The hull's
answer depends on the declared domain. It must not be reported as exact
volt-watt behavior.

Run the external paths in the isolated optional environment:

```sh
julia --project=. scripts/instantiate_pinned.jl
julia scripts/formulations/setup.jl /tmp/pol-lowering-env  # new, empty directory
POL_FORMULATION_RESULTS=/tmp/control-lowering.toml \
  julia --project=/tmp/pol-lowering-env scripts/formulations/control_lowering_optional.jl
```

That script explicitly configures:

```julia
using HiGHS, PiecewiseLinearOpt, MadNLP
using MathOptComplements, NLPModelsJuMP, CCOpt

smooth = FormulationMethod("softplus / MadNLP",SoftplusFormulation(.01),
    MadNLP.Optimizer;configure! = set_silent,options=(tol=1e-9,))
exact = FormulationMethod("exact / HiGHS",ExactPWLGraph(method=:Logarithmic),
    HiGHS.Optimizer;configure! = set_silent)
hull = FormulationMethod("hull / HiGHS",PWLConvexHull(),
    HiGHS.Optimizer;configure! = set_silent)
mpcc = FormulationMethod("MPCC / CCOpt",ComplementarityGraph(scale=1.),
    CCOpt.Optimizer;
    configure! = m -> (MathOptComplements.Bridges.add_all_bridges(m);set_silent(m)),
    options=(tol=1e-9,))
```

`configure!` runs **before** optimizer attachment so bridges can be installed.
The pinned CCOpt interface constructs its MPCC representation and calls CCOpt's
homotopy solver when complementarity is present. PowerOptLab calls `optimize!`
once; it adds no local continuation, retry or re-solve loop. Inspect CCOpt's
[interface and options](https://github.com/madsuite-org/CCOpt.jl) for its own
algorithm configuration. An approximate MPCC result is not automatically an
exact graph point or a stationarity certificate. The optional tutorial retains
raw statuses, candidates and physical residuals separately.

## Configuration: which choice belongs where?

| Setting | Meaning and location |
|:--|:--|
| Curves, `sensing` | Physical law and measurement policy in `VoltVarWattIntent`; sensing is `:worst_phase`, `:average_voltage`, or `:positive_sequence` |
| `volt_watt_basis` | `:available` or `:rated`, applied by the full policy |
| `conflict_policy`, `conflict_epsilon` | Worst-phase firmware conflict rule and its dimensionless blend width, including exact replay; inactive for the other sensing policies |
| `worst_phase_watt_guard` | Guard choice for positive-sequence sensing; inactive for the other policies |
| Per-curve representation and width | `VoltVarWattEncoding`; each present curve needs an explicit representation; smooth widths are physical volts |
| `extrema_epsilon` | Physical-voltage width for the existing algebraic phase extrema/guard in smooth full-policy lowering; default 0.05 V; not used by an already sensed curve port |
| `domain` | Required finite voltage domain for a curve port; a physical modeling restriction, distinct from the intent's flat-tail knots |
| `input_scale`, `output_scale` | Positive physical units per working-coordinate unit; used by curve-port/model lowering |
| `ComplementarityGraph(scale=...)` | Positive physical voltage per hinge-variable unit; independent of `input_scale`. Default `nothing` uses the input scale |
| `ExactPWLGraph(method=...)` | PiecewiseLinearOpt's graph encoding method; default `:Logarithmic`; backend support still matters |
| `reuse` | Share the identical curve block at one variable in one model; default `true` |
| `FormulationMethod` | Optimizer constructor, solver-specific `options`, bridge/configuration callback and metadata |
| Experiment configuration | Caller-defined starts, scenario data and objectives. The inverter adapter accepts `s_base` and `selection_objective=:loss` or `:zero`; custom cases expose additional settings |
| Controller and plant settings | Sequence control, voltage floors, current allocation, limiter priorities, ratings and losses remain in their existing typed controller/plant objects |

Changing scales changes the numerical coordinates, not the physical smoothing
width. For a smooth equality `y=f(sᵢ x)/sₒ`, the derivative multipliers are
`sᵢ/sₒ` and `sᵢ²/sₒ`. Very small smoothing widths can still produce large curvature.
Choose width from a physical approximation budget, then choose sensible variable
and residual scales, solver tolerances and starts. A tolerance in normalized
solver coordinates is not directly a voltage, watt or graph-error guarantee.
See [physical error and conditioning](error_budgets.md).

## Lower efficiently, and measure the right costs

There are two different kinds of reuse:

* **Operator reuse** avoids registering the same scalar function and derivative
  callbacks repeatedly. Different voltage variables can share an operator but
  still need distinct output equations.
* **Curve-block reuse** avoids duplicating those output variables and equations
  when the same curve is requested at the **same variable**, with the same
  domain, encoding, scales and lowering target. `formulate_control_curve!`
  returns the existing handle. Set `reuse=false` for a distinct block.

```@example control_compilation
shared = ControlLoweringExample.reuse_demo(ports=10)
fresh = ControlLoweringExample.reuse_demo(ports=10,reuse=false)
@assert shared.shared_curve_blocks == 1
@assert shared.pwl_operators == fresh.pwl_operators == 1
@assert shared.variables == 2 && fresh.variables == 11
(shared=shared,fresh=fresh)
```

The shared model has one input and one curve output; the fresh model has one
input and ten outputs. These are structural counts, not a solve-speed benchmark.
The default cache applies only to `VariableRef` inputs. General/mutable
expressions get fresh blocks. Caches belong to their model; do not transplant
handles to another model or manually delete cached rows and expect automatic
repair. Build a fresh model after such structural changes.

JuMP does not automatically eliminate common nonlinear subexpressions. Reusing
a Julia variable holding `@expression` does not by itself guarantee that its
nonlinear computation is evaluated only once. For a costly repeated sensed
quantity, introduce one physically scaled auxiliary variable and its defining
row, then feed that variable to the curve ports. The full controller already
uses normalized auxiliary definitions to limit expression growth. Every added
auxiliary also enlarges the solver's linear system, so use sharing where there
is actual repetition. See JuMP's
[common-subexpression discussion](https://jump.dev/JuMP.jl/stable/manual/nonlinear/#Common-subexpressions)
and [performance guidance](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/performance_tips/).

Prepare immutable intent and encodings outside scenario loops; lower a full
policy once per configuration. Within a model, reuse ports and operator
registrations. The prepared curve stores slopes and nonzero hinge coefficients,
so PowerOptLab's derivative callbacks do not rebuild them at each evaluation.
Softplus still delegates numeric evaluation to BMOPFTools, whose helper retains
its own preparation/validation cost. Legacy `PiecewiseLinearLaw` stores mutable
vectors, so its numeric path is not blindly cached. This is useful bounded reuse,
not a claim of allocation-free evaluation or an optimizing compiler for arbitrary
control expression graphs.

For timing, warm Julia on a disposable representative model first, then measure
fresh-model construction separately from `optimize!`. Report model dimensions,
allocations, solver iterations, physical errors and environment versions along
with elapsed times. `run_formulation_experiment` records construction and solve
times, but first-use compilation may be included. `lowering_statistics` reports
PowerOptLab cache sizes and model dimensions; it does not count upstream
BMOPFTools registrations or time MOI bridging separately. Changes to curve data,
widths or representation require new lowering; this frontend does not implement
parameter updates across existing models. The experiment runner deliberately
builds fresh models so warm starts and mutable solver state do not silently
cross configurations.

Continue with [the full physical IBR example](controllers.md), or use the
[analytic feeder comparison](index.md) to study feedback and multiple equilibria.
For the mathematical distinction between graph formulations and relaxations,
see Huchette and Vielma's
[piecewise-linear formulation study](https://arxiv.org/abs/1708.00050).
