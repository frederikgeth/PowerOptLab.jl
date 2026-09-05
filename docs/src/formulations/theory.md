# Function formulations: one curve, different mathematical models

A canonical bounded continuous piecewise-linear (PWL) function can have several
numerical representations. This page derives their mathematical meaning and
compares them on an analytically solvable circuit. For choosing an interface or
a worked example, start with [the overview](index.md).

| Representation | Meaning | Implementation |
|:--|:--|:--|
| `SoftplusFormulation(width)` | Smooth surrogate graph, C∞ | BMOPFTools telescoping softplus |
| `LocalC2Formulation(width)` | Smooth surrogate graph, C2 | Compact quartic hinge patches |
| `AlgebraicFormulation(width)` | Smooth surrogate graph, C∞ | Square-root hinge; existing IBR selector family |
| `ComplementarityGraph()` | Exact graph, subject to numerical MPCC tolerances | Native JuMP/MOI complementarity; external solver |
| `ExactPWLGraph()` | Exact bounded segment graph, subject to numerical MIP tolerances | Optional PiecewiseLinearOpt extension |
| `PWLConvexHull()` | Continuous convex hull of the bounded graph | Vertex convex combinations |

A solver can solve each representation successfully while answering different
questions. In particular, a hull point need not be a realizable controller state.
Primitive hulls composed with network constraints generally yield an outer
relaxation of the full feasible set. They do not certify a safe DOE.

## Canonical data and physical units

`PWLFunction` copies finite, strictly increasing breakpoints and corresponding
values into immutable tuples. The first and last breakpoints define its closed
physical domain. Public numeric evaluation rejects points outside this domain.
Smooth construction uses the curve's flat extension, including its endpoint kinks.
Unit symbols are labels; the caller supplies coordinate conversions explicitly.

```@example pwl
using PowerOptLab, JuMP, Ipopt
curve = PWLFunction([220.,240.,250.,270.], [20.,20.,0.,0.];
    input_unit=:V, output_unit=:A)
@assert primitive_value(curve,245.) == 10.
primitive_value(curve,245.)
```

A hinge, min/max with a constant, clamp or deadband can be specified by its knots
on a finite domain. A binary max can be composed as `b + max(0,a-b)`, provided the
physical domain of `a-b` is known. Arbitrary multivariate surfaces, jumps, and
stateful hysteresis are not supported by this scalar API. A continuous function cannot
approximate a discontinuous jump with arbitrarily small uniform error across the
jump. A switching model must specify its boundary and memory semantics.

## Smoothing guarantees and derivatives

For the flat-extended curve, write

```math
f(x)=f(x_1)+\sum_j c_j(x-k_j)_+,
```

where the coefficients are slope changes, including both endpoints. Zero slope
changes are removed. BMOPFTools evaluates the softplus telescoping sum, with
stable `log1pexp`. The compact hinge uses

```math
h_\delta(x)=\begin{cases}
0,&x\leq-\delta,\\
\delta(3/16+t/2+3t^2/8-t^4/16),\quad t=x/\delta,&|x|<\delta,\\
x,&x\geq\delta.
\end{cases}
```

Its value, first derivative, and second derivative match at both joins. It is
monotone and convex, with maximum positive hinge error `3δ/16` and maximum
curvature `3/(4δ)`. It equals the original hinge outside the patch. Patch bands
may overlap: the error bounds still hold, but more of the curve is then modified.

For softplus the maximum positive hinge error is `ε*log(2)` and the maximum
hinge curvature is `1/(4ε)`. For either representation, let `B` be its maximum
hinge error. Then

```math
B\sum_{c_j<0}c_j\ \leq\ f_{\rm smooth}(x)-f(x)\ \leq\ B\sum_{c_j>0}c_j.
```

Mixed coefficient signs preclude inferring a one-sided bound for the complete
curve from the one-sided hinge bound. These are real-arithmetic function bounds;
floating-point evaluation, optimization residuals, and model errors are additional.
They do not certify equilibrium error in a general AC network.

```@example pwl
budget = 0.01 # A
soft = smoothing_for_error(curve,SoftplusFormulation,budget)
local_c2 = smoothing_for_error(curve,LocalC2Formulation,budget)
for rep in (soft,local_c2)
    contract = formulation_contract(curve,rep)
    @assert isapprox(contract.error_upper,budget)
    @assert isapprox(contract.error_lower,-budget)
end
primitive_derivatives(curve,245.,local_c2)
```

The local formulation preserves sufficiently distant plateaus exactly. Shrinking
either width increases curvature; compare formulations at equal physical error
budgets. Registered scalar JuMP operators supply analytic gradients and Hessians,
including patch joins. Numeric local-C2 evaluation supports ForwardDiff; numeric
softplus delegates to BMOPFTools' Float64 oracle. On a BMOPFTools staged context,
`formulate_pwl!` reuses its public expression builder and registration cache.

## An executable electrical reference

Consider a real-current inverter exporting through a 1 Ω resistor to a source at
230 V. Its static control law is the curve above. The equations are

```math
V=230+I,\qquad I=f(V),\qquad P=VI.
```

This is a scalar resistive electrical model, not an approximation claimed to
represent a general unbalanced AC feeder. On the sloped segment, `I=500-2V`, so
`V=243⅓ V` and `I=13⅓ A`. Segment enumeration verifies that this is the only
canonical equilibrium. The linear objective maximizes export current; uniqueness
makes its choice irrelevant to the exact controller equilibrium.

```@example pwl
function solve_reference(rep)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model,"tol",1e-9)
    @variable(model,v,start=245/230)
    h = formulate_pwl!(model,curve,v,rep;
        input_scale=230.,output_scale=20.)
    # Physical voltage is 230*v and physical current is 20*h.output.
    @constraint(model,v == 1 + (20/230)*h.output)
    @objective(model,Max,h.output)
    optimize!(model)
    @assert termination_status(model) == MOI.LOCALLY_SOLVED
    return audit_pwl(h)
end
for rep in (soft,local_c2)
    a = solve_reference(rep)
    @assert abs(a.input-230-a.output) < 1e-6
    @assert abs(a.output-40/3) <= budget+1e-6
    @assert abs(a.surrogate_equation_error) < 1e-6
end
```

For a non-increasing continuous controller `f` and `R≥0`, define
`F(I)=I-f(Vs+RI)`. For any `I₂>I₁` with voltages in the declared domain,

```math
F(I_2)-F(I_1)=(I_2-I_1)+f(V_s+RI_1)-f(V_s+RI_2)\geq I_2-I_1.
```

Thus `F` is strongly monotone with modulus at least one, including at PWL kinks.
An equilibrium, **if it exists in the domain**, is unique. If the surrogate has
uniform law error at most `e`, and both an exact equilibrium `I*` and a surrogate
equilibrium `Iδ` exist in the domain, then

```math
|I_\delta-I^*|\leq |F(I_\delta)-F(I^*)|=|F(I_\delta)|\leq e.
```

There is no resistance-dependent amplification factor in this scalar model.
This sensitivity result does not imply convergence of the iteration
`Iₖ₊₁=f(Vs+RIₖ)`: its slope magnitude can exceed one.
For an inexact surrogate solve with current-equation residual bounded by `η`
and an exact voltage relation, the last bound becomes `e+η`. A voltage-equation
residual requires additional accounting, for example using a Lipschitz bound on
`f`. In a general AC network, the coupled residual map need not be strongly
monotone; smoothness alone does not restore this property.

The hull admits convex combinations of graph vertices without requiring adjacent
segments. Its maximum-current state instead has `I=16 A, V=246 V`. The canonical
law commands only `8 A` at that voltage:

```@example pwl
relaxed = solve_reference(PWLConvexHull())
@assert relaxed.semantics == :outer_relaxation
@assert isapprox(relaxed.output,16.;atol=1e-6)
@assert isapprox(relaxed.exact_graph_error,8.;atol=1e-6)
relaxed
```

This witness is electrically consistent with the resistor equation and within
current limits, yet violates the controller. It demonstrates why a physical
binding check must include the prescribed control law.

## Multiple equilibria and initialization

Strong monotonicity fails for an increasing segment with sufficient feedback.
This synthetic tent law has two canonical equilibria; it is an analytic test
function, not a prescribed inverter characteristic:

```@example multiroot
using PowerOptLab,JuMP,Ipopt
tent=PWLFunction([0.,1.,2.],[0.,2.,0.])
roots=resistive_equilibria(tent,0.,1.)
@assert isapprox([p.current for p in roots.points],[0.,4/3])
roots
```

The zero-current root is at an endpoint kink. Compact smoothing makes the command
positive there and removes this root; starting near zero does not guarantee a
nearby smooth equilibrium. To illustrate initialization without conflating it
with a changed root set, shift the tent downward by 0.2 A. Its roots are now
strictly inside affine segments, at `I=0.2 A` and `I=19/15 A`; local C2 leaves
both neighborhoods unchanged.

```@example multiroot
shifted=PWLFunction([0.,1.,2.],[-.2,1.8,-.2];input_unit=:V,output_unit=:A)
case=resistive_control_case(shifted;source_voltage=0.,resistance=1.)
method=FormulationMethod("C2",LocalC2Formulation(.01),Ipopt.Optimizer;
    configure! = set_silent,options=(tol=1e-10,bound_relax_factor=0.))
rows=run_formulation_experiment([case],[method];
    configurations=[(start_input=s,objective=:zero) for s in (.15,1.5)],on_error=:throw)
@assert all(r->r["strict_solver_success"],rows)
@assert isapprox([r["metrics"].current_A for r in rows],[.2,19/15];atol=1e-7)
[(r["configuration"].start_input,r["metrics"].current_A,
  only(r["observations"]).exact_graph_error) for r in rows]
```

The two runs solve a feasibility problem with the same equations and different
starts. Their successful local outcomes are not uniqueness or global-selection
certificates. MPCC and mixed-integer encodings represent different mechanisms for
choosing graph segments; neither removes the need to specify the physical
selection rule or optimization objective. Compare root sets and returned
candidates as well as termination flags.

## Exact graphs and external complementarity solvers

Load `PiecewiseLinearOpt` to activate `ExactPWLGraph`; its default logarithmic
formulation can be solved with HiGHS. We delegate segment encoding to that package.
The continuous hull itself needs only affine vertex equations and JuMP.

For each interior hinge, complementarity introduces `p,n >= 0`,
`p-n = x-k`, and `p ⟂ n`. The hinge value is `p`. Endpoint hinges are affine on
the bounded domain and are eliminated. Classical constraint qualifications fail
for the standard NLP reformulation of complementarity; a successful ordinary NLP
solve is not a general MPCC optimality certificate.

The optional runner loads `MathOptComplements`, `NLPModelsJuMP`, and `CCOpt`, adds
its bridges, and selects `CCOpt.Optimizer`. CCOpt owns the relaxation/homotopy
algorithm. PowerOptLab supplies no continuation schedule or solve retries.

With `ComplementarityGraph(scale=s)`, stored hinge variables are normalized:
`p-n=(V-k)/s` and the physical hinge is `s*p`. This scale is independent of the
model's input coordinate base. `audit_pwl` converts minimum and product residuals
back to physical input units and their square. The default `scale=nothing`
preserves input-base normalization; changing the scale requires fresh evidence,
not reinterpretation of a previous solver tolerance as a physical tolerance.

```sh
julia --project=. scripts/instantiate_pinned.jl
julia scripts/formulations/setup.jl /tmp/pol-pwl-env
julia --project=/tmp/pol-pwl-env scripts/formulations/optional_tests.jl
```

The destination must be empty. The setup preserves the BMOPFTools revision and
pins a compatible set of optional packages in an isolated environment. Set `POL_FORMULATION_RESULTS` to a TOML output
path to record package versions and every external-backend run.

Run `julia --project=. scripts/formulations/comparison.jl` for the three core
representations. Both core and optional tests compare plateaus, kinks, and the
ramp over independent starts and physical coordinate scales. Each run reports
solver status separately from electrical, domain, current-limit, exact-controller,
and complementarity residuals. MPCC stationarity is explicitly recorded as not
independently assessed; neither a small primal residual nor an MOI success flag
is promoted to a stationarity or network safety certificate.

## Interpreting external-backend evidence

An exact encoding does not guarantee an accurate numerical candidate. Near a
biactive complementarity point, a small product residual can coexist with larger
individual hinge errors. Product residuals scale quadratically with physical
input units. Report both physical minimum and product residuals, the canonical
controller discrepancy, and the raw termination status.

The optional scripts export every run, including unsuccessful candidates. They
check exact-graph accuracy and complementarity integration; a passing integration
check is not a reliability claim for every MPCC solve. Reproduce the matrix with
`POL_FORMULATION_RESULTS` set, and inspect its package versions and options before
comparing results. A historical reference dataset and reproduction instructions
are available in [the experiment scripts](https://github.com/frederikgeth/PowerOptLab.jl/tree/main/scripts/formulations).

## Scope and scientific references

The layer covers continuous scalar PWL representations, purpose-specific norms
and positive-domain roots, affine error intervals and local sensitivity estimates.
Configurable cases and opt-in controller families support further experiments.
A general nonlinear error-propagation calculus, discontinuous/hysteretic controller
semantics and automatic graph reformulation of whole AC controllers are not implemented. See [the configurable runner](experiments.md) and
[controller integration](controllers.md) for the expanded toolkit.

- Chen and Mangasarian (1996), [smoothing functions for complementarity](https://doi.org/10.1007/BF00249052): construction via integrated probability densities.
- Chen (2012), [nonsmooth nonconvex smoothing](https://doi.org/10.1007/s10107-012-0569-0): gradient consistency and limiting stationarity require assumptions beyond function approximation.
- Huchette and Vielma, [PWL formulations](https://arxiv.org/abs/1708.00050): exact graph encodings and their relaxation strength; implemented by PiecewiseLinearOpt.
- Nurkanović, Pozharskiy and Diehl (2024), [MPCC methods and control benchmarks](https://arxiv.org/abs/2312.11022): stationarity distinctions and limits of transferring benchmark performance.
- Pozharskiy et al. (2026 preprint), [CCOpt](https://arxiv.org/abs/2604.18726): specialized relaxation/penalty and active-set methods built on MadNLP.

See [references and related software](references.md) for the scientific context
and the responsibilities of each external library.
