# Function formulations: one curve, different mathematical models

This experimental layer separates a canonical bounded continuous PWL function
from its numerical representation. It is a staging area for reusable BMOPFTools
function primitives. The existing inverter controllers retain their established
defaults and can now select a smoothing family independently for each curve.
Start with [a configurable experiment](experiments.md), then use this page for
the mathematical models and the original illustrative backend comparison.

| Representation | Meaning | Implementation |
|:--|:--|:--|
| `SoftplusFormulation(width)` | Smooth surrogate graph, C∞ | BMOPFTools telescoping softplus |
| `LocalC2Formulation(width)` | Smooth surrogate graph, C2 | Compact quartic hinge patches |
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
stateful hysteresis are outside this first API. A continuous function cannot
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

For this particular nonincreasing scalar controller and nonnegative resistance,
the uniform current-law error bounds equilibrium current error: if the current
increases, voltage increases, and the exact controller command cannot increase.
This negative-feedback argument does not generalize automatically to AC systems.

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
pins experimental packages, including DiffOpt 0.6.2: the baseline local environment's
0.6.1 caps MathOptInterface below the version required by MathOptComplements 0.1.1.
Dependency changes are isolated. Set `POL_FORMULATION_RESULTS` to a TOML output
path to record package versions and every external-backend run.

Run `julia --project=. scripts/formulations/comparison.jl` for the three core
representations. Both core and optional tests compare plateaus, kinks, and the
ramp over independent starts and physical coordinate scales. Each run reports
solver status separately from electrical, domain, current-limit, exact-controller,
and complementarity residuals. MPCC stationarity is explicitly recorded as not
independently assessed; neither a small primal residual nor an MOI success flag
is promoted to a stationarity or network safety certificate.

## Measured external-backend limitations

On the initial macOS / Julia 1.12.6 matrix with the pinned optional packages:

| Representation | Strict solver success | Canonical equations within 1e-6 physical units |
|:--|--:|--:|
| Exact graph / HiGHS | 20/20 | 20/20 |
| Complementarity / CCOpt 0.1.0 | 0/20 | 8/20 |
| Hull / HiGHS | 1/1 | 0/1 (expected relaxation witness) |

All CCOpt runs return `ALMOST_LOCALLY_SOLVED`; none is promoted to strict success.
The largest exact-controller discrepancy is about 0.0376 A, exceeding the smooth
comparison's 0.01 A budget. The complementarity encoding is mathematically exact;
these numerical candidates need not be. Product-based complementarity residuals
scale quadratically with physical coordinate scales and can allow much larger
individual hinge errors near biactive points. Changing variable coordinates does
not change the exact graph, but it can change finite-tolerance solver accuracy.

The optional CI job checks exact-graph accuracy and CCOpt integration, and uploads
the complete characterization matrix. It does **not** assert convergence of every
CCOpt case. Its pass status must not be read as an MPCC reliability claim. Solver
configuration, physical complementarity normalization, and stationarity validation
remain explicit research work. The raw snapshot accompanies the comparison
scripts as `scripts/formulations/reference_results.toml`; rerunning with a results
path produces fresh evidence, including unsuccessful statuses.

## Scope and scientific references

The layer covers continuous scalar PWL representations, purpose-specific norms
and positive-domain roots, affine error intervals and local sensitivity estimates.
Configurable cases and opt-in controller families support further experiments.
A general nonlinear error-propagation calculus, discontinuous/hysteretic controller
semantics and automatic graph reformulation of whole AC controllers remain outside
this delivery. See [the configurable runner](experiments.md) and
[controller integration](controllers.md) for the expanded toolkit.

- Chen and Mangasarian (1996), [smoothing functions for complementarity](https://doi.org/10.1007/BF00249052): construction via integrated probability densities.
- Chen (2012), [nonsmooth nonconvex smoothing](https://doi.org/10.1007/s10107-012-0569-0): gradient consistency and limiting stationarity require assumptions beyond function approximation.
- Huchette and Vielma, [PWL formulations](https://arxiv.org/abs/1708.00050): exact graph encodings and their relaxation strength; implemented by PiecewiseLinearOpt.
- Nurkanović, Pozharskiy and Diehl (2024), [MPCC methods and control benchmarks](https://arxiv.org/abs/2312.11022): stationarity distinctions and limits of transferring benchmark performance.
- Pozharskiy et al. (2026 preprint), [CCOpt](https://arxiv.org/abs/2604.18726): specialized relaxation/penalty and active-set methods built on MadNLP.

See the [implementation plan](plan.md) for the staging and migration boundary.
