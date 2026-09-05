# Physical error, composition and conditioning

A useful comparison distinguishes four questions: how much the function changed,
how accurately the encoded equations were solved, how far the equilibrium moved,
and whether the physical operating constraints hold. A bound answering the first
question does not automatically answer the other three.

## From output budgets to input widths

For a flat-extended PWL curve, `f(x)=f₀+Σ cⱼ(x-kⱼ)₊`. Suppose each smoothed hinge
has signed error in `[a,b]`. Then a real-arithmetic interval for the complete
curve error is

```math
\sum_j\min(c_j a,c_j b)\leq \widetilde f(x)-f(x)
\leq\sum_j\max(c_j a,c_j b).
```

This follows by multiplying each interval by its signed slope change and adding.
It ignores correlation between hinge errors, so it can overestimate the attainable
error. `formulation_contract` performs this accounting, including endpoint slope
changes. For the built-in homogeneous hinge families, `smoothing_for_error` inverts
the bound to select a width. It controls absolute error in **curve output units**.
A volt-watt fraction budget of `1e-3` is not `1e-3 W`: multiplication by an available
power of 12 kW can contribute up to 12 W error before other controller operations.

```@example budgets
using PowerOptLab, LinearAlgebra
curve = PWLFunction([220.,240.,250.,270.],[1.,1.,.2,.2];
    input_unit=:V,output_unit=:pu)
r = smoothing_for_error(curve,LocalC2Formulation,1e-3)
c = formulation_contract(curve,r)
power_error = affine_error_bound([12_000.],[(c.error_lower,c.error_upper)])
@assert isapprox(power_error.upper,12.)
(r.width,power_error)
```

`affine_error_bound` applies to fixed coefficients and compatible output units.
If the coefficient is itself uncertain, state-dependent, or another surrogate,
its error and range require additional analysis. Nonlinear composition needs
bounds such as a Lipschitz constant over a declared domain; no automatic global
composition certificate is supplied. Clipping a surrogate afterward changes its
regularity and should not be assumed harmless for an NLP solver.

## Norms: direction follows the purpose

Let `q=‖z‖₂` and `ε>0`, in the same physical units as the components. Define

```math
q_-(z)=\sqrt{z^\top z+\epsilon^2}-\epsilon,\qquad
q_+(z)=\sqrt{z^\top z+\epsilon^2}.
```

Both are infinitely differentiable. Directly from `q ≤ sqrt(q²+ε²) ≤ q+ε`,
`-ε ≤ q₋-q ≤ 0` and `0 ≤ q₊-q ≤ ε`. Their Hessian spectral norm is at most
`1/ε`; reducing bias therefore increases maximum curvature. The lower version
preserves zero and is often useful inside a controller. The upper version bounds
a magnitude conservatively for an upper-limit inequality, provided the component
vector itself is the physical quantity of interest.

| Purpose | Appropriate construction | Consequence |
|:--|:--|:--|
| Smooth magnitude signal with zero output at zero | `MagnitudeApproximation(ε)` | Underestimates by at most ε |
| Conservative approximate upper magnitude constraint | `direction=:upper` | May exclude exact feasible states |
| Exact circular capability `‖i‖≤Imax`, `Imax≥0` | `sum(i_k^2)≤Imax^2` | No root or magnitude approximation needed |
| Exact root away from zero | `positive_root_expression(...; lower_bound=L)` | Adds the explicit domain restriction `radicand≥L>0` |

Enforcing `q₋≤Imax` alone can admit up to ε magnitude excess; smoothing is not a
substitute for the exact physical cap. Conversely `q₊≤Imax` can be infeasible even
at zero if ε exceeds Imax. `positive_root_expression` does not cure a root whose
true domain includes zero; it deliberately changes the allowed domain.

```@example budgets
lower = MagnitudeApproximation(.01;unit=:A)
upper = MagnitudeApproximation(.01;direction=:upper,unit=:A)
@assert magnitude_value([0.,0.],lower)==0.
@assert magnitude_value([3.,4.],lower) <= 5. <= magnitude_value([3.,4.],upper)
magnitude_contract(lower)
```

`magnitude_expression` accepts a JuMP model or staged context. Components use
working coordinates; `component_scale` converts them to physical quantities and
`output_scale` converts the returned magnitude to the desired model coordinate.
The staged lower construction reuses BMOPFTools' shifted-root helper. Bounds are
real-arithmetic statements; cancellation and solver residuals still need checking.

## From a function perturbation to an equilibrium shift

For `F(x)=0`, replacing a component by a surrogate adds a residual perturbation
`e(x)`. At an isolated regular equilibrium, a first-order estimate is

```math
\Delta x\approx-J_F(x)^{-1}e(x).
```

For the scalar resistor example, `F(V)=V-Vs-R*f(V)` and `J=1-R*f′(V)`.
If the controller increases current by `ΔI`, its contribution to the equation
perturbation is `-R*ΔI`. On a segment of slope `-2 A/V` with `R=1 Ω`, a `.03 A`
perturbation therefore predicts a `.01 V` increase:

```@example budgets
estimate = local_error_response([3.;;],[-.03])
@assert isapprox(only(estimate.delta),.01)
@assert local_error_response([0.;;],[-.03]).delta === nothing
estimate
```

This diagnostic uses the supplied dense, square Jacobian and flags poor numerical
conditioning. A scalar nonzero Jacobian has condition number one even when its
absolute inverse is large: inspect the response magnitude and scaling as well.
Near a kink there may be no single classical derivative; near a fold the inverse
may cease to exist. A local estimate does not prove an exact-law equilibrium
exists nearby, is unique, or satisfies limits. Whole-network evaluation should
report exact physical residuals and use an independently specified equilibrium or
verification method when its scientific question requires one.

## Scientific context and attribution

The smooth-hinge approach belongs to established smoothing methods; see
[Chen and Mangasarian (1996), *A class of smoothing functions for nonlinear and
mixed complementarity problems*](https://doi.org/10.1007/BF00249052).
The compact quartic patch and the bounds above are derived explicitly in this
documentation and [the formulation guide](index.md); no claim of a novel algorithm
or a universal reliability advantage is made.

Exact MIP PWL graphs are implemented by the external library described by
[Huchette and Vielma, *Nonconvex piecewise linear functions: Advanced formulations
and simple modeling tools*](https://arxiv.org/abs/1708.00050). Exact segment semantics
do not make a nonlinear AC model globally solvable by a linear MIP solver.

For MPCCs, constraint qualification and stationarity require specialized care.
[Nurkanović, Pozharskiy and Diehl](https://arxiv.org/abs/2312.11022) review these
concepts and show why benchmark success does not automatically generalize across
problem classes. [CCOpt](https://github.com/madsuite-org/CCOpt.jl) supplies the
external algorithm; this layer supplies the encoding and physical diagnostics.
Primitive hulls give useful outer relaxations, but a feasible hull point alone is
not evidence of a realizable controller or a safe operating envelope.
