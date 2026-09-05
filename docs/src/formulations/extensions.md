# Implement a smoothing family

The extension interface separates a mathematical hinge approximation from model
construction. Implement an immutable subtype of `AbstractPWLSmoothing` and three
methods: value, first/second derivative, and a contract. The existing PWL
composition, JuMP scalar operator and physical error accounting then use it.

As an example, consider the classical square-root hinge

```math
h_\delta(x)=\frac{x+\sqrt{x^2+\delta^2}}{2},\qquad\delta>0.
```

Its derivatives are `(1+x/sqrt(x²+δ²))/2` and
`δ²/(2(x²+δ²)^(3/2))`. Subtracting `max(0,x)` gives
`(sqrt(x²+δ²)-|x|)/2`, which lies in `[0,δ/2]`. The maximum second derivative is
`1/(2δ)`. Unlike the compact patch, this family changes the hinge everywhere.
This family is available as `AlgebraicFormulation`. The following small
reimplementation is pedagogical: it shows every extension method a researcher
would supply for a different family.

```@example extension
using PowerOptLab, JuMP, Ipopt
struct RootHinge <: AbstractPWLSmoothing
    width::Float64
    function RootHinge(width)
        isfinite(width) && width>0 || throw(ArgumentError("positive finite width required"))
        new(Float64(width))
    end
end
# Rationalize the negative branch to avoid subtracting similar magnitudes.
function PowerOptLab.hinge_value(x,r::RootHinge)
    d = r.width
    root = hypot(x,d)
    return x>=0 ? x/2+root/2 : (d/2)*(d/(root-x))
end
function PowerOptLab.hinge_derivatives(x,r::RootHinge)
    root = hypot(x,r.width)
    return ((1+x/root)/2, (r.width/root)^2/(2root))
end
PowerOptLab.hinge_contract(r::RootHinge) = (
    error_lower=0.,error_upper=r.width/2,regularity=:C_infinity,
    width=r.width,second_derivative_bound=1/(2r.width))
curve = PWLFunction([220.,240.,250.,270.],[20.,20.,0.,0.];
    input_unit=:V,output_unit=:A)
rep = RootHinge(.01)
@assert isapprox(hinge_value(0.,rep),.005)
@assert isapprox(hinge_derivatives(0.,rep)[2],50.)
formulation_contract(curve,rep)
```

These declarations are supplied by the author; the library does not prove them.
The hinge contract must provide all five fields shown above. `width` is in input
units, signed errors in hinge-output units, and `second_derivative_bound` bounds
absolute curvature. For general signed hinge bounds, the PWL contract propagates
both ends through each slope change. `regularity` describes actual global
smoothness, including joins. A custom family must supply at least two continuous
derivatives for the default Hessian-based NLP operator to be appropriate.

```@example extension
case = resistive_control_case(curve;source_voltage=230.,resistance=1.)
method = FormulationMethod("root hinge",rep,Ipopt.Optimizer;
    options=(tol=1e-9,),configure! = set_silent,
    metadata=(width_V=rep.width,optimizer="Ipopt"))
row = only(run_formulation_experiment([case],[method];on_error=:throw))
@assert row["strict_solver_success"]
@assert abs(only(row["observations"]).surrogate_equation_error)<1e-6
nothing # hide
```

If you want physical-budget selection, implement a specialized
`smoothing_for_error(curve, ::Type{YourFamily}, budget)` method or use a method
callback. The generic helper intentionally does not assume every new family's
error scales linearly with a parameter called width. You may also select different
families/widths for different curves via a named-tuple representation interpreted
by your own case builder.

Validate representative values, join behavior, independent derivative comparisons,
and model integration. Distinguish analytically proved properties from numerical
checks: sampling alone does not establish a uniform error bound or smoothness at
every join. Document the domain, units and any unverified contract assumptions.

## Graphs, jumps and state

`AbstractPWLSmoothing` is specifically a continuous scalar hinge interface. It
cannot represent a set-valued convex hull, an exact jump, or memory. A continuous
surrogate of a jump of height `J` has worst-case error at least `J/2` when both
one-sided limits are included; reducing the transition width does not remove that
uniform-error obstruction. Choose boundary values, one-sided semantics, and any
hysteresis state explicitly in a custom `FormulationCase`. Use native JuMP graph,
integer or complementarity constraints and an appropriate external solver.

A graph formulation can be added through specialized construction/contract
methods; it should not pretend to be a smooth hinge. The existing exact bounded
PWL extension delegates to PiecewiseLinearOpt instead of maintaining a second
implementation of its graph encodings. Reusing that implementation also keeps
encoding options aligned with the external package.
