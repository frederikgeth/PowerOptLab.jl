"""
    PWLFunction(breakpoints, values; input_unit=:unitless, output_unit=:unitless)

Immutable continuous scalar piecewise-linear function on the finite closed domain
`(first(breakpoints), last(breakpoints))`. Coordinates and values are physical;
unit symbols are metadata, not automatic unit conversions. Outside-domain public
numeric evaluation is rejected. Smoothing uses the flat extension of this curve.
"""
struct PWLFunction
    breakpoints::Tuple{Vararg{Float64}}
    values::Tuple{Vararg{Float64}}
    input_unit::Symbol
    output_unit::Symbol
    function PWLFunction(breakpoints, values; input_unit::Symbol=:unitless,
                         output_unit::Symbol=:unitless)
        xs, ys = Tuple(Float64.(breakpoints)), Tuple(Float64.(values))
        length(xs) == length(ys) >= 2 || throw(ArgumentError("Need at least two paired knots"))
        all(isfinite, xs) && all(isfinite, ys) || throw(ArgumentError("Curve data must be finite"))
        steps = diff(collect(xs))
        all(isfinite,steps) && all(>(0),steps) || throw(ArgumentError("Knot intervals must be finite and positive"))
        slopes = diff(collect(ys)) ./ diff(collect(xs))
        all(isfinite, slopes) && all(isfinite, diff([0.; slopes; 0.])) ||
            throw(ArgumentError("Curve slopes and slope changes must be finite"))
        new(xs, ys, input_unit, output_unit)
    end
end

"""Supertype of PWL graph, surrogate and relaxation representations."""
abstract type AbstractPWLFormulation end

"""
    AbstractPWLSmoothing <: AbstractPWLFormulation

Extension interface for a smooth hinge family. Implement `hinge_value`,
`hinge_derivatives`, and `hinge_contract`. Use immutable configuration objects.
The scalar JuMP operator and complete PWL error accounting are then supplied here.
"""
abstract type AbstractPWLSmoothing <: AbstractPWLFormulation end

function _positive_width(x)
    value = Float64(x)
    isfinite(value) && value > 0 ||
        throw(ArgumentError("Width must be finite and positive in physical input units"))
    return value
end

"""`SoftplusFormulation(width)` uses BMOPFTools' telescoping softplus with physical input width."""
struct SoftplusFormulation <: AbstractPWLSmoothing
    width::Float64
    SoftplusFormulation(width::Real) = new(_positive_width(width))
end

"""`LocalC2Formulation(width)` replaces each hinge in a band of physical half-width `width` by a C2 quartic."""
struct LocalC2Formulation <: AbstractPWLSmoothing
    width::Float64
    LocalC2Formulation(width::Real) = new(_positive_width(width))
end

"""`AlgebraicFormulation(width)` uses the C∞ hinge `(x + sqrt(x² + width²))/2`, with physical input width. This is the established algebraic family used by inverter-control selectors."""
struct AlgebraicFormulation <: AbstractPWLSmoothing
    width::Float64
    AlgebraicFormulation(width::Real) = new(_positive_width(width))
end

"""`ComplementarityGraph(; scale=nothing)` encodes an exact graph with MOI complementarity. `scale` normalizes hinge variables in physical input units independently of network coordinates; `nothing` retains legacy input-scale normalization."""
struct ComplementarityGraph <: AbstractPWLFormulation
    scale::Union{Nothing,Float64}
    function ComplementarityGraph(; scale=nothing)
        new(isnothing(scale) ? nothing : _positive_width(scale))
    end
end

"""Continuous convex hull of a bounded PWL graph. This is an outer relaxation, not a controller equation."""
struct PWLConvexHull <: AbstractPWLFormulation end

"""`ExactPWLGraph(; method=:Logarithmic)` delegates the exact bounded segment graph to the optional PiecewiseLinearOpt extension."""
struct ExactPWLGraph <: AbstractPWLFormulation
    method::Symbol
end
ExactPWLGraph(; method::Symbol=:Logarithmic) = ExactPWLGraph(method)

function _pwl_hinges(f::PWLFunction)
    slopes = diff(collect(f.values)) ./ diff(collect(f.breakpoints))
    changes = diff([0.; slopes; 0.])
    return [(c, k) for (c, k) in zip(changes, f.breakpoints) if !iszero(c)]
end

function _pwl_exact(f::PWLFunction, x)
    xs, ys = f.breakpoints, f.values
    x <= first(xs) && return first(ys) + zero(x)
    x >= last(xs) && return last(ys) + zero(x)
    i = findlast(k -> k <= x, xs)
    return ys[i] + (ys[i+1]-ys[i]) * ((x-xs[i])/(xs[i+1]-xs[i]))
end

# Integrating the compact nonnegative kernel 3(1-t^2)/4 twice yields this
# monotone, convex C2 hinge; its maximum error is 3*width/16 at the kink.
function _c2_hinge(x, width)
    x <= -width && return zero(x)
    x >= width && return x
    t = x/width
    return width * (3/16 + t/2 + 3*t^2/8 - t^4/16)
end
function _c2_hinge_derivatives(x, width)
    x <= -width && return (zero(x), zero(x))
    x >= width && return (one(x), zero(x))
    t = x/width
    return (1/2 + 3*t/4 - t^3/4, 3*(1-t^2)/(4*width))
end
function _softplus_hinge_derivatives(x, width)
    z = x/width
    e = exp(-abs(z))
    first = z >= 0 ? 1/(1+e) : e/(1+e)
    return first, e/((1+e)^2*width)
end

_pwl_smooth(f, x, r::SoftplusFormulation) = BMOPFTools.piecewise_linear_value(
    x, collect(f.breakpoints), collect(f.values); epsilon=r.width)
_pwl_smooth(f, x, r::AbstractPWLSmoothing) = first(f.values) + sum(
    c*hinge_value(x-k,r) for (c,k) in _pwl_hinges(f); init=zero(x))

"""`hinge_value(x, smoothing)` evaluates the smooth approximation to max(0,x). Extension method for new families."""
function hinge_value end
hinge_value(x,r::LocalC2Formulation) = _c2_hinge(x,r.width)
function hinge_value(x,r::AlgebraicFormulation)
    root = hypot(x,r.width)
    # Rationalize the negative tail to avoid subtractive cancellation.
    return x>=0 ? x/2+root/2 : (r.width/2)*(r.width/(root-x))
end
# The PWL softplus value path delegates directly to BMOPFTools; this scalar
# identity is stable and useful for independently inspecting its hinge family.
function hinge_value(x,r::SoftplusFormulation)
    z = x/r.width
    return max(x,zero(x)) + r.width*log1p(exp(-abs(z)))
end

"""`hinge_derivatives(x, smoothing)` returns analytic first and second derivatives. Extension method for new families."""
function hinge_derivatives end
hinge_derivatives(x,r::LocalC2Formulation) = _c2_hinge_derivatives(x,r.width)
hinge_derivatives(x,r::SoftplusFormulation) = _softplus_hinge_derivatives(x,r.width)
function hinge_derivatives(x,r::AlgebraicFormulation)
    root = hypot(x,r.width)
    return ((1+x/root)/2,(r.width/root)^2/(2root))
end

"""
    hinge_contract(smoothing)

Declare a hinge's signed error interval, regularity, physical width and maximum
absolute second derivative. Custom declarations are researcher-supplied claims,
not automatically proved. The built-in bounds hold in real arithmetic.
"""
function hinge_contract end
hinge_contract(r::SoftplusFormulation) = (error_lower=0.,error_upper=r.width*log(2.),
    regularity=:C_infinity,width=r.width,second_derivative_bound=1/(4r.width))
hinge_contract(r::LocalC2Formulation) = (error_lower=0.,error_upper=3r.width/16,
    regularity=:C2,width=r.width,second_derivative_bound=3/(4r.width))
hinge_contract(r::AlgebraicFormulation) = (error_lower=0.,error_upper=r.width/2,
    regularity=:C_infinity,width=r.width,second_derivative_bound=1/(2r.width))

"""
    primitive_value(curve, x[, formulation])

Evaluate in physical units on the declared domain. The default is the canonical
PWL function; a smooth formulation evaluates its surrogate. `domain_policy=:flat_extension`
explicitly permits evaluating the flat-extended curve (or its smoothing) outside
the declared domain, as needed for controller tails. A convex hull has no
single-valued evaluator. Numeric softplus delegates to BMOPFTools (Float64);
use `primitive_derivatives` for its analytic derivatives.
"""
function primitive_value(f::PWLFunction, x::Real; domain_policy=:error)
    _check_pwl_domain(f, x, domain_policy)
    return _pwl_exact(f, x)
end
function primitive_value(f::PWLFunction, x::Real,
                         r::AbstractPWLSmoothing; domain_policy=:error)
    _check_pwl_domain(f, x, domain_policy)
    return _pwl_smooth(f, x, r)
end
primitive_value(f::PWLFunction, x::Real, ::Union{ComplementarityGraph,ExactPWLGraph}) =
    primitive_value(f, x)
_pwl_primal(x) = x
_pwl_primal(x::ForwardDiff.Dual) = _pwl_primal(ForwardDiff.value(x))
function _check_pwl_domain(f, x, policy=:error)
    policy in (:error,:flat_extension) || throw(ArgumentError("Unknown domain policy $policy"))
    value = _pwl_primal(x)
    isfinite(value) && (policy == :flat_extension || first(f.breakpoints) <= value <= last(f.breakpoints)) ||
        throw(DomainError(x, "Input is outside the declared physical domain"))
end

"""`primitive_derivatives(curve, x, smooth_formulation)` returns physical first and second derivatives, including at patch joins."""
function primitive_derivatives(f::PWLFunction, x::Real,
                               r::AbstractPWLSmoothing; domain_policy=:error)
    _check_pwl_domain(f, x, domain_policy)
    return _pwl_derivatives(f, x, r)
end
function _pwl_derivatives(f, x, r)
    d1, d2 = zero(x), zero(x)
    for (c,k) in _pwl_hinges(f)
        h1, h2 = hinge_derivatives(x-k,r)
        d1 += c*h1
        d2 += c*h2
    end
    return d1, d2
end

"""
    formulation_contract(curve, formulation)

Inspect domain, units, semantics, regularity and a signed physical output-error
interval (`surrogate - exact`). Hull error bounds are `nothing`: its feasible
points are not a surrogate function. Smooth bounds are analytic real-arithmetic
bounds; floating-point evaluation and solver feasibility errors are additional.
"""
function formulation_contract(f::PWLFunction, r::AbstractPWLFormulation)
    r isa Union{AbstractPWLSmoothing,PWLConvexHull,ExactPWLGraph,ComplementarityGraph} ||
        throw(ArgumentError("Define formulation_contract for $(typeof(r))"))
    smooth = r isa AbstractPWLSmoothing
    hull = r isa PWLConvexHull
    hc = smooth ? hinge_contract(r) : nothing
    lower, upper = if smooth
        coeffs = first.(_pwl_hinges(f))
        (sum(min(c*hc.error_lower,c*hc.error_upper) for c in coeffs; init=0.),
         sum(max(c*hc.error_lower,c*hc.error_upper) for c in coeffs; init=0.))
    elseif hull
        (nothing,nothing)
    else
        (0.,0.)
    end
    return (domain=(first(f.breakpoints),last(f.breakpoints)),
        input_unit=f.input_unit,output_unit=f.output_unit,
        semantics=hull ? :outer_relaxation : smooth ? :surrogate_graph : :exact_graph,
        regularity=smooth ? hc.regularity : :not_applicable,
        error_lower=lower,error_upper=upper,width=smooth ? hc.width : nothing)
end

"""
    smoothing_for_error(curve, family, budget)

Choose a physical input width from a positive absolute output-error budget using
signed slope-change sums. This controls function approximation, not solver or
network error. A constant curve needs no smoothing and requires an explicit width
if a smoothing object is still desired. Custom families can extend this method.
"""
function smoothing_for_error(f::PWLFunction, family::Type{<:AbstractPWLSmoothing}, budget::Real)
    family in (SoftplusFormulation,LocalC2Formulation,AlgebraicFormulation) ||
        throw(ArgumentError("Custom families must define their error-to-width mapping"))
    error = _positive_width(budget)
    contract = formulation_contract(f,family(1.))
    magnitude = max(abs(contract.error_lower),abs(contract.error_upper))
    magnitude > 0 || throw(ArgumentError("Constant curves need no smoothing; choose a width explicitly"))
    return family(error/magnitude)
end
