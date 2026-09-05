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
        all(>(0), diff(collect(xs))) || throw(ArgumentError("Knots must be strictly increasing"))
        slopes = diff(collect(ys)) ./ diff(collect(xs))
        all(isfinite, slopes) && all(isfinite, diff([0.; slopes; 0.])) ||
            throw(ArgumentError("Curve slopes and slope changes must be finite"))
        new(xs, ys, input_unit, output_unit)
    end
end

abstract type AbstractPWLFormulation end

function _positive_width(x)
    value = Float64(x)
    isfinite(value) && value > 0 ||
        throw(ArgumentError("Width must be finite and positive in physical input units"))
    return value
end

"""`SoftplusFormulation(width)` uses BMOPFTools' telescoping softplus with physical input width."""
struct SoftplusFormulation <: AbstractPWLFormulation
    width::Float64
    SoftplusFormulation(width::Real) = new(_positive_width(width))
end

"""`LocalC2Formulation(width)` replaces each hinge in a band of physical half-width `width` by a C2 quartic."""
struct LocalC2Formulation <: AbstractPWLFormulation
    width::Float64
    LocalC2Formulation(width::Real) = new(_positive_width(width))
end

"""Exact PWL graph encoded by nonnegative hinge pairs and native MOI complementarity; requires a suitable external solver."""
struct ComplementarityGraph <: AbstractPWLFormulation end

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
_pwl_smooth(f, x, r::LocalC2Formulation) = first(f.values) + sum(
    c*_c2_hinge(x-k, r.width) for (c,k) in _pwl_hinges(f); init=zero(x))

"""
    primitive_value(curve, x[, formulation])

Evaluate in physical units on the declared domain. The default is the canonical
PWL function; a smooth formulation evaluates its surrogate. A convex hull has no
single-valued evaluator. Numeric softplus delegates to BMOPFTools (Float64);
use `primitive_derivatives` for its analytic derivatives.
"""
function primitive_value(f::PWLFunction, x::Real)
    _check_pwl_domain(f, x)
    return _pwl_exact(f, x)
end
function primitive_value(f::PWLFunction, x::Real,
                         r::Union{SoftplusFormulation,LocalC2Formulation})
    _check_pwl_domain(f, x)
    return _pwl_smooth(f, x, r)
end
primitive_value(f::PWLFunction, x::Real, ::Union{ComplementarityGraph,ExactPWLGraph}) =
    primitive_value(f, x)
_pwl_primal(x) = x
_pwl_primal(x::ForwardDiff.Dual) = _pwl_primal(ForwardDiff.value(x))
function _check_pwl_domain(f, x)
    value = _pwl_primal(x)
    isfinite(value) && first(f.breakpoints) <= value <= last(f.breakpoints) ||
        throw(DomainError(x, "Input is outside the declared physical domain"))
end

"""`primitive_derivatives(curve, x, smooth_formulation)` returns physical first and second derivatives, including at patch joins."""
function primitive_derivatives(f::PWLFunction, x::Real,
                               r::Union{SoftplusFormulation,LocalC2Formulation})
    _check_pwl_domain(f, x)
    return _pwl_derivatives(f, x, r)
end
function _pwl_derivatives(f, x, r)
    d1, d2 = zero(x), zero(x)
    for (c,k) in _pwl_hinges(f)
        h1, h2 = r isa LocalC2Formulation ? _c2_hinge_derivatives(x-k, r.width) :
            _softplus_hinge_derivatives(x-k, r.width)
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
    smooth = r isa Union{SoftplusFormulation,LocalC2Formulation}
    hull = r isa PWLConvexHull
    factor = r isa SoftplusFormulation ? log(2.) : 3/16
    lower, upper = if smooth
        coeffs = first.(_pwl_hinges(f))
        (r.width*factor*sum(min(0.,c) for c in coeffs; init=0.),
         r.width*factor*sum(max(0.,c) for c in coeffs; init=0.))
    elseif hull
        (nothing, nothing)
    else
        (0., 0.)
    end
    return (domain=(first(f.breakpoints),last(f.breakpoints)),
        input_unit=f.input_unit, output_unit=f.output_unit,
        semantics=hull ? :outer_relaxation : smooth ? :surrogate_graph : :exact_graph,
        regularity=r isa SoftplusFormulation ? :C_infinity : r isa LocalC2Formulation ? :C2 : :piecewise_linear,
        error_lower=lower, error_upper=upper,
        width=smooth ? r.width : nothing)
end
