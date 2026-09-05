"""
    MagnitudeApproximation(epsilon; direction=:lower, unit=:unitless)

Smooth Euclidean magnitude with a physical absolute error budget. `:lower` uses
sqrt(sum(x²)+epsilon²)-epsilon and preserves zero; `:upper` omits the shift and
bounds the exact magnitude from above. Direction is mathematical, not a complete
certificate that a composed network inequality is conservative.
"""
struct MagnitudeApproximation
    epsilon::Float64
    direction::Symbol
    unit::Symbol
    function MagnitudeApproximation(epsilon::Real;direction::Symbol=:lower,unit::Symbol=:unitless)
        direction in (:lower,:upper) || throw(ArgumentError("Choose :lower or :upper"))
        new(_positive_width(epsilon),direction,unit)
    end
end

"""`magnitude_value(components, approximation)` evaluates a smooth magnitude in physical units; components must be real."""
function magnitude_value(components,r::MagnitudeApproximation)
    isempty(components) && throw(ArgumentError("Provide at least one component"))
    all(x -> x isa Real && isfinite(x),components) || throw(ArgumentError("Components must be finite real numbers"))
    upper = norm([components...;r.epsilon])
    return r.direction == :lower ? upper-r.epsilon : upper
end

"""`magnitude_contract(approximation)` reports signed magnitude error and the physical Hessian spectral-norm bound 1/epsilon."""
magnitude_contract(r::MagnitudeApproximation) = (
    error_lower=r.direction == :lower ? -r.epsilon : 0.,
    error_upper=r.direction == :lower ? 0. : r.epsilon,
    unit=r.unit,regularity=:C_infinity,hessian_bound=1/r.epsilon,
    preserves_zero=r.direction == :lower)

"""
    magnitude_expression(model_or_context, components, approximation;
                         component_scale=1, output_scale=1)

Build the magnitude from working-coordinate real components. Physical components
are `component_scale * components`; the returned expression is in units of
`physical magnitude / output_scale`. On staged contexts reuse BMOPFTools.smooth_norm.
The optional `name` labels the upstream differentiability annotation.
This adds no capability constraint: squared physical inequalities remain available.
"""
function magnitude_expression(target,components,r::MagnitudeApproximation;
                              component_scale::Real=1,output_scale::Real=1,
                              name::AbstractString="")
    isempty(components) && throw(ArgumentError("Provide at least one component"))
    si,so = _pwl_scale(component_scale),_pwl_scale(output_scale)
    if target isa JuMP.Model
        root = sqrt(sum(c^2 for c in components)+(r.epsilon/si)^2)
        return (r.direction == :lower ? root-r.epsilon/si : root)*(si/so)
    end
    lower = BMOPFTools.smooth_norm(target,collect(components);scale=1.,
        eps_rel=r.epsilon/si,name=name)
    return lower*(si/so)+(r.direction == :upper ? r.epsilon/so : 0.)
end

"""
    positive_root_expression(model_or_context, radicand; lower_bound)

Add the strictly positive domain `radicand >= lower_bound` and return its exact
square root. Both radicand and lower bound use the same working units. This changes
the allowed domain explicitly; it is not a smoothing of a root at zero.
"""
function positive_root_expression(target,radicand;lower_bound::Real)
    lower = _positive_width(lower_bound)
    model = target isa JuMP.Model ? target : BMOPFTools.opf_model(target)
    @constraint(model,radicand >= lower)
    return sqrt(radicand)
end

"""
    affine_error_bound(weights, intervals)

Propagate signed real-arithmetic error intervals through a fixed affine
combination. Correlations are ignored, so the interval may be conservative. All
weighted errors must have compatible output units, which the caller supplies.
"""
function affine_error_bound(weights,intervals)
    length(weights)==length(intervals) || throw(DimensionMismatch("Weights and intervals must match"))
    lo,hi = 0.,0.
    for (w,(a,b)) in zip(weights,intervals)
        all(isfinite,(w,a,b)) && a<=b || throw(ArgumentError("Need finite weights and ordered intervals"))
        lo += min(w*a,w*b)
        hi += max(w*a,w*b)
    end
    return (lower=lo,upper=hi)
end

"""
    local_error_response(jacobian, residual_perturbation; condition_limit=1e10)

Compute the first-order estimate `delta = -J \\ perturbation` for a square
linearization. Return no estimate for singular/poorly conditioned matrices. This
is a local diagnostic, not a nonlinear error bound or a uniqueness certificate.
Condition numbers depend on the supplied row/state scaling. The diagnostic uses a
dense matrix and is intended for small examples or reduced sensitivity systems.
"""
function local_error_response(jacobian::AbstractMatrix,perturbation::AbstractVector;
                              condition_limit::Real=1e10)
    n,m = size(jacobian)
    n==m==length(perturbation) && n>0 || throw(DimensionMismatch("Need a nonempty square compatible system"))
    condition_limit>=1 || throw(ArgumentError("Condition limit must be at least one"))
    J,p = Matrix{Float64}(jacobian),Vector{Float64}(perturbation)
    all(isfinite,J) && all(isfinite,p) || throw(ArgumentError("Linearization must be finite"))
    condition = cond(J)
    usable = isfinite(condition) && condition<=condition_limit
    return (status=usable ? :estimated : :ill_conditioned,
        delta=usable ? -(J\p) : nothing,condition_number=condition,
        interpretation=:first_order_estimate)
end
