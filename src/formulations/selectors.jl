"""
    hinge_expression(model_or_context, input, smoothing; input_scale=1, output_scale=1)

Build a global smooth positive-part expression. Physical input is
`input_scale * input`; output is physical hinge value divided by `output_scale`.
Custom families reuse their analytic value/derivative interface. Unlike a bounded
PWL graph this primitive imposes no domain constraint.
"""
function hinge_expression(target,x,r::AbstractPWLSmoothing;input_scale::Real=1,output_scale::Real=1)
    si,so = _pwl_scale(input_scale),_pwl_scale(output_scale)
    x isa Real && return hinge_value(x*si,r)/so
    r isa AlgebraicFormulation && return (x+sqrt(x^2+(r.width/si)^2))*(si/(2so))
    model = target isa JuMP.Model ? target : BMOPFTools.opf_model(target)
    cache = get!(model.ext,:PowerOptLabHingeOperators) do
        Dict{Any,Any}()
    end
    op = get!(cache,(r,si,so)) do
        add_nonlinear_operator(model,1,
            x -> hinge_value(x*si,r)/so,
            x -> hinge_derivatives(x*si,r)[1]*si/so,
            x -> hinge_derivatives(x*si,r)[2]*si^2/so;
            name=gensym(:pol_hinge))
    end
    return op(x)
end

function _selector_kind(kind,r)
    kind in (:max,:min,:nonnegative_min) || throw(ArgumentError("Unknown selector kind $kind"))
    if kind == :nonnegative_min
        h = hinge_contract(r)
        h.error_lower>=0 && isfinite(h.error_upper) && h.error_upper>=h.error_lower &&
            hinge_value(0.,r)>0 || throw(ArgumentError(
                "Nonnegative minimum needs an upper hinge with finite error and a positive value at zero"))
    end
end

# Preserve the established symmetric arithmetic graph of IBR selectors. Generic
# hinge composition is supplied separately for other smooth families.
function _algebraic_selector(a,b,width,kind)
    upper = (a+b+sqrt((a-b)^2+width^2))/2
    return kind == :max ? upper : kind == :min ?
        (a+b-sqrt((a-b)^2+width^2))/2 : a*b/upper
end

"""
    selector_value(a, b, smoothing; kind=:max)

Numeric smooth binary maximum/minimum in physical units. `kind=:nonnegative_min`
uses `a*b/smooth_max(a,b)`, requires `a,b>=0`, preserves
zero and lies below the exact minimum. Ordinary smooth minima may be negative
at zero and should not be substituted blindly into capability scale composition.
"""
function selector_value(a::Real,b::Real,r::AbstractPWLSmoothing;kind::Symbol=:max)
    _selector_kind(kind,r)
    all(isfinite,(a,b)) || throw(ArgumentError("Selector inputs must be finite"))
    kind == :nonnegative_min && min(a,b)<0 && throw(DomainError((a,b),"Need nonnegative inputs"))
    r isa AlgebraicFormulation && return _algebraic_selector(a,b,r.width,kind)
    h = hinge_value(a-b,r)
    return kind == :max ? b+h : kind == :min ? a-h : a*b/(b+h)
end

"""
    selector_expression(model_or_context, a, b, smoothing;
                        kind=:max, input_scale=1, output_scale=1)

Symbolic counterpart of `selector_value`. Inputs use working coordinates; both
share `input_scale`, and the returned output is divided by `output_scale`.
For `:nonnegative_min`, the caller must ensure nonnegative inputs in the model.
No domain inequalities or extra variables are introduced by this function.
"""
function selector_expression(target,a,b,r::AbstractPWLSmoothing;kind::Symbol=:max,
                             input_scale::Real=1,output_scale::Real=1)
    _selector_kind(kind,r)
    si,so = _pwl_scale(input_scale),_pwl_scale(output_scale)
    a isa Real && b isa Real && return selector_value(a*si,b*si,r;kind)/so
    r isa AlgebraicFormulation && return _algebraic_selector(a,b,r.width/si,kind)*(si/so)
    h = hinge_expression(target,a-b,r;input_scale=si,output_scale=so)
    upper = b*(si/so)+h
    return kind == :max ? upper : kind == :min ? a*(si/so)-h :
        (a*(si/so))*(b*(si/so))/upper
end

"""
    selector_contract(smoothing; kind=:max)

Signed real-arithmetic error interval relative to the selected hard maximum or
minimum, with its domain assumptions. Bounds for custom families inherit the
author's hinge contract. A nonnegative minimum requires an upper hinge with finite
error bound `B` and a positive value at zero; its deficit is at most `B`.
"""
function selector_contract(r::AbstractPWLSmoothing;kind::Symbol=:max)
    _selector_kind(kind,r)
    h = hinge_contract(r)
    lo,hi = kind == :max ? (h.error_lower,h.error_upper) :
        kind == :min ? (-h.error_upper,-h.error_lower) : (-h.error_upper,0.)
    return (error_lower=lo,error_upper=hi,regularity=h.regularity,width=h.width,
        domain=kind == :nonnegative_min ? :nonnegative_inputs : :real_inputs,
        preserves_zero=kind == :nonnegative_min)
end

"""
    symmetric_clip_value(value, limit, smoothing)

Smooth clipping to `[-limit, limit]`, `limit>=0`, using the symmetric sum
`smooth_min(value,limit)+smooth_max(value,-limit)-value`. For built-in hinges
with error `[0,B]`, its signed error at a fixed limit is in `[-B,B]`.
"""
function symmetric_clip_value(x::Real,limit::Real,r::AbstractPWLSmoothing)
    isfinite(limit) && limit>=0 || throw(DomainError(limit,"Need a finite nonnegative limit"))
    return selector_value(x,limit,r;kind=:min)+selector_value(x,-limit,r;kind=:max)-x
end

"""
    symmetric_clip_expression(model_or_context, value, limit, smoothing;
                              input_scale=1, output_scale=1)

Symbolic counterpart of `symmetric_clip_value`. The caller must ensure a
nonnegative limit. This expression adds no exact physical capability constraint.
"""
function symmetric_clip_expression(target,x,limit,r::AbstractPWLSmoothing;
                                  input_scale::Real=1,output_scale::Real=1)
    si,so = _pwl_scale(input_scale),_pwl_scale(output_scale)
    limit isa Real && limit<0 && throw(DomainError(limit,"Need a nonnegative limit"))
    return selector_expression(target,x,limit,r;kind=:min,input_scale=si,output_scale=so)+
        selector_expression(target,x,-limit,r;kind=:max,input_scale=si,output_scale=so)-x*(si/so)
end
