"""
    PWLFormulationHandle

A built representation: canonical `curve`, `formulation`, working-coordinate
`input` and `output`, coordinate scales, physical voltage/input `domain`, and
complementarity `pairs`. The domain can differ from the curve's outer knots for
smooth control ports, without changing the smoothed curve. Inspect its
semantics with `formulation_contract(handle.curve, handle.formulation)` and its
candidate values with `audit_pwl`. Handles do not certify a network solution.
"""
struct PWLFormulationHandle
    model::JuMP.Model
    curve::PWLFunction
    formulation::AbstractPWLFormulation
    input::Any
    output::Any
    input_scale::Float64
    output_scale::Float64
    complementarity_scale::Float64
    pairs::Vector{Tuple{JuMP.VariableRef,JuMP.VariableRef}}
    domain::Tuple{Float64,Float64}
end

function _pwl_scale(x)
    value = Float64(x)
    isfinite(value) && value > 0 || throw(ArgumentError("Coordinate scales must be finite and positive"))
    return value
end

# Register one analytic univariate operator per curve/formulation/coordinate
# system. Numeric softplus values come from BMOPFTools. Explicit derivatives
# avoid passing Dual numbers through its Float64 numeric API.
function _pwl_operator(m, f, r, si, so)
    cache = get!(m.ext, :PowerOptLabPWLOperators) do
        Dict{Any,Any}()
    end
    key = (f.breakpoints, f.values, r, si, so)
    return get!(cache,key) do
        # Capture vector conversion once. The upstream softplus oracle still
        # owns its validation/evaluation; derivative coefficients are prepared
        # in the immutable curve instead of reconstructed at every callback.
        xs,ys = collect(f.breakpoints),collect(f.values)
        evaluate = r isa SoftplusFormulation ?
            (x -> BMOPFTools.piecewise_linear_value(x*si,xs,ys;epsilon=r.width)/so) :
            (x -> _pwl_smooth(f,x*si,r)/so)
        JuMP.add_nonlinear_operator(m, 1,
            evaluate,
            x -> _pwl_derivatives(f,x*si,r)[1]*si/so,
            x -> _pwl_derivatives(f,x*si,r)[2]*si^2/so;
            name=gensym(:pol_pwl))
    end
end

"""
    formulate_pwl!(model_or_context, curve, input, formulation;
                   input_scale=1, output_scale=1)

Build a PWL representation and constrain its input to the declared domain. A
physical input is `input_scale * input`; a physical output is
`output_scale * handle.output`. Both scales must be finite and positive.

A BMOPFTools staged context uses its public softplus expression builder. A plain
JuMP model uses cached scalar operators with analytic first/second derivatives.
Complementarity requires external solver support; the exact segment graph requires
`using PiecewiseLinearOpt`. No optimizer, homotopy, or retry is selected here.
"""
function formulate_pwl!(m::JuMP.Model, f::PWLFunction, x, r::AbstractPWLFormulation;
                        input_scale::Real=1, output_scale::Real=1)
    return _formulate_pwl!(m,nothing,f,x,r,input_scale,output_scale)
end
function formulate_pwl!(ctx, f::PWLFunction, x, r::AbstractPWLFormulation;
                        input_scale::Real=1, output_scale::Real=1)
    return _formulate_pwl!(BMOPFTools.opf_model(ctx),ctx,f,x,r,input_scale,output_scale)
end
function _formulate_pwl!(m,ctx,f,x,r,input_scale,output_scale)
    si,so = _pwl_scale(input_scale),_pwl_scale(output_scale)
    cs = r isa ComplementarityGraph && r.scale !== nothing ? r.scale : si
    # Fail before mutating the model if an optional backend is absent.
    r isa ExactPWLGraph && _require_pwl_graph_extension()
    @constraint(m, first(f.breakpoints)/si <= x <= last(f.breakpoints)/si)
    pairs = Tuple{JuMP.VariableRef,JuMP.VariableRef}[]
    if all(==(first(f.values)),f.values)
        y = first(f.values)/so
    elseif r isa AbstractPWLSmoothing
        expression = smooth_pwl_expression(isnothing(ctx) ? m : ctx,f,x,r;
            input_scale=si,output_scale=so)
        y = @variable(m, base_name="pwl_output")
        @constraint(m, y == expression)
    elseif r isa ComplementarityGraph
        expression = first(f.values)/so
        for (c,k) in _pwl_hinges(f)
            # Endpoint hinges are affine on the closed domain. Eliminating
            # them avoids unnecessary complementarity at the endpoints.
            if k == first(f.breakpoints)
                expression += c*si/so*x-c*k/so
            elseif k != last(f.breakpoints)
                positive = @variable(m, lower_bound=0, base_name="pwl_positive")
                negative = @variable(m, lower_bound=0, base_name="pwl_negative")
                @constraint(m, positive-negative == (si/cs)*x-k/cs)
                @constraint(m, [positive,negative] in JuMP.MOI.Complements(2))
                push!(pairs,(positive,negative))
                expression += c*cs/so*positive
            end
        end
        y = @variable(m, base_name="pwl_output")
        @constraint(m, y == expression)
    elseif r isa PWLConvexHull
        weights = @variable(m, [1:length(f.breakpoints)], lower_bound=0, base_name="pwl_weight")
        @constraint(m, sum(weights) == 1)
        @constraint(m, x == sum(f.breakpoints[i]/si*weights[i] for i in eachindex(weights)))
        y = @expression(m, sum(f.values[i]/so*weights[i] for i in eachindex(weights)))
    elseif r isa ExactPWLGraph
        y = _exact_pwl_graph!(m,x,collect(f.breakpoints)./si,collect(f.values)./so,r)
    else
        throw(ArgumentError("Unsupported PWL formulation $(typeof(r))"))
    end
    return PWLFormulationHandle(m,f,r,x,y,si,so,cs,pairs,
        (first(f.breakpoints),last(f.breakpoints)))
end
function _require_pwl_graph_extension()
    Base.get_extension(@__MODULE__, :PowerOptLabPiecewiseLinearOptExt) === nothing &&
        throw(ArgumentError("ExactPWLGraph requires loading the optional PiecewiseLinearOpt package"))
end
function _exact_pwl_graph! end

"""
    audit_pwl(handle)

Inspect a solver candidate in physical coordinates, independently of solver status.
Reports domain violation, signed exact-graph error, smooth-equation error when
applicable, and complementarity minimum/product residuals in physical input units
and their square. A hull candidate can satisfy its encoding while missing the graph.
No status is promoted and no physical safety or stationarity certificate is issued.
"""
function audit_pwl(h::PWLFormulationHandle)
    JuMP.has_values(h.model) || throw(ArgumentError("Model has no candidate values"))
    x,y = JuMP.value(h.input)*h.input_scale,JuMP.value(h.output)*h.output_scale
    all(isfinite,(x,y)) || throw(ArgumentError("Candidate input/output is nonfinite"))
    lo,hi = h.domain
    domain_violation = max(lo-x,x-hi,0.)
    # The flat extension is used only for residual reporting outside the domain;
    # the contract's error guarantee applies only when the domain is satisfied.
    exact_error = y-_pwl_exact(h.curve,x)
    smooth = h.formulation isa AbstractPWLSmoothing
    surrogate_error = smooth ? y-_pwl_smooth(h.curve,x,h.formulation) : nothing
    values = [(JuMP.value(a)*h.complementarity_scale,JuMP.value(b)*h.complementarity_scale) for (a,b) in h.pairs]
    minimum_residual = maximum((abs(min(a,b)) for (a,b) in values); init=0.)
    product_residual = maximum((abs(a*b) for (a,b) in values); init=0.)
    return (input=x, output=y, domain_violation=domain_violation,
        exact_graph_error=exact_error, surrogate_equation_error=surrogate_error,
        complementarity_minimum=minimum_residual, complementarity_product=product_residual,
        semantics=formulation_contract(h.curve,h.formulation).semantics,
        termination_status=JuMP.termination_status(h.model))
end


"""
    smooth_pwl_expression(model_or_context, curve, input, smoothing;
                          input_scale=1, output_scale=1)

Return a smooth expression of the canonical flat extension, without imposing a
bounded domain or adding an output variable. Suitable for controller curves whose
tails remain active beyond their breakpoints. Use `formulate_pwl!` when the domain
must be enforced. Staged softplus construction delegates to BMOPFTools.
"""
function smooth_pwl_expression(target,f::PWLFunction,x,r::AbstractPWLSmoothing;
                               input_scale::Real=1,output_scale::Real=1)
    si,so = _pwl_scale(input_scale),_pwl_scale(output_scale)
    all(==(first(f.values)),f.values) && return first(f.values)/so
    if !(target isa JuMP.Model) && r isa SoftplusFormulation
        return BMOPFTools.opf_piecewise_linear_expression(target,x,
            collect(f.breakpoints)./si,collect(f.values)./so;epsilon=r.width/si)
    end
    m = target isa JuMP.Model ? target : BMOPFTools.opf_model(target)
    return _pwl_operator(m,f,r,si,so)(x)
end
