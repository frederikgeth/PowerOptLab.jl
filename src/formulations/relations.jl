"""
    PWLRelationPlan

Inspectable, model-independent lowering decision for a PWL relation on a physical
`domain`. Fields include `relation` (`:equal`, `:upper`, `:lower`), `shape`,
`strategy`, `semantics`, `reason`/`reason_code`, `conservative`, physical
`output_shift`, affine `lines`, and remaining nonlinear smooth
`active_hinges` (empty for affine or graph strategies).
Plans describe scalar relations, not convexity or solvability of an entire model.
"""
struct PWLRelationPlan
    curve::PWLFunction
    domain::Tuple{Float64,Float64}
    relation::Symbol
    formulation::Union{Symbol,AbstractPWLFormulation}
    shape::Symbol
    strategy::Symbol
    semantics::Symbol
    reason::String
    reason_code::Symbol
    conservative::Bool
    output_shift::Float64
    lines::Tuple
    active_hinges::Tuple
end

function _relation_domain(domain)
    length(domain)==2 || throw(ArgumentError("Domain needs two physical endpoints"))
    lo,hi = Float64.(domain)
    all(isfinite,(lo,hi)) && lo<=hi || throw(ArgumentError("Need a finite, nonempty domain"))
    return (lo,hi)
end

function _relation_segments(f,domain)
    lo,hi = domain
    lo==hi && return (((0.,_pwl_exact(f,lo)),),:affine)
    lines = Tuple{Float64,Float64}[]
    # Exact rational comparisons of the supplied Float64 data avoid declaring
    # convexity merely because small curvature was rounded to a tolerance.
    slopes = Rational{BigInt}[]
    if lo<first(f.breakpoints)
        push!(lines,(0.,first(f.values))); push!(slopes,0)
    end
    for i in eachindex(f.slopes)
        a,b = f.breakpoints[i:i+1]
        a<hi && b>lo || continue
        s = f.slopes[i]
        push!(lines,(s,f.values[i]-s*a))
        push!(slopes,(Rational{BigInt}(f.values[i+1])-Rational{BigInt}(f.values[i]))/
            (Rational{BigInt}(b)-Rational{BigInt}(a)))
    end
    if hi>last(f.breakpoints)
        push!(lines,(0.,last(f.values))); push!(slopes,0)
    end
    changes = diff(slopes)
    shape = all(iszero,changes) ? :affine : all(<=(0),changes) ? :concave :
        all(>=(0),changes) ? :convex : :neither
    return Tuple(unique(lines)),shape
end

"""
    plan_pwl_relation(curve, domain; relation=:equal, formulation=:auto, specialize=true, conservative=false)

Plan `y == f(x)`, `y <= f(x)` (`:upper`), or `y >= f(x)` (`:lower`) in
physical coordinates. `:auto` selects only proven exact affine or polyhedral
relations; otherwise the plan is `:unresolved` and building requires an explicit
formulation. Concave upper bounds and convex lower bounds use supporting lines.
No near-convexity tolerance or inference from the objective is used.

Explicit smooth formulations preserve the chosen surrogate. Local C2 hinges are
removed exactly outside their patches; softplus/algebraic tails are retained.
An explicit graph/hull representation may also simplify when the relation is
exactly polyhedral. `specialize=false` retains an explicit representation for
comparison. A hull that is not exact remains labeled `:outer_relaxation`.

For an explicit smooth upper/lower relation, `conservative=true` shifts the
surrogate by `-error_upper` / `-error_lower`, respectively. This is an inner
approximation of the canonical scalar relation under the declared real-arithmetic
error contract. It can conflict with other constraints (e.g. nonnegative power
at a zero cap); no clipping or fallback changes that model. Equality and nonsmooth
requests reject this option. Finite solver residuals remain additional errors.
"""
function plan_pwl_relation(f::PWLFunction,domain;relation::Symbol=:equal,
        formulation=:auto,specialize::Bool=true,conservative::Bool=false)
    relation in (:equal,:upper,:lower) || throw(ArgumentError("Choose :equal, :upper or :lower"))
    formulation isa AbstractPWLFormulation || formulation === :auto ||
        throw(ArgumentError("Supply :auto or a PWL formulation"))
    formulation === :auto && !specialize && throw(ArgumentError("Disabling specialization needs an explicit formulation"))
    conservative && (!(formulation isa AbstractPWLSmoothing) || relation==:equal) &&
        throw(ArgumentError("Conservative correction needs an explicit smooth upper or lower relation"))
    bounds = _relation_domain(domain)
    lines,shape = _relation_segments(f,bounds)
    r = formulation
    strategy,semantics,reason,active = :unresolved,:exact_relation,"An explicit encoding is required",()
    if r isa AbstractPWLSmoothing
        semantics = :smooth_surrogate
        if bounds[1]==bounds[2]
            lines = ((0.,_pwl_smooth(f,bounds[1],r)),)
            strategy,reason = :affine,"Fixed input: evaluate the selected surrogate"
        elseif specialize && r isa LocalC2Formulation
            slope,intercept = 0.,first(f.values)
            remaining = Tuple{Float64,Float64}[]
            for (c,k) in f.hinges
                if bounds[1]>=k+r.width
                    slope += c; intercept -= c*k
                elseif bounds[2]>k-r.width
                    push!(remaining,(c,k))
                end
            end
            lines,active = ((slope,intercept),),Tuple(remaining)
            strategy = isempty(active) ? :affine : :local_c2
            reason = "Remove only hinges whose entire compact patch is outside the domain"
        elseif all(==(first(f.values)),f.values)
            lines = ((0.,first(f.values)),)
            strategy,reason = :affine,"Constant function is unchanged by smoothing"
        else
            strategy,reason = :smooth,"Preserve the full selected smooth function"
            active = f.hinges
        end
    elseif specialize && (shape==:affine || (relation==:upper && shape==:concave) ||
                           (relation==:lower && shape==:convex))
        strategy = shape==:affine ? :affine : :supporting_lines
        shape==:affine && (lines=(first(lines),))
        reason = "Exact $relation relation for a $shape curve on the enforced domain"
    elseif bounds[1]==bounds[2]
        strategy,reason = :affine,"Fixed input: evaluate the canonical curve"
    elseif r isa Union{ExactPWLGraph,ComplementarityGraph,PWLConvexHull}
        strategy = :graph
        hull_exact = shape==:affine || (relation==:upper && shape==:concave) ||
            (relation==:lower && shape==:convex)
        semantics = r isa PWLConvexHull && !hull_exact ? :outer_relaxation : :exact_relation
        reason = "Retain the explicitly selected bounded graph representation"
    end
    reason_code = if strategy==:unresolved
        :explicit_encoding_required
    elseif r isa AbstractPWLSmoothing
        bounds[1]==bounds[2] ? :fixed_surrogate :
            specialize && r isa LocalC2Formulation ? :compact_patch_specialization :
            strategy==:affine ? :constant_surrogate : :full_surrogate
    elseif strategy==:graph
        :explicit_graph
    elseif bounds[1]==bounds[2]
        :fixed_canonical
    else
        :exact_polyhedral
    end
    shift = 0.
    if conservative
        contract = formulation_contract(f,r)
        all(isfinite,(contract.error_lower,contract.error_upper)) &&
            contract.error_lower<=contract.error_upper ||
            throw(ArgumentError("Conservative correction needs finite ordered signed error bounds"))
        shift = -(relation==:upper ? contract.error_upper : contract.error_lower)
        semantics = :inner_approximation
    end
    return PWLRelationPlan(f,bounds,relation,r,shape,strategy,semantics,reason,
        reason_code,conservative,shift,lines,active)
end


function Base.show(io::IO,p::PWLRelationPlan)
    print(io,"PWLRelationPlan(:",p.relation,", ",p.shape," on ",p.domain,
        " → ",p.strategy,", ",length(p.lines)," lines, ",p.semantics)
    p.conservative && print(io,", output shift=",p.output_shift)
    print(io,")")
end
Base.show(io::IO,::MIME"text/plain",p::PWLRelationPlan) = show(io,p)

_built_pwl_formulation(p::PWLRelationPlan) = p.strategy==:graph ?
    (p.formulation isa PWLConvexHull ? :convex_hull :
     p.formulation isa ExactPWLGraph ? :exact_pwl_graph : :complementarity_graph) :
    p.strategy==:supporting_lines ? :linear_inequalities : p.strategy

"""A built contextual relation. Inspect `plan`, `rows`, optional `graph`, and use `audit_pwl_relation` on candidates."""
struct PWLRelationHandle
    model::JuMP.Model
    plan::PWLRelationPlan
    input::Any
    output::Any
    input_scale::Float64
    output_scale::Float64
    rows::Vector{Any}
    graph::Union{Nothing,PWLFormulationHandle}
end

function _declared_relation_domain(x,domain,si)
    lo,hi = domain === nothing ? (-Inf,Inf) : _relation_domain(domain)
    if x isa VariableRef
        if is_fixed(x)
            lo,hi = max(lo,fix_value(x)*si),min(hi,fix_value(x)*si)
        else
            has_lower_bound(x) && (lo=max(lo,lower_bound(x)*si))
            has_upper_bound(x) && (hi=min(hi,upper_bound(x)*si))
        end
    end
    return _relation_domain((lo,hi))
end

function _relation_row!(m,y,expression,relation)
    relation == :equal && return @constraint(m,y == expression)
    relation == :upper && return @constraint(m,y <= expression)
    return @constraint(m,y >= expression)
end

function _restrict_relation_curve(f,domain)
    lo,hi = domain
    knots = [lo; [k for k in f.breakpoints if lo<k<hi]; hi]
    PWLFunction(knots,[_pwl_exact(f,k) for k in knots];input_unit=f.input_unit,output_unit=f.output_unit)
end

"""
    formulate_pwl_relation!(target, curve, input, output; domain=nothing,
        relation=:equal, formulation=:auto, specialize=true, conservative=false,
        input_scale=1, output_scale=1, reuse=true)

Build a contextual scalar relation. Intersect an explicit physical `domain` with
finite declared bounds (or a fixed value) of a `VariableRef` input. General
expressions require an explicit domain: no network bound propagation is performed.
Both endpoints must be finite. Always enforce the effective domain in added rows,
so later loosening variable bounds does not invalidate the specialization.

`:auto` refuses nonpolyhedral relations before adding rows; choose a smooth,
exact-graph, complementarity or hull formulation explicitly. The input and output
are working coordinates with positive physical scales. A shared cache is keyed
by both variables, effective domain, relation, formulation, specialization,
conservative correction and
scales. A law used at another variable or as another relation can lower differently.
No optimizer or objective is selected. Return a `PWLRelationHandle`.
"""
function formulate_pwl_relation!(target,f::PWLFunction,x,y;domain=nothing,
        relation::Symbol=:equal,formulation=:auto,specialize::Bool=true,conservative::Bool=false,
        input_scale::Real=1,output_scale::Real=1,reuse::Bool=true)
    m = target isa JuMP.Model ? target : BMOPFTools.opf_model(target)
    si,so = _pwl_scale(input_scale),_pwl_scale(output_scale)
    for v in (x,y)
        v isa VariableRef && owner_model(v)!==m && throw(ArgumentError("Variable belongs to another model"))
    end
    bounds = _declared_relation_domain(x,domain,si)
    x = x isa Union{AffExpr,QuadExpr} ? copy(x) : x
    y = y isa Union{AffExpr,QuadExpr} ? copy(y) : y
    build = function()
        # Physical analysis is reusable across distinct occurrences/scales;
        # instantiated rows below remain specific to their model variables.
        plans = get!(m.ext,:PowerOptLabPWLRelationPlans) do
            Dict{Any,PWLRelationPlan}()
        end
        plan = get!(plans,(f.breakpoints,f.values,f.input_unit,f.output_unit,
                bounds,relation,formulation,specialize,conservative)) do
            plan_pwl_relation(f,bounds;relation,formulation,specialize,conservative)
        end
        plan.strategy == :unresolved && throw(UnsupportedFormulation(plan.reason))
        plan.strategy == :graph && formulation isa ExactPWLGraph && _require_pwl_graph_extension()
        # Graph lowering already adds its domain row; do not stamp it twice.
        rows = Any[]
        plan.strategy==:graph || push!(rows,@constraint(m,bounds[1]/si <= x <= bounds[2]/si))
        graph = nothing
        if plan.strategy in (:affine,:supporting_lines)
            for (a,b) in plan.lines
                push!(rows,_relation_row!(m,y,(a*si/so)*x+(b+plan.output_shift)/so,relation))
            end
        elseif plan.strategy == :local_c2
            a,b = only(plan.lines)
            expression = (a*si/so)*x+b/so
            for (c,k) in plan.active_hinges
                expression += c*si/so*hinge_expression(target,x-k/si,
                    LocalC2Formulation(formulation.width/si))
            end
            push!(rows,_relation_row!(m,y,expression+plan.output_shift/so,relation))
        elseif plan.strategy == :smooth
            expression = smooth_pwl_expression(target,f,x,formulation;input_scale=si,output_scale=so)
            push!(rows,_relation_row!(m,y,expression+plan.output_shift/so,relation))
        else
            graph = formulate_pwl!(target,_restrict_relation_curve(f,bounds),x,formulation;
                input_scale=si,output_scale=so)
            push!(rows,_relation_row!(m,y,graph.output,relation))
        end
        PWLRelationHandle(m,plan,x,y,si,so,rows,graph)
    end
    reuse && x isa VariableRef && y isa VariableRef || return build()
    cache = get!(m.ext,:PowerOptLabPWLRelations) do
        Dict{Any,PWLRelationHandle}()
    end
    get!(build,cache,(target,f.breakpoints,f.values,f.input_unit,f.output_unit,
        bounds,x,y,relation,formulation,specialize,conservative,si,so))
end

"""
    formulate_control_relation!(target, intent, role, voltage, fraction; kwargs...)

Lower one `VoltVarWattIntent` curve as an equality or bound on an existing output.
`role` is `:volt_watt` or `:volt_var`. Keywords are those of
`formulate_pwl_relation!`, including `relation`, `formulation` and optional domain.
No sensing, power-base multiplication, availability or plant constraints are added.
"""
function formulate_control_relation!(target,intent::VoltVarWattIntent,role::Symbol,x,y;kwargs...)
    role in (:volt_watt,:volt_var) || throw(ArgumentError("Choose :volt_watt or :volt_var"))
    f = getproperty(intent,role)
    f === nothing && throw(ArgumentError("Intent has no $role curve"))
    formulate_pwl_relation!(target,f,x,y;kwargs...)
end

"""
    audit_pwl_relation(handle)

Report signed `output - canonical(input)`, one-sided/equality relation violation,
selected-surrogate violation (including any conservative output shift) when
relevant, and physical domain violation. A slack
upper limit is not a tracking error. These are candidate checks, not solve or
stationarity certificates; a hull's relation violation may be positive. The
optional `graph_audit` separately checks the lifted auxiliary graph, including
complementarity residuals; it does not treat that auxiliary as dispatched power.
"""
function audit_pwl_relation(h::PWLRelationHandle)
    has_values(h.model) || throw(ArgumentError("Model has no candidate values"))
    x,y = value(h.input)*h.input_scale,value(h.output)*h.output_scale
    all(isfinite,(x,y)) || throw(ArgumentError("Candidate input/output is nonfinite"))
    residual = y-_pwl_exact(h.plan.curve,x)
    violation(r) = h.plan.relation==:equal ? abs(r) : h.plan.relation==:upper ? max(r,0.) : max(-r,0.)
    r = h.plan.formulation
    surrogate = r isa AbstractPWLSmoothing ?
        violation(y-(_pwl_smooth(h.plan.curve,x,r)+h.plan.output_shift)) : nothing
    lo,hi = h.plan.domain
    return (input=x,output=y,canonical_residual=residual,
        canonical_violation=violation(residual),surrogate_violation=surrogate,
        domain_violation=max(lo-x,x-hi,0.),strategy=h.plan.strategy,semantics=h.plan.semantics,
        conservative=h.plan.conservative,output_shift=h.plan.output_shift,
        graph_audit=h.graph===nothing ? nothing : audit_pwl(h.graph))
end

function _formulation_observation(h::PWLRelationHandle)
    f,r = h.plan.curve,h.plan.formulation
    # A projected hull relation can be exact even though the auxiliary graph is
    # relaxed. Never overwrite relation semantics with the graph's contract.
    merge(audit_pwl_relation(h),
        (observation_kind=:relation,relation=h.plan.relation,
         domain=h.plan.domain,shape=h.plan.shape,reason=h.plan.reason,reason_code=h.plan.reason_code,
         requested_formulation=r===:auto ? "auto" : string(typeof(r)),
         built_formulation=_built_pwl_formulation(h.plan),
         active_hinges=length(h.plan.active_hinges),
         formulation_type=r===:auto ? "auto" : string(typeof(r)),
         approximation_contract=r isa AbstractPWLSmoothing ? formulation_contract(f,r) : nothing,
         graph_method=r isa ExactPWLGraph ? r.method : nothing,
         complementarity_scale=h.graph===nothing ? nothing : h.graph.complementarity_scale,
         input_scale=h.input_scale,output_scale=h.output_scale,
         curve=(breakpoints=f.breakpoints,values=f.values,
                input_unit=f.input_unit,output_unit=f.output_unit)))
end
