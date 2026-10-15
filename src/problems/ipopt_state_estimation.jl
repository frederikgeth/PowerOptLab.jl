const _SEMOI = JuMP.MOI

# Optimise scaled voltage increments, retaining the common evaluator in SI.
struct _SEIpoptEvaluator <: _SEMOI.AbstractNLPEvaluator
    structure::SEStructure
    parameters::SEParameters
    origin::Vector{Float64}
    levels::Vector{Float64}
    constraint_scales::Vector{Float64}
    constraint_rows::Vector{Int}
end
_SEMOI.features_available(::_SEIpoptEvaluator) = [:Grad, :Jac]
function _SEMOI.initialize(::_SEIpoptEvaluator, features::Vector{Symbol})
    all(f -> f in (:Grad, :Jac), features) || throw(ArgumentError("shared Ipopt evaluator provides first derivatives only; use limited-memory Hessians"))
    nothing
end
_se_ipopt_state(d, y) = d.origin + d.levels .* y
function _SEMOI.eval_objective(d::_SEIpoptEvaluator, y)
    e = _se_trial_evaluation(d.structure, d.parameters, _se_ipopt_state(d, y))
    e === nothing ? NaN : .5sum(abs2, e.residual)
end
function _SEMOI.eval_objective_gradient(d::_SEIpoptEvaluator, g, y)
    x = _se_ipopt_state(d, y)
    e = _se_trial_evaluation(d.structure, d.parameters, x)
    jac = e === nothing ? nothing : _se_jacobians(d.structure, d.parameters, x)
    if jac === nothing
        fill!(g, NaN)
    else
        g .= d.levels .* (jac[1]' * e.residual)
    end
    nothing
end
function _SEMOI.eval_constraint(d::_SEIpoptEvaluator, c, y)
    e = _se_trial_evaluation(d.structure, d.parameters, _se_ipopt_state(d, y))
    e === nothing ? fill!(c, NaN) : (c .= e.constraints[d.constraint_rows] ./ d.constraint_scales)
    nothing
end
# Deliberately complete structure: derivatives can vanish at the initial
# state and become nonzero later (devices and user functions). Do not infer
# symbolic sparsity by sampling one numerical Jacobian.
_SEMOI.jacobian_structure(d::_SEIpoptEvaluator) =
    [(i, j) for j in eachindex(d.origin) for i in eachindex(d.constraint_scales)]
function _SEMOI.eval_constraint_jacobian(d::_SEIpoptEvaluator, J, y)
    jac = _se_jacobians(d.structure, d.parameters, _se_ipopt_state(d, y))
    if jac === nothing
        fill!(J, NaN)
    else
        C = jac[2]; k = 1
        for j in eachindex(d.origin), i in eachindex(d.constraint_scales)
            J[k] = C[d.constraint_rows[i], j] * d.levels[j] / d.constraint_scales[i]; k += 1
        end
    end
    nothing
end

# A Jacobian sampled at one state cannot establish an equation identity.
# Build c(x) = A*x + b + B*j_power(x): A/b include passive physics and
# constant I/Y laws; B retains an independent complex current for every
# nonzero constant-power branch, even when its derivative happens to vanish.
# A dependency of this complete coefficient matrix holds at every state.
function _se_ipopt_constraint_signature(s, p)
    n=length(s.nodes);q=length(s.constraint_pattern)+size(s.return_current_pattern,1)
    selector=zeros(Float64,2q,2n)
    for (k,i) in enumerate(s.constraint_pattern)
        selector[2k-1,i]=1.;selector[2k,n+i]=1.
    end
    offset=length(s.constraint_pattern)
    for (r,i,v) in zip(findnz(s.return_current_pattern)...)
        selector[2(offset+r)-1,i]=v;selector[2(offset+r),n+i]=v
    end
    linear=deepcopy(p)
    fill!(linear.device_powers,0im)
    # The zero-power part is affine and defined at zero voltage. Alpha=0
    # avoids applying the physical nonzero-voltage domain to this decomposition.
    linear.continuation_alpha=0.
    linear.regularization_voltage=1.
    x=zeros(size(s.voltage_state_jacobian,2))
    vr,vi,ir,ii=_se_parts(s,linear,x)
    dr,di,Jr,Ji=_device_parts(s,linear,vr,vi)
    A=selector*(s.current_state_jacobian+vcat(Jr,Ji))
    b=selector*vcat(ir+dr,ii+di)
    nonlinear=Tuple{Int,Int}[]
    for device in s.device_pattern
        device.kind in (:constant_power,:zip) || continue
        for (j,k) in enumerate(device.parameter_range)
            iszero(p.device_powers[k]) || push!(nonlinear,(device.positive[j],device.negative[j]))
        end
    end
    B=zeros(Float64,2n,2length(nonlinear))
    for (k,(i,j)) in enumerate(nonlinear)
        B[i,2k-1]=1.;B[n+i,2k]=1.
        if j!=0
            B[j,2k-1]=-1.;B[n+j,2k]=-1.
        end
    end
    hcat(A,b,selector*B)
end

function _se_ipopt_constraint_basis(s,p)
    _se_constraint_count(s)==0 && return Int[],zeros(0,0)
    signature=_se_ipopt_constraint_signature(s,p)
    m=size(signature,1)
    m==0 && return Int[],zeros(0,0)
    # Equilibrate units before pivoting. This uses a near-roundoff cutoff,
    # independent of the much looser user observability/rank tolerance.
    for j in axes(signature,2)
        level=norm(signature[:,j]);level>0 && (signature[:,j]./=level)
    end
    row_levels=ones(m)
    for i in axes(signature,1)
        level=norm(signature[i,:])
        if level>0
            row_levels[i]=level;signature[i,:]./=level
        end
    end
    factor=qr(transpose(signature),ColumnNorm())
    diagonal=abs.(diag(factor.R))
    tolerance=8eps(Float64)*max(size(signature)...)
    rank=count(>(tolerance*maximum(diagonal;init=1.)),diagonal)
    rows=sort!(factor.p[1:rank])
    length(rows)==m && return rows,Matrix{Float64}(I,m,m)
    # Check the complete coefficient identity, not just the rank estimate.
    # If numerical cancellation prevents that check, retain the full system.
    expansion=isempty(rows) ? zeros(m,0) :
        collect(transpose(transpose(signature[rows,:]) \ transpose(signature)))
    reconstruction=expansion*signature[rows,:]
    for row in axes(signature,1)
        norm(signature[row,:]-reconstruction[row,:],Inf)<=tolerance ||
            return collect(1:m),Matrix{Float64}(I,m,m)
    end
    # Selected rows are identities exactly. Undo row normalization so this
    # map also bounds physical-current amplification (e.g. turns ratios).
    expansion[rows,:].=Matrix{Float64}(I,length(rows),length(rows))
    expansion=expansion.*row_levels./transpose(row_levels[rows])
    rows,expansion
end
_se_ipopt_constraint_rows(s,p)=first(_se_ipopt_constraint_basis(s,p))

"""
Outcome of `solve_ipopt_state_estimator`. `state` and `evaluation` use the
same SI model as the dense and sparse solvers. Optimizer statuses are retained
separately from checked physical feasibility, stationarity, curvature and rank.
Use `solve_status(result).publishable`, not an Ipopt termination code alone.
`constraint_rows` identifies the physical equations sent to Ipopt after removing
verified coefficient identities. `evaluation.constraints` always retains every
physical equation; diagnostics report both counts.
"""
struct IpoptStateEstimationResult <: AbstractSolveResult
    status::Symbol
    state::Vector{Float64}
    evaluation::SEEvaluation
    constraint_rank::Int
    tangent_dimension::Int
    observable_dimension::Int
    termination_status::String
    primal_status::String
    constraint_rows::Vector{Int}
end
# Preserve the original result-construction signature for callers retaining
# results produced before structural equation reduction was introduced.
IpoptStateEstimationResult(status,state,evaluation,constraint_rank,tangent_dimension,
                          observable_dimension,termination_status,primal_status) =
    IpoptStateEstimationResult(status,state,evaluation,constraint_rank,tangent_dimension,
        observable_dimension,termination_status,primal_status,collect(eachindex(evaluation.constraints)))
function solve_status(r::IpoptStateEstimationResult)
    ok = r.status in (:converged_unique, :converged_underobserved)
    _result_solve_status(string(r.status), ok; primal_status=r.primal_status)
end
solve_diagnostics(r::IpoptStateEstimationResult) =
    (constraint_rank=r.constraint_rank, tangent_dimension=r.tangent_dimension,
     observable_dimension=r.observable_dimension, objective=.5sum(abs2, r.evaluation.residual),
     optimizer_termination=r.termination_status, optimizer_primal=r.primal_status,
     physical_constraint_count=length(r.evaluation.constraints),
     optimizer_constraint_count=length(r.constraint_rows), optimizer_constraint_rows=copy(r.constraint_rows))

"""
    solve_ipopt_state_estimator(structure, parameters, x0; kwargs...)
    solve_state_estimation(structure, parameters, x0; kwargs...)

Classic Ipopt nonlinear-programming solve over the **same compiled model** as
the dense and sparse constrained estimators. All measurement functions,
transformers, grounding, exact devices, priors and profiled customer powers
therefore have identical semantics. This does not invoke the legacy BMOPF OPF
wrapper or impose its complete P/Q coverage contract.

Uses analytic/automatic first derivatives, voltage/current scaling and Ipopt's
limited-memory Hessian approximation. Constraint Jacobian storage is currently
dense; this backend is a reference/compatibility path, not a scaling claim.
State-independent dependencies are checked using passive/device affine
coefficients and independent nonlinear branch-current incidences. Only a
verified basis is sent to Ipopt; final checks retain all original equations.
An initially rank-deficient Jacobian alone never permits dropping an equation.
`solver_options` accepts Ipopt attribute pairs or a NamedTuple; exact Hessians
are not provided. Final physical feasibility and reduced stationarity are
checked independently of the optimizer's status. Internal homotopy stages
remain nonpublishable.
"""
function solve_ipopt_state_estimator(s::SEStructure, p::SEParameters, x0::AbstractVector{<:Real};
        max_iterations::Integer=500, constraint_tolerance::Real=1e-8,
        constraint_rtol::Real=1e-11, optimality_tolerance::Real=1e-6,
        rank_rtol::Real=sqrt(eps(Float64)), verbose::Bool=false, solver_options=())
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    for v in (constraint_tolerance, constraint_rtol, optimality_tolerance, rank_rtol)
        isfinite(Float64(v)) && v >= 0 || throw(ArgumentError("tolerances must be finite and nonnegative"))
    end
    options = collect(solver_options isa NamedTuple ? pairs(solver_options) : solver_options)
    for (k, v) in options
        string(k) == "hessian_approximation" && v != "limited-memory" &&
            throw(ArgumentError("shared Ipopt backend requires hessian_approximation=limited-memory"))
        string(k) in ("dependency_detector", "dependency_detection_with_rhs") &&
            throw(ArgumentError("unverified constraint dropping is not supported; only state-independent coefficient identities may reduce equations"))
    end
    x = Float64.(x0)
    active_rows = Int[]
    e = _se_trial_evaluation(s, p, x)
    result(status, state, evaluation, rc=0, td=0, od=0, term="NOT_RUN", primal="NO_SOLUTION") =
        IpoptStateEstimationResult(status, state, evaluation, rc, td, od, term, primal, copy(active_rows))
    e === nothing && return result(:invalid_initial_domain, x, _se_invalid_evaluation(s, p))
    jac = _se_jacobians(s, p, x)
    jac === nothing && return result(:undefined_derivative, x, e)
    ctol = _se_constraint_tolerance(constraint_tolerance, constraint_rtol, _se_current_scale(s, p))
    if isempty(x)
        status = norm(e.constraints) <= ctol ? :converged_unique : :constraint_not_satisfied
        if status === :converged_unique
            status = _se_checked_status(s,p,x,e,jac[1],jac[2],zeros(0,0),status;rank_rtol)
        end
        return result(status, x, e, 0, 0, 0, "OPTIMAL", "FEASIBLE_POINT")
    end
    levels = 1 ./ _se_voltage_scale(s, p, x)
    # Mixed voltage, current and transformer-power channels can differ by many
    # orders of magnitude even after voltage-base scaling. Equilibrate observed
    # columns by the whitened Jacobian, retaining voltage scales in directions
    # with no initial information. This is numerical scaling, not a prior.
    for j in eachindex(levels)
        information = norm(jac[1][:, j])
        information > 0 && (levels[j] = min(levels[j], 1 / information))
    end
    active_rows, expansion = _se_ipopt_constraint_basis(s,p)
    scales = [max(norm(jac[2][i, :] .* levels), 1.) for i in active_rows]
    evaluator = _SEIpoptEvaluator(s, p, x, levels, scales, active_rows)
    model = Model(Ipopt.Optimizer)
    verbose || set_silent(model)
    @variable(model, y[1:length(x)], start=0.)
    set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
    # Mixed winding P/Q channels have strongly coupled voltage curvature.
    # Ipopt's six-vector default can repeatedly forget it even on small
    # circuits; a bounded larger history remains linear in the state count.
    set_optimizer_attribute(model, "limited_memory_max_history", 30)
    set_optimizer_attribute(model, "nlp_scaling_method", "none")
    set_optimizer_attribute(model, "max_iter", Int(max_iterations))
    set_optimizer_attribute(model, "tol", max(Float64(optimality_tolerance) / 100, 1e-12))
    # Ipopt checks scaled constraint components. Translate the physical 2-norm
    # budget to those components; its default feasibility tolerance can accept
    # a raw current error above our requested tolerance after row scaling.
    # If |c_active ./ scales| <= t componentwise, the complete physical
    # constraint norm is bounded by norm(abs(expansion)*scales)*t.
    amplification = norm(abs.(expansion)*scales)
    scaled_ctol = max(ctol / (2max(amplification,1.)), 1e-12)
    set_optimizer_attribute(model, "constr_viol_tol", scaled_ctol)
    set_optimizer_attribute(model, "acceptable_constr_viol_tol", scaled_ctol)
    _set_solver_options!(model, options)
    bounds = [_SEMOI.NLPBoundsPair(0., 0.) for _ in scales]
    _SEMOI.set(model, _SEMOI.NLPBlock(), _SEMOI.NLPBlockData(bounds, evaluator, true))
    set_objective_sense(model, _SEMOI.MIN_SENSE)
    optimize!(model)
    outcome = _solve_outcome(model)
    term, primal = string(outcome.termination_status), string(outcome.primal_status)
    if !outcome.has_primal
        return result(:optimizer_failed, x, e, 0, 0, 0, term, primal)
    end
    xf = _se_ipopt_state(evaluator, value.(y))
    ef = _se_trial_evaluation(s, p, xf)
    ef === nothing && return result(:invalid_final_domain, xf, _se_invalid_evaluation(s, p), 0, 0, 0, term, primal)
    jf = _se_jacobians(s, p, xf)
    jf === nothing && return result(:undefined_derivative, xf, ef, 0, 0, 0, term, primal)
    H, C = jf
    rc, Z = _se_rank_nullspace(C; rtol=rank_rtol)
    od, _ = _se_rank_nullspace(H * Z; rtol=rank_rtol, reference_scale=norm(H))
    td = size(Z, 2)
    status = !_publishable(outcome) ? :optimizer_failed :
             norm(ef.constraints) > ctol ? :constraint_not_satisfied :
             _se_stationarity(H, Z' * (H' * ef.residual)) > optimality_tolerance ? :optimality_not_satisfied :
             od == td ? :converged_unique : :converged_underobserved
    if status in (:converged_unique, :converged_underobserved)
        status = _se_checked_status(s, p, xf, ef, H, C, Z, status; rank_rtol)
    end
    result(status, xf, ef, rc, td, od, term, primal)
end

solve_state_estimation(s::SEStructure, p::SEParameters, x0::AbstractVector{<:Real}; kwargs...) =
    solve_ipopt_state_estimator(s, p, x0; kwargs...)
