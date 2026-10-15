const _SEMOI = JuMP.MOI

# Optimise scaled voltage increments, retaining the common evaluator in SI.
struct _SEIpoptEvaluator <: _SEMOI.AbstractNLPEvaluator
    structure::SEStructure
    parameters::SEParameters
    origin::Vector{Float64}
    levels::Vector{Float64}
    constraint_scales::Vector{Float64}
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
    e === nothing ? fill!(c, NaN) : (c .= e.constraints ./ d.constraint_scales)
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
            J[k] = C[i, j] * d.levels[j] / d.constraint_scales[i]; k += 1
        end
    end
    nothing
end

"""
Outcome of `solve_ipopt_state_estimator`. `state` and `evaluation` use the
same SI model as the dense and sparse solvers. Optimizer statuses are retained
separately from checked physical feasibility, stationarity, curvature and rank.
Use `solve_status(result).publishable`, not an Ipopt termination code alone.
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
end
function solve_status(r::IpoptStateEstimationResult)
    ok = r.status in (:converged_unique, :converged_underobserved)
    _result_solve_status(string(r.status), ok; primal_status=r.primal_status)
end
solve_diagnostics(r::IpoptStateEstimationResult) =
    (constraint_rank=r.constraint_rank, tangent_dimension=r.tangent_dimension,
     observable_dimension=r.observable_dimension, objective=.5sum(abs2, r.evaluation.residual),
     optimizer_termination=r.termination_status, optimizer_primal=r.primal_status)

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
            throw(ArgumentError("constraint dropping is not supported; exact equations must remain in the model"))
    end
    x = Float64.(x0)
    e = _se_trial_evaluation(s, p, x)
    result(status, state, evaluation, rc=0, td=0, od=0, term="NOT_RUN", primal="NO_SOLUTION") =
        IpoptStateEstimationResult(status, state, evaluation, rc, td, od, term, primal)
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
    scales = [max(norm(jac[2][i, :] .* levels), 1.) for i in axes(jac[2], 1)]
    evaluator = _SEIpoptEvaluator(s, p, x, levels, scales)
    model = Model(Ipopt.Optimizer)
    verbose || set_silent(model)
    @variable(model, y[1:length(x)], start=0.)
    set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
    set_optimizer_attribute(model, "nlp_scaling_method", "none")
    set_optimizer_attribute(model, "max_iter", Int(max_iterations))
    set_optimizer_attribute(model, "tol", max(Float64(optimality_tolerance) / 100, 1e-12))
    # Ipopt checks scaled constraint components. Translate the physical 2-norm
    # budget to those components; its default feasibility tolerance can accept
    # a raw current error above our requested tolerance after row scaling.
    scaled_ctol = max(ctol / (2sqrt(max(1,length(scales))) * maximum(scales;init=1.)), 1e-12)
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
    od, _ = _se_rank_nullspace(H * Z; rtol=rank_rtol)
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
