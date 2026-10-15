const _SEMOI = JuMP.MOI

# Derivative patterns come from compiled dependencies, never sampled zeros.
include("state_estimation_patterns.jl")

mutable struct _SEIpoptEvaluator <: _SEMOI.AbstractNLPEvaluator
    structure::SEStructure
    parameters::SEParameters
    origin::Vector{Float64}
    levels::Vector{Float64}
    constraint_scales::Vector{Float64}
    constraint_rows::Vector{Int}
    jacobian_pattern::Vector{Tuple{Int,Int}}
    hessian_pattern::Vector{Tuple{Int,Int}}
    cached_y::Vector{Float64}
    cached_evaluation::Union{Nothing,SEEvaluation}
    cached_residual_jacobian::Union{Nothing,SparseMatrixCSC{Float64,Int}}
    cached_constraint_jacobian::Union{Nothing,SparseMatrixCSC{Float64,Int}}
    model_evaluations::Int
    residual_jacobian_evaluations::Int
    constraint_jacobian_evaluations::Int
end

function _SEIpoptEvaluator(s, p, origin, levels, scales, rows)
    H, C = _se_symbolic_jacobians(s, p)
    ci, cj, _ = findnz(C[rows, :])
    hi, hj, _ = findnz(tril(sparse(H' * H)))
    _SEIpoptEvaluator(
        s,
        p,
        origin,
        levels,
        scales,
        rows,
        collect(zip(ci, cj)),
        collect(zip(hi, hj)),
        fill(NaN, length(origin)),
        nothing,
        nothing,
        nothing,
        0,
        0,
        0,
    )
end

_SEMOI.features_available(::_SEIpoptEvaluator) = [:Grad, :Jac, :Hess]
function _SEMOI.initialize(::_SEIpoptEvaluator, features::Vector{Symbol})
    all(f->f in (:Grad, :Jac, :Hess), features) ||
        throw(ArgumentError("unsupported derivative feature"))
    nothing
end
_se_ipopt_state(d, y) = d.origin + d.levels .* y

function _se_ipopt_cache!(d, y)
    if !isequal(d.cached_y, y)
        copyto!(d.cached_y, y)
        d.model_evaluations+=1
        d.cached_evaluation = _se_trial_evaluation(d.structure, d.parameters, _se_ipopt_state(d, y))
        d.cached_residual_jacobian = nothing
        d.cached_constraint_jacobian = nothing
    end
    d.cached_evaluation
end
function _se_ipopt_H!(d, y)
    _se_ipopt_cache!(d, y) === nothing && return nothing
    if d.cached_residual_jacobian === nothing
        d.residual_jacobian_evaluations+=1
        d.cached_residual_jacobian = try
            residual_jacobian(d.structure, d.parameters, _se_ipopt_state(d, y))
        catch error
            error isa DomainError || rethrow()
            nothing
        end
    end
    d.cached_residual_jacobian
end
function _se_ipopt_C!(d, y)
    _se_ipopt_cache!(d, y) === nothing && return nothing
    if d.cached_constraint_jacobian === nothing
        d.constraint_jacobian_evaluations+=1
        d.cached_constraint_jacobian = try
            constraint_jacobian(d.structure, d.parameters, _se_ipopt_state(d, y))
        catch error
            error isa DomainError || rethrow()
            nothing
        end
    end
    d.cached_constraint_jacobian
end
function _SEMOI.eval_objective(d::_SEIpoptEvaluator, y)
    e = _se_ipopt_cache!(d, y)
    e === nothing ? NaN : 0.5sum(abs2, e.residual)
end
function _SEMOI.eval_objective_gradient(d::_SEIpoptEvaluator, g, y)
    H = _se_ipopt_H!(d, y)
    H === nothing ? fill!(g, NaN) : (g .= d.levels .* (H' * d.cached_evaluation.residual))
    nothing
end
function _SEMOI.eval_constraint(d::_SEIpoptEvaluator, c, y)
    e = _se_ipopt_cache!(d, y)
    e === nothing ? fill!(c, NaN) : (c .= e.constraints[d.constraint_rows] ./ d.constraint_scales)
    nothing
end
_SEMOI.jacobian_structure(d::_SEIpoptEvaluator) = d.jacobian_pattern
function _SEMOI.eval_constraint_jacobian(d::_SEIpoptEvaluator, J, y)
    C = _se_ipopt_C!(d, y)
    if C === nothing
        fill!(J, NaN)
    else
        for (k, (i, j)) in enumerate(d.jacobian_pattern)
            J[k] = C[d.constraint_rows[i], j] * d.levels[j] / d.constraint_scales[i]
        end
    end
    nothing
end
_SEMOI.hessian_lagrangian_structure(d::_SEIpoptEvaluator) = d.hessian_pattern
function _SEMOI.eval_hessian_lagrangian(
    d::_SEIpoptEvaluator,
    values,
    y,
    objective_factor,
    multipliers,
)
    H = _se_ipopt_H!(d, y)
    if H === nothing
        fill!(values, NaN)
    else
        scaled = H * spdiagm(0=>d.levels)
        curvature = scaled' * scaled
        for (k, (i, j)) in enumerate(d.hessian_pattern)
            values[k] = objective_factor * curvature[i, j]
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
    n=length(s.nodes)
    q=length(s.constraint_pattern)+size(s.return_current_pattern, 1)
    selector=spzeros(Float64, 2q, 2n)
    for (k, i) in enumerate(s.constraint_pattern)
        selector[2k-1, i]=1.0
        selector[2k, n+i]=1.0
    end
    offset=length(s.constraint_pattern)
    for (r, i, v) in zip(findnz(s.return_current_pattern)...)
        selector[2(offset+r)-1, i]=v
        selector[2(offset+r), n+i]=v
    end
    linear=deepcopy(p)
    fill!(linear.device_powers, 0im)
    for device in s.device_pattern
        device.kind===:constant_magnitude_current &&
            (linear.device_currents[device.parameter_range].=0im)
    end
    # The zero-power part is affine and defined at zero voltage. Alpha=0
    # avoids applying the physical nonzero-voltage domain to this decomposition.
    linear.continuation_alpha=0.0
    linear.regularization_voltage=1.0
    x=zeros(size(s.voltage_state_jacobian, 2))
    vr, vi, ir, ii=_se_parts(s, linear, x)
    dr, di, Jr, Ji=_device_parts(s, linear, vr, vi)
    operator=s.current_state_jacobian+vcat(Jr, Ji)
    if !isempty(s.injection_state_map)
        nv=2length(s.free_state_map)
        nu=size(s.injection_incidence, 2)
        operator[1:n, (nv+1):(nv+nu)]-=s.injection_incidence
        operator[(n+1):2n, (nv+nu+1):(nv+2nu)]-=s.injection_incidence
    end
    A=selector*operator
    b=selector*vcat(ir+dr, ii+di)
    nonlinear=Tuple{Int,Int}[]
    for device in s.device_pattern
        device.kind in (:constant_power, :zip, :constant_magnitude_current) || continue
        for (j, k) in enumerate(device.parameter_range)
            nonlinear_current=device.kind===:constant_magnitude_current ?
                              !iszero(p.device_currents[k]) : !iszero(p.device_powers[k])
            nonlinear_current && push!(nonlinear, (device.positive[j], device.negative[j]))
        end
    end
    B=spzeros(Float64, 2n, 2length(nonlinear))
    for (k, (i, j)) in enumerate(nonlinear)
        B[i, 2k-1]=1.0
        B[n+i, 2k]=1.0
        if j!=0
            B[j, 2k-1]=-1.0
            B[n+j, 2k]=-1.0
        end
    end
    hcat(A, sparse(reshape(b, :, 1)), selector*B)
end

function _se_ipopt_constraint_basis_compute(s, p)
    _se_constraint_count(s)==0 && return Int[], zeros(0, 0)
    signature=_se_ipopt_constraint_signature(s, p)
    m=size(signature, 1)
    m==0 && return Int[], zeros(0, 0)
    # Equilibrate units before pivoting. This uses a near-roundoff cutoff,
    # independent of the much looser user observability/rank tolerance.
    for j in axes(signature, 2)
        level=norm(signature[:, j])
        level>0 && (signature[:, j]./=level)
    end
    row_levels=ones(m)
    for i in axes(signature, 1)
        level=norm(signature[i, :])
        if level>0
            row_levels[i]=level
            signature[i, :]./=level
        end
    end
    factor=qr(sparse(transpose(signature)); tol = 8eps(Float64)*max(size(signature)...))
    diagonal=abs.(diag(factor.R))
    tolerance=8eps(Float64)*max(size(signature)...)
    rank=count(>(tolerance*maximum(diagonal; init = 1.0)), diagonal)
    rows=sort!(factor.pcol[1:rank])
    length(rows)==m && return rows, spdiagm(0=>ones(m))
    # Check the complete coefficient identity, not just the rank estimate.
    # If numerical cancellation prevents that check, retain the full system.
    expansion=if isempty(rows)
        zeros(m, 0)
    else
        selected=qr(sparse(transpose(signature[rows, :])); tol = tolerance)
        collect(transpose(selected \ collect(Float64, transpose(signature))))
    end
    reconstruction=expansion*signature[rows, :]
    for row in axes(signature, 1)
        norm(signature[row, :]-reconstruction[row, :], Inf)<=tolerance ||
            return collect(1:m), spdiagm(0=>ones(m))
    end
    # Selected rows are identities exactly. Undo row normalization so this
    # map also bounds physical-current amplification (e.g. turns ratios).
    expansion[rows, :].=Matrix{Float64}(I, length(rows), length(rows))
    expansion=expansion .* row_levels ./ transpose(row_levels[rows])
    rows, expansion
end
function _se_ipopt_constraint_basis(s, p)
    # Reading/sigma/prior updates leave this coefficient proof unchanged.
    # Physical source or device changes require a fresh proof, including b.
    key=(
        Tuple(p.fixed_voltages),
        Tuple(p.device_powers),
        Tuple(p.device_currents),
        Tuple(p.device_admittances),
    )
    lock(s.constraint_basis_lock) do
        if !haskey(s.constraint_basis_cache, key)
            empty!(s.constraint_basis_cache)
            s.constraint_basis_cache[key]=_se_ipopt_constraint_basis_compute(s, p)
        end
        s.constraint_basis_cache[key]
    end
end
_se_ipopt_constraint_rows(s, p) = first(_se_ipopt_constraint_basis(s, p))

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
    statistics::NamedTuple
end
function solve_status(r::IpoptStateEstimationResult)
    ok = r.status in (:converged_unique, :converged_underobserved)
    _result_solve_status(string(r.status), ok; primal_status = r.primal_status)
end
solve_diagnostics(r::IpoptStateEstimationResult) = merge(
    r.statistics,
    (
        constraint_rank = r.constraint_rank,
        tangent_dimension = r.tangent_dimension,
        observable_dimension = r.observable_dimension,
        objective = sum(abs2, r.evaluation.residual),
        optimizer_termination = r.termination_status,
        optimizer_primal = r.primal_status,
        physical_constraint_count = length(r.evaluation.constraints),
        optimizer_constraint_count = length(r.constraint_rows),
        optimizer_constraint_rows = copy(r.constraint_rows),
    ),
)

"""
    solve_ipopt_state_estimator(structure, parameters, x0; kwargs...)

Classic Ipopt nonlinear-programming solve over the **same compiled model** as
the dense and sparse constrained estimators. All measurement functions,
transformers, grounding, exact devices, priors and profiled customer powers
therefore have identical semantics. This does not invoke the legacy BMOPF OPF
wrapper or impose its complete P/Q coverage contract.

Uses analytic/automatic first derivatives, voltage/current scaling, compiled
sparse derivative patterns and Gauss–Newton curvature `H' * H`. Residual and
nonlinear constraint second derivatives are omitted. Evaluations and individual
Jacobians are cached at each iterate. Final rank checks remain dense; tested
feeder sizes do not establish general scaling.
State-independent dependencies are checked using passive/device affine
coefficients and independent nonlinear branch-current incidences. Only a
verified basis is sent to Ipopt; final checks retain all original equations.
An initially rank-deficient Jacobian alone never permits dropping an equation.
`solver_options` accepts Ipopt attribute pairs or a NamedTuple; the supplied
curvature approximation uses Ipopt's `hessian_approximation="exact"` callback mode. Final physical feasibility and reduced stationarity are
checked independently of the optimizer's status. Internal homotopy stages
remain nonpublishable.
"""
function _solve_ipopt_state_estimator(
    s::SEStructure,
    p::SEParameters,
    x0::AbstractVector{<:Real};
    max_iterations::Integer = 500,
    constraint_tolerance::Real = 1e-8,
    constraint_rtol::Real = 1e-11,
    optimality_tolerance::Real = 1e-6,
    rank_rtol::Real = sqrt(eps(Float64)),
    verbose::Bool = false,
    solver_options = (),
)
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    for v in (constraint_tolerance, constraint_rtol, optimality_tolerance, rank_rtol)
        isfinite(Float64(v)) && v >= 0 ||
            throw(ArgumentError("tolerances must be finite and nonnegative"))
    end
    options = collect(solver_options isa NamedTuple ? pairs(solver_options) : solver_options)
    for (k, v) in options
        string(k) == "hessian_approximation" &&
            v != "exact" &&
            throw(
                ArgumentError(
                    "shared Ipopt backend supplies Gauss-Newton curvature; hessian_approximation must be exact",
                ),
            )
        string(k) in ("dependency_detector", "dependency_detection_with_rhs") && throw(
            ArgumentError(
                "unverified constraint dropping is not supported; only state-independent coefficient identities may reduce equations",
            ),
        )
    end
    x = Float64.(x0)
    active_rows = Int[]
    e = _se_trial_evaluation(s, p, x)
    evaluator=nothing
    function result(
        status,
        state,
        evaluation,
        rc = 0,
        td = 0,
        od = 0,
        term = "NOT_RUN",
        primal = "NO_SOLUTION",
    )
        stats=evaluator===nothing ? NamedTuple() :
              (
            optimizer_jacobian_entries = length(evaluator.jacobian_pattern),
            optimizer_hessian_entries = length(evaluator.hessian_pattern),
            model_evaluations = evaluator.model_evaluations,
            residual_jacobian_evaluations = evaluator.residual_jacobian_evaluations,
            constraint_jacobian_evaluations = evaluator.constraint_jacobian_evaluations,
        )
        IpoptStateEstimationResult(
            status,
            state,
            evaluation,
            rc,
            td,
            od,
            term,
            primal,
            copy(active_rows),
            stats,
        )
    end
    e === nothing && return result(:invalid_initial_domain, x, _se_invalid_evaluation(s, p))
    jac = _se_jacobians(s, p, x)
    jac === nothing && return result(:undefined_derivative, x, e)
    ctol = _se_constraint_tolerance(constraint_tolerance, constraint_rtol, _se_current_scale(s, p))
    if isempty(x)
        status = norm(e.constraints) <= ctol ? :converged_unique : :constraint_not_satisfied
        if status === :converged_unique
            status = _se_checked_status(s, p, x, e, jac[1], jac[2], zeros(0, 0), status; rank_rtol)
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
    active_rows, expansion = _se_ipopt_constraint_basis(s, p)
    scales = [max(norm(jac[2][i, :] .* levels), 1.0) for i in active_rows]
    evaluator = _SEIpoptEvaluator(s, p, x, levels, scales, active_rows)
    model = Model(Ipopt.Optimizer)
    verbose || set_silent(model)
    @variable(model, y[1:length(x)], start=0.0)
    set_optimizer_attribute(model, "hessian_approximation", "exact")
    set_optimizer_attribute(model, "nlp_scaling_method", "none")
    set_optimizer_attribute(model, "max_iter", Int(max_iterations))
    set_optimizer_attribute(model, "tol", max(Float64(optimality_tolerance) / 10, 1e-12))
    # Ipopt checks scaled constraint components. Translate the physical 2-norm
    # budget to those components; its default feasibility tolerance can accept
    # a raw current error above our requested tolerance after row scaling.
    # If |c_active ./ scales| <= t componentwise, the complete physical
    # constraint norm is bounded by norm(abs(expansion)*scales)*t.
    amplification = norm(abs.(expansion)*scales)
    scaled_ctol = max(ctol / (2max(amplification, 1.0)), 1e-12)
    set_optimizer_attribute(model, "constr_viol_tol", scaled_ctol)
    set_optimizer_attribute(model, "acceptable_constr_viol_tol", scaled_ctol)
    _set_solver_options!(model, options)
    bounds = [_SEMOI.NLPBoundsPair(0.0, 0.0) for _ in scales]
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
    ef === nothing && return result(
        :invalid_final_domain,
        xf,
        _se_invalid_evaluation(s, p),
        0,
        0,
        0,
        term,
        primal,
    )
    jf = _se_jacobians(s, p, xf)
    jf === nothing && return result(:undefined_derivative, xf, ef, 0, 0, 0, term, primal)
    H, C = jf
    rc, Z = _se_rank_nullspace(C; rtol = rank_rtol)
    od, _ = _se_rank_nullspace(H * Z; rtol = rank_rtol, reference_scale = norm(H))
    td = size(Z, 2)
    status =
        norm(ef.constraints) > ctol ? :constraint_not_satisfied :
        _se_stationarity(H, Z' * (H' * ef.residual)) > optimality_tolerance ?
        :optimality_not_satisfied : od == td ? :converged_unique : :converged_underobserved
    if status in (:converged_unique, :converged_underobserved)
        status = _se_checked_status(s, p, xf, ef, H, C, Z, status; rank_rtol)
    end
    result(status, xf, ef, rc, td, od, term, primal)
end

solve_state_estimation(s::SEStructure, p::SEParameters, x0::AbstractVector{<:Real}; kwargs...) =
    solve_sparse_state_estimator(s, p, x0; kwargs...)

function solve_ipopt_state_estimator(
    s::SEStructure,
    p::SEParameters,
    x0::AbstractVector{<:Real};
    kwargs...,
)
    allowed=(
        :max_iterations,
        :constraint_tolerance,
        :constraint_rtol,
        :optimality_tolerance,
        :rank_rtol,
        :verbose,
        :solver_options,
    )
    invalid=setdiff(collect(keys(kwargs)), collect(allowed))
    isempty(invalid) || throw(
        ArgumentError(
            "Unsupported Ipopt options $invalid; trust-region options belong to dense/sparse solvers",
        ),
    )
    r=_solve_ipopt_state_estimator(s, p, x0; kwargs...)
    IpoptStateEstimationResult(
        r.status,
        r.state,
        r.evaluation,
        r.constraint_rank,
        r.tangent_dimension,
        r.observable_dimension,
        r.termination_status,
        r.primal_status,
        r.constraint_rows,
        merge(r.statistics, _se_statistics(s, p, r)),
    )
end

@doc (@doc _solve_ipopt_state_estimator) solve_ipopt_state_estimator
