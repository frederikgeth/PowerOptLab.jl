# Inverse Carson reconstruction from diagonal sequence data.
#
# Discrete construction candidates are enumerated outside the NLP. Each
# candidate has a small, smooth, scaled continuous problem for geometry and
# conductor temperature. The forward physics comes from BMOPFTools' public
# primitive overhead-line constants kernel; PowerOptLab owns sequence
# projection, fitting, multistart, identifiability, and ambiguity reporting.

const _IC_GEOMETRIES = (:horizontal_3, :triangle_3, :horizontal_4,
                        :neutral_under_4)

_ic_parameter_names(g::Symbol) =
    g == :horizontal_3 ? (:half_span, :height, :temperature) :
    g == :triangle_3 ? (:half_span, :height, :temperature) :
    g == :horizontal_4 ? (:inner_offset, :outer_gap, :height, :temperature) :
    g == :neutral_under_4 ? (:half_span, :phase_height, :neutral_drop,
                              :temperature) :
    throw(ArgumentError("unsupported overhead geometry :$g"))

_ic_conductor_count(g::Symbol) = g in (:horizontal_3, :triangle_3) ? 3 : 4

function _ic_validate_candidate_domain(geometry, radius, lo, hi, angle)
    if geometry in (:horizontal_3, :triangle_3)
        lo[1] > max(radius[1] + radius[2], radius[2] + radius[3]) ||
            throw(ArgumentError("minimum half-span permits overlapping conductors"))
        lo[2] > maximum(radius) ||
            throw(ArgumentError("minimum height must keep conductors above ground"))
        geometry == :triangle_3 && !(0 <= angle < pi / 2) &&
            throw(ArgumentError("triangle angle must lie in [0, π/2)"))
    elseif geometry == :horizontal_4
        lo[2] > max(radius[1] + radius[2], radius[3] + radius[4]) ||
            throw(ArgumentError("minimum outer gap permits overlapping conductors"))
        2lo[1] > radius[2] + radius[3] ||
            throw(ArgumentError("minimum inner offset permits overlapping conductors"))
        lo[3] > maximum(radius) ||
            throw(ArgumentError("minimum height must keep conductors above ground"))
    else
        lo[1] > max(radius[1] + radius[2], radius[2] + radius[3]) ||
            throw(ArgumentError("minimum half-span permits overlapping phase conductors"))
        lo[3] > radius[2] + radius[4] ||
            throw(ArgumentError("minimum neutral drop permits overlapping conductors"))
        lo[2] - hi[3] > radius[4] ||
            throw(ArgumentError("bounds permit the neutral to reach ground"))
    end
end

"""
    SequenceLineObservation(; z0, z1, frequency, ...)

Diagonal zero- and positive-sequence line data for [`solve_inverse_carson`](@ref).
`z0` and `z1` are complex impedances. Optional `b0` and `b1` are real shunt
susceptances. Values are converted to SI per metre internally.

# Keywords
- `z_units=:ohm_per_km` — `:ohm_per_m` or `:ohm_per_km`.
- `b_units=:micro_siemens_per_km` — `:siemens_per_m`,
  `:siemens_per_km`, or `:micro_siemens_per_km`.
- `sigma` — standard deviations in the same input units, ordered
  `(R0, X0, R1, X1[, B0, B1])`. If omitted, a descriptive 1% tolerance is
  used with an absolute floor in the declared input units; candidate scores
  then have no formal statistical interpretation.
- `covariance` — full positive-definite covariance matrix in the same ordered
  input units. It is mutually exclusive with `sigma`; correlations affect the
  Mahalanobis objective while `standardized_residual` remains marginal.
- `frequency` [Hz] and `earth_resistivity=100.0` [Ω·m].

The observation is assumed to come from a three-phase matrix, with any circuit
neutral Kron-reduced before the symmetrical-component transform.
"""
struct SequenceLineObservation
    z0::ComplexF64
    z1::ComplexF64
    b0::Union{Nothing,Float64}
    b1::Union{Nothing,Float64}
    sigma::Vector{Float64}
    covariance::Matrix{Float64}
    covariance_cholesky::Matrix{Float64}
    frequency::Float64
    earth_resistivity::Float64
end

function SequenceLineObservation(; z0, z1, b0=nothing, b1=nothing,
                                 frequency,
                                 earth_resistivity::Real=100.0,
                                 z_units::Symbol=:ohm_per_km,
                                 b_units::Symbol=:micro_siemens_per_km,
                                 sigma=nothing,
                                 covariance=nothing)
    zfactor = z_units == :ohm_per_m ? 1.0 :
              z_units == :ohm_per_km ? 1e-3 :
              throw(ArgumentError("z_units must be :ohm_per_m or :ohm_per_km"))
    bfactor = b_units == :siemens_per_m ? 1.0 :
              b_units == :siemens_per_km ? 1e-3 :
              b_units == :micro_siemens_per_km ? 1e-9 :
              throw(ArgumentError("unsupported b_units :$b_units"))
    (b0 === nothing) == (b1 === nothing) ||
        throw(ArgumentError("b0 and b1 must either both be supplied or both omitted"))
    frequency > 0 || throw(ArgumentError("frequency must be > 0 Hz"))
    earth_resistivity > 0 ||
        throw(ArgumentError("earth_resistivity must be > 0 Ω·m"))

    zin = Float64[real(z0), imag(z0), real(z1), imag(z1)]
    yin = b0 === nothing ? zin : [zin; Float64(b0); Float64(b1)]
    factors = b0 === nothing ? fill(zfactor, 4) :
              [fill(zfactor, 4); bfactor; bfactor]
    sigma !== nothing && covariance !== nothing &&
        throw(ArgumentError("supply either sigma or covariance, not both"))
    covariance_in = if covariance === nothing
        sigin = sigma === nothing ?
            [max(0.01 * abs(v), i <= 4 ? 1e-6 : 1e-3)
             for (i, v) in enumerate(yin)] : Float64.(collect(sigma))
        length(sigin) == length(yin) ||
            throw(DimensionMismatch("sigma must have $(length(yin)) entries"))
        all(>(0), sigin) || throw(ArgumentError("all sigma values must be > 0"))
        Diagonal(sigin .^ 2)
    else
        cov = Float64.(collect(covariance))
        size(cov) == (length(yin), length(yin)) ||
            throw(DimensionMismatch("covariance must be $(length(yin))×$(length(yin))"))
        all(isfinite, cov) || throw(ArgumentError("covariance must be finite"))
        isapprox(cov, transpose(cov); rtol=1e-10, atol=0.0) ||
            throw(ArgumentError("covariance must be symmetric"))
        cov
    end
    covariance_si = factors .* covariance_in .* transpose(factors)
    covariance_si = 0.5 .* (covariance_si .+ transpose(covariance_si))
    factor = try
        cholesky(Symmetric(covariance_si))
    catch err
        err isa PosDefException || rethrow()
        throw(ArgumentError("covariance must be positive definite"))
    end
    sigma_si = sqrt.(diag(covariance_si))

    SequenceLineObservation(ComplexF64(z0) * zfactor,
        ComplexF64(z1) * zfactor,
        b0 === nothing ? nothing : Float64(b0) * bfactor,
        b1 === nothing ? nothing : Float64(b1) * bfactor,
        sigma_si, covariance_si, collect(factor.L),
        Float64(frequency), Float64(earth_resistivity))
end

"""
    OverheadCarsonCandidate(; id, geometry, r_ac_ref, gmr, radius, lower, upper, ...)

One discrete overhead construction considered by [`solve_inverse_carson`](@ref).
Conductor arrays are ordered `(a,b,c[,n])` and use SI units: `r_ac_ref` [Ω/m],
`gmr`, `radius`, and `cap_radius` [m]. Separate phase and neutral conductors are
supported.

The continuous parameter order is returned in each fit and depends on geometry:

- `:horizontal_3`: `(half_span, height, temperature)`
- `:triangle_3`: `(half_span, height, temperature)`; `angle` is fixed
- `:horizontal_4`: `(inner_offset, outer_gap, height, temperature)`
- `:neutral_under_4`: `(half_span, phase_height, neutral_drop, temperature)`

Geometry values are metres and temperature is °C. `lower`, `upper`, and
`initial` follow that order. `alpha_20` may be a scalar or one value per
conductor; resistance is corrected from `temperature_ref` using IEC's linear
temperature relation.
"""
struct OverheadCarsonCandidate
    id::String
    geometry::Symbol
    r_ac_ref::Vector{Float64}
    gmr::Vector{Float64}
    radius::Vector{Float64}
    cap_radius::Vector{Float64}
    temperature_ref::Float64
    alpha_20::Vector{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
    initial::Vector{Float64}
    angle::Float64
    metadata::Dict{String,Any}
end

function OverheadCarsonCandidate(; id, geometry::Symbol, r_ac_ref, gmr, radius,
                                 lower, upper, initial=nothing,
                                 cap_radius=radius,
                                 temperature_ref::Real=20.0,
                                 alpha_20=0.004,
                                 angle::Real=21.67pi / 180,
                                 metadata=Dict{String,Any}())
    geometry in _IC_GEOMETRIES ||
        throw(ArgumentError("unsupported overhead geometry :$geometry"))
    n = _ic_conductor_count(geometry)
    arrays = Float64.(collect(r_ac_ref)), Float64.(collect(gmr)),
             Float64.(collect(radius)), Float64.(collect(cap_radius))
    all(length(v) == n for v in arrays) ||
        throw(DimensionMismatch("geometry :$geometry requires $n conductors"))
    all(>(0), arrays[1]) || throw(ArgumentError("r_ac_ref must be positive"))
    all(>(0), arrays[2]) || throw(ArgumentError("gmr must be positive"))
    all(>(0), arrays[3]) || throw(ArgumentError("radius must be positive"))
    all(>(0), arrays[4]) || throw(ArgumentError("cap_radius must be positive"))
    all(arrays[2] .<= arrays[3] .* (1 + 1e-9)) ||
        throw(ArgumentError("gmr must not exceed radius"))

    names = _ic_parameter_names(geometry)
    lo, hi = Float64.(collect(lower)), Float64.(collect(upper))
    length(lo) == length(names) == length(hi) ||
        throw(DimensionMismatch("lower and upper must follow parameter order $names"))
    all(lo .< hi) || throw(ArgumentError("every lower bound must be below upper"))
    init = initial === nothing ? (lo .+ hi) ./ 2 : Float64.(collect(initial))
    length(init) == length(names) ||
        throw(DimensionMismatch("initial must follow parameter order $names"))
    all((lo .<= init) .& (init .<= hi)) ||
        throw(ArgumentError("initial must lie within bounds"))
    _ic_validate_candidate_domain(geometry, arrays[3], lo, hi, Float64(angle))

    alpha = alpha_20 isa Real ? fill(Float64(alpha_20), n) :
            Float64.(collect(alpha_20))
    length(alpha) == n || throw(DimensionMismatch("alpha_20 must have $n entries"))
    all(>=(0), alpha) || throw(ArgumentError("alpha_20 must be nonnegative"))
    all(1 + alpha[i] * (Float64(temperature_ref) - 20) > 0 for i in 1:n) ||
        throw(ArgumentError("temperature_ref and alpha_20 imply nonpositive resistance"))
    all(1 + alpha[i] * (lo[end] - 20) > 0 for i in 1:n) ||
        throw(ArgumentError("lower temperature bound implies nonpositive resistance"))

    OverheadCarsonCandidate(string(id), geometry, arrays[1], arrays[2],
        arrays[3], arrays[4], Float64(temperature_ref), alpha, lo, hi, init,
        Float64(angle), Dict{String,Any}(metadata))
end

function _ic_normal_quantile(p::Real)
    0 < p < 1 || throw(ArgumentError("normal quantile requires 0 < p < 1"))
    # Acklam's rational approximation; absolute error is below 1.2e-9.
    a = (-3.969683028665376e1, 2.209460984245205e2,
         -2.759285104469687e2, 1.383577518672690e2,
         -3.066479806614716e1, 2.506628277459239)
    b = (-5.447609879822406e1, 1.615858368580409e2,
         -1.556989798598866e2, 6.680131188771972e1,
         -1.328068155288572e1)
    c = (-7.784894002430293e-3, -3.223964580411365e-1,
         -2.400758277161838, -2.549732539343734,
          4.374664141464968, 2.938163982698783)
    d = (7.784695709041462e-3, 3.224671290700398e-1,
         2.445134137142996, 3.754408661907416)
    plow = 0.02425
    if p < plow
        q = sqrt(-2log(p))
        return (((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
               ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1)
    elseif p > 1 - plow
        return -_ic_normal_quantile(1 - p)
    end
    q = p - 0.5
    r = q^2
    (((((a[1] * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * r + a[6]) * q /
        (((((b[1] * r + b[2]) * r + b[3]) * r + b[4]) * r + b[5]) * r + 1)
end

"Connected profile-likelihood interval for one inverse-Carson parameter."
struct InverseCarsonProfileInterval
    parameter::Symbol
    estimate::Float64
    lower::Float64
    upper::Float64
    lower_status::Symbol
    upper_status::Symbol
    confidence_level::Float64
    delta_objective::Float64
end

"One fitted candidate returned by [`solve_inverse_carson`](@ref)."
struct InverseCarsonFit
    candidate_id::String
    geometry::Symbol
    termination_status::String
    objective::Float64
    max_standardized_residual::Float64
    compatible::Bool
    parameter_names::Vector{Symbol}
    parameters::Vector{Float64}
    predicted::Vector{Float64}
    standardized_residual::Vector{Float64}
    Z_primitive::Matrix{ComplexF64}
    C_primitive::Matrix{Float64}
    Z_sequence::Matrix{ComplexF64}
    B_sequence::Union{Nothing,Matrix{ComplexF64}}
    jacobian_singular_values::Vector{Float64}
    jacobian_rank::Int
    local_parameter_covariance::Union{Nothing,Matrix{Float64}}
    local_confidence_intervals::Union{Nothing,Matrix{Float64}}
    confidence_level::Float64
    local_solutions::Vector{Vector{Float64}}
    frequency::Float64
    earth_resistivity::Float64
end

"Result of [`solve_inverse_carson`](@ref), retaining all candidate fits."
struct InverseCarsonResult
    fits::Vector{InverseCarsonFit}
    compatible_candidates::Vector{String}
    warnings::Vector{String}
end

_ic_observed(obs::SequenceLineObservation) = obs.b0 === nothing ?
    Float64[real(obs.z0), imag(obs.z0), real(obs.z1), imag(obs.z1)] :
    Float64[real(obs.z0), imag(obs.z0), real(obs.z1), imag(obs.z1),
            obs.b0, obs.b1]

function _ic_coordinates(c::OverheadCarsonCandidate, p)
    g = c.geometry
    if g == :horizontal_3
        u, h = p[1], p[2]
        return [-u, zero(u), u], fill(h, 3)
    elseif g == :triangle_3
        u, h = p[1], p[2]
        return [-u, zero(u), u], [h, h + u * tan(c.angle), h]
    elseif g == :horizontal_4
        inner, gap, h = p[1], p[2], p[3]
        outer = inner + gap
        return [-outer, -inner, inner, outer], fill(h, 4)
    end
    u, h, drop = p[1], p[2], p[3]
    [-u, zero(u), u, zero(u)], [h, h, h, h - drop]
end

function _ic_physical(c::OverheadCarsonCandidate, u)
    c.lower .+ (c.upper .- c.lower) .* u
end

function _ic_forward(c::OverheadCarsonCandidate, p,
                     frequency::Real, earth_resistivity::Real)
    x, y = _ic_coordinates(c, p)
    temperature = p[end]
    rac = [c.r_ac_ref[i] *
           (1 + c.alpha_20[i] * (temperature - 20)) /
           (1 + c.alpha_20[i] * (c.temperature_ref - 20))
           for i in eachindex(c.r_ac_ref)]
    BMOPFTools.overhead_line_constants(rac, c.gmr, c.radius, x, y;
        cap_radius=c.cap_radius, frequency=frequency,
        earth_model="modified_carson", earth_resistivity=earth_resistivity)
end

function _ic_phase_matrices(Z, C)
    if size(Z, 1) == 3
        return Z, C
    end
    Zabc = Z[1:3, 1:3] -
           Z[1:3, 4:4] * (Z[4:4, 4:4] \ Z[4:4, 1:3])
    # Grounding the eliminated circuit neutral fixes Vn=0; the phase shunt
    # relation is therefore the phase principal block of the nodal C matrix.
    Zabc, C[1:3, 1:3]
end

function _ic_sequence_matrices(Z, C, frequency)
    Zabc, Cabc = _ic_phase_matrices(Z, C)
    a = cis(2pi / 3)
    A = [1.0 1.0 1.0; 1.0 a^2 a; 1.0 a a^2]
    Ainv = inv(A)
    Z012 = Ainv * Zabc * A
    B012 = Ainv * ((2pi * frequency) .* Cabc) * A
    Z012, B012
end

function _ic_predict(c::OverheadCarsonCandidate, obs::SequenceLineObservation, u)
    p = _ic_physical(c, u)
    constants = _ic_forward(c, p, obs.frequency, obs.earth_resistivity)
    Z012, B012 = _ic_sequence_matrices(constants.Z, constants.C, obs.frequency)
    out = [real(Z012[1, 1]), imag(Z012[1, 1]),
           real(Z012[2, 2]), imag(Z012[2, 2])]
    if obs.b0 !== nothing
        append!(out, (real(B012[1, 1]), real(B012[2, 2])))
    end
    out
end

function _ic_residual(c, obs, u)
    LowerTriangular(obs.covariance_cholesky) \
        (_ic_predict(c, obs, u) .- _ic_observed(obs))
end

_ic_standardized_residual(c, obs, u) =
    (_ic_predict(c, obs, u) .- _ic_observed(obs)) ./ obs.sigma

_ic_objective(c, obs, u) = sum(abs2, _ic_residual(c, obs, u))

_ic_success_status(status) = status in
    (JuMP.MOI.OPTIMAL, JuMP.MOI.LOCALLY_SOLVED,
     JuMP.MOI.ALMOST_OPTIMAL, JuMP.MOI.ALMOST_LOCALLY_SOLVED)

_ic_usable_profile_solution(model, status) = _ic_success_status(status) ||
    (status == JuMP.MOI.SLOW_PROGRESS && JuMP.has_values(model) &&
     JuMP.primal_status(model) in
        (JuMP.MOI.FEASIBLE_POINT, JuMP.MOI.NEARLY_FEASIBLE_POINT))

function _ic_configure_model!(model, verbose::Bool, solver_options)
    verbose || JuMP.set_silent(model)
    # Ipopt relaxes variable bounds by default. That can expose the physical
    # forward kernel to u just outside the domain validated by the candidate
    # constructor (for example, slightly overlapping conductors). Keep its
    # trial points inside the exact JuMP bounds. Do not send these raw
    # Ipopt-specific attributes to another configurable optimizer.
    if JuMP.unsafe_backend(model) isa Ipopt.Optimizer
        JuMP.set_attribute(model, "hessian_approximation", "limited-memory")
        JuMP.set_attribute(model, "bound_relax_factor", 0.0)
    end
    options = solver_options isa NamedTuple ? pairs(solver_options) : solver_options
    for (name, value) in options
        JuMP.set_attribute(model, string(name), value)
    end
    nothing
end

function _ic_starts(c::OverheadCarsonCandidate, count::Int)
    count >= 1 || throw(ArgumentError("starts must be at least 1"))
    n = length(c.lower)
    initial = (c.initial .- c.lower) ./ (c.upper .- c.lower)
    points = Vector{Vector{Float64}}([initial])
    roots = sqrt.([2.0, 3.0, 5.0, 7.0, 11.0, 13.0][1:n])
    for k in 1:count-1
        push!(points, [0.05 + 0.90 * mod(0.5 + k * roots[j], 1.0) for j in 1:n])
    end
    points
end

function _ic_failed_fit(c, obs, status, confidence_level)
    InverseCarsonFit(c.id, c.geometry, status, Inf, Inf, false,
        collect(_ic_parameter_names(c.geometry)), Float64[], Float64[], Float64[],
        zeros(ComplexF64, 0, 0), zeros(0, 0), zeros(ComplexF64, 0, 0), nothing,
        Float64[], 0, nothing, nothing, confidence_level,
        Vector{Vector{Float64}}(), obs.frequency, obs.earth_resistivity)
end

function _ic_fit_candidate(c::OverheadCarsonCandidate,
                           obs::SequenceLineObservation;
                           starts::Int, acceptance_sigma::Float64,
                           rank_tolerance::Float64, confidence_level::Float64,
                           optimizer,
                           verbose::Bool, solver_options)
    n = length(c.lower)
    model = JuMP.Model(optimizer)
    _ic_configure_model!(model, verbose, solver_options)
    u = JuMP.@variable(model, 0 <= u[1:n] <= 1)
    objective(args::T...) where {T<:Real} =
        _ic_objective(c, obs, collect(args))
    op = JuMP.add_nonlinear_operator(model, n, objective;
                                     name=:inverse_carson_objective)
    JuMP.@objective(model, Min, op(u...))

    local_u = Vector{Vector{Float64}}()
    local_obj = Float64[]
    local_status = String[]
    statuses = String[]
    for start in _ic_starts(c, starts)
        JuMP.set_start_value.(u, start)
        JuMP.optimize!(model)
        status = JuMP.termination_status(model)
        push!(statuses, string(status))
        _ic_success_status(status) || continue
        sol = JuMP.value.(u)
        obj = _ic_objective(c, obs, sol)
        duplicate = findfirst(old -> norm(sol - old) <= 1e-6, local_u)
        if duplicate === nothing
            push!(local_u, sol); push!(local_obj, obj); push!(local_status, string(status))
        elseif obj < local_obj[duplicate]
            local_u[duplicate] = sol
            local_obj[duplicate] = obj
            local_status[duplicate] = string(status)
        end
    end
    isempty(local_u) && return _ic_failed_fit(c, obs,
        join(unique(statuses), ","), confidence_level)

    best = argmin(local_obj)
    ubest = local_u[best]
    pbest = _ic_physical(c, ubest)
    pred = Float64.(_ic_predict(c, obs, ubest))
    residual = _ic_standardized_residual(c, obs, ubest)
    constants = _ic_forward(c, pbest, obs.frequency, obs.earth_resistivity)
    Z012, B012 = _ic_sequence_matrices(constants.Z, constants.C, obs.frequency)
    J = ForwardDiff.jacobian(v -> _ic_residual(c, obs, v), ubest)
    jacobian_svd = svd(J)
    singular = jacobian_svd.S
    cutoff = isempty(singular) ? Inf : rank_tolerance * maximum(singular)
    numerical_rank = count(>(cutoff), singular)
    parameter_covariance = if numerical_rank == n
        scale = Diagonal(c.upper .- c.lower)
        covariance_u = jacobian_svd.V * Diagonal(1 ./ singular .^ 2) *
                       transpose(jacobian_svd.V)
        collect(scale * covariance_u * scale)
    else
        nothing
    end
    local_intervals = if parameter_covariance === nothing
        nothing
    else
        zcrit = _ic_normal_quantile((1 + confidence_level) / 2)
        se = sqrt.(max.(diag(parameter_covariance), 0.0))
        hcat(max.(c.lower, pbest .- zcrit .* se),
             min.(c.upper, pbest .+ zcrit .* se))
    end
    ordered_local = local_u[sortperm(local_obj)]

    InverseCarsonFit(c.id, c.geometry, local_status[best], local_obj[best],
        maximum(abs, residual), maximum(abs, residual) <= acceptance_sigma,
        collect(_ic_parameter_names(c.geometry)), Float64.(pbest), pred,
        Float64.(residual), ComplexF64.(constants.Z), Float64.(constants.C),
        ComplexF64.(Z012), obs.b0 === nothing ? nothing : ComplexF64.(B012),
        Float64.(singular), numerical_rank, parameter_covariance,
        local_intervals, confidence_level,
        [Float64.(_ic_physical(c, v)) for v in ordered_local],
        obs.frequency, obs.earth_resistivity)
end

"""
    solve_inverse_carson(observation, candidates; kwargs...) -> InverseCarsonResult

Fit every discrete overhead construction candidate to diagonal sequence data.
Each candidate is solved independently as a bound-constrained smooth NLP with
deterministic multistart. Fits are returned in increasing weighted-residual
order; ambiguity is retained rather than collapsed to a single winner.

# Keywords
- `starts=16` — deterministic starts per candidate.
- `acceptance_sigma=3.0` — a candidate is compatible only when every
  standardized residual is within this threshold.
- `rank_tolerance=1e-6` — relative threshold for the singular values of the
  standardized prediction Jacobian.
- `confidence_level=0.95` — level for the local linearized parameter intervals.
- `optimizer=Ipopt.Optimizer`, `verbose=false`; `solver_options` accepts an
  iterable of name-value pairs or a named tuple. Ipopt receives safe defaults
  `hessian_approximation="limited-memory"` and `bound_relax_factor=0.0` before
  user options are applied. Other optimizers receive no Ipopt-specific options.

The solver uses modified Carson only. Returned `Z_primitive` and `C_primitive`
retain an explicit circuit neutral; Kron reduction is used solely to compare
against the sequence observation.
"""
function solve_inverse_carson(obs::SequenceLineObservation,
                              candidates::AbstractVector{<:OverheadCarsonCandidate};
                              starts::Int=16,
                              acceptance_sigma::Real=3.0,
                              rank_tolerance::Real=1e-6,
                              confidence_level::Real=0.95,
                              optimizer=Ipopt.Optimizer,
                              verbose::Bool=false,
                              solver_options=())
    isempty(candidates) && throw(ArgumentError("at least one candidate is required"))
    length(unique(c.id for c in candidates)) == length(candidates) ||
        throw(ArgumentError("candidate ids must be unique"))
    acceptance_sigma > 0 || throw(ArgumentError("acceptance_sigma must be > 0"))
    rank_tolerance > 0 || throw(ArgumentError("rank_tolerance must be > 0"))
    0 < confidence_level < 1 ||
        throw(ArgumentError("confidence_level must lie strictly between 0 and 1"))

    fits = [_ic_fit_candidate(c, obs; starts=starts,
                acceptance_sigma=Float64(acceptance_sigma),
                rank_tolerance=Float64(rank_tolerance), optimizer=optimizer,
                confidence_level=Float64(confidence_level),
                verbose=verbose, solver_options=solver_options)
            for c in candidates]
    sort!(fits, by=f -> f.objective)
    compatible = [f.candidate_id for f in fits if f.compatible]
    warnings = String[]
    isempty(compatible) && push!(warnings,
        "no candidate fits within the stated measurement uncertainty")
    length(compatible) > 1 && push!(warnings,
        "multiple candidates are compatible; construction is ambiguous")
    for fit in fits
        fit.compatible && fit.jacobian_rank < length(fit.parameters) &&
            push!(warnings,
                "candidate $(fit.candidate_id) is locally rank-deficient")
    end
    InverseCarsonResult(fits, compatible, unique(warnings))
end

function _ic_profile_endpoint(c, obs, best_u, index, target, threshold;
                              points, bisection_steps, optimizer, verbose,
                              solver_options)
    n = length(best_u)
    model = JuMP.Model(optimizer)
    _ic_configure_model!(model, verbose, solver_options)
    u = JuMP.@variable(model, 0 <= u[1:n] <= 1)
    objective(args::T...) where {T<:Real} =
        _ic_objective(c, obs, collect(args))
    op = JuMP.add_nonlinear_operator(model, n, objective;
                                     name=:inverse_carson_profile_objective)
    JuMP.@objective(model, Min, op(u...))

    function evaluate(fixed_value, start)
        JuMP.fix(u[index], fixed_value; force=true)
        JuMP.set_start_value.(u, start)
        JuMP.set_start_value(u[index], fixed_value)
        JuMP.optimize!(model)
        status = JuMP.termination_status(model)
        _ic_usable_profile_solution(model, status) || return (Inf, start, false)
        solution = JuMP.value.(u)
        (_ic_objective(c, obs, solution), solution, true)
    end

    inside = best_u[index]
    inside_solution = copy(best_u)
    # Log spacing resolves both very tight metrology-driven intervals and
    # weakly identified profiles that extend to the candidate bounds.
    fractions = exp.(range(log(1e-8), 0.0; length=points))
    for fraction in fractions
        fixed_value = best_u[index] + fraction * (target - best_u[index])
        value, solution, ok = evaluate(fixed_value, inside_solution)
        if !ok || value > threshold
            outside = fixed_value
            crossed = ok
            for _ in 1:bisection_steps
                midpoint = (inside + outside) / 2
                midvalue, midsolution, midok = evaluate(midpoint, inside_solution)
                if !midok
                    outside = midpoint
                elseif midvalue <= threshold
                    inside = midpoint
                    inside_solution = midsolution
                else
                    outside = midpoint
                    crossed = true
                end
            end
            return crossed ? (inside, :threshold) : (target, :failed)
        end
        inside = fixed_value
        inside_solution = solution
    end
    (target, :bound)
end

"""
    profile_inverse_carson(fit, candidate, observation; kwargs...)

Compute connected one-parameter profile-likelihood confidence intervals around
an inverse-Carson fit. Each parameter is fixed successively while all remaining
parameters are reoptimized. An endpoint status of `:threshold` means the
chi-square threshold was crossed, `:bound` means the candidate bound was reached,
and `:failed` conservatively returns the bound because the profile solve failed.

The default threshold is the one-degree-of-freedom chi-square quantile implied
by `fit.confidence_level`. `points=12` traces each side before bisection;
`bisection_steps=16` refines the first crossing. These are local connected
profiles, not a guarantee that disconnected feasible regions do not exist.
"""
function profile_inverse_carson(fit::InverseCarsonFit,
                                c::OverheadCarsonCandidate,
                                obs::SequenceLineObservation;
                                confidence_level::Real=fit.confidence_level,
                                points::Int=12,
                                bisection_steps::Int=16,
                                optimizer=Ipopt.Optimizer,
                                verbose::Bool=false,
                                solver_options=())
    fit.candidate_id == c.id ||
        throw(ArgumentError("fit and candidate ids do not match"))
    isempty(fit.parameters) && throw(ArgumentError("cannot profile a failed fit"))
    isapprox(fit.frequency, obs.frequency) ||
        throw(ArgumentError("fit and observation frequencies do not match"))
    isapprox(fit.earth_resistivity, obs.earth_resistivity) ||
        throw(ArgumentError("fit and observation earth resistivities do not match"))
    0 < confidence_level < 1 ||
        throw(ArgumentError("confidence_level must lie strictly between 0 and 1"))
    points >= 2 || throw(ArgumentError("points must be at least 2"))
    bisection_steps >= 1 ||
        throw(ArgumentError("bisection_steps must be at least 1"))

    best_u = (fit.parameters .- c.lower) ./ (c.upper .- c.lower)
    zcrit = _ic_normal_quantile((1 + confidence_level) / 2)
    delta = zcrit^2
    threshold = fit.objective + delta
    intervals = InverseCarsonProfileInterval[]
    for (i, name) in enumerate(fit.parameter_names)
        lower_u, lower_status = _ic_profile_endpoint(c, obs, best_u, i, 0.0,
            threshold; points=points, bisection_steps=bisection_steps,
            optimizer=optimizer, verbose=verbose, solver_options=solver_options)
        upper_u, upper_status = _ic_profile_endpoint(c, obs, best_u, i, 1.0,
            threshold; points=points, bisection_steps=bisection_steps,
            optimizer=optimizer, verbose=verbose, solver_options=solver_options)
        scale = c.upper[i] - c.lower[i]
        push!(intervals, InverseCarsonProfileInterval(name, fit.parameters[i],
            c.lower[i] + scale * lower_u, c.lower[i] + scale * upper_u,
            lower_status, upper_status, Float64(confidence_level), delta))
    end
    intervals
end

"""
    materialize_inverse_carson(fit, candidate) -> Dict

Create BMOPF-ready `wire_data` and `line_geometry` blocks for a successful
inverse fit. The returned dictionary is not inserted into a network and no
linecode is compiled automatically.
"""
function materialize_inverse_carson(fit::InverseCarsonFit,
                                    c::OverheadCarsonCandidate)
    fit.candidate_id == c.id ||
        throw(ArgumentError("fit and candidate ids do not match"))
    isempty(fit.parameters) && throw(ArgumentError("cannot materialize a failed fit"))
    x, y = _ic_coordinates(c, fit.parameters)
    terminals = size(fit.Z_primitive, 1) == 3 ? ["a", "b", "c"] :
                ["a", "b", "c", "n"]
    wire_data = Dict{String,Any}()
    conductors = Any[]
    for i in eachindex(terminals)
        wid = "$(c.id)_wire_$i"
        wire_data[wid] = Dict{String,Any}(
            "kind" => "overhead", "r_ac" => c.r_ac_ref[i],
            "gmr" => c.gmr[i], "radius" => c.radius[i],
            "cap_radius" => c.cap_radius[i],
            "temperature_ref" => c.temperature_ref,
            "alpha_20" => c.alpha_20[i])
        push!(conductors, Dict{String,Any}(
            "wire_data" => wid, "x" => Float64(x[i]), "y" => Float64(y[i]),
            "terminal" => terminals[i]))
    end
    Dict{String,Any}(
        "wire_data" => wire_data,
        "line_geometry" => Dict{String,Any}(c.id => Dict{String,Any}(
            "frequency" => fit.frequency,
            "earth_model" => "modified_carson",
            "earth_resistivity" => fit.earth_resistivity,
            "temperature" => fit.parameters[end],
            "conductors" => conductors)))
end
