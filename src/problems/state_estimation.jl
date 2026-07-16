# Weighted-least-squares (WLS) state estimation.
#
# A DIFFERENT problem specification over the same network physics: given noisy
# measurements of an energised network, find the bus voltage state that best fits
# them. It reuses the BMOPFTools device model but with (a) no operational bounds,
# (b) free injection currents at measured buses instead of fixed loads, and (c) a
# measurement-residual objective instead of generation cost — all via
# `build_opf_model(add_objective=false, model_hook! = …)`.
#
# The estimation network is a CONTRACT, not a convention. A net carrying loads,
# generators, IBRs, or operational limits does not describe an estimation problem
# — those devices impose fixed injections and bounds that silently bias or
# constrain the estimate. `solve_state_estimation` therefore rejects such a net
# up front (`allow_operational=true` overrides with a warning). Buses are either
# measured (a `:pinj`+`:qinj` pair gives them a free injection) or must be
# declared `zero_injection`; absence of telemetry is NOT read as zero injection.
#
# Assumptions (documented, not all enforceable):
#   * independent, zero-mean measurement errors with diagonal covariance
#     (weight `1/σ²`); the estimator is a maximum-likelihood fit under this model;
#   * network topology and line/transformer parameters are exact;
#   * the voltage source magnitude/angle is a hard boundary (the angle reference);
#   * missing injections are zero ONLY where explicitly declared;
#   * the fit is nonconvex — Ipopt returns a LOCAL stationary point, not a
#     certified global/unique one (see the observability diagnostic below).

using LinearAlgebra: svdvals

# Bus operational-limit fields that must not appear in an estimation net.
const _SE_BUS_LIMIT_FIELDS = ("v_min", "v_max", "vpn_min", "vpn_max",
                              "vpp_min", "vpp_max", "vpos_min", "vpos_max",
                              "v_dc_min", "v_dc_max")
# Branch operational-limit fields.
const _SE_BRANCH_LIMIT_FIELDS = ("i_max", "s_max", "s_rating")
# Device blocks that inject current / impose a dispatch (not passive physics).
const _SE_INJECTING_BLOCKS = ("load", "generator", "ibr", "dc_source")

"""
    Measurement(; kind, bus, value, sigma, terminal="1", reference=missing)

A single scalar measurement for [`solve_state_estimation`](@ref). `value` and
`sigma` are SI (volts for `:vmag`, watts for `:pinj`, vars for `:qinj`); `sigma`
is the measurement standard deviation (WLS weight `1/sigma²`), and must be
finite and strictly positive.

- `kind::Symbol` — `:vr`, `:vi`, or `:vmag` (a rectangular component or the
  magnitude of the voltage across `(bus, terminal)`→`reference`), `:pinj`
  (active power injected into the network at that terminal pair), or `:qinj`
  (reactive power injection).
- `bus::String`, `terminal::String="1"` — the measured phase conductor.
- `reference` — the return terminal the quantity is referenced to. `missing`
  (the default) inherits the solve's `neutral`; a `String` names an explicit
  return terminal on the same bus; `nothing` references terminal-to-ground.
  Voltage and power measurements at a bus therefore share one reference — a
  smart-meter reading is phase-to-neutral for both.

Construction validates `kind`, finiteness, `sigma > 0`, and non-empty
identifiers; it does not check that the identifiers exist in a particular net.
"""
struct Measurement <: AbstractMeasurement
    kind::Symbol
    bus::String
    value::Float64
    sigma::Float64
    terminal::String
    reference::Union{String,Nothing,Missing}

    function Measurement(; kind::Symbol, bus::String, value::Real, sigma::Real,
                         terminal::String="1",
                         reference::Union{String,Nothing,Missing}=missing)
        kind in (:vr, :vi, :vmag, :pinj, :qinj) ||
            throw(ArgumentError("unknown measurement kind :$(kind); expected :vr, :vi, :vmag, :pinj, or :qinj"))
        _validate_measurement_scalar(kind, value, sigma)
        isempty(bus) && throw(ArgumentError("measurement bus must be non-empty"))
        isempty(terminal) && throw(ArgumentError("measurement terminal must be non-empty"))
        reference isa String && isempty(reference) &&
            throw(ArgumentError("measurement reference terminal must be non-empty"))
        new(kind, bus, value, Float64(sigma), terminal, reference)
    end
end

# Resolved return terminal for a measurement given the solve-level neutral:
# `missing` inherits the neutral, an explicit value is used as-is.
_resolve_ref(m::Measurement, neutral) = m.reference === missing ? neutral : m.reference

"""
    StateEstimationResult

Result of [`solve_state_estimation`](@ref).

# Fields
- `termination_status::String`, `primal_status::String` — the solver's
  termination status and primal-point status. Trust the estimate only when
  `primal_status == "FEASIBLE_POINT"`.
- `objective::Float64` — the optimal weighted-residual sum `∑ (z−h)²/σ²`
  (`NaN` if no feasible point was found).
- `bus::Dict{String,Any}` — the estimated SI bus voltages (`vr`, `vi`, `vm`,
  `va` per terminal). **`NaN` throughout when `primal_status` is not
  `FEASIBLE_POINT`** — an unconverged solver iterate is not an estimate.
- `residuals::Vector{NamedTuple}` — per input measurement, in order:
  `(kind, bus, terminal, reference, measured, estimated, residual, standardized)`
  with `residual = measured − estimated` and `standardized = residual/σ`.
  `standardized` is the σ-normalised RAW residual, **not** the classical
  leverage-adjusted normalised residual (`rᴺ = r / √(Sᵢᵢ)`) used for bad-data
  identification; it is a scale-free residual, not a χ²/rᴺ test statistic.
- `observability::NamedTuple` — a LOCAL numerical identifiability diagnostic
  `(observable, n_states, rank, redundancy, min_singular, cond)` from the
  measurement Jacobian at the returned point (see [`solve_state_estimation`](@ref)).
"""
struct StateEstimationResult <: AbstractSolveResult
    termination_status::String
    primal_status::String
    objective::Float64
    bus::Dict{String,Any}
    residuals::Vector{NamedTuple}
    observability::NamedTuple
    solve::SolveStatus
end

solve_status(result::StateEstimationResult) = result.solve

solve_diagnostics(result::StateEstimationResult) =
    (objective=result.objective, observability=result.observability,
     residual_count=length(result.residuals))

_vscale(ctx, bus) = ctx.bases === nothing ? 1.0 : ctx.bases.v_base[bus]

# ── network contract ────────────────────────────────────────────────────────

# Reject (or, if allowed, warn about) any content that makes `net` an operational
# model rather than a pure estimation physics model: injecting devices and
# operational limits. Returns nothing; throws unless `allow_operational`.
function _reject_operational(net::Dict{String,Any}, allow_operational::Bool)
    offenders = String[]
    for blk in _SE_INJECTING_BLOCKS
        d = get(net, blk, nothing)
        if d isa AbstractDict && !isempty(d)
            push!(offenders, "$(length(d)) $(blk) device(s) impose fixed injections/dispatch")
        end
    end
    for (bid, bus) in get(net, "bus", Dict())
        bus isa AbstractDict || continue
        for f in _SE_BUS_LIMIT_FIELDS
            haskey(bus, f) && push!(offenders, "bus '$bid' carries operational limit '$f'")
        end
    end
    for blk in ("line", "switch", "linecode", "transformer", "dc_branch")
        d = get(net, blk, nothing)
        d isa AbstractDict || continue
        for (id, el) in d
            el isa AbstractDict || continue
            for f in _SE_BRANCH_LIMIT_FIELDS
                haskey(el, f) && push!(offenders, "$blk '$id' carries operational limit '$f'")
            end
        end
    end
    isempty(offenders) && return nothing
    msg = "the network is not a pure estimation model; it carries operational " *
          "content that would bias or constrain the estimate:\n  - " *
          join(offenders, "\n  - ") *
          "\nRemove these (loads/generators/IBRs become measured or zero-injection " *
          "buses; strip voltage/thermal limits), or pass `allow_operational=true` " *
          "to estimate against this model anyway."
    allow_operational || throw(ArgumentError(msg))
    @warn "solve_state_estimation: estimating against an operational network " *
          "(allow_operational=true); the estimate reflects its devices/limits.\n" *
          join(offenders, "\n")
    return nothing
end

# Per-bus nominal |V| (SI volts) for seeding, covering every bus: a :vmag reading
# at the bus if one exists, else a representative source magnitude fallback.
function _nominal_v(net, measurements)
    src_mags = Float64[]
    for (_, vs) in get(net, "voltage_source", Dict())
        append!(src_mags, Float64.(get(vs, "v_magnitude", Float64[])))
    end
    fallback = !isempty(src_mags) ? maximum(src_mags) :
               (any(m.kind == :vmag for m in measurements) ?
                maximum(m.value for m in measurements if m.kind == :vmag) : 1.0)
    vnom = Dict{String,Float64}(b => fallback for b in keys(get(net, "bus", Dict())))
    for m in measurements
        m.kind == :vmag && (vnom[m.bus] = m.value)
    end
    vnom
end

# ── injection coverage / P–Q pairing ────────────────────────────────────────

# Normalise a `zero_injection` argument (bus ids and/or (bus, terminal) tuples)
# into a set of (bus, terminal) pairs, expanding a bare bus id over its phase
# terminals (all terminal_names except `neutral` and grounded terminals).
function _zero_injection_set(net, zero_injection, neutral)
    zi = Set{Tuple{String,String}}()
    buses = get(net, "bus", Dict())
    for z in zero_injection
        if z isa Tuple{String,String} || (z isa Tuple && length(z) == 2)
            push!(zi, (String(z[1]), String(z[2])))
        elseif z isa AbstractString
            b = String(z)
            bus = get(buses, b, nothing)
            bus === nothing && throw(ArgumentError("zero_injection bus '$b' not in net"))
            for t in _phase_terminals(bus, neutral)
                push!(zi, (b, t))
            end
        else
            throw(ArgumentError("zero_injection entries must be a bus id or a (bus, terminal) tuple; got $(z)"))
        end
    end
    zi
end

# Phase terminals of a bus: declared terminals minus the neutral and any
# perfectly-grounded terminals.
function _phase_terminals(bus::AbstractDict, neutral)
    terms = String.(get(bus, "terminal_names", String[]))
    grounded = Set(String.(get(bus, "perfectly_grounded_terminals", String[])))
    [t for t in terms if t != neutral && !(t in grounded)]
end

# Buses that host a voltage source (their fixed terminals are the boundary, not
# unknowns to be covered by measurements).
function _source_buses(net)
    Set(String(get(vs, "bus", "")) for (_, vs) in get(net, "voltage_source", Dict()))
end

# Enforce the injection-coverage contract: every non-source phase terminal is
# either a measured injection (both :pinj and :qinj present) or declared
# zero-injection, and the two are mutually exclusive. Throws on any violation.
function _check_coverage(net, measurements, neutral, zi::Set{Tuple{String,String}})
    has_p = Set{Tuple{String,String}}()
    has_q = Set{Tuple{String,String}}()
    for m in measurements
        m.kind == :pinj && push!(has_p, (m.bus, m.terminal))
        m.kind == :qinj && push!(has_q, (m.bus, m.terminal))
    end

    problems = String[]
    # P–Q pairing: an injection is a complex quantity; a lone P or Q is ill-posed.
    for k in union(has_p, has_q)
        k in has_p || push!(problems, "bus '$(k[1])' terminal '$(k[2])' has a :qinj but no :pinj")
        k in has_q || push!(problems, "bus '$(k[1])' terminal '$(k[2])' has a :pinj but no :qinj")
    end
    # A declared zero-injection terminal must not also be measured.
    for k in intersect(zi, union(has_p, has_q))
        push!(problems, "bus '$(k[1])' terminal '$(k[2])' is declared zero_injection but also has an injection measurement")
    end
    # Coverage: no silent zero-injection assumption for un-declared, un-measured buses.
    srcs = _source_buses(net)
    measured = intersect(has_p, has_q)
    for (bid, bus) in get(net, "bus", Dict())
        bid in srcs && continue
        bus isa AbstractDict || continue
        for t in _phase_terminals(bus, neutral)
            k = (bid, t)
            (k in measured || k in zi) && continue
            push!(problems, "bus '$bid' terminal '$t' has no injection measurement and is not declared zero_injection " *
                            "(absence of telemetry is not evidence of zero injection)")
        end
    end
    isempty(problems) && return nothing
    throw(ArgumentError("injection specification is incomplete:\n  - " * join(problems, "\n  - ") *
                        "\nDeclare zero-injection buses via `zero_injection=[...]` and pair every :pinj with a :qinj."))
end

"""
    solve_state_estimation(net, measurements; kwargs...) -> StateEstimationResult

Estimate the network state of `net` — a *physics-only* BMOPFTools net (buses,
lines, transformers, shunts, and a voltage source; **no** loads, generators,
IBRs, or operational limits) — that best fits `measurements` in a
weighted-least-squares sense.

The network is treated as a contract and validated up front (see
`allow_operational`). Every non-source phase terminal must be either a measured
injection (a `:pinj`+`:qinj` pair) or declared in `zero_injection`; an
un-declared, un-measured bus is an error, never a silent zero injection.

# Keywords
- `neutral="n"` — default return terminal for injection measurements, free
  injections, and voltage references. Pass `nothing` if phase terminals are
  referenced directly to ground.
- `zero_injection=String[]` — buses (or `(bus, terminal)` pairs) known to carry
  no injection. A bare bus id expands over its phase terminals.
- `allow_operational=false` — when `true`, downgrade the network-contract check
  to a warning instead of an error (estimate against loads/limits deliberately).
- `check_observability=true` — compute a local identifiability diagnostic.
- `per_unit=true`, `s_base=1e6` — engine unit handling; measurements stay SI.
- `optimizer=Ipopt.Optimizer`, `verbose=false`, `solver_options=()`.

# Observability
The returned `observability` NamedTuple reports a LOCAL numerical check: the
Jacobian of the measurement + zero-injection equations with respect to the
rectangular node voltages is formed at the returned point (reusing
`ybus_passive`), and `observable = rank == n_states`. It reports `redundancy`
(surplus equations), the smallest singular value, and the condition number. This
detects local rank deficiency / critically-weak measurement sets; it is not a
global uniqueness proof, and solver convergence alone never establishes one.

# Returns
A [`StateEstimationResult`](@ref); its `bus` voltages are `NaN` unless
`primal_status == "FEASIBLE_POINT"`.
"""
function solve_state_estimation(net::Dict{String,Any}, measurements::AbstractVector;
                                neutral::Union{String,Nothing}="n",
                                zero_injection=String[],
                                allow_operational::Bool=false,
                                check_observability::Bool=true,
                                per_unit::Bool=true,
                                s_base::Float64=1e6,
                                optimizer=Ipopt.Optimizer,
                                verbose::Bool=false,
                                solver_options=())
    isempty(measurements) && throw(ArgumentError("no measurements supplied"))
    all(m -> m isa Measurement, measurements) ||
        throw(ArgumentError("measurements must be a collection of `Measurement`"))
    all(m -> m.kind in (:vmag, :pinj, :qinj), measurements) ||
        throw(ArgumentError("solve_state_estimation currently supports :vmag, :pinj, and :qinj; " *
                            "use the compiled constrained-NLLS evaluator for :vr/:vi measurements"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError("s_base must be finite and > 0"))

    _reject_operational(net, allow_operational)
    zi = _zero_injection_set(net, zero_injection, neutral)
    _check_coverage(net, measurements, neutral, zi)

    # Seed injection starts from the measured P/Q at each (bus, terminal).
    pmeas = Dict{Tuple{String,String},Float64}()
    qmeas = Dict{Tuple{String,String},Float64}()
    for m in measurements
        m.kind == :pinj && (pmeas[(m.bus, m.terminal)] = m.value)
        m.kind == :qinj && (qmeas[(m.bus, m.terminal)] = m.value)
    end
    vnom = _nominal_v(net, measurements)   # per-bus nominal |V| (SI volts)

    # (measurement, SI-valued h-expression) pairs, filled by the hook for residuals.
    probes = Vector{Tuple{Measurement,Any}}()

    function wls!(ctx)
        m = ctx.model
        vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
        sb = _sbase(ctx)

        # A free injection current at each measured (bus, phase), added to KCL so
        # those buses' voltages stay free to fit. Referenced to the measurement's
        # own return terminal. Seeded from the measured (P, Q).
        inj_r = Dict{Tuple{String,String},Any}()
        inj_i = Dict{Tuple{String,String},Any}()
        for meas in measurements
            (meas.kind in (:pinj, :qinj)) || continue
            key = (meas.bus, meas.terminal)
            haskey(inj_r, key) && continue
            ref = _resolve_ref(meas, neutral)
            # Start: I ≈ conj(S)/conj(V). In MODEL units V_model ≈ v_nom_SI/_vscale
            # (≈1 in per-unit, ≈v_nom in SI), so this seeds consistently in both.
            p0 = get(pmeas, key, 0.0); q0 = get(qmeas, key, 0.0)
            vmod = max(get(vnom, meas.bus, 1.0) / _vscale(ctx, meas.bus), 1e-9)
            cr0 = (p0 / sb) / vmod; ci0 = -(q0 / sb) / vmod
            cr = JuMP.@variable(m, base_name = "seinj_r_$(meas.bus)_$(meas.terminal)", start = cr0)
            ci = JuMP.@variable(m, base_name = "seinj_i_$(meas.bus)_$(meas.terminal)", start = ci0)
            inj_r[key] = cr; inj_i[key] = ci
            JuMP.add_to_expression!(ctx.kcl_r[(meas.bus, meas.terminal)], cr)
            JuMP.add_to_expression!(ctx.kcl_i[(meas.bus, meas.terminal)], ci)
            if ref !== nothing
                JuMP.add_to_expression!(ctx.kcl_r[(meas.bus, ref)], -cr)
                JuMP.add_to_expression!(ctx.kcl_i[(meas.bus, ref)], -ci)
            end
        end

        obj = zero(JuMP.QuadExpr)
        for meas in measurements
            w = 1.0 / meas.sigma^2
            b = meas.bus; t = meas.terminal
            ref = _resolve_ref(meas, neutral)
            # Phase-to-reference voltage drop — the SAME reference for |V| and P/Q.
            if ref === nothing
                dvr = vr[(b,t)]; dvi = vi[(b,t)]
            else
                dvr = JuMP.@expression(m, vr[(b,t)] - vr[(b,ref)])
                dvi = JuMP.@expression(m, vi[(b,t)] - vi[(b,ref)])
            end
            if meas.kind == :vmag
                vb = _vscale(ctx, b)
                # Auxiliary |V| ≥ 0 with |V|² = dvr²+dvi², seeded from the reading
                # (a zero start sits on the degenerate |V|²=… Jacobian).
                vm = JuMP.@variable(m, base_name = "sevm_$(b)_$(t)", lower_bound = 0.0,
                                    start = meas.value / vb)
                JuMP.@constraint(m, vm^2 == dvr^2 + dvi^2)
                h_si = JuMP.@expression(m, vm * vb)                        # → volts
                obj += w * (h_si - meas.value)^2
                push!(probes, (meas, h_si))
            else  # :pinj / :qinj
                cr = inj_r[(b,t)]; ci = inj_i[(b,t)]
                p_or_q = meas.kind == :pinj ?
                    JuMP.@expression(m, dvr*cr + dvi*ci) :
                    JuMP.@expression(m, dvi*cr - dvr*ci)
                h_si = JuMP.@expression(m, p_or_q * sb)                    # → W / var
                obj += w * (h_si - meas.value)^2
                push!(probes, (meas, h_si))
            end
        end
        JuMP.@objective(m, Min, obj)
    end

    ctx = build_opf_model(net; per_unit=per_unit, s_base=s_base,
                          add_objective=false, model_hook! = wls!,
                          optimizer=optimizer, verbose=verbose)
    _set_solver_options!(ctx.model, solver_options)
    enforce_kcl!(ctx)
    JuMP.optimize!(ctx.model)

    outcome = _solve_outcome(ctx.model)
    status = string(outcome.termination_status)
    pstatus = string(outcome.primal_status)
    solved = _publishable(outcome)
    obj = solved ? JuMP.objective_value(ctx.model) : NaN

    # Residuals (SI) from the probed h-expressions while the model is still live.
    residuals = NamedTuple[]
    for (meas, h) in probes
        est = solved ? JuMP.value(h) : NaN
        r = meas.value - est
        push!(residuals, (kind=meas.kind, bus=meas.bus, terminal=meas.terminal,
                          reference=_resolve_ref(meas, neutral),
                          measured=meas.value, estimated=est,
                          residual=r, standardized=r / meas.sigma))
    end

    result = _extract_result(ctx, outcome)
    est_bus = result["bus"]

    # Observability from the returned operating point (SI, via ybus_passive),
    # BEFORE we potentially NaN-out an unconverged estimate.
    obsv = check_observability ?
        _observability(net, measurements, neutral, zi, est_bus) :
        (observable=missing, n_states=0, rank=0, redundancy=0,
         min_singular=NaN, cond=NaN)
    if check_observability && obsv.observable === false
        @warn "solve_state_estimation: state is locally UNOBSERVABLE " *
              "(rank $(obsv.rank) < $(obsv.n_states) states); the estimate is not " *
              "unique. Add measurements or declare more zero-injection buses."
    end

    if !solved
        # An unconverged iterate is not an estimate — do not publish it as one.
        for (_, bt) in est_bus, (_, v) in bt
            for k in keys(v); v[k] = NaN; end
        end
    end

    return StateEstimationResult(status, pstatus, obj, est_bus, residuals, obsv,
                                 SolveStatus(outcome))
end

# ── observability diagnostic ─────────────────────────────────────────────────

# Local numerical identifiability: form H = ∂[measurements; zero-injection]/∂x
# at the returned voltages, where x is the rectangular voltage of every
# non-source node (source terminals are the fixed boundary). Reuses the passive
# system admittance `ybus_passive` (I = Y·V, current into the network) rather
# than re-deriving the network. `rank(H) == n_states` ⇔ locally observable.
function _observability(net, measurements, neutral, zi::Set{Tuple{String,String}},
                        est_bus)
    yb = ybus_passive(net)
    nodes = yb.nodes                                  # Vector{(bus, term)}
    N = length(nodes)
    posof = Dict(nd => k for (k, nd) in enumerate(nodes))
    # Non-broadcast real/imag: these return SparseMatrixCSC{Float64}. The
    # broadcast form `Float64.(real.(Y))` infers eltype Any on Julia 1.10.
    Yr = real(yb.Y); Yi = imag(yb.Y)

    # Fixed (source-boundary) node voltages: v_magnitude∠v_angle per source terminal.
    fixed = Dict{Tuple{String,String},ComplexF64}()
    for (_, vs) in get(net, "voltage_source", Dict())
        b = String(get(vs, "bus", "")); tmap = String.(get(vs, "terminal_map", String[]))
        mag = Float64.(get(vs, "v_magnitude", Float64[]))
        ang = Float64.(get(vs, "v_angle", zeros(length(tmap))))
        for (k, t) in enumerate(tmap)
            k <= length(mag) || continue
            fixed[(b, t)] = mag[k] * cis(k <= length(ang) ? ang[k] : 0.0)
        end
    end

    # State layout: every node that is not source-fixed contributes (vr, vi).
    state_nodes = [nd for nd in nodes if !haskey(fixed, nd)]
    ns = length(state_nodes)
    spos = Dict(nd => k for (k, nd) in enumerate(state_nodes))
    n_states = 2 * ns

    # Voltage of a node at state x (returns (re, im), 0 for ground / off-grid).
    _vnode(nd, x) = haskey(spos, nd) ? (x[spos[nd]], x[ns + spos[nd]]) :
                    haskey(fixed, nd) ? (real(fixed[nd]), imag(fixed[nd])) :
                    (zero(eltype(x)), zero(eltype(x)))

    # Zero-injection phase nodes (declared, not measured, in the Ybus node set).
    measured = Set((m.bus, m.terminal) for m in measurements if m.kind in (:pinj, :qinj))
    zi_nodes = [nd for nd in nodes if nd in zi && !(nd in measured) && !haskey(fixed, nd)]

    # Row-scaled residual map h(x): measurement rows σ-normalised; zero-injection
    # rows scaled by a representative voltage weight so ranks are comparable.
    wz = isempty(measurements) ? 1.0 : 1.0 / minimum(m.sigma for m in measurements)
    function hvec(x)
        T = eltype(x)
        Vr = Vector{T}(undef, N); Vi = Vector{T}(undef, N)
        for k in 1:N
            re, im = _vnode(nodes[k], x); Vr[k] = re; Vi[k] = im
        end
        Ir = Yr * Vr - Yi * Vi
        Ii = Yr * Vi + Yi * Vr
        rows = T[]
        for meas in measurements
            b = meas.bus; t = meas.terminal; ref = _resolve_ref(meas, neutral)
            vtr, vti = _vnode((b, t), x)
            rr, ri = ref === nothing ? (zero(T), zero(T)) : _vnode((b, ref), x)
            dvr = vtr - rr; dvi = vti - ri
            kt = get(posof, (b, t), 0)
            w = 1.0 / meas.sigma
            if meas.kind == :vmag
                push!(rows, w * sqrt(dvr^2 + dvi^2))
            elseif kt != 0 && meas.kind == :pinj
                push!(rows, w * (dvr * Ir[kt] + dvi * Ii[kt]))
            elseif kt != 0  # :qinj
                push!(rows, w * (dvi * Ir[kt] - dvr * Ii[kt]))
            else
                push!(rows, zero(T))
            end
        end
        for nd in zi_nodes
            k = posof[nd]
            push!(rows, wz * Ir[k]); push!(rows, wz * Ii[k])
        end
        rows
    end

    # Evaluate at the returned voltages (falls back to source-flat if missing).
    x0 = zeros(Float64, n_states)
    for (nd, k) in spos
        b, t = nd
        vr = get(get(get(est_bus, b, Dict()), t, Dict()), "vr", NaN)
        vi = get(get(get(est_bus, b, Dict()), t, Dict()), "vi", NaN)
        x0[k]      = isfinite(vr) ? vr : _flat_re(fixed, b)
        x0[ns + k] = isfinite(vi) ? vi : 0.0
    end

    n_states == 0 && return (observable=true, n_states=0, rank=0, redundancy=0,
                             min_singular=Inf, cond=1.0)
    H = _fd_jacobian(hvec, x0)
    sv = svdvals(H)
    smax = isempty(sv) ? 0.0 : maximum(sv)
    tol = smax * max(size(H)...) * eps(Float64)
    rk = count(>(tol), sv)
    smin = length(sv) >= n_states ? sv[n_states] : 0.0
    condn = smin > 0 ? smax / smin : Inf
    (observable = rk == n_states,
     n_states = n_states,
     rank = rk,
     redundancy = size(H, 1) - n_states,
     min_singular = smin,
     cond = condn)
end

# Central-difference Jacobian of `f: Rⁿ → Rᵐ` at `x` (steps scaled to |x| so it
# is well-conditioned across SI volt magnitudes). Used only for the rank check.
function _fd_jacobian(f, x::Vector{Float64})
    m = length(f(x)); n = length(x)
    J = Matrix{Float64}(undef, m, n)
    for j in 1:n
        h = max(1e-6, 1e-6 * abs(x[j]))
        xp = copy(x); xm = copy(x); xp[j] += h; xm[j] -= h
        J[:, j] = (f(xp) .- f(xm)) ./ (2h)
    end
    J
end

# Any source magnitude as a flat fallback voltage for a bus with no estimate.
function _flat_re(fixed, bus)
    for ((b, _), v) in fixed
        b == bus && return abs(v)
    end
    isempty(fixed) ? 1.0 : abs(first(values(fixed)))
end
