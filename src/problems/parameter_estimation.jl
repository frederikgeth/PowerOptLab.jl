# Network parameter estimation (calibration).
#
# A THIRD problem specification over the same network physics — the shared-parameter
# dual of state estimation, and a multi-time-step extension of it. State estimation
# fixes the parameters and estimates the per-snapshot state; calibration fixes
# nothing about the uncertain elements and estimates the shared, time-invariant
# parameters (line lengths, transformer tap ratios) that best reproduce many
# snapshots of smart-meter data at once.
#
# This mirrors the formulation of Vanin, Geth, Heidari & Van Hertem, "Distribution
# System State and Impedance Estimation Augmented with Carson's Equations"
# (arXiv:2506.04949): a weighted measurement-residual objective over a time series,
# the multiconductor IVR (rectangular current–voltage) power flow as the physics,
# time-invariant impedances written as nominal·length (their Eq. 12), and the
# smart-meter measurement set as noisy (P, Q, |V|) triples per user (their Eqs.
# 4–6). Identifiability needs multiple time steps: a single snapshot cannot
# separate a long line from a heavy load, whereas diverse loads over T snapshots
# over-determine the shared parameters.
#
# ── How the uncertain elements are stamped ──────────────────────────────────
# The two element classes are handled DIFFERENTLY, because the BMOPFTools engine
# exposes them differently:
#
#   * Transformer taps use the engine's NATIVE free-tap variable: a transformer
#     carrying `tap_min`/`tap_max` gets a `:tap` decision variable that the engine
#     threads through its (per-unit-correct, base-referred) winding constraints.
#     We keep the transformer in every snapshot net and couple its tap variable
#     equal across snapshots. Nothing is re-derived — the engine's transformer
#     physics (leakage, losses, base referral) is reused unchanged.
#
#   * Line lengths have NO native free variable, so the uncertain lines are OMITTED
#     from the physics nets and re-stamped here with a shared free length ℓ:
#         V_f − V_t = (r₀+jx₀)·ℓ·I         (per-unit: r₀,x₀ scaled by z_base)
#     This is the one genuine "replace Ohm's law" workaround (see the dev notes in
#     BMOPFTools docs/src/dev/opf_engine.md and the open issue on a free line
#     impedance/length variable).
#
# Both the ℓ·I product and the P/Q measurement projections are bilinear, so the
# calibration is a smooth NLP (Ipopt), consistent with the rest of this package.
#
# Prototype scope: single-phase elements (one terminal per end — add one CalibLine
# per phase for a multi-phase feeder). Per-unit and SI are both supported and give
# identical results.

"""
    CalibLine(; id, bus_from, bus_to, r_per_length, x_per_length, ...)

An uncertain line whose **length** is a free parameter estimated by
[`solve_parameter_estimation`](@ref). The series impedance is
`(r_per_length + j·x_per_length)·length` [Ω]; only `length` is estimated (the
nominal per-length impedance is treated as known, following the construction-code
model `Z = Zⁿᵒᵐ·ℓ` of Vanin et al.).

The uncertain lines must be **omitted** from the physics nets — they are the
unknowns, stamped by this function. Single-phase: it connects `(bus_from, terminal)`
to `(bus_to, terminal)`; add one `CalibLine` per phase for a multi-phase line.

# Keywords
- `id::String` — label for reporting.
- `bus_from::String`, `bus_to::String`, `terminal::String="1"` — the endpoints.
- `r_per_length::Float64`, `x_per_length::Float64` — per-unit-length series R/X [Ω].
- `length_init::Float64=1.0` — starting guess handed to the solver.
- `length_min::Float64=0.1`, `length_max::Float64=10.0` — bounds on the estimate.
"""
Base.@kwdef struct CalibLine
    id::String
    bus_from::String
    bus_to::String
    terminal::String = "1"
    r_per_length::Float64
    x_per_length::Float64
    length_init::Float64 = 1.0
    length_min::Float64 = 0.1
    length_max::Float64 = 10.0
end

"""
    CalibTap(; id, tap_min=0.9, tap_max=1.1)

An uncertain transformer/regulator **tap** estimated by
[`solve_parameter_estimation`](@ref). Unlike [`CalibLine`](@ref), the transformer
stays **in** the physics nets: this just names the transformer `id` and the tap
bounds. The function sets `tap_min`/`tap_max` on that transformer so the engine
creates its native free-tap variable, then couples the tap equal across all
snapshots and reports the estimate.

The reported estimate is the tap **multiplier** `τ` on the nominal ratio
`N₀ = v_nom_from / v_nom_to` (so `τ = 1` is nominal); the effective turns ratio is
`N = N₀·τ`.

# Keywords
- `id::String` — transformer id present in every snapshot net.
- `tap_min::Float64=0.9`, `tap_max::Float64=1.1` — bounds on the multiplier `τ`.
"""
Base.@kwdef struct CalibTap
    id::String
    tap_min::Float64 = 0.9
    tap_max::Float64 = 1.1
end

"""
    ParameterEstimationResult

Result of [`solve_parameter_estimation`](@ref).

# Fields
- `termination_status::String`, `objective::Float64` — solver status and the
  optimal weighted-residual value.
- `line_length::Dict{String,Float64}` — estimated length per `CalibLine` id.
- `tap::Dict{String,Float64}` — estimated tap multiplier `τ` per `CalibTap` id.
- `residual_rms::Float64` — root-mean-square voltage-measurement residual [V], a
  quick goodness-of-fit / noise-floor check.
- `snapshots::Vector{Dict{String,Any}}` — per-snapshot `extract_result` (SI), the
  fitted network state at the estimated parameters.
"""
struct ParameterEstimationResult
    termination_status::String
    objective::Float64
    line_length::Dict{String,Float64}
    tap::Dict{String,Float64}
    residual_rms::Float64
    snapshots::Vector{Dict{String,Any}}
end

# Per-bus voltage base (volts) — 1.0 in SI mode, v_base[bus] in per-unit mode.
_pe_vbase(ctx, bus) = ctx.bases === nothing ? 1.0 : ctx.bases.v_base[bus]
# Per-bus impedance base (Ω) — 1.0 in SI mode.
_pe_zbase(ctx, bus) = ctx.bases === nothing ? 1.0 : ctx.bases.z_base[bus]

# Stamp a variable-length line: KVL with the free length ℓ (per-unit-scaled R/X)
# and injection of its current into both endpoints' KCL.
function _stamp_calib_line!(ctx, ln::CalibLine, ell, t)
    m = ctx.model; vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
    f = ln.bus_from; g = ln.bus_to; c = ln.terminal
    zb = _pe_zbase(ctx, f)                 # line: same base both ends (no ratio)
    r = ln.r_per_length / zb; x = ln.x_per_length / zb
    cr = JuMP.@variable(m, base_name = "calib_cr_$(ln.id)_$t")
    ci = JuMP.@variable(m, base_name = "calib_ci_$(ln.id)_$t")
    JuMP.@constraint(m, vr[(f,c)] - vr[(g,c)] == r*ell*cr - x*ell*ci)
    JuMP.@constraint(m, vi[(f,c)] - vi[(g,c)] == r*ell*ci + x*ell*cr)
    JuMP.add_to_expression!(ctx.kcl_r[(f,c)], -cr); JuMP.add_to_expression!(ctx.kcl_i[(f,c)], -ci)
    JuMP.add_to_expression!(ctx.kcl_r[(g,c)],  cr); JuMP.add_to_expression!(ctx.kcl_i[(g,c)],  ci)
end

# Locate a transformer id across subtype blocks; return (subtype_entries, N0) and
# stamp the tap bounds into a deep-copied net so the engine frees the tap.
function _prepare_tap_nets(nets, taps)
    isempty(taps) && return nets, Dict{String,Float64}()
    n0 = Dict{String,Float64}()
    prepared = map(nets) do net
        net = deepcopy(net)
        xdict = get(net, "transformer", nothing)
        xdict isa Dict || throw(ArgumentError("nets have no transformer for CalibTap"))
        for tp in taps
            found = false
            for (_, entries) in xdict
                entries isa Dict || continue
                haskey(entries, tp.id) || continue
                e = entries[tp.id]
                e["tap_min"] = tp.tap_min; e["tap_max"] = tp.tap_max
                vf = Float64(get(e, "v_nom_from", 1.0)); vt = Float64(get(e, "v_nom_to", 1.0))
                n0[tp.id] = iszero(vt) ? 1.0 : vf/vt
                found = true; break
            end
            found || throw(ArgumentError("CalibTap id '$(tp.id)' not found among transformers"))
        end
        net
    end
    return prepared, n0
end

"""
    solve_parameter_estimation(nets, measurements; lines, taps, kwargs...)
        -> ParameterEstimationResult

Calibrate uncertain **line lengths** and **transformer tap ratios** from
smart-meter time series. The uncertain elements' parameters are shared,
time-invariant unknowns; each snapshot supplies noisy `(P, Q, |V|)` meter readings
that are fit in a weighted-least-squares (or robust weighted-least-absolute-value)
sense across the whole horizon.

The meters' loads are **not** baked in as exact injections — each measured user
gets a free injection current fit to its noisy `P`/`Q` readings (as in state
estimation), and the phase-to-neutral voltage magnitude is fit to `|V|`. Buses
with no injection measurement are zero-injection.

# Arguments
- `nets::AbstractVector` — `T` per-snapshot **physics** nets (`parse_bmopf`
  output): the known source, the known lines, and the transformers named by
  `taps` (kept, so the engine frees their taps). The uncertain `lines` are
  **omitted** — they are the unknowns. No load objects are needed; injections come
  from the measurements.
- `measurements::AbstractVector` — parallel to `nets`; `measurements[t]` is a
  `Vector{Measurement}` of that snapshot's meter readings. `:vmag` (SI volts),
  `:pinj`, `:qinj` (SI W/var, injection into the network — negative for a load)
  are all treated as noisy, weighted by their `sigma`.

# Keywords
- `lines::AbstractVector=CalibLine[]`, `taps::AbstractVector=CalibTap[]` — the
  uncertain elements to estimate (at least one required).
- `neutral="n"` — return terminal for the phase-to-neutral projections; pass
  `nothing` if phase terminals are referenced directly to ground.
- `objective=:wls` — `:wls` (weighted least squares, smooth) or `:wlav` (weighted
  least absolute value, better bad-data rejection; the choice in Vanin et al.).
- `per_unit=true`, `s_base=1e6` — engine unit handling; measurements stay SI.
- `optimizer=Ipopt.Optimizer`, `verbose=false`, `solver_options=()`.

# Returns
A [`ParameterEstimationResult`](@ref) with the estimated lengths and tap
multipliers, the RMS voltage residual, and the per-snapshot fitted state.
"""
function solve_parameter_estimation(nets::AbstractVector, measurements::AbstractVector;
                                    lines::AbstractVector = CalibLine[],
                                    taps::AbstractVector = CalibTap[],
                                    neutral::Union{String,Nothing} = "n",
                                    objective::Symbol = :wls,
                                    per_unit::Bool = true,
                                    s_base::Float64 = 1e6,
                                    optimizer = Ipopt.Optimizer,
                                    verbose::Bool = false,
                                    solver_options = ())
    T = length(nets)
    T >= 1 || throw(ArgumentError("need at least one snapshot"))
    length(measurements) == T ||
        throw(ArgumentError("measurements must be parallel to nets (got $(length(measurements)) vs $T)"))
    (isempty(lines) && isempty(taps)) &&
        throw(ArgumentError("nothing to estimate: supply at least one CalibLine or CalibTap"))
    objective in (:wls, :wlav) || throw(ArgumentError("objective must be :wls or :wlav, got :$objective"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError("s_base must be finite and > 0"))
    allunique([[l.id for l in lines]; [t.id for t in taps]]) ||
        throw(ArgumentError("CalibLine/CalibTap ids must be unique"))

    nets, n0 = _prepare_tap_nets(nets, taps)

    model = JuMP.Model(optimizer)
    verbose || JuMP.set_silent(model)
    _set_solver_options!(model, solver_options)

    # Shared unknown: one free length per uncertain line (dimensionless, SI).
    ell = Dict(l.id => JuMP.@variable(model, lower_bound = l.length_min,
                   upper_bound = l.length_max, start = l.length_init,
                   base_name = "len_$(l.id)") for l in lines)

    # (measured value, sigma, SI h-expression) triples for the residual objective.
    probes = Vector{Tuple{Float64,Float64,Any}}()
    vprobes = Vector{Tuple{Float64,Any}}()      # (measured |V|, |V|-expr) for RMS

    stamp(t) = ctx -> begin
        m = ctx.model; vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
        sb = ctx.bases === nothing ? 1.0 : ctx.bases.s_base
        for l in lines; _stamp_calib_line!(ctx, l, ell[l.id], t); end

        # Free injection currents at each measured (bus, phase) carrying a P/Q
        # measurement, added to KCL so those buses' voltages stay free to fit.
        inj = Dict{Tuple{String,String},Tuple{Any,Any}}()
        for meas in measurements[t]
            meas.kind in (:pinj, :qinj) || continue
            key = (meas.bus, meas.terminal); haskey(inj, key) && continue
            cr = JuMP.@variable(m, base_name = "peinj_r_$(meas.bus)_$(meas.terminal)_$t")
            ci = JuMP.@variable(m, base_name = "peinj_i_$(meas.bus)_$(meas.terminal)_$t")
            inj[key] = (cr, ci)
            JuMP.add_to_expression!(ctx.kcl_r[(meas.bus, meas.terminal)], cr)
            JuMP.add_to_expression!(ctx.kcl_i[(meas.bus, meas.terminal)], ci)
            if neutral !== nothing
                JuMP.add_to_expression!(ctx.kcl_r[(meas.bus, neutral)], -cr)
                JuMP.add_to_expression!(ctx.kcl_i[(meas.bus, neutral)], -ci)
            end
        end

        for meas in measurements[t]
            b = meas.bus; c = meas.terminal
            if neutral === nothing
                dvr = vr[(b,c)]; dvi = vi[(b,c)]
            else
                dvr = JuMP.@expression(m, vr[(b,c)] - vr[(b,neutral)])
                dvi = JuMP.@expression(m, vi[(b,c)] - vi[(b,neutral)])
            end
            if meas.kind == :vmag
                vm = JuMP.@variable(m, lower_bound = 0.0, start = meas.value / _pe_vbase(ctx,b),
                                    base_name = "pevm_$(b)_$(c)_$t")
                JuMP.@constraint(m, vm^2 == dvr^2 + dvi^2)
                h = JuMP.@expression(m, vm * _pe_vbase(ctx, b))       # → volts
                push!(probes, (meas.value, meas.sigma, h)); push!(vprobes, (meas.value, h))
            elseif meas.kind in (:pinj, :qinj)
                cr, ci = inj[(b,c)]
                pq = meas.kind == :pinj ?
                    JuMP.@expression(m, dvr*cr + dvi*ci) :
                    JuMP.@expression(m, dvi*cr - dvr*ci)
                h = JuMP.@expression(m, pq * sb)                       # → W / var
                push!(probes, (meas.value, meas.sigma, h))
            else
                throw(ArgumentError("unknown measurement kind :$(meas.kind)"))
            end
        end
    end

    ctxs = [build_opf_model(nets[t]; model = model, per_unit = per_unit, s_base = s_base,
                            add_objective = false, model_hook! = stamp(t))
            for t in 1:T]

    # Couple each transformer's native free tap equal across all snapshots.
    for tp in taps
        haskey(ctxs[1].vars, :tap) && haskey(ctxs[1].vars[:tap], tp.id) ||
            throw(ArgumentError("transformer '$(tp.id)' did not expose a free tap; check tap bounds"))
        for t in 2:T
            JuMP.@constraint(model, ctxs[t].vars[:tap][tp.id] == ctxs[1].vars[:tap][tp.id])
        end
    end

    isempty(probes) && throw(ArgumentError("no measurements supplied"))
    if objective == :wls
        JuMP.@objective(model, Min, sum(((z - h)/s)^2 for (z, s, h) in probes))
    else                                    # weighted least absolute value
        obj = JuMP.AffExpr(0.0)
        for (z, s, h) in probes
            rho = JuMP.@variable(model, lower_bound = 0.0)
            JuMP.@constraint(model, rho >= (h - z)/s)
            JuMP.@constraint(model, rho >= (z - h)/s)
            JuMP.add_to_expression!(obj, rho)
        end
        JuMP.@objective(model, Min, obj)
    end

    foreach(enforce_kcl!, ctxs)
    JuMP.optimize!(model)

    outcome = _solve_outcome(model)
    status = string(outcome.termination_status)
    solved = _publishable(outcome)
    obj = solved ? JuMP.objective_value(model) : NaN

    line_length = Dict(id => (solved ? JuMP.value(v) : NaN) for (id, v) in ell)
    tap = Dict(tp.id => (solved ? JuMP.value(ctxs[1].vars[:tap][tp.id]) / get(n0, tp.id, 1.0) : NaN)
               for tp in taps)

    residual_rms = if solved && !isempty(vprobes)
        sqrt(sum((z - JuMP.value(h))^2 for (z, h) in vprobes) / length(vprobes))
    else
        NaN
    end

    snapshots = [_extract_result(ctxs[t], outcome) for t in 1:T]
    return ParameterEstimationResult(status, obj, line_length, tap, residual_rms, snapshots)
end
