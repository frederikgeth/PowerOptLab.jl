# Network parameter estimation (calibration).
#
# Distribution GIS records are often wrong: line lengths are approximate and
# transformer/regulator tap positions drift or are mislogged. This is a THIRD
# problem specification over the same network physics — a sibling of state
# estimation. State estimation fixes the parameters and estimates the per-snapshot
# STATE; calibration inverts that: it fixes nothing about the uncertain elements
# and instead estimates the SHARED, TIME-INVARIANT parameters (line lengths, tap
# ratios) that best reproduce many snapshots of smart-meter data at once.
#
# Identifiability needs multiple time steps: a single snapshot has too few
# equations to separate a long line from a heavy load. With T snapshots of diverse
# loads the shared parameters are over-determined and the fit converges. Each
# snapshot contributes the network's KVL/KCL physics; the smart meters supply the
# known loads (as fixed injections in each snapshot net) and the measured voltage
# magnitudes (the WLS targets). The unknowns — one length per uncertain line, one
# tap multiplier per uncertain transformer — are SHARED variables on a single
# JuMP model spanning every snapshot, exactly the multi-period build pattern.
#
# The uncertain elements are the unknowns, so they are NOT part of the per-snapshot
# physics nets (which carry the known source, known lines, and the metered loads);
# they are stamped into the shared model by the hook with variable parameters:
#   * variable-length line  V_f − V_t = (r₀+jx₀)·ℓ·I   (ℓ the free length)
#   * variable-ratio ideal transformer  V_f = a·V_t,  I_s = a·I_p,  a = a₀·τ
# Both introduce bilinear terms (ℓ·I, τ·V), so the calibration is a smooth NLP
# (Ipopt), like every other nonconvex piece of this package.
#
# Prototype scope: single-phase (one terminal per element end; add one uncertain
# element per phase for a multi-phase feeder), lossless ideal-transformer taps,
# and SI units throughout (the WLS residual is in volts, and the line/tap physics
# read exactly as written — per-unit transformer base referral is deliberately out
# of scope for this didactic example).

"""
    CalibLine(; id, bus_from, bus_to, r_per_length, x_per_length, ...)

An uncertain line whose **length** is a free parameter estimated by
[`solve_parameter_estimation`](@ref). The series impedance is
`(r_per_length + j·x_per_length)·length` [Ω]; only `length` is estimated.

Single-phase: it connects `(bus_from, terminal)` to `(bus_to, terminal)`. For a
multi-phase line add one `CalibLine` per phase conductor.

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
    CalibTap(; id, bus_from, bus_to, ratio_nom=1.0, ...)

An uncertain (ideal, lossless) transformer/regulator whose **tap** is a free
parameter estimated by [`solve_parameter_estimation`](@ref). The turns ratio is
`a = ratio_nom · τ` with `V_from = a·V_to`; the tap multiplier `τ` is estimated
(so `τ=1` is the nominal ratio, `τ=1.05` a +5 % boost).

`bus_from` is the primary (regulated/HV) side, `bus_to` the secondary. Single-phase
(one `terminal` per side).

# Keywords
- `id::String` — label for reporting.
- `bus_from::String`, `bus_to::String`, `terminal::String="1"` — primary/secondary.
- `ratio_nom::Float64=1.0` — nominal primary/secondary turns ratio `a₀`.
- `tap_init::Float64=1.0` — starting guess for `τ`.
- `tap_min::Float64=0.9`, `tap_max::Float64=1.1` — bounds on `τ`.
"""
Base.@kwdef struct CalibTap
    id::String
    bus_from::String
    bus_to::String
    terminal::String = "1"
    ratio_nom::Float64 = 1.0
    tap_init::Float64 = 1.0
    tap_min::Float64 = 0.9
    tap_max::Float64 = 1.1
end

"""
    ParameterEstimationResult

Result of [`solve_parameter_estimation`](@ref).

# Fields
- `termination_status::String`, `objective::Float64` — solver status and the
  optimal weighted-residual sum `∑ (z−|V|)²/σ²` over all snapshots.
- `line_length::Dict{String,Float64}` — estimated length per `CalibLine` id.
- `tap::Dict{String,Float64}` — estimated tap multiplier `τ` per `CalibTap` id.
- `residual_rms::Float64` — root-mean-square voltage residual across all
  measurements [V] (a quick goodness-of-fit / noise-floor check).
- `snapshots::Vector{Dict{String,Any}}` — per-snapshot `extract_result` (SI),
  the fitted network state at the estimated parameters.
"""
struct ParameterEstimationResult
    termination_status::String
    objective::Float64
    line_length::Dict{String,Float64}
    tap::Dict{String,Float64}
    residual_rms::Float64
    snapshots::Vector{Dict{String,Any}}
end

# Stamp a variable-length line into the shared model: KVL with the free length ℓ
# and injection of its current into both endpoints' KCL.
function _stamp_calib_line!(ctx, ln::CalibLine, ell, t)
    m = ctx.model; vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
    f = ln.bus_from; g = ln.bus_to; c = ln.terminal
    cr = JuMP.@variable(m, base_name = "calib_cr_$(ln.id)_$t")
    ci = JuMP.@variable(m, base_name = "calib_ci_$(ln.id)_$t")
    # V_f − V_t = (r₀+jx₀)·ℓ·(cr+jci)   (bilinear in ℓ·c)
    JuMP.@constraint(m, vr[(f,c)] - vr[(g,c)] == ln.r_per_length*ell*cr - ln.x_per_length*ell*ci)
    JuMP.@constraint(m, vi[(f,c)] - vi[(g,c)] == ln.r_per_length*ell*ci + ln.x_per_length*ell*cr)
    JuMP.add_to_expression!(ctx.kcl_r[(f,c)], -cr); JuMP.add_to_expression!(ctx.kcl_i[(f,c)], -ci)
    JuMP.add_to_expression!(ctx.kcl_r[(g,c)],  cr); JuMP.add_to_expression!(ctx.kcl_i[(g,c)],  ci)
end

# Stamp a variable-ratio ideal transformer: V_f = a·V_t, I_s = a·I_p, a = a₀·τ.
# Primary current I_p leaves bus_from; secondary current I_s enters bus_to.
function _stamp_calib_tap!(ctx, tr::CalibTap, tau, t)
    m = ctx.model; vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
    f = tr.bus_from; g = tr.bus_to; c = tr.terminal
    a = JuMP.@expression(m, tr.ratio_nom * tau)
    ipr = JuMP.@variable(m, base_name = "calib_ipr_$(tr.id)_$t")
    ipi = JuMP.@variable(m, base_name = "calib_ipi_$(tr.id)_$t")
    isr = JuMP.@variable(m, base_name = "calib_isr_$(tr.id)_$t")
    isi = JuMP.@variable(m, base_name = "calib_isi_$(tr.id)_$t")
    JuMP.@constraint(m, vr[(f,c)] == a * vr[(g,c)])   # V_from = a·V_to
    JuMP.@constraint(m, vi[(f,c)] == a * vi[(g,c)])
    JuMP.@constraint(m, isr == a * ipr)               # ampere-turns: I_s = a·I_p
    JuMP.@constraint(m, isi == a * ipi)
    JuMP.add_to_expression!(ctx.kcl_r[(f,c)], -ipr); JuMP.add_to_expression!(ctx.kcl_i[(f,c)], -ipi)
    JuMP.add_to_expression!(ctx.kcl_r[(g,c)],  isr); JuMP.add_to_expression!(ctx.kcl_i[(g,c)],  isi)
end

"""
    solve_parameter_estimation(nets, measurements; lines, taps, kwargs...)
        -> ParameterEstimationResult

Calibrate uncertain **line lengths** and **transformer tap ratios** from
smart-meter data spread over multiple time steps. The uncertain elements'
parameters are shared, time-invariant unknowns; each snapshot supplies the known
loads (as the injections already baked into its net) and the measured voltage
magnitudes to fit.

Every snapshot is built into one shared JuMP model; the `lines` and `taps` are
stamped into it with a single free length / tap variable each (reused across all
snapshots), and the objective is the combined weighted-least-squares voltage
residual. Diverse loads across snapshots are what make the parameters
identifiable — a single snapshot generally cannot separate a long line from a
heavy load.

# Arguments
- `nets::AbstractVector` — `T` per-snapshot physics nets (`parse_bmopf` output),
  each carrying the known source, the **known** lines/transformers, and that
  snapshot's metered loads, but **not** the uncertain elements in `lines`/`taps`
  (those are the unknowns, stamped by this function).
- `measurements::AbstractVector` — parallel to `nets`; `measurements[t]` is a
  `Vector{Measurement}` of that snapshot's meter readings. Only `:vmag` kind is
  used (voltage magnitude, SI volts, with its `sigma` as the WLS weight).

# Keywords
- `lines::AbstractVector=CalibLine[]`, `taps::AbstractVector=CalibTap[]` — the
  uncertain elements to estimate (at least one required).
- `optimizer=Ipopt.Optimizer`, `verbose=false`, `solver_options=()`.

# Returns
A [`ParameterEstimationResult`](@ref) with the estimated lengths and tap
multipliers, the RMS voltage residual, and the per-snapshot fitted state.
"""
function solve_parameter_estimation(nets::AbstractVector, measurements::AbstractVector;
                                    lines::AbstractVector = CalibLine[],
                                    taps::AbstractVector = CalibTap[],
                                    optimizer = Ipopt.Optimizer,
                                    verbose::Bool = false,
                                    solver_options = ())
    T = length(nets)
    T >= 1 || throw(ArgumentError("need at least one snapshot"))
    length(measurements) == T ||
        throw(ArgumentError("measurements must be parallel to nets (got $(length(measurements)) vs $T)"))
    (isempty(lines) && isempty(taps)) &&
        throw(ArgumentError("nothing to estimate: supply at least one CalibLine or CalibTap"))
    allunique([[l.id for l in lines]; [t.id for t in taps]]) ||
        throw(ArgumentError("CalibLine/CalibTap ids must be unique"))

    model = JuMP.Model(optimizer)
    verbose || JuMP.set_silent(model)
    for (name, value) in solver_options
        JuMP.set_attribute(model, string(name), value)
    end

    # Shared, time-invariant unknowns: one length per line, one tap per transformer.
    ell = Dict(l.id => JuMP.@variable(model, lower_bound = l.length_min,
                   upper_bound = l.length_max, start = l.length_init,
                   base_name = "len_$(l.id)") for l in lines)
    tau = Dict(t.id => JuMP.@variable(model, lower_bound = t.tap_min,
                   upper_bound = t.tap_max, start = t.tap_init,
                   base_name = "tap_$(t.id)") for t in taps)

    # (measured value, |V|-expression) pairs collected across all snapshots for WLS.
    probes = Vector{Tuple{Float64,Float64,Any}}()   # (z, sigma, |V|)

    stamp(t) = ctx -> begin
        m = ctx.model; vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
        for l in lines; _stamp_calib_line!(ctx, l, ell[l.id], t); end
        for tr in taps; _stamp_calib_tap!(ctx, tr, tau[tr.id], t); end
        for meas in measurements[t]
            meas.kind == :vmag ||
                throw(ArgumentError("parameter estimation uses :vmag measurements; got :$(meas.kind)"))
            b = meas.bus; c = meas.terminal
            vm = JuMP.@variable(m, lower_bound = 0.0, start = meas.value,
                                base_name = "pevm_$(b)_$(c)_$t")
            JuMP.@constraint(m, vm^2 == vr[(b,c)]^2 + vi[(b,c)]^2)
            push!(probes, (meas.value, meas.sigma, vm))
        end
    end

    ctxs = [build_opf_model(nets[t]; model = model, per_unit = false,
                            add_objective = false, model_hook! = stamp(t))
            for t in 1:T]

    isempty(probes) && throw(ArgumentError("no :vmag measurements supplied"))
    JuMP.@objective(model, Min, sum(((z - vm) / s)^2 for (z, s, vm) in probes))

    foreach(enforce_kcl!, ctxs)
    JuMP.optimize!(model)

    status = string(JuMP.termination_status(model))
    solved = JuMP.primal_status(model) == JuMP.MOI.FEASIBLE_POINT
    obj = solved ? JuMP.objective_value(model) : NaN

    line_length = Dict(id => (solved ? JuMP.value(v) : NaN) for (id, v) in ell)
    tap = Dict(id => (solved ? JuMP.value(v) : NaN) for (id, v) in tau)

    # RMS voltage residual (SI) as a plain goodness-of-fit summary.
    if solved
        sse = sum((z - JuMP.value(vm))^2 for (z, _, vm) in probes)
        residual_rms = sqrt(sse / length(probes))
    else
        residual_rms = NaN
    end

    snapshots = [extract_result(ctxs[t]) for t in 1:T]
    return ParameterEstimationResult(status, obj, line_length, tap, residual_rms, snapshots)
end
