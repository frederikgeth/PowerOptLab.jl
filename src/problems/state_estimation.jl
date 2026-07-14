# Weighted-least-squares (WLS) state estimation.
#
# A DIFFERENT problem specification over the same network physics: given noisy
# measurements of an energised network, find the bus voltage state that best fits
# them. It reuses the BMOPFTools device model but with (a) no operational bounds
# (a bounds-free net yields a pure physics model with free voltages), (b) free
# injection currents at measured buses instead of fixed loads, and (c) a
# measurement-residual objective instead of generation cost — all via
# `build_opf_model(add_objective=false, model_hook! = …)`.
#
# Buses without an injection measurement are treated as zero-injection (their KCL
# residual is pinned to zero — the classic zero-injection pseudo-measurement).
# Give a bus a :pinj/:qinj measurement to make its injection a free estimated
# quantity.

"""
    Measurement(; kind, bus, value, sigma, terminal="1")

A single scalar measurement for [`solve_state_estimation`](@ref). `value` and
`sigma` are SI (volts for `:vmag`, watts for `:pinj`, vars for `:qinj`); `sigma`
is the measurement standard deviation (WLS weight `1/sigma²`).

- `kind::Symbol` — `:vmag` (voltage magnitude at `(bus, terminal)`), `:pinj`
  (active power injected into the network at `(bus, terminal)`), or `:qinj`
  (reactive power injection).
- `bus::String`, `terminal::String="1"` — where the quantity is measured
  (`terminal` is the phase conductor; injections are referenced to `neutral`,
  see [`solve_state_estimation`](@ref)).
"""
Base.@kwdef struct Measurement
    kind::Symbol
    bus::String
    value::Float64
    sigma::Float64
    terminal::String = "1"
end

"""
    StateEstimationResult

Result of [`solve_state_estimation`](@ref).

# Fields
- `termination_status::String`, `objective::Float64` — solver status and the
  optimal weighted-residual sum `∑ (z−h)²/σ²`.
- `bus::Dict{String,Any}` — the estimated SI bus voltages (`vr`, `vi`, `vm`,
  `va` per terminal), i.e. the BMOPFTools `result["bus"]`.
- `residuals::Vector{NamedTuple}` — per input measurement, in order:
  `(kind, bus, terminal, measured, estimated, residual, normalized)` with
  `residual = measured − estimated` and `normalized = residual/σ` (SI).
"""
struct StateEstimationResult
    termination_status::String
    objective::Float64
    bus::Dict{String,Any}
    residuals::Vector{NamedTuple}
end

_vscale(ctx, bus) = ctx.bases === nothing ? 1.0 : ctx.bases.v_base[bus]

"""
    solve_state_estimation(net, measurements; kwargs...) -> StateEstimationResult

Estimate the network state of `net` (a physics-only BMOPFTools net: buses, lines,
transformers, and a voltage source, ideally carrying no operational limits) that
best fits `measurements` in a weighted-least-squares sense.

# Keywords
- `neutral="n"` — return terminal for injection measurements and free injections.
  Pass `nothing` if phase terminals are referenced directly to ground.
- `per_unit=true`, `s_base=1e6` — engine unit handling; measurements stay SI.
- `optimizer=Ipopt.Optimizer`, `verbose=false`, `solver_options=()`.

# Returns
A [`StateEstimationResult`](@ref) with the estimated SI bus voltages and the
per-measurement residuals.
"""
function solve_state_estimation(net::Dict{String,Any}, measurements::AbstractVector;
                                neutral::Union{String,Nothing}="n",
                                per_unit::Bool=true,
                                s_base::Float64=1e6,
                                optimizer=Ipopt.Optimizer,
                                verbose::Bool=false,
                                solver_options=())
    isempty(measurements) && throw(ArgumentError("no measurements supplied"))

    # (measurement, SI-valued h-expression) pairs, filled by the hook for residuals.
    probes = Vector{Tuple{Measurement,Any}}()

    function wls!(ctx)
        m = ctx.model
        vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
        sb = _sbase(ctx)

        # Free injection currents at each (bus, phase) carrying an injection
        # measurement, added to KCL so those buses' voltages stay free to fit.
        inj_r = Dict{Tuple{String,String},Any}()
        inj_i = Dict{Tuple{String,String},Any}()
        for meas in measurements
            (meas.kind in (:pinj, :qinj)) || continue
            key = (meas.bus, meas.terminal)
            haskey(inj_r, key) && continue
            cr = JuMP.@variable(m, base_name = "seinj_r_$(meas.bus)_$(meas.terminal)")
            ci = JuMP.@variable(m, base_name = "seinj_i_$(meas.bus)_$(meas.terminal)")
            inj_r[key] = cr; inj_i[key] = ci
            JuMP.add_to_expression!(ctx.kcl_r[(meas.bus, meas.terminal)],  cr)
            JuMP.add_to_expression!(ctx.kcl_i[(meas.bus, meas.terminal)],  ci)
            if neutral !== nothing
                JuMP.add_to_expression!(ctx.kcl_r[(meas.bus, neutral)], -cr)
                JuMP.add_to_expression!(ctx.kcl_i[(meas.bus, neutral)], -ci)
            end
        end

        obj = zero(JuMP.QuadExpr)
        for meas in measurements
            w = 1.0 / meas.sigma^2
            b = meas.bus; t = meas.terminal
            if meas.kind == :vmag
                # Auxiliary |V| ≥ 0 with |V|² = vr²+vi² keeps the objective quadratic.
                vm = JuMP.@variable(m, base_name = "sevm_$(b)_$(t)", lower_bound = 0.0)
                JuMP.@constraint(m, vm^2 == vr[(b,t)]^2 + vi[(b,t)]^2)
                h_si = JuMP.@expression(m, vm * _vscale(ctx, b))          # → volts
                obj += w * (h_si - meas.value)^2
                push!(probes, (meas, h_si))
            elseif meas.kind in (:pinj, :qinj)
                cr = inj_r[(b,t)]; ci = inj_i[(b,t)]
                if neutral === nothing
                    dvr = vr[(b,t)]; dvi = vi[(b,t)]
                else
                    dvr = JuMP.@expression(m, vr[(b,t)] - vr[(b,neutral)])
                    dvi = JuMP.@expression(m, vi[(b,t)] - vi[(b,neutral)])
                end
                p_or_q = meas.kind == :pinj ?
                    JuMP.@expression(m, dvr*cr + dvi*ci) :
                    JuMP.@expression(m, dvi*cr - dvr*ci)
                h_si = JuMP.@expression(m, p_or_q * sb)                   # → W / var
                obj += w * (h_si - meas.value)^2
                push!(probes, (meas, h_si))
            else
                throw(ArgumentError("unknown measurement kind :$(meas.kind)"))
            end
        end
        JuMP.@objective(m, Min, obj)
    end

    ctx = build_opf_model(net; per_unit=per_unit, s_base=s_base,
                          add_objective=false, model_hook! = wls!,
                          optimizer=optimizer, verbose=verbose)
    for (name, value) in solver_options
        JuMP.set_attribute(ctx.model, string(name), value)
    end
    enforce_kcl!(ctx)
    JuMP.optimize!(ctx.model)

    status = string(JuMP.termination_status(ctx.model))
    solved = JuMP.primal_status(ctx.model) == JuMP.MOI.FEASIBLE_POINT
    obj = solved ? JuMP.objective_value(ctx.model) : NaN

    # Residuals (SI) from the probed h-expressions while the model is still live.
    residuals = NamedTuple[]
    for (meas, h) in probes
        est = solved ? JuMP.value(h) : NaN
        r = meas.value - est
        push!(residuals, (kind=meas.kind, bus=meas.bus, terminal=meas.terminal,
                          measured=meas.value, estimated=est,
                          residual=r, normalized=r / meas.sigma))
    end

    result = extract_result(ctx)
    return StateEstimationResult(status, obj, result["bus"], residuals)
end
