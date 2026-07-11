# Prototype advanced inverter model.
#
# A more detailed inverter-based-resource (IBR) than the BMOPFTools engine's
# built-in current-injection IBR, following the BMOPFTools "IBR Model Extensions"
# design doc (docs/ibr_model_extensions.md). The engine's IBR is a bounded
# current source at the point of connection (POC); this prototype adds the
# design's core structural idea — an explicit INTERNAL AC NODE behind the
# converter — and the five feature phases layered on it:
#
#   POC bus ──[filter r+jx (+ grid shunt b)]── internal node ──[converter]── DC
#     network sets V                            EMF lives here          losses/ripple
#
#   Phase 0  Output filter        series r+jx from POC to internal node (+ optional
#                                  grid-side shunt susceptance).
#   Phase 1  Internal EMF bounds   |V_int| box, or a DC-link modulation limit
#                                  |V_int| ≤ modulation_max·V_dc/√3.
#   Phase 2  Grid-forming          balanced 120° internal EMF with a bounded
#                                  magnitude decision variable v_gfm.
#   Phase 3  Converter losses      P_dc = P_ac + P_loss, non-branching, with
#                                  P_loss = p_loss_fixed + a_loss·|I| + c_loss·|I|².
#   Phase 4  Double-frequency      |Σ_k V_int_k·I_k|² ≤ p_ripple_max² (2ω ripple).
#
# Built on the BMOPFTools staged API through a `model_hook!`; it does not modify
# the engine. Runs in SI (per_unit=false) so device parameters — ohms, volts,
# amps, watts, siemens — enter directly; the per-unit base coupling the design
# doc flags (AC v_base vs DC v_dc_base) is an engine-integration concern the
# prototype avoids by staying in SI.

const _SQRT3 = sqrt(3.0)

"""
    AdvancedInverter(; id, bus, s_max, kwargs...)

A prototype inverter with an internal AC node behind an output filter. All
parameters are SI. Features are opt-in: with only `id`, `bus`, `s_max` (and a
zero filter) it is a plain grid-following converter at the POC.

# Required
- `id::String`, `bus::String` — identifier and point-of-connection bus.
- `s_max::Float64` — converter apparent-power rating (VA), applied on the
  converter side (internal-node voltage × current), per the design doc.

# Connection
- `phase_terminals=["1"]`, `neutral="n"` — phase conductor(s) and return terminal
  (`nothing` ⇒ referenced to ground). Grid-forming 120° requires three phases.

# Phase 0 — output filter
- `r_filter=0.0`, `x_filter=0.0` — series filter impedance per phase (Ω).
- `b_filter_shunt=0.0` — grid-side (POC) shunt susceptance (S).

# Phase 1 — internal EMF bounds (choose at most one mechanism)
- `v_int_min`, `v_int_max` — per-phase EMF magnitude box (V).
- `modulation_max`, `v_dc` — DC-link modulation cap `|V_int| ≤ modulation_max·v_dc/√3` (V).

# Phase 2 — grid-forming
- `grid_forming=false` — enforce a balanced positive-sequence internal EMF; the
  magnitude `v_gfm ∈ [v_int_min, v_int_max]` is a decision variable.

# Phase 3 — converter losses
- `p_loss_fixed=0.0` (W), `a_loss=0.0` (W/A), `c_loss=0.0` (W/A²).

# Phase 4 — double-frequency ripple
- `p_ripple_max` — bound on the 2ω power-ripple magnitude (VA); mainly bites for
  single-phase / unbalanced operation.

- `i_max=nothing` — optional per-conductor current limit (A).
"""
Base.@kwdef struct AdvancedInverter
    id::String
    bus::String
    phase_terminals::Vector{String} = ["1"]
    neutral::Union{String,Nothing} = "n"
    s_max::Float64
    i_max::Union{Float64,Nothing} = nothing
    r_filter::Float64 = 0.0
    x_filter::Float64 = 0.0
    b_filter_shunt::Float64 = 0.0
    v_int_min::Union{Float64,Nothing} = nothing
    v_int_max::Union{Float64,Nothing} = nothing
    modulation_max::Union{Float64,Nothing} = nothing
    v_dc::Union{Float64,Nothing} = nothing
    grid_forming::Bool = false
    p_loss_fixed::Float64 = 0.0
    a_loss::Float64 = 0.0
    c_loss::Float64 = 0.0
    p_ripple_max::Union{Float64,Nothing} = nothing
end

"""
    InverterResult

Result of [`solve_advanced_inverter`](@ref). Powers SI (W / var / VA), voltages V,
currents A.

# Fields
- `termination_status::String`
- `p_poc`, `q_poc` — active/reactive power injected at the POC (grid side).
- `p_conv`, `q_conv` — converter-side power (at the internal node).
- `p_loss`, `p_dc` — converter loss and DC-link power (`p_dc = p_conv + p_loss`).
- `v_int_mag::Vector{Float64}` — internal EMF magnitude per phase.
- `i_mag::Vector{Float64}` — converter current magnitude per phase.
- `ripple` — 2ω power-ripple magnitude `|Σ V_int·I|`.
- `bus::Dict{String,Any}` — the BMOPFTools `result["bus"]` (POC voltages, …).
"""
struct InverterResult
    termination_status::String
    p_poc::Float64
    q_poc::Float64
    p_conv::Float64
    q_conv::Float64
    p_loss::Float64
    p_dc::Float64
    v_int_mag::Vector{Float64}
    i_mag::Vector{Float64}
    ripple::Float64
    bus::Dict{String,Any}
end

# Handles the hook publishes for post-solve reporting.
struct _InvHandles
    p_poc; q_poc; p_conv; q_conv; p_loss; p_dc
    vrint::Vector{Any}; viint::Vector{Any}
    cri::Vector{Any}; cii::Vector{Any}
    ripple_re; ripple_im
end

_start_or(v, default) = (s = JuMP.start_value(v); s === nothing ? default : s)

# Stamp the inverter's internal-node model into `ctx`. Returns `_InvHandles`.
function _stamp_inverter!(ctx, inv::AdvancedInverter)
    m = ctx.model
    vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
    bus = inv.bus; phases = inv.phase_terminals; neutral = inv.neutral
    r = inv.r_filter; x = inv.x_filter
    nph = length(phases)

    vrint = Vector{Any}(undef, nph); viint = Vector{Any}(undef, nph)
    cri   = Vector{Any}(undef, nph); cii   = Vector{Any}(undef, nph)
    imag_aux = Vector{Any}(undef, nph)

    Pconv = zero(JuMP.QuadExpr); Qconv = zero(JuMP.QuadExpr)
    Ppoc  = zero(JuMP.QuadExpr); Qpoc  = zero(JuMP.QuadExpr)
    rip_re = zero(JuMP.QuadExpr); rip_im = zero(JuMP.QuadExpr)
    isq_sum = zero(JuMP.QuadExpr)     # Σ |I|²  (for c_loss)
    imag_sum = zero(JuMP.AffExpr)     # Σ |I|   (for a_loss)

    for (k, ph) in enumerate(phases)
        # POC phase-to-neutral voltage.
        if neutral === nothing
            dvr = vr[(bus, ph)]; dvi = vi[(bus, ph)]
            vseed = _start_or(vr[(bus, ph)], 0.0)
        else
            dvr = JuMP.@expression(m, vr[(bus, ph)] - vr[(bus, neutral)])
            dvi = JuMP.@expression(m, vi[(bus, ph)] - vi[(bus, neutral)])
            vseed = _start_or(vr[(bus, ph)], 0.0) - _start_or(vr[(bus, neutral)], 0.0)
        end

        vrint[k] = JuMP.@variable(m, base_name = "vrint_$(inv.id)_$(ph)")
        viint[k] = JuMP.@variable(m, base_name = "viint_$(inv.id)_$(ph)")
        cri[k]   = JuMP.@variable(m, base_name = "cri_$(inv.id)_$(ph)")
        cii[k]   = JuMP.@variable(m, base_name = "cii_$(inv.id)_$(ph)")
        # Seed the internal node near the POC voltage (physical, high-V root).
        JuMP.set_start_value(vrint[k], vseed)
        JuMP.set_start_value(viint[k], 0.0)

        # Phase 0: filter Ohm's law  V_int − V_poc = (r+jx)·I.
        JuMP.@constraint(m, vrint[k] - dvr == r*cri[k] - x*cii[k])
        JuMP.@constraint(m, viint[k] - dvi == r*cii[k] + x*cri[k])

        # Converter-side and POC-side power (same series current I).
        Pconv += JuMP.@expression(m, vrint[k]*cri[k] + viint[k]*cii[k])
        Qconv += JuMP.@expression(m, viint[k]*cri[k] - vrint[k]*cii[k])
        Ppoc  += JuMP.@expression(m, dvr*cri[k] + dvi*cii[k])
        Qpoc  += JuMP.@expression(m, dvi*cri[k] - dvr*cii[k])

        # Phase 4 accumulators: non-conjugated Σ V_int·I (converter side).
        rip_re += JuMP.@expression(m, vrint[k]*cri[k] - viint[k]*cii[k])
        rip_im += JuMP.@expression(m, vrint[k]*cii[k] + viint[k]*cri[k])

        # Inject the converter current into the POC bus KCL.
        JuMP.add_to_expression!(ctx.kcl_r[(bus, ph)], cri[k])
        JuMP.add_to_expression!(ctx.kcl_i[(bus, ph)], cii[k])
        if neutral !== nothing
            JuMP.add_to_expression!(ctx.kcl_r[(bus, neutral)], -cri[k])
            JuMP.add_to_expression!(ctx.kcl_i[(bus, neutral)], -cii[k])
        end

        # |I|² and a |I| ≥ √(cri²+cii²) auxiliary (tight whenever loss/limits bite).
        isq = JuMP.@expression(m, cri[k]^2 + cii[k]^2)
        isq_sum += isq
        if inv.i_max !== nothing
            JuMP.@constraint(m, isq <= inv.i_max^2)
        end
        im = JuMP.@variable(m, base_name = "imag_$(inv.id)_$(ph)", lower_bound = 0.0)
        JuMP.@constraint(m, im^2 >= isq)
        imag_aux[k] = im
        imag_sum += im

        # Phase 1: internal EMF magnitude bounds (skipped under grid-forming,
        # which pins the magnitude via v_gfm below).
        if !inv.grid_forming
            vmag2 = JuMP.@expression(m, vrint[k]^2 + viint[k]^2)
            if inv.modulation_max !== nothing && inv.v_dc !== nothing
                cap = inv.modulation_max * inv.v_dc / _SQRT3
                JuMP.@constraint(m, vmag2 <= cap^2)
            else
                inv.v_int_max !== nothing && JuMP.@constraint(m, vmag2 <= inv.v_int_max^2)
                inv.v_int_min !== nothing && JuMP.@constraint(m, vmag2 >= inv.v_int_min^2)
            end
        end
    end

    # Phase 0: optional grid-side shunt susceptance at the POC (phase-to-neutral).
    if inv.b_filter_shunt != 0.0
        b = inv.b_filter_shunt
        for ph in phases
            if neutral === nothing
                vrp = vr[(bus, ph)]; vip = vi[(bus, ph)]
            else
                vrp = JuMP.@expression(m, vr[(bus, ph)] - vr[(bus, neutral)])
                vip = JuMP.@expression(m, vi[(bus, ph)] - vi[(bus, neutral)])
            end
            # Capacitor y=jb injects −(jb)·V into the bus (current leaves into cap):
            #   into-bus real += b·vi ,  imag −= b·vr   (engine shunt sign).
            JuMP.add_to_expression!(ctx.kcl_r[(bus, ph)],  b, vip)
            JuMP.add_to_expression!(ctx.kcl_i[(bus, ph)], -b, vrp)
            if neutral !== nothing
                JuMP.add_to_expression!(ctx.kcl_r[(bus, neutral)], -b, vip)
                JuMP.add_to_expression!(ctx.kcl_i[(bus, neutral)],  b, vrp)
            end
        end
    end

    # Phase 2: grid-forming balanced 120° internal EMF (three-phase only).
    if inv.grid_forming && nph == 3
        h = _SQRT3 / 2
        # V_b = a²·V_a , V_c = a·V_a with a = −½ + j√3/2.
        JuMP.@constraint(m, vrint[2] == -0.5*vrint[1] + h*viint[1])
        JuMP.@constraint(m, viint[2] == -h*vrint[1]  - 0.5*viint[1])
        JuMP.@constraint(m, vrint[3] == -0.5*vrint[1] - h*viint[1])
        JuMP.@constraint(m, viint[3] ==  h*vrint[1]  - 0.5*viint[1])
    end
    if inv.grid_forming
        vgfm = JuMP.@variable(m, base_name = "vgfm_$(inv.id)", lower_bound = 0.0)
        inv.v_int_min !== nothing && JuMP.@constraint(m, vgfm >= inv.v_int_min)
        vmax_gfm = inv.v_int_max
        if inv.modulation_max !== nothing && inv.v_dc !== nothing
            cap = inv.modulation_max * inv.v_dc / _SQRT3
            vmax_gfm = vmax_gfm === nothing ? cap : min(vmax_gfm, cap)
        end
        vmax_gfm !== nothing && JuMP.@constraint(m, vgfm <= vmax_gfm)
        JuMP.@constraint(m, vrint[1]^2 + viint[1]^2 == vgfm^2)   # |V_a| = v_gfm
        JuMP.set_start_value(vgfm, _start_or(vrint[1], inv.v_int_max === nothing ? 1.0 : inv.v_int_max))
    end

    # Converter apparent-power circle P_conv²+Q_conv² ≤ s_max² (auxiliaries keep
    # it a quadratic, not quartic, constraint).
    pv = JuMP.@variable(m, base_name = "pconv_$(inv.id)")
    qv = JuMP.@variable(m, base_name = "qconv_$(inv.id)")
    JuMP.@constraint(m, pv == Pconv)
    JuMP.@constraint(m, qv == Qconv)
    JuMP.@constraint(m, pv^2 + qv^2 <= inv.s_max^2)

    # Phase 3: converter loss and DC-link power.
    P_loss = JuMP.@expression(m,
        inv.p_loss_fixed + inv.a_loss*imag_sum + inv.c_loss*isq_sum)
    P_dc = JuMP.@expression(m, Pconv + P_loss)

    # Phase 4: 2ω ripple magnitude bound (converter side).
    if inv.p_ripple_max !== nothing
        rr = JuMP.@variable(m, base_name = "rip_re_$(inv.id)")
        ri = JuMP.@variable(m, base_name = "rip_im_$(inv.id)")
        JuMP.@constraint(m, rr == rip_re)
        JuMP.@constraint(m, ri == rip_im)
        JuMP.@constraint(m, rr^2 + ri^2 <= inv.p_ripple_max^2)
    end

    return _InvHandles(Ppoc, Qpoc, Pconv, Qconv, P_loss, P_dc,
                       vrint, viint, cri, cii, rip_re, rip_im)
end

"""
    solve_advanced_inverter(net, inverter; objective=:max_export, kwargs...)
        -> InverterResult

Stamp `inverter` into `net` and solve, demonstrating the prototype internal-node
model. The network supplies the surrounding grid (a voltage source, lines, any
loads); the inverter injects at its POC bus.

# Objective
- `:max_export` — maximise the active power delivered to the grid at the POC
  (watch the converter rating, filter, EMF/modulation and ripple limits bind).
- `:min_loss` — minimise converter loss subject to `p_set` active-power delivery;
  pass `p_set` (W).

# Keywords
- `p_set=nothing` — required active-power target (W) for `:min_loss`.
- `q_set=nothing` — optional reactive-power constraint at the POC (var).
- `per_unit=false`, `s_base=1e6`, `optimizer=Ipopt.Optimizer`, `verbose=false`,
  `solver_options=()`.
"""
function solve_advanced_inverter(net::Dict{String,Any}, inverter::AdvancedInverter;
                                 objective::Symbol=:max_export,
                                 p_set::Union{Float64,Nothing}=nothing,
                                 q_set::Union{Float64,Nothing}=nothing,
                                 per_unit::Bool=false,
                                 s_base::Float64=1e6,
                                 optimizer=Ipopt.Optimizer,
                                 verbose::Bool=false,
                                 solver_options=())
    objective in (:max_export, :min_loss) ||
        throw(ArgumentError("objective must be :max_export or :min_loss, got :$objective"))
    objective == :min_loss && p_set === nothing &&
        throw(ArgumentError(":min_loss requires a p_set (W) active-power target"))

    handles = Ref{_InvHandles}()
    hook! = ctx -> begin
        h = _stamp_inverter!(ctx, inverter)
        handles[] = h
        q_set === nothing || JuMP.@constraint(ctx.model, h.q_poc == q_set)
        if objective == :max_export
            JuMP.@objective(ctx.model, Max, h.p_poc)
        else
            JuMP.@constraint(ctx.model, h.p_poc == p_set)
            JuMP.@objective(ctx.model, Min, h.p_loss)
        end
    end

    ctx = build_opf_model(net; per_unit=per_unit, s_base=s_base,
                          add_objective=false, model_hook! = hook!,
                          optimizer=optimizer, verbose=verbose)
    for (name, value) in solver_options
        JuMP.set_attribute(ctx.model, string(name), value)
    end
    enforce_kcl!(ctx)
    JuMP.optimize!(ctx.model)

    status = string(JuMP.termination_status(ctx.model))
    solved = JuMP.primal_status(ctx.model) == JuMP.MOI.FEASIBLE_POINT
    val(e) = solved ? JuMP.value(e) : NaN

    h = handles[]
    nph = length(inverter.phase_terminals)
    vint = [solved ? hypot(JuMP.value(h.vrint[k]), JuMP.value(h.viint[k])) : NaN for k in 1:nph]
    imag = [solved ? hypot(JuMP.value(h.cri[k]),   JuMP.value(h.cii[k]))   : NaN for k in 1:nph]
    ripple = solved ? hypot(JuMP.value(h.ripple_re), JuMP.value(h.ripple_im)) : NaN

    result = extract_result(ctx)
    return InverterResult(status, val(h.p_poc), val(h.q_poc), val(h.p_conv),
                          val(h.q_conv), val(h.p_loss), val(h.p_dc),
                          vint, imag, ripple, result["bus"])
end
