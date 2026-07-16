# Prototype advanced inverter model.
#
# A more detailed inverter-based-resource (IBR) than the BMOPFTools engine's
# built-in current-injection IBR. It combines the "IBR Model Extensions" design
# doc's internal-AC-node structure with exact three-phase converter
# feasible-region models drawn from ongoing research (fundamental-frequency
# switching-polytope voltage feasibility with a 2ω bus-ripple derating):
#
#   POC bus ──[filter r+jx (+ grid shunt b)]── internal node ──[converter]── DC
#     network sets V                            EMF lives here          losses/ripple
#
# Layered, opt-in features:
#   - Output filter (series r+jx + optional grid-side shunt).
#   - Internal EMF magnitude box, or the exact switching-polytope DC-utilisation
#     limit of a chosen three-phase topology (see below).
#   - Grid-forming balanced 120° internal EMF (bounded magnitude decision var).
#   - Converter losses: non-branching P_dc = P_conv + P_loss.
#   - Double-frequency (2ω) power ripple, either a standalone cap (single-phase)
#     or a bus-ripple phasor that endogenously derates the DC rails (three-phase).
#
# Three-phase topologies (`topology`), with EXACT time-sampled switching-polytope
# voltage feasibility on the converter output U_x = internal-node voltage V_int:
#   :THREE_LEG  3-leg 3-wire  — pairwise line-to-line rails, no neutral current.
#   :FOUR_LEG   4-leg 4-wire  — pairwise + per-phase rails, neutral via the 4th leg.
#   :SPLIT_DC   split-cap 4-wire — half-bus per-phase rails + midpoint ripple.
# `:SINGLE_PHASE` (default) keeps the simpler internal-EMF model.
#
# Built on the BMOPFTools staged API through a `model_hook!`; it does not modify
# the engine. Works in SI (per_unit=false) or per-unit (per_unit=true): device
# parameters are SI and scaled to model units via `ctx.bases` (Vdc/Cdc/In_max
# stay SI, the AC↔DC coupling scales through the POC bus's v_base/i_base/s_base).

const _SQRT3 = sqrt(3.0)
# Positive start for the current-magnitude aux im (im² = |I|², im ≥ 0), to keep
# interior-point iterates off the I = 0 point where that equality's Jacobian is
# degenerate (per-unit current).
const _IMAG_START = 1e-3
const _PAIRS_IDX = ((1, 2), (2, 3), (3, 1))
_sample_grid(N::Int) = [2pi * (k - 1) / N for k in 1:N]

"""
    AdvancedInverter(; id, bus, s_max, kwargs...)

A prototype inverter with an internal AC node behind an output filter. All
parameters are SI. Features are opt-in: with only `id`, `bus`, `s_max` (and a
zero filter) it is a plain grid-following converter at the POC.

# Required
- `id::String`, `bus::String` — identifier and point-of-connection bus.
- `s_max::Float64` — converter apparent-power rating (VA), applied on the
  converter side (internal-node voltage × current).

# Connection
- `phase_terminals=["1"]`, `neutral="n"` — phase conductor(s) and return terminal
  (`nothing` ⇒ referenced to ground). The three-phase topologies and grid-forming
  require three phases.

# Topology (three-phase switching-polytope DC-utilisation model)
- `topology=:SINGLE_PHASE` — one of `:SINGLE_PHASE` (internal-EMF model),
  `:THREE_LEG`, `:FOUR_LEG`, `:SPLIT_DC`. The three-phase topologies apply the
  exact time-sampled switching-polytope voltage feasibility on the internal-node
  voltage and require `v_dc` and `c_dc` (and `In_max` for the 4-wire ones).
- `v_dc` — DC-link voltage (V). `c_dc` — DC-link capacitance (F; per half for split).
- `m_max=1.0` — modulation index limit (dimensionless).
- `In_max` — neutral current limit (A): split = cap ripple rating, 4-leg = leg rating.
- `dv2_max` — optional cap on the 2ω bus-ripple amplitude (V).
- `n_samples=36` — time-sampling grid for the exact feasibility (gap ~(π/N)²).
- `f=50.0` — fundamental frequency (Hz).

# Output filter
- `r_filter=0.0`, `x_filter=0.0` — series filter impedance per phase (Ω).
- `b_filter_shunt=0.0` — grid-side (POC) shunt susceptance (S).

# Internal EMF box / single-phase modulation
- `v_int_min`, `v_int_max` — per-phase EMF magnitude box (V; applies to all topologies).
- `modulation_max` — SINGLE_PHASE only: DC-link cap `|V_int| ≤ modulation_max·v_dc/√3`.

# Grid-forming
- `grid_forming=false` — balanced positive-sequence internal EMF; the magnitude
  `v_gfm ∈ [v_int_min, v_int_max]` is a decision variable (composes with a topology).

# Converter losses
- `p_loss_fixed=0.0` (W), `a_loss=0.0` (W/A), `c_loss=0.0` (W/A²).

# Double-frequency ripple / current
- `p_ripple_max` — SINGLE_PHASE only: bound on the 2ω power-ripple magnitude (VA).
- `i_max=nothing` — optional per-conductor current limit (A).
"""
Base.@kwdef struct AdvancedInverter <: AbstractDevice
    id::String
    bus::String
    phase_terminals::Vector{String} = ["1"]
    neutral::Union{String,Nothing} = "n"
    s_max::Float64
    i_max::Union{Float64,Nothing} = nothing
    topology::Symbol = :SINGLE_PHASE
    v_dc::Union{Float64,Nothing} = nothing
    c_dc::Union{Float64,Nothing} = nothing
    m_max::Float64 = 1.0
    In_max::Union{Float64,Nothing} = nothing
    dv2_max::Union{Float64,Nothing} = nothing
    n_samples::Int = 36
    f::Float64 = 50.0
    r_filter::Float64 = 0.0
    x_filter::Float64 = 0.0
    b_filter_shunt::Float64 = 0.0
    v_int_min::Union{Float64,Nothing} = nothing
    v_int_max::Union{Float64,Nothing} = nothing
    modulation_max::Union{Float64,Nothing} = nothing
    grid_forming::Bool = false
    p_loss_fixed::Float64 = 0.0
    a_loss::Float64 = 0.0
    c_loss::Float64 = 0.0
    p_ripple_max::Union{Float64,Nothing} = nothing
end

const _THREE_PHASE_TOPOLOGIES = (:THREE_LEG, :FOUR_LEG, :SPLIT_DC)

function _validate_inverter(inv::AdvancedInverter, nets=())
    _validate_connection(inv.id, inv.bus, inv.phase_terminals, inv.neutral, nets)
    inv.topology in (:SINGLE_PHASE, _THREE_PHASE_TOPOLOGIES...) ||
        throw(ArgumentError("unknown topology :$(inv.topology)"))

    isfinite(inv.s_max) && inv.s_max > 0 || throw(ArgumentError(
        "inverter '$(inv.id)' s_max must be finite and > 0"))
    for (name, value) in (("r_filter", inv.r_filter), ("x_filter", inv.x_filter),
                          ("b_filter_shunt", inv.b_filter_shunt),
                          ("m_max", inv.m_max), ("f", inv.f),
                          ("p_loss_fixed", inv.p_loss_fixed),
                          ("a_loss", inv.a_loss), ("c_loss", inv.c_loss))
        isfinite(value) || throw(ArgumentError(
            "inverter '$(inv.id)' $name must be finite"))
    end
    inv.r_filter >= 0 || throw(ArgumentError(
        "inverter '$(inv.id)' r_filter must be >= 0"))
    inv.m_max > 0 || throw(ArgumentError(
        "inverter '$(inv.id)' m_max must be > 0"))
    inv.f > 0 || throw(ArgumentError(
        "inverter '$(inv.id)' f must be > 0"))
    all(x -> x >= 0, (inv.p_loss_fixed, inv.a_loss, inv.c_loss)) ||
        throw(ArgumentError("inverter '$(inv.id)' loss coefficients must be >= 0"))

    for (name, value) in (("i_max", inv.i_max), ("v_dc", inv.v_dc),
                          ("c_dc", inv.c_dc), ("In_max", inv.In_max),
                          ("dv2_max", inv.dv2_max), ("v_int_max", inv.v_int_max),
                          ("modulation_max", inv.modulation_max),
                          ("p_ripple_max", inv.p_ripple_max))
        value === nothing || (isfinite(value) && value > 0) || throw(ArgumentError(
            "inverter '$(inv.id)' $name must be finite and > 0 when supplied"))
    end
    if inv.v_int_min !== nothing
        isfinite(inv.v_int_min) && inv.v_int_min >= 0 || throw(ArgumentError(
            "inverter '$(inv.id)' v_int_min must be finite and >= 0"))
    end
    inv.v_int_min !== nothing && inv.v_int_max !== nothing &&
        inv.v_int_min > inv.v_int_max && throw(ArgumentError(
            "inverter '$(inv.id)' requires v_int_min <= v_int_max"))

    is_3ph = inv.topology in _THREE_PHASE_TOPOLOGIES
    if is_3ph
        length(inv.phase_terminals) == 3 || throw(ArgumentError(
            "topology :$(inv.topology) requires 3 phase_terminals"))
        (inv.v_dc !== nothing && inv.c_dc !== nothing) || throw(ArgumentError(
            "topology :$(inv.topology) requires v_dc and c_dc"))
        inv.n_samples >= 4 || throw(ArgumentError(
            "topology :$(inv.topology) requires n_samples >= 4"))
        inv.modulation_max === nothing || throw(ArgumentError(
            "inverter '$(inv.id)' modulation_max applies only to :SINGLE_PHASE"))
        inv.p_ripple_max === nothing || throw(ArgumentError(
            "inverter '$(inv.id)' p_ripple_max applies only to :SINGLE_PHASE"))
        (inv.topology in (:FOUR_LEG, :SPLIT_DC) && inv.In_max === nothing) &&
            throw(ArgumentError("topology :$(inv.topology) requires In_max"))
    elseif inv.modulation_max !== nothing && inv.v_dc === nothing
        throw(ArgumentError(
            "inverter '$(inv.id)' modulation_max requires v_dc"))
    end
    inv.grid_forming && length(inv.phase_terminals) != 3 && throw(ArgumentError(
        "grid-forming inverter '$(inv.id)' requires 3 phase_terminals"))
    return nothing
end

function validate_device(inv::AdvancedInverter, nets=(); periods::Integer=length(nets))
    periods == length(nets) || throw(ArgumentError(
        "period count $periods does not match $(length(nets)) network snapshots"))
    _validate_inverter(inv, nets)
end

"""
    InverterResult

Result of [`solve_advanced_inverter`](@ref). Powers SI (W / var / VA), voltages V,
currents A.

# Fields
- `termination_status::String`, `topology::Symbol`
- `p_poc`, `q_poc` — active/reactive power injected at the POC (grid side).
- `p_conv`, `q_conv` — converter-side power (at the internal node).
- `p_loss`, `p_dc` — converter loss and DC-link power (`p_dc = p_conv + p_loss`).
- `v_int_mag::Vector{Float64}` — internal EMF magnitude per phase.
- `i_mag::Vector{Float64}` — converter current magnitude per phase.
- `i_neutral` — neutral current magnitude `|I_n|` (0 for 3-wire / single-phase).
- `ripple` — 2ω power-ripple magnitude `|Σ V_int·I|` (VA).
- `dv2` — 2ω bus-ripple amplitude (V; three-phase topologies, else 0).
- `bus::Dict{String,Any}` — the BMOPFTools `result["bus"]` (POC voltages, …).
"""
struct InverterResult <: AbstractSolveResult
    termination_status::String
    topology::Symbol
    p_poc::Float64
    q_poc::Float64
    p_conv::Float64
    q_conv::Float64
    p_loss::Float64
    p_dc::Float64
    v_int_mag::Vector{Float64}
    i_mag::Vector{Float64}
    i_neutral::Float64
    ripple::Float64
    dv2::Float64
    bus::Dict{String,Any}
    solve::SolveStatus
end

solve_status(result::InverterResult) = result.solve

solve_diagnostics(result::InverterResult) =
    (topology=result.topology, p_loss=result.p_loss, ripple=result.ripple,
     neutral_current=result.i_neutral)

# Handles the hook publishes for post-solve reporting. Power/current/voltage
# expressions are in MODEL units; `sb`/`vb`/`ib` convert them back to SI.
struct _InvHandles
    p_poc; q_poc; p_conv; q_conv; p_loss; p_dc
    vrint::Vector{Any}; viint::Vector{Any}
    cri::Vector{Any}; cii::Vector{Any}
    ripple_re; ripple_im
    in_re; in_im          # SI neutral-current components (or nothing)
    dre; dim              # SI 2ω ripple-phasor components (or nothing)
    sb::Float64; vb::Float64; ib::Float64
end

_start_or(v, default) = (s = JuMP.start_value(v); s === nothing ? default : s)

# Stamp the inverter's internal-node model into `ctx`. Returns `_InvHandles`.
function _stamp_inverter!(ctx, inv::AdvancedInverter)
    m = ctx.model
    vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
    bus = inv.bus; phases = inv.phase_terminals; neutral = inv.neutral
    nph = length(phases)

    # Per-unit base factors (1.0 in SI mode), looked up per bus like the engine.
    bs = ctx.bases
    sb = bs === nothing ? 1.0 : bs.s_base
    vb = bs === nothing ? 1.0 : get(bs.v_base, bus, 1.0)
    ib = bs === nothing ? 1.0 : get(bs.i_base, bus, 1.0)
    zb = bs === nothing ? 1.0 : get(bs.z_base, bus, 1.0)

    r = inv.r_filter / zb           # Ω → per-unit impedance
    x = inv.x_filter / zb

    is_3ph_topo = inv.topology in _THREE_PHASE_TOPOLOGIES

    vrint = Vector{Any}(undef, nph); viint = Vector{Any}(undef, nph)
    cri   = Vector{Any}(undef, nph); cii   = Vector{Any}(undef, nph)

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

        # Output filter Ohm's law  V_int − V_poc = (r+jx)·I  (per-unit impedance).
        JuMP.@constraint(m, vrint[k] - dvr == r*cri[k] - x*cii[k])
        JuMP.@constraint(m, viint[k] - dvi == r*cii[k] + x*cri[k])

        # Converter-side and POC-side power (same series current I), model units.
        Pconv += JuMP.@expression(m, vrint[k]*cri[k] + viint[k]*cii[k])
        Qconv += JuMP.@expression(m, viint[k]*cri[k] - vrint[k]*cii[k])
        Ppoc  += JuMP.@expression(m, dvr*cri[k] + dvi*cii[k])
        Qpoc  += JuMP.@expression(m, dvi*cri[k] - dvr*cii[k])

        # Non-conjugated Σ V_int·I (converter side) — the 2ω oscillating power S2.
        rip_re += JuMP.@expression(m, vrint[k]*cri[k] - viint[k]*cii[k])
        rip_im += JuMP.@expression(m, vrint[k]*cii[k] + viint[k]*cri[k])

        # Inject the converter current into the POC bus KCL.
        JuMP.add_to_expression!(ctx.kcl_r[(bus, ph)], cri[k])
        JuMP.add_to_expression!(ctx.kcl_i[(bus, ph)], cii[k])
        if neutral !== nothing
            JuMP.add_to_expression!(ctx.kcl_r[(bus, neutral)], -cri[k])
            JuMP.add_to_expression!(ctx.kcl_i[(bus, neutral)], -cii[k])
        end

        # |I|² (for the quadratic loss term and the current cap) and the current
        # magnitude |I| = √(cri²+cii²) via an auxiliary `im ≥ 0`. The magnitude only
        # feeds the LINEAR loss term a_loss·|I|, and how tightly `im` is pinned to
        # |I| matters only there:
        #   • a_loss ≠ 0 → the implicit square-root EQUALITY `im² = |I|²`. An epigraph
        #     `im² ≥ isq` is tight only when the objective necessarily minimises im,
        #     which FAILS for charging / negative-price objectives that would then
        #     inflate im to manufacture losses and reverse the requested dispatch.
        #     (Degenerate Jacobian only at I = 0, ∂(im²−isq)/∂im = 2·im → 0, so im
        #     gets a small positive start.)
        #   • a_loss == 0 → the epigraph `im² ≥ isq`. im is unused in the objective,
        #     so the loose form is harmless and better conditioned for the solve.
        isq = JuMP.@expression(m, cri[k]^2 + cii[k]^2)
        isq_sum += isq
        if inv.i_max !== nothing
            JuMP.@constraint(m, isq <= (inv.i_max / ib)^2)
        end
        im = JuMP.@variable(m, base_name = "imag_$(inv.id)_$(ph)", lower_bound = 0.0)
        if inv.a_loss != 0.0
            JuMP.set_start_value(im, _IMAG_START)   # keep iterates off the I=0 degeneracy
            JuMP.@constraint(m, im^2 == isq)        # exact: no loss inflation
        else
            JuMP.@constraint(m, im^2 >= isq)        # epigraph: unused ⇒ harmless slack
        end
        imag_sum += im

        # Internal EMF magnitude bounds (skipped under grid-forming, which pins the
        # magnitude via v_gfm below). The v_int box applies to every topology; the
        # scalar modulation cap only to SINGLE_PHASE (three-phase topologies use the
        # exact switching-polytope instead).
        if !inv.grid_forming
            vmag2 = JuMP.@expression(m, vrint[k]^2 + viint[k]^2)
            inv.v_int_max !== nothing && JuMP.@constraint(m, vmag2 <= (inv.v_int_max / vb)^2)
            inv.v_int_min !== nothing && JuMP.@constraint(m, vmag2 >= (inv.v_int_min / vb)^2)
            if !is_3ph_topo && inv.modulation_max !== nothing && inv.v_dc !== nothing
                cap = inv.modulation_max * inv.v_dc / _SQRT3
                JuMP.@constraint(m, vmag2 <= (cap / vb)^2)
            end
        end
    end

    # Optional grid-side shunt susceptance at the POC (phase-to-neutral). Its
    # reactive injection is part of the total grid-side POC exchange.
    if inv.b_filter_shunt != 0.0
        b = inv.b_filter_shunt * zb          # S → per-unit susceptance
        for ph in phases
            if neutral === nothing
                vrp = vr[(bus, ph)]; vip = vi[(bus, ph)]
            else
                vrp = JuMP.@expression(m, vr[(bus, ph)] - vr[(bus, neutral)])
                vip = JuMP.@expression(m, vi[(bus, ph)] - vi[(bus, neutral)])
            end
            JuMP.add_to_expression!(ctx.kcl_r[(bus, ph)],  b, vip)
            JuMP.add_to_expression!(ctx.kcl_i[(bus, ph)], -b, vrp)
            Qpoc += JuMP.@expression(m, b * (vrp^2 + vip^2))
            if neutral !== nothing
                JuMP.add_to_expression!(ctx.kcl_r[(bus, neutral)], -b, vip)
                JuMP.add_to_expression!(ctx.kcl_i[(bus, neutral)],  b, vrp)
            end
        end
    end

    # Grid-forming balanced 120° internal EMF (three-phase only), magnitude v_gfm.
    if inv.grid_forming && nph == 3
        h = _SQRT3 / 2                       # a = −½ + j√3/2 ; V_b=a²V_a, V_c=aV_a
        JuMP.@constraint(m, vrint[2] == -0.5*vrint[1] + h*viint[1])
        JuMP.@constraint(m, viint[2] == -h*vrint[1]  - 0.5*viint[1])
        JuMP.@constraint(m, vrint[3] == -0.5*vrint[1] - h*viint[1])
        JuMP.@constraint(m, viint[3] ==  h*vrint[1]  - 0.5*viint[1])
    end
    if inv.grid_forming
        vgfm = JuMP.@variable(m, base_name = "vgfm_$(inv.id)", lower_bound = 0.0)
        inv.v_int_min !== nothing && JuMP.@constraint(m, vgfm >= inv.v_int_min / vb)
        vmax_gfm = inv.v_int_max === nothing ? nothing : inv.v_int_max / vb
        if !is_3ph_topo && inv.modulation_max !== nothing && inv.v_dc !== nothing
            cap = inv.modulation_max * inv.v_dc / _SQRT3 / vb
            vmax_gfm = vmax_gfm === nothing ? cap : min(vmax_gfm, cap)
        end
        vmax_gfm !== nothing && JuMP.@constraint(m, vgfm <= vmax_gfm)
        JuMP.@constraint(m, vrint[1]^2 + viint[1]^2 == vgfm^2)
        JuMP.set_start_value(vgfm, _start_or(vrint[1], 1.0))
    end

    # ── Three-phase switching-polytope feasibility (exact, time-sampled) ────────
    in_re = nothing; in_im = nothing; dre = nothing; dim = nothing
    if is_3ph_topo
        w = 2pi * inv.f
        v_dc = inv.v_dc; c_dc = inv.c_dc
        # Converter output (pole) voltages and 2ω oscillating power, in SI.
        UreSI = [JuMP.@expression(m, vrint[k]*vb) for k in 1:3]
        UimSI = [JuMP.@expression(m, viint[k]*vb) for k in 1:3]
        S2re_SI = JuMP.@expression(m, rip_re * sb)
        S2im_SI = JuMP.@expression(m, rip_im * sb)

        # 2ω bus-ripple phasor D = j·S2/(2ω·C_eq·Vdc); split link has C_eq=C/2.
        denom = inv.topology == :SPLIT_DC ? (w * c_dc * v_dc) : (2w * c_dc * v_dc)
        dre = JuMP.@variable(m, base_name = "dre_$(inv.id)")
        dim = JuMP.@variable(m, base_name = "dim_$(inv.id)")
        JuMP.@constraint(m, dre == -S2im_SI / denom)
        JuMP.@constraint(m, dim ==  S2re_SI / denom)
        inv.dv2_max !== nothing && JuMP.@constraint(m, dre^2 + dim^2 <= inv.dv2_max^2)

        # Split-link midpoint fundamental ripple phasor N (RMS), merged into W_x.
        if inv.topology == :SPLIT_DC
            NreSI = JuMP.@expression(m,  sum(cii[k] for k in 1:3)*ib / (2w * c_dc))
            NimSI = JuMP.@expression(m, -sum(cri[k] for k in 1:3)*ib / (2w * c_dc))
        end

        # Sampled instantaneous rail constraints (linear per sample in dre/dim).
        for th in _sample_grid(inv.n_samples)
            c1, s1 = cos(th), sin(th)
            c2, s2 = cos(2th), sin(2th)
            rail = JuMP.@expression(m, inv.m_max * (v_dc + dre*c2 - dim*s2))
            if inv.topology == :THREE_LEG || inv.topology == :FOUR_LEG
                for (a, b) in _PAIRS_IDX, sg in (1.0, -1.0)
                    JuMP.@constraint(m,
                        sg*sqrt(2)*((UreSI[a]-UreSI[b])*c1 - (UimSI[a]-UimSI[b])*s1) <= rail)
                end
            end
            if inv.topology == :FOUR_LEG
                for k in 1:3, sg in (1.0, -1.0)
                    JuMP.@constraint(m, sg*sqrt(2)*(UreSI[k]*c1 - UimSI[k]*s1) <= rail)
                end
            end
            if inv.topology == :SPLIT_DC
                railh = JuMP.@expression(m, (inv.m_max/2) * (v_dc + dre*c2 - dim*s2))
                for k in 1:3, sg in (1.0, -1.0)
                    JuMP.@constraint(m,
                        sg*sqrt(2)*((UreSI[k]+NreSI)*c1 - (UimSI[k]+NimSI)*s1) <= railh)
                end
            end
        end

        # Neutral current: 3-wire carries none; 4-wire limited by In_max.
        if inv.topology == :THREE_LEG
            JuMP.@constraint(m, sum(cri[k] for k in 1:3) == 0)
            JuMP.@constraint(m, sum(cii[k] for k in 1:3) == 0)
        else
            sumcr = JuMP.@expression(m, sum(cri[k] for k in 1:3))
            sumci = JuMP.@expression(m, sum(cii[k] for k in 1:3))
            JuMP.@constraint(m, sumcr^2 + sumci^2 <= (inv.In_max / ib)^2)
        end
        # Neutral current (SI), for reporting: I_n = −Σ I_phase.
        in_re = JuMP.@expression(m, -sum(cri[k] for k in 1:3) * ib)
        in_im = JuMP.@expression(m, -sum(cii[k] for k in 1:3) * ib)
    end

    # Converter apparent-power circle P_conv²+Q_conv² ≤ s_max² (aux keeps it quadratic).
    pv = JuMP.@variable(m, base_name = "pconv_$(inv.id)")
    qv = JuMP.@variable(m, base_name = "qconv_$(inv.id)")
    JuMP.@constraint(m, pv == Pconv)
    JuMP.@constraint(m, qv == Qconv)
    JuMP.@constraint(m, pv^2 + qv^2 <= (inv.s_max / sb)^2)

    # Converter loss and DC-link power (model units): each coeff scaled to per-unit.
    P_loss = JuMP.@expression(m,
        inv.p_loss_fixed/sb + (inv.a_loss/vb)*imag_sum + (inv.c_loss*sb/vb^2)*isq_sum)
    P_dc = JuMP.@expression(m, Pconv + P_loss)

    # Single-phase standalone 2ω ripple magnitude bound (VA → per-unit power).
    if inv.p_ripple_max !== nothing && !is_3ph_topo
        rr = JuMP.@variable(m, base_name = "rip_re_$(inv.id)")
        ri = JuMP.@variable(m, base_name = "rip_im_$(inv.id)")
        JuMP.@constraint(m, rr == rip_re)
        JuMP.@constraint(m, ri == rip_im)
        JuMP.@constraint(m, rr^2 + ri^2 <= (inv.p_ripple_max / sb)^2)
    end

    return _InvHandles(Ppoc, Qpoc, Pconv, Qconv, P_loss, P_dc,
                       vrint, viint, cri, cii, rip_re, rip_im,
                       in_re, in_im, dre, dim, sb, vb, ib)
end

stamp_device!(ctx, inverter::AdvancedInverter; kwargs...) =
    _stamp_inverter!(ctx, inverter)

link_device!(model, inverter::AdvancedInverter, handles, sb, grid::TimeGrid) = nothing

function extract_device(inverter::AdvancedInverter, h::_InvHandles,
                        status::SolveStatus)
    solved = status.publishable
    sb = h.sb; vb = h.vb; ib = h.ib
    scaled(e, scale) = solved ? JuMP.value(e) * scale : NaN
    mag(a, b, scale) = solved ? hypot(JuMP.value(a), JuMP.value(b)) * scale : NaN
    nph = length(inverter.phase_terminals)
    return (
        p_poc=scaled(h.p_poc, sb), q_poc=scaled(h.q_poc, sb),
        p_conv=scaled(h.p_conv, sb), q_conv=scaled(h.q_conv, sb),
        p_loss=scaled(h.p_loss, sb), p_dc=scaled(h.p_dc, sb),
        v_int_mag=[mag(h.vrint[k], h.viint[k], vb) for k in 1:nph],
        i_mag=[mag(h.cri[k], h.cii[k], ib) for k in 1:nph],
        i_neutral=h.in_re === nothing ? 0.0 :
            (solved ? hypot(JuMP.value(h.in_re), JuMP.value(h.in_im)) : NaN),
        ripple=mag(h.ripple_re, h.ripple_im, sb),
        dv2=h.dre === nothing ? 0.0 :
            (solved ? hypot(JuMP.value(h.dre), JuMP.value(h.dim)) : NaN),
    )
end

"""
    solve_advanced_inverter(net, inverter; objective=:max_export, kwargs...)
        -> InverterResult

Stamp `inverter` into `net` and solve, demonstrating the prototype internal-node
model. The network supplies the surrounding grid (a voltage source, lines, any
loads); the inverter injects at its POC bus.

# Objective
- `:max_export` — maximise the active power delivered to the grid at the POC.
- `:min_loss` — minimise converter loss subject to `p_set` active-power delivery.

# Keywords
- `p_set=nothing` — required active-power target (W) for `:min_loss`.
- `q_set=nothing` — optional reactive-power constraint at the POC (var).
- `per_unit=false`, `s_base=1e6`, `optimizer=Ipopt.Optimizer`, `verbose=false`,
  `solver_options=()`. Per-unit results are returned in SI regardless.
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
    p_set === nothing || isfinite(p_set) || throw(ArgumentError("p_set must be finite"))
    q_set === nothing || isfinite(q_set) || throw(ArgumentError("q_set must be finite"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError("s_base must be finite and > 0"))
    validate_device(inverter, (net,); periods=1)

    handles = Ref{_InvHandles}()
    hook! = ctx -> begin
        h = stamp_device!(ctx, inverter)
        handles[] = h
        q_set === nothing || JuMP.@constraint(ctx.model, h.q_poc == q_set / h.sb)
        if objective == :max_export
            JuMP.@objective(ctx.model, Max, h.p_poc)
        else
            JuMP.@constraint(ctx.model, h.p_poc == p_set / h.sb)
            JuMP.@objective(ctx.model, Min, h.p_loss)
        end
    end

    ctx = build_opf_model(net; per_unit=per_unit, s_base=s_base,
                          add_objective=false, model_hook! = hook!,
                          optimizer=optimizer, verbose=verbose)
    _set_solver_options!(ctx.model, solver_options)
    enforce_kcl!(ctx)
    JuMP.optimize!(ctx.model)

    outcome = _solve_outcome(ctx.model)
    status = string(outcome.termination_status)
    h = handles[]
    device_result = extract_device(inverter, h, SolveStatus(outcome))

    result = _extract_result(ctx, outcome)
    return InverterResult(status, inverter.topology,
                          device_result.p_poc, device_result.q_poc,
                          device_result.p_conv, device_result.q_conv,
                          device_result.p_loss, device_result.p_dc,
                          device_result.v_int_mag, device_result.i_mag,
                          device_result.i_neutral, device_result.ripple,
                          device_result.dv2, result["bus"], SolveStatus(outcome))
end
