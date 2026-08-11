# Prototype advanced inverter model.
#
# A more detailed inverter-based-resource (IBR) than the BMOPFTools engine's
# built-in current-injection IBR. It combines the "IBR Model Extensions" design
# doc's internal-AC-node structure with sampled three-phase converter
# feasible-region models drawn from ongoing research (fundamental-frequency
# switching-polytope voltage feasibility with a 2ω bus-ripple derating):
#
#   converter ──[Zc]── filter midpoint ──[Zg]── POC bus
#                         │
#                      Rd + Cf
#                         │
#                      neutral
#   The midpoint/grid arm are optional; omitting them recovers the original
#   one-series-primitive internal-node model.
#
# Layered, opt-in features:
#   - Output filter (reduced series primitive or explicit damped LCL + POC shunt).
#   - Internal EMF magnitude box, or a sampled switching-polytope DC-utilisation
#     limit for a chosen three-phase topology (see below).
#   - Grid-forming balanced 120° internal EMF (bounded magnitude decision var).
#   - Converter losses: non-branching P_dc = P_conv + P_loss.
#   - Double-frequency (2ω) power ripple, either a standalone cap (single-phase)
#     or a bus-ripple phasor that endogenously derates the DC rails (three-phase).
#
# Three-phase topologies (`topology`), with time-sampled switching-polytope
# voltage feasibility on the converter output U_x = internal-node voltage V_int:
#   :THREE_LEG  3-leg 3-wire  — pairwise line-to-line rails, no neutral current.
#   :FOUR_LEG   4-leg 4-wire  — pairwise + per-phase rails, neutral via the 4th leg.
#   :SPLIT_DC   split-cap 4-wire — half-bus per-phase rails + midpoint ripple.
# `:SINGLE_PHASE` (default) keeps the simpler internal-EMF model.
#
# Built on the BMOPFTools staged API through a `model_hook!`; it does not modify
# the engine. Works in SI (per_unit=false) or per-unit (per_unit=true): device
# parameters are SI and scaled to model units via `opf_bases(ctx)` (Vdc/Cdc/In_max
# stay SI, the AC↔DC coupling scales through the POC bus's v_base/i_base/s_base).

const _SQRT3 = sqrt(3.0)
const _HALF_SQRT3 = _SQRT3 / 2
# Physical smoothing scale for the current-linear loss norm. The shifted smooth
# norm sqrt(|I|² + ε²) - ε is zero at I=0 and differs from |I| by less than ε.
# It removes the singular Jacobian of im²=|I|² at zero current without creating
# an epigraph that an unrelated objective could inflate.
const _IMAG_EPS_SI = 1e-3  # A
const _PAIRS_IDX = ((1, 2), (2, 3), (3, 1))
_sample_grid(N::Int) = [2pi * (k - 1) / N for k in 1:N]

# Components of the unconjugated phasor sum sum(V_x I_x). Kept as a small,
# numeric-generic helper so the source expression and literature-derived unit
# tests use the same sign/phase convention (RMS fundamental phasors).
function _ripple_components(vre, vim, ire, iim)
    idx = eachindex(vre, vim, ire, iim)
    return (sum(vre[k]*ire[k] - vim[k]*iim[k] for k in idx),
            sum(vre[k]*iim[k] + vim[k]*ire[k] for k in idx))
end

"""
    solve_advanced_inverter(net, inverter; objective=:max_export, kwargs...)
        -> InverterResult

Solve the experimental circuit-aware inverter model. When
`inverter.pwm_strategy != :NONE`, the solver alternates between the smooth NLP
and a carrier-level switching audit. The capacitor reserve is increased until
it covers the quadrature sum of the modeled carrier ripple and the residual
manual `i_sw` allowance. This is a conservative sequential closure, not a
mixed-integer switched-converter optimisation.

PWM closure keywords are `pwm_tolerance=1e-3` (relative current tolerance),
`pwm_max_iterations=8`, and `pwm_reserve_factor=1.01`. Inspect
`result.pwm_reserve_margin`, `result.pwm_modulation_margin`, and
`result.pwm_iterations` for every publishable PWM-enabled solve.
"""
function solve_advanced_inverter(net::Dict{String,Any}, inverter::AdvancedInverter;
                                 objective::Symbol=:max_export,
                                 p_set::Union{Float64,Nothing}=nothing,
                                 q_set::Union{Float64,Nothing}=nothing,
                                 per_unit::Bool=false,
                                 s_base::Float64=1e6,
                                 optimizer=Ipopt.Optimizer,
                                 verbose::Bool=false,
                                 solver_options=(),
                                 pwm_tolerance::Float64=1e-3,
                                 pwm_max_iterations::Int=8,
                                 pwm_reserve_factor::Float64=1.01)
    validate_device(inverter, (net,); periods=1)
    isfinite(pwm_tolerance) && pwm_tolerance >= 0 || throw(ArgumentError(
        "pwm_tolerance must be finite and >= 0"))
    pwm_max_iterations >= 1 || throw(ArgumentError(
        "pwm_max_iterations must be >= 1"))
    isfinite(pwm_reserve_factor) && pwm_reserve_factor >= 1 || throw(ArgumentError(
        "pwm_reserve_factor must be finite and >= 1"))

    common = (; objective, p_set, q_set, per_unit, s_base, optimizer,
               verbose, solver_options)
    if inverter.pwm_strategy == :NONE
        return _solve_advanced_inverter_once(net, inverter; common...)
    end

    reserve = inverter.i_sw
    last_result = nothing
    for iteration in 1:pwm_max_iterations
        working = _with_pwm_reserve(inverter, reserve)
        result = _solve_advanced_inverter_once(net, working; common...,
            pwm_audit_inverter=inverter, pwm_iterations=iteration)
        last_result = result
        solve_status(result).publishable || return result
        isfinite(result.i_cap_switching) || return result
        target = hypot(inverter.i_sw, result.i_cap_switching)
        tolerance = pwm_tolerance * max(target, 1.0)
        result.i_cap_switching_reserved + tolerance >= target && return result
        reserve = max(reserve, pwm_reserve_factor * target)
    end
    return last_result
end

# Fortescue current components for the terminal order a-b-c. Each component is
# returned as a (real, imaginary) pair so the helper remains generic over JuMP
# expressions as well as ordinary numbers.
function _sequence_components(ire, iim)
    length(ire) == length(iim) == 3 || throw(ArgumentError(
        "symmetrical components require exactly three phase currents"))
    ar, br, cr = ire; ai, bi, ci = iim; h = _HALF_SQRT3
    i0 = ((ar + br + cr) / 3, (ai + bi + ci) / 3)
    i1 = ((ar - 0.5br - h*bi - 0.5cr + h*ci) / 3,
          (ai + h*br - 0.5bi - h*cr - 0.5ci) / 3)
    i2 = ((ar - 0.5br + h*bi - 0.5cr - h*ci) / 3,
          (ai - h*br - 0.5bi + h*cr - 0.5ci) / 3)
    return i0, i1, i2
end

# Primitive conductor-domain impedance. Conductor order is the phase-terminal
# order followed by the physical neutral, when present. The returned phase
# drops are phase-to-neutral drops, so the neutral-conductor primitive drop is
# subtracted after applying Z to [I_phase; I_neutral], with I_neutral=-sum(I_phase).
function _series_filter_matrices(inv, nph::Int, zb::Real; side::Symbol=:converter)
    side in (:converter, :grid) || throw(ArgumentError("unknown filter side :$side"))
    has_neutral = inv.neutral !== nothing
    ncond = nph + has_neutral
    if side == :converter
        rp, xp = inv.r_filter, inv.x_filter
        rn, xn = inv.r_filter_neutral, inv.x_filter_neutral
        rmat, xmat = inv.r_filter_matrix, inv.x_filter_matrix
    else
        rp, xp = inv.r_filter_grid, inv.x_filter_grid
        rn, xn = inv.r_filter_grid_neutral, inv.x_filter_grid_neutral
        rmat, xmat = inv.r_filter_grid_matrix, inv.x_filter_grid_matrix
    end
    matrix_mode = rmat !== nothing || xmat !== nothing
    if matrix_mode
        R = rmat === nothing ? zeros(ncond, ncond) : rmat
        X = xmat === nothing ? zeros(ncond, ncond) : xmat
        return R ./ zb, X ./ zb
    end
    R = zeros(ncond, ncond); X = zeros(ncond, ncond)
    for k in 1:nph
        R[k, k] = rp / zb
        X[k, k] = xp / zb
    end
    if has_neutral
        R[end, end] = rn / zb
        X[end, end] = xn / zb
    end
    return R, X
end

_lcl_active(inv) = inv.c_filter_mid != 0.0 || inv.r_filter_damping != 0.0 ||
    inv.r_filter_grid != 0.0 || inv.x_filter_grid != 0.0 ||
    inv.r_filter_grid_neutral != 0.0 || inv.x_filter_grid_neutral != 0.0 ||
    inv.r_filter_grid_matrix !== nothing || inv.x_filter_grid_matrix !== nothing

# Fundamental-frequency admittance of a capacitor C in series with damping R.
# Positive B means a capacitive shunt current I=(G+jB)V flowing phase-to-neutral.
function _mid_filter_admittance(inv)
    inv.c_filter_mid == 0.0 && return (0.0, 0.0)
    xc = 1 / (2pi * inv.f * inv.c_filter_mid)
    den = inv.r_filter_damping^2 + xc^2
    return inv.r_filter_damping / den, xc / den
end

# Undamped scalar LCL estimate. Matrix-valued or incomplete cases have no unique
# scalar resonance and deliberately report NaN rather than a misleading number.
function _filter_resonance_hz(inv)
    (!_lcl_active(inv) || inv.c_filter_mid <= 0 || inv.x_filter <= 0 ||
     inv.x_filter_grid <= 0 || inv.r_filter_matrix !== nothing ||
     inv.x_filter_matrix !== nothing || inv.r_filter_grid_matrix !== nothing ||
     inv.x_filter_grid_matrix !== nothing) && return NaN
    l1 = inv.x_filter / (2pi * inv.f)
    l2 = inv.x_filter_grid / (2pi * inv.f)
    return sqrt((l1 + l2) / (l1 * l2 * inv.c_filter_mid)) / (2pi)
end

function _filter_voltage_drops(R, X, ire, iim, nph::Int, has_neutral::Bool)
    ncond = nph + has_neutral
    length(ire) == length(iim) == nph || throw(DimensionMismatch(
        "phase-current vectors must have length $nph"))
    size(R) == size(X) == (ncond, ncond) || throw(DimensionMismatch(
        "filter matrices must both be $ncond×$ncond"))
    cre = has_neutral ? [ire; -sum(ire)] : collect(ire)
    cim = has_neutral ? [iim; -sum(iim)] : collect(iim)
    dr = [sum(R[k, j]*cre[j] - X[k, j]*cim[j] for j in 1:ncond)
          for k in 1:ncond]
    di = [sum(R[k, j]*cim[j] + X[k, j]*cre[j] for j in 1:ncond)
          for k in 1:ncond]
    if has_neutral
        return [dr[k] - dr[end] for k in 1:nph],
               [di[k] - di[end] for k in 1:nph]
    end
    return dr, di
end

# A dense post-solve audit of the continuous-time switching inequalities. This
# is deliberately separate from the NLP's sample grid: it detects between-grid
# violations but is not presented as an analytic global certificate.
function _switching_margin(inv, ure, uim, dre, dim, nre, nim, nmean=0.0;
                           n_eval::Int=max(3600, 100*inv.n_samples))
    inv.topology in _THREE_PHASE_TOPOLOGIES || return 0.0
    length(ure) == length(uim) == 3 || throw(DimensionMismatch(
        "switching audit requires three voltage phasors"))
    n_eval >= 4 || throw(ArgumentError("switching audit requires n_eval >= 4"))
    dvr = dre === nothing ? 0.0 : dre
    dvi = dim === nothing ? 0.0 : dim
    nr = nre === nothing ? 0.0 : nre
    ni = nim === nothing ? 0.0 : nim
    margin = Inf
    for th in _sample_grid(n_eval)
        c1, s1 = cos(th), sin(th)
        c2, s2 = cos(2th), sin(2th)
        rail = inv.m_max * (inv.v_dc + dvr*c2 - dvi*s2)
        if inv.topology in (:THREE_LEG, :FOUR_LEG)
            for (a, b) in _PAIRS_IDX
                vll = sqrt(2)*((ure[a]-ure[b])*c1 - (uim[a]-uim[b])*s1)
                margin = min(margin, rail - abs(vll))
            end
        end
        if inv.topology == :FOUR_LEG
            for k in 1:3
                vpn = sqrt(2)*(ure[k]*c1 - uim[k]*s1)
                margin = min(margin, rail - abs(vpn))
            end
        elseif inv.topology == :SPLIT_DC
            railh = rail / 2
            for k in 1:3
                vpn = sqrt(2)*((ure[k]+nr)*c1 - (uim[k]+ni)*s1) + nmean
                margin = min(margin, railh - abs(vpn))
            end
        end
    end
    return margin
end

# Post-solve carrier-level PWM audit. Fundamental voltages/currents are frozen
# over each switching period, as in Mandrioli et al. (Energies 2021, 14, 1434).
# The shared triangular carrier is integrated explicitly, which preserves the
# switching-function correlation that independent-leg RMS approximations lose.
function _pwm_ripple_audit(inv, ure, uim, ire, iim, dre=0.0, dim=0.0,
                           nre=0.0, nim=0.0, nmean=0.0)
    inv.pwm_strategy == :NONE && return (
        i_rms=0.0, dv_rms=0.0, dv_pp=0.0, modulation_margin=NaN)
    length(ure) == length(uim) == length(ire) == length(iim) == 3 ||
        throw(DimensionMismatch("PWM ripple audit requires three voltage and current phasors"))

    nf = inv.pwm_fundamental_samples
    nc = inv.pwm_carrier_samples
    fsw = inv.f_sw
    ceq = if inv.topology == :SPLIT_DC
        cu, cl = _split_capacitances(inv)
        cu * cl / (cu + cl)
    else
        inv.c_dc
    end
    i_rms_sq = 0.0
    dv_rms_sq = 0.0
    dv_pp = 0.0
    modulation_margin = Inf
    ihat = zeros(nc)
    vhat = zeros(nc)

    for th in _sample_grid(nf)
        c1, s1 = cos(th), sin(th)
        c2, s2 = cos(2th), sin(2th)
        rail = inv.v_dc + dre*c2 - dim*s2
        rail > 0 || return (i_rms=NaN, dv_rms=NaN, dv_pp=NaN,
                            modulation_margin=-Inf)
        uph = [sqrt(2)*(ure[k]*c1 - uim[k]*s1) for k in 1:3]
        iph = [sqrt(2)*(ire[k]*c1 - iim[k]*s1) for k in 1:3]

        commands, currents = if inv.topology == :SPLIT_DC
            nr = sqrt(2)*(nre*c1 - nim*s1) + nmean
            uph .+ nr, iph
        elseif inv.topology == :FOUR_LEG
            base = [uph; 0.0]
            gamma = inv.pwm_strategy == :CENTERED ?
                -(maximum(base) + minimum(base))/2 : 0.0
            base .+ gamma, [iph; -sum(iph)]
        else # :THREE_LEG
            gamma = inv.pwm_strategy == :CENTERED ?
                -(maximum(uph) + minimum(uph))/2 : 0.0
            uph .+ gamma, iph
        end

        modulation_margin = min(modulation_margin,
            rail/2 - maximum(abs, commands))
        refs = commands ./ rail                 # carrier spans [-1/2, +1/2]
        duties = 0.5 .+ refs
        if minimum(duties) < -1e-9 || maximum(duties) > 1 + 1e-9
            return (i_rms=NaN, dv_rms=NaN, dv_pp=NaN,
                    modulation_margin=modulation_margin)
        end
        iavg = sum(duties .* currents)
        for j in 1:nc
            tau = (j - 0.5) / nc
            carrier = tau < 0.5 ? -0.5 + 2tau : 1.5 - 2tau
            ihat[j] = sum((refs[k] >= carrier ? 1.0 : 0.0) * currents[k]
                          for k in eachindex(refs)) - iavg
        end

        # Remove the finite-grid mean before integration. In the continuous
        # carrier limit this correction is zero; here it prevents numerical
        # charge drift from contaminating the voltage-ripple RMS.
        ihat .-= sum(ihat) / nc
        i_rms_sq += sum(abs2, ihat) / nc
        q = 0.0
        dt = 1 / (fsw * nc)
        for j in 1:nc
            q -= ihat[j] * dt
            vhat[j] = q / ceq
        end
        vhat .-= sum(vhat) / nc
        dv_rms_sq += sum(abs2, vhat) / nc
        dv_pp = max(dv_pp, maximum(vhat) - minimum(vhat))
    end
    return (i_rms=sqrt(i_rms_sq / nf), dv_rms=sqrt(dv_rms_sq / nf),
            dv_pp=dv_pp, modulation_margin=modulation_margin)
end

"""
    AdvancedInverter(; id, bus, s_max, kwargs...)

An experimental inverter with an internal AC node behind an output filter. All
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
  time-sampled switching-polytope voltage feasibility on the internal-node
  voltage and require `v_dc` and `c_dc`. The 4-wire topologies additionally need
  their neutral current bounded: `:FOUR_LEG` requires `In_max` (the 4th-leg
  rating), `:SPLIT_DC` requires `In_max` **or** `i_cap_max` (its neutral flows
  through the capacitors, so a bank rating bounds it on its own).
- `v_dc` — DC-link voltage (V). `c_dc` — DC-link capacitance (F; per half for split).
- `c_dc_upper`, `c_dc_lower` — optional split-link half-bank capacitances (F).
  Each defaults to `c_dc`; supplying either exposes unequal series capacitance,
  unequal neutral-current sharing, and a mean midpoint offset.
- `m_max=1.0` — utilisation factor applied to the ideal two-level switching
  hull (dimensionless, `0 < m_max ≤ 1`).
- `In_max` — neutral current limit (A). For `:FOUR_LEG` this is the 4th leg's
  device (thermal) rating. For `:SPLIT_DC` it is the half-bank capacitor ripple
  rating **net of the 2ω allocation**: the same capacitors carry both the
  fundamental neutral current and the 2ω bus current, which sit at different
  frequencies and so combine in RMS. For a symmetric link,
  `In_max = 2·√(I_half_rated² − I_2ω,rms² − I_sw,rms²)` with
  `I_2ω,rms = |S̃|/(√2·v_dc)`. Passing the raw bank rating double-counts the
  capacitors and overstates neutral capability. Prefer `i_cap_max`, which makes
  that allocation endogenous instead of asking you to pre-compute it.
- `i_cap_max` — capacitor RMS ripple-current rating (A), three-phase topologies
  only: **per half-bank** for `:SPLIT_DC`, the whole DC-link bank otherwise.
  Supplying it stamps the simultaneous thermal allocation on both split
  half-banks. In the symmetric, equal-weight case this is
  `(|I_n|/2)² + I_2ω,rms² + i_sw² ≤ i_cap_max²`; unequal banks use their
  capacitance-proportional neutral shares. The neutral term is present only for
  `:SPLIT_DC`; the 4-leg's neutral flows through its 4th leg, not the caps. Thus
  the split between capacitor-rating and neutral-current limits is decided by the
  solve at the operating point rather than fixed as a nameplate. Composes with
  `In_max` — whichever binds, binds. For `:SPLIT_DC` it can also be supplied
  *instead of* `In_max`, since it bounds `|I_n|` on its own.
  **Rating convention:** because electrolytic ESR falls with frequency, the
  equal-weight sum above is only conservative if `i_cap_max` is referred to the
  *lowest* frequency in it — the 50/60 Hz neutral term. Refer the datasheet
  rating to that frequency using the capacitor manufacturer's own multiplier;
  passing a raw 100/120 Hz rating understates 50/60 Hz heating.
- `i_cap_upper_max`, `i_cap_lower_max` — optional per-half split-link thermal-
  equivalent current ratings (A). They compose with the common `i_cap_max`.
- `cap_thermal_weights=(1,1,1)` — squared-current weights for the fundamental
  neutral, 2ω, and switching components, respectively. Use ESR ratios referred
  to the rating frequency; the default reproduces the original RMS sum.
- `esr_dc` — monolithic-bank reference ESR (Ω). `esr_dc_upper` and
  `esr_dc_lower` are the split half-bank values. Together with
  `cap_thermal_weights`, these produce capacitor dissipation and add it to
  `p_dc`; zero defaults preserve the lossless-link model.
- `q_mid_balance_max` — optional split-link charge-transfer authority (C). The
  actuator can shift mean midpoint voltage by `2q/(C_upper+C_lower)`.
- `v_mid_mean_max` — optional magnitude limit on mean midpoint offset (V).
- `i_sw` — optional constant switching-frequency RMS allowance (A) reserved out
  of a supplied common or half-bank capacitor rating (use when the
  electrolytics, not a parallel film cap, carry the `f_sw` component).
- `pwm_strategy=:NONE` — optional carrier-level DC-link ripple model. `:SPWM`
  uses zero common-mode injection; `:CENTERED` uses centered carrier PWM (the
  carrier-based SVM equivalent). The convenience solver closes the resulting
  capacitor-current reserve by conservative outer iteration; direct device
  stamping continues to require a fixed `i_sw`.
- `f_sw` — switching frequency (Hz), required when `pwm_strategy != :NONE`.
  `pwm_fundamental_samples` and `pwm_carrier_samples` set the two numerical
  quadrature grids used by the post-solve switching-state audit.
- `dv2_max` — optional cap on the 2ω bus-ripple amplitude (V).
- `dv_mid_max` — `:SPLIT_DC` only: optional cap on the RMS fundamental
  midpoint-to-ideal-midpoint voltage ripple (V).
- `i_zero_max`, `i_positive_max`, `i_negative_max` — optional RMS symmetrical-
  component current limits (A). These are useful for grid-code limits and for
  separating zero-sequence neutral stress from negative-sequence 2ω stress.
- `n_samples=36` — time-sampling grid for an outer approximation of continuous-
  time voltage feasibility. Increase it for boundary-sensitive studies.
- `f=50.0` — fundamental frequency (Hz).

# Output filter
- `r_filter=0.0`, `x_filter=0.0` — series filter impedance per phase (Ω).
- `r_filter_neutral=0.0`, `x_filter_neutral=0.0` — neutral-conductor series
  impedance (Ω). Its common voltage drop couples all phase-to-neutral KVL
  equations; it requires a physical `neutral` terminal.
- `r_filter_matrix=nothing`, `x_filter_matrix=nothing` — optional primitive
  conductor-domain series-impedance matrices (Ω), ordered as `phase_terminals`
  followed by `neutral` when present. They capture unequal conductors and mutual
  coupling. Supplying either matrix supersedes the scalar impedances; mixing the
  two parameterisations is rejected. Matrices must be finite and symmetric, and
  the resistance matrix must be positive semidefinite.
- `r_filter_grid`, `x_filter_grid` and their `_neutral` / `_matrix` variants —
  optional grid-side series arm. Supplying any grid-side arm or `c_filter_mid`
  activates an explicit LCL midpoint; the original filter becomes the
  converter-side arm.
- `c_filter_mid=0.0` — per-phase midpoint filter capacitance (F).
  `r_filter_damping=0.0` is its optional series damping resistance (Ω).
- `i_grid_max=nothing` — optional grid-side arm current limit (A), distinct from
  converter-side `i_max` when the midpoint capacitor carries current.
- `b_filter_shunt=0.0` — additional grid-side (POC) shunt susceptance (S), kept
  separate from the LCL midpoint capacitor for backward compatibility.

# Internal EMF box / single-phase modulation
- `v_int_min`, `v_int_max` — per-phase EMF magnitude box (V; applies to all topologies).
- `modulation_max` — SINGLE_PHASE only: legacy scalar DC-link convention
  `|V_int| ≤ modulation_max·v_dc/√3`. This is not a bridge-specific
  full-/half-bridge switching model.

# Grid-forming
- `grid_forming=false` — balanced positive-sequence internal EMF; the magnitude
  `v_gfm ∈ [v_int_min, v_int_max]` is a decision variable (composes with a topology).

# Converter losses
- `p_loss_fixed=0.0` (W), `a_loss=0.0` (W/A), `c_loss=0.0` (W/A²).

# Double-frequency ripple / current
- `p_ripple_max` — SINGLE_PHASE only: bound on the 2ω power-ripple amplitude
  (VA); it does not by itself calculate capacitor voltage or RMS current.
- `i_max=nothing` — optional per-phase-conductor current limit (A). Strongly
  recommended for every three-phase topology: the aggregate `s_max` circle does
  not bound phase-redistribution currents whose per-phase powers cancel.
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
    c_dc_upper::Union{Float64,Nothing} = nothing
    c_dc_lower::Union{Float64,Nothing} = nothing
    m_max::Float64 = 1.0
    In_max::Union{Float64,Nothing} = nothing
    i_cap_max::Union{Float64,Nothing} = nothing
    i_cap_upper_max::Union{Float64,Nothing} = nothing
    i_cap_lower_max::Union{Float64,Nothing} = nothing
    cap_thermal_weights::NTuple{3,Float64} = (1.0, 1.0, 1.0)
    esr_dc::Float64 = 0.0
    esr_dc_upper::Float64 = 0.0
    esr_dc_lower::Float64 = 0.0
    i_sw::Float64 = 0.0
    pwm_strategy::Symbol = :NONE
    f_sw::Union{Float64,Nothing} = nothing
    pwm_fundamental_samples::Int = 360
    pwm_carrier_samples::Int = 256
    dv2_max::Union{Float64,Nothing} = nothing
    dv_mid_max::Union{Float64,Nothing} = nothing
    q_mid_balance_max::Union{Float64,Nothing} = nothing
    v_mid_mean_max::Union{Float64,Nothing} = nothing
    i_zero_max::Union{Float64,Nothing} = nothing
    i_positive_max::Union{Float64,Nothing} = nothing
    i_negative_max::Union{Float64,Nothing} = nothing
    n_samples::Int = 36
    f::Float64 = 50.0
    r_filter::Float64 = 0.0
    x_filter::Float64 = 0.0
    r_filter_neutral::Float64 = 0.0
    x_filter_neutral::Float64 = 0.0
    r_filter_matrix::Union{Matrix{Float64},Nothing} = nothing
    x_filter_matrix::Union{Matrix{Float64},Nothing} = nothing
    r_filter_grid::Float64 = 0.0
    x_filter_grid::Float64 = 0.0
    r_filter_grid_neutral::Float64 = 0.0
    x_filter_grid_neutral::Float64 = 0.0
    r_filter_grid_matrix::Union{Matrix{Float64},Nothing} = nothing
    x_filter_grid_matrix::Union{Matrix{Float64},Nothing} = nothing
    c_filter_mid::Float64 = 0.0
    r_filter_damping::Float64 = 0.0
    i_grid_max::Union{Float64,Nothing} = nothing
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

function _with_pwm_reserve(inv::AdvancedInverter, reserve::Float64)
    names = fieldnames(AdvancedInverter)
    values = Tuple(getfield(inv, name) for name in names)
    parameters = NamedTuple{names}(values)
    return AdvancedInverter(; merge(parameters,
        (i_sw=reserve, pwm_strategy=:NONE))...)
end

const _THREE_PHASE_TOPOLOGIES = (:THREE_LEG, :FOUR_LEG, :SPLIT_DC)

_split_capacitances(inv) =
    (something(inv.c_dc_upper, inv.c_dc), something(inv.c_dc_lower, inv.c_dc))

function _validate_inverter(inv::AdvancedInverter, nets=())
    _validate_connection(inv.id, inv.bus, inv.phase_terminals, inv.neutral, nets)
    inv.topology in (:SINGLE_PHASE, _THREE_PHASE_TOPOLOGIES...) ||
        throw(ArgumentError("unknown topology :$(inv.topology)"))
    inv.pwm_strategy in (:NONE, :SPWM, :CENTERED) || throw(ArgumentError(
        "inverter '$(inv.id)' pwm_strategy must be :NONE, :SPWM, or :CENTERED"))

    isfinite(inv.s_max) && inv.s_max > 0 || throw(ArgumentError(
        "inverter '$(inv.id)' s_max must be finite and > 0"))
    for (name, value) in (("r_filter", inv.r_filter), ("x_filter", inv.x_filter),
                          ("r_filter_neutral", inv.r_filter_neutral),
                          ("x_filter_neutral", inv.x_filter_neutral),
                          ("r_filter_grid", inv.r_filter_grid),
                          ("x_filter_grid", inv.x_filter_grid),
                          ("r_filter_grid_neutral", inv.r_filter_grid_neutral),
                          ("x_filter_grid_neutral", inv.x_filter_grid_neutral),
                          ("c_filter_mid", inv.c_filter_mid),
                          ("r_filter_damping", inv.r_filter_damping),
                          ("esr_dc", inv.esr_dc),
                          ("esr_dc_upper", inv.esr_dc_upper),
                          ("esr_dc_lower", inv.esr_dc_lower),
                          ("b_filter_shunt", inv.b_filter_shunt),
                          ("m_max", inv.m_max), ("f", inv.f),
                          ("p_loss_fixed", inv.p_loss_fixed),
                          ("a_loss", inv.a_loss), ("c_loss", inv.c_loss))
        isfinite(value) || throw(ArgumentError(
            "inverter '$(inv.id)' $name must be finite"))
    end
    inv.r_filter >= 0 || throw(ArgumentError(
        "inverter '$(inv.id)' r_filter must be >= 0"))
    inv.r_filter_neutral >= 0 || throw(ArgumentError(
        "inverter '$(inv.id)' r_filter_neutral must be >= 0"))
    inv.r_filter_grid >= 0 || throw(ArgumentError(
        "inverter '$(inv.id)' r_filter_grid must be >= 0"))
    inv.r_filter_grid_neutral >= 0 || throw(ArgumentError(
        "inverter '$(inv.id)' r_filter_grid_neutral must be >= 0"))
    inv.c_filter_mid >= 0 || throw(ArgumentError(
        "inverter '$(inv.id)' c_filter_mid must be >= 0"))
    inv.r_filter_damping >= 0 || throw(ArgumentError(
        "inverter '$(inv.id)' r_filter_damping must be >= 0"))
    inv.r_filter_damping > 0 && inv.c_filter_mid == 0 && throw(ArgumentError(
        "inverter '$(inv.id)' r_filter_damping requires c_filter_mid > 0"))
    all(x -> x >= 0, (inv.esr_dc, inv.esr_dc_upper, inv.esr_dc_lower)) ||
        throw(ArgumentError("inverter '$(inv.id)' DC-link ESR values must be >= 0"))
    all(isfinite, inv.cap_thermal_weights) && all(x -> x >= 0, inv.cap_thermal_weights) ||
        throw(ArgumentError("inverter '$(inv.id)' cap_thermal_weights must be finite and >= 0"))

    ncond = length(inv.phase_terminals) + (inv.neutral === nothing ? 0 : 1)
    function validate_primitive(label, rp, xp, rn, xn, rmat, xmat)
        matrix_mode = rmat !== nothing || xmat !== nothing
        matrix_mode || return
        any(x -> !iszero(x), (rp, xp, rn, xn)) && throw(ArgumentError(
            "inverter '$(inv.id)' $label matrix impedance cannot be mixed with " *
            "scalar phase/neutral impedances"))
        for (name, M) in (("r_$(label)_matrix", rmat), ("x_$(label)_matrix", xmat))
            M === nothing && continue
            size(M) == (ncond, ncond) || throw(ArgumentError(
                "inverter '$(inv.id)' $name must be $ncond×$ncond in conductor order"))
            all(isfinite, M) || throw(ArgumentError(
                "inverter '$(inv.id)' $name must contain only finite values"))
            isapprox(M, transpose(M); rtol=1e-10, atol=1e-12) || throw(ArgumentError(
                "inverter '$(inv.id)' $name must be symmetric"))
        end
        if rmat !== nothing
            R = Symmetric(rmat)
            tol = 1e-10 * max(1.0, opnorm(rmat, Inf))
            eigmin(R) >= -tol || throw(ArgumentError(
                "inverter '$(inv.id)' r_$(label)_matrix must be positive semidefinite"))
        end
    end
    validate_primitive("filter", inv.r_filter, inv.x_filter,
        inv.r_filter_neutral, inv.x_filter_neutral,
        inv.r_filter_matrix, inv.x_filter_matrix)
    validate_primitive("filter_grid", inv.r_filter_grid, inv.x_filter_grid,
        inv.r_filter_grid_neutral, inv.x_filter_grid_neutral,
        inv.r_filter_grid_matrix, inv.x_filter_grid_matrix)
    0 < inv.m_max <= 1 || throw(ArgumentError(
        "inverter '$(inv.id)' m_max must satisfy 0 < m_max <= 1"))
    inv.f > 0 || throw(ArgumentError(
        "inverter '$(inv.id)' f must be > 0"))
    all(x -> x >= 0, (inv.p_loss_fixed, inv.a_loss, inv.c_loss)) ||
        throw(ArgumentError("inverter '$(inv.id)' loss coefficients must be >= 0"))

    for (name, value) in (("i_max", inv.i_max), ("i_grid_max", inv.i_grid_max),
                          ("v_dc", inv.v_dc),
                          ("c_dc", inv.c_dc), ("c_dc_upper", inv.c_dc_upper),
                          ("c_dc_lower", inv.c_dc_lower), ("In_max", inv.In_max),
                          ("i_cap_max", inv.i_cap_max),
                          ("i_cap_upper_max", inv.i_cap_upper_max),
                          ("i_cap_lower_max", inv.i_cap_lower_max),
                          ("f_sw", inv.f_sw),
                          ("dv2_max", inv.dv2_max), ("dv_mid_max", inv.dv_mid_max),
                          ("q_mid_balance_max", inv.q_mid_balance_max),
                          ("v_mid_mean_max", inv.v_mid_mean_max),
                          ("i_zero_max", inv.i_zero_max),
                          ("i_positive_max", inv.i_positive_max),
                          ("i_negative_max", inv.i_negative_max),
                          ("v_int_max", inv.v_int_max),
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
    sequence_limits = (inv.i_zero_max, inv.i_positive_max, inv.i_negative_max)
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
        if inv.pwm_strategy != :NONE
            inv.f_sw !== nothing || throw(ArgumentError(
                "inverter '$(inv.id)' pwm_strategy requires f_sw"))
            inv.f_sw > inv.f || throw(ArgumentError(
                "inverter '$(inv.id)' f_sw must exceed fundamental frequency f"))
            inv.pwm_fundamental_samples >= 12 || throw(ArgumentError(
                "inverter '$(inv.id)' pwm_fundamental_samples must be >= 12"))
            inv.pwm_carrier_samples >= 16 && iseven(inv.pwm_carrier_samples) ||
                throw(ArgumentError("inverter '$(inv.id)' pwm_carrier_samples " *
                                    "must be even and >= 16"))
            inv.topology == :SPLIT_DC && inv.pwm_strategy == :CENTERED &&
                throw(ArgumentError("inverter '$(inv.id)' :SPLIT_DC supports " *
                                    "pwm_strategy=:SPWM only; its midpoint is fixed"))
        end
        # The 4-leg's neutral current flows through its fourth leg, whose device
        # rating is independent of anything on the DC side — so In_max is always
        # required. The split link's neutral flows through the capacitors, so a
        # bank rating (i_cap_max) bounds it on its own; either knob suffices.
        if inv.topology == :FOUR_LEG && inv.In_max === nothing
            throw(ArgumentError("topology :FOUR_LEG requires In_max (the 4th-leg rating)"))
        elseif inv.topology == :SPLIT_DC && inv.In_max === nothing &&
               inv.i_cap_max === nothing && inv.i_cap_upper_max === nothing &&
               inv.i_cap_lower_max === nothing
            throw(ArgumentError("topology :SPLIT_DC requires In_max or a capacitor-current rating " *
                                "to bound the neutral current through the capacitors"))
        end
        if inv.topology == :SPLIT_DC && inv.v_mid_mean_max !== nothing &&
           inv.q_mid_balance_max === nothing
            cu, cl = _split_capacitances(inv)
            natural_mean = inv.v_dc * (cu - cl) / (2(cu + cl))
            abs(natural_mean) <= inv.v_mid_mean_max || throw(ArgumentError(
                "inverter '$(inv.id)' natural mean midpoint offset exceeds " *
                "v_mid_mean_max; supply q_mid_balance_max or relax the limit"))
        end
    elseif inv.modulation_max !== nothing && inv.v_dc === nothing
        throw(ArgumentError(
            "inverter '$(inv.id)' modulation_max requires v_dc"))
    end
    if !is_3ph && inv.pwm_strategy != :NONE
        throw(ArgumentError("inverter '$(inv.id)' pwm_strategy requires a " *
                            "three-phase topology"))
    end
    if !is_3ph && any(x -> x !== nothing, sequence_limits)
        throw(ArgumentError("inverter '$(inv.id)' symmetrical-component current " *
                            "limits require a three-phase topology"))
    end
    if inv.dv_mid_max !== nothing && inv.topology != :SPLIT_DC
        throw(ArgumentError("inverter '$(inv.id)' dv_mid_max applies only to :SPLIT_DC"))
    end
    split_only = (inv.c_dc_upper, inv.c_dc_lower, inv.i_cap_upper_max,
                  inv.i_cap_lower_max, inv.q_mid_balance_max, inv.v_mid_mean_max)
    if inv.topology != :SPLIT_DC && (any(x -> x !== nothing, split_only) ||
                                     inv.esr_dc_upper != 0.0 || inv.esr_dc_lower != 0.0)
        throw(ArgumentError("inverter '$(inv.id)' upper/lower DC-link parameters " *
                            "apply only to :SPLIT_DC"))
    end
    if inv.topology == :SPLIT_DC && inv.esr_dc != 0.0
        throw(ArgumentError("inverter '$(inv.id)' esr_dc is for a monolithic link; " *
                            "use esr_dc_upper/esr_dc_lower for :SPLIT_DC"))
    end
    if inv.neutral === nothing &&
       (inv.r_filter_neutral != 0.0 || inv.x_filter_neutral != 0.0 ||
        inv.r_filter_grid_neutral != 0.0 || inv.x_filter_grid_neutral != 0.0)
        throw(ArgumentError("inverter '$(inv.id)' neutral filter impedance requires " *
                            "a physical neutral terminal"))
    end
    if !is_3ph && (inv.i_cap_max !== nothing || inv.esr_dc != 0.0 ||
                   inv.cap_thermal_weights != (1.0, 1.0, 1.0))
        throw(ArgumentError("inverter '$(inv.id)' DC-bank thermal parameters apply " *
                            "only to the three-phase topologies"))
    end
    isfinite(inv.i_sw) && inv.i_sw >= 0 || throw(ArgumentError(
        "inverter '$(inv.id)' i_sw must be finite and >= 0"))
    cap_ratings = (inv.i_cap_max, inv.i_cap_upper_max, inv.i_cap_lower_max)
    if inv.i_sw > 0 && all(x -> x === nothing, cap_ratings)
        throw(ArgumentError("inverter '$(inv.id)' i_sw is an allowance reserved out " *
                            "of a capacitor rating; supply one or leave i_sw = 0"))
    end
    thermal_sw = sqrt(inv.cap_thermal_weights[3]) * inv.i_sw
    if any(rating -> rating !== nothing && thermal_sw >= rating, cap_ratings)
        throw(ArgumentError("inverter '$(inv.id)' weighted i_sw must be smaller than " *
                            "every supplied capacitor-current rating"))
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
- `p_loss`, `p_cap_loss`, `p_dc` — semiconductor loss, capacitor ESR loss,
  and DC-link power (`p_dc = p_conv + p_loss + p_cap_loss`).
- `p_filter_loss` — total real loss in both series arms and midpoint damping (W).
- `v_int_mag::Vector{Float64}` — internal EMF magnitude per phase.
- `v_filter_mag` — LCL midpoint phase-to-neutral voltage magnitudes (V); equals
  the POC voltage in the reduced single-series-arm model.
- `i_mag`, `i_grid_mag`, `i_filter_shunt_mag` — converter-side, grid-side, and
  midpoint shunt-branch RMS current magnitudes per phase (A).
- `i_neutral` — neutral current magnitude `|I_n|` (0 for 3-wire / single-phase).
- `i_zero`, `i_positive`, `i_negative` — RMS symmetrical-component current
  magnitudes (A; zero for `:SINGLE_PHASE`). For a four-wire connection,
  `i_neutral = 3i_zero` up to numerical tolerance.
- `ripple` — 2ω power-ripple magnitude `|Σ V_int·I|` (VA).
- `dv2` — 2ω bus-ripple amplitude (V; three-phase topologies, else 0).
- `dv_mid` — split-link midpoint fundamental-ripple RMS magnitude (V;
  `:SPLIT_DC`, else 0).
- `v_mid_mean`, `q_mid_balance` — signed mean midpoint shift (V) used by the
  switching hull and balancing charge transfer (C; split link only).
- `i_cap` — **total** capacitor RMS ripple current (A; three-phase topologies,
  else 0): the 2ω component, neutral share, and `i_sw`. For an asymmetric split
  link it is the larger half-bank RMS current.
- `i_cap_thermal` — maximum thermally equivalent current after applying
  `cap_thermal_weights`; this is the quantity constrained by current ratings.
- `i_cap_upper`, `i_cap_lower`, `i_cap_thermal_upper`,
  `i_cap_thermal_lower` — split half-bank physical and thermal-equivalent
  currents (A; zero for monolithic links).
- `i_cap_switching` — carrier-model prediction of DC-link switching-current RMS
  (A); `i_cap_switching_reserved` is the conservative value actually allocated
  by the NLP. `pwm_reserve_margin` is reserved minus the quadrature combination
  of modeled ripple and the user's residual `i_sw` allowance.
- `dv_switching_rms`, `dv_switching_pp` — carrier-model total DC-bus switching
  voltage ripple (V RMS and maximum local peak-to-peak). These are distinct from
  the low-frequency `dv2` amplitude.
- `pwm_modulation_margin` — minimum carrier-reference headroom (V); a negative
  value or `NaN` means the selected PWM strategy cannot realize the solved
  fundamental reference. `pwm_iterations` reports outer reserve solves.
- `switching_margin` — minimum rail headroom (V) on a dense, independent
  post-solve time grid. Positive values pass that audit; negative values reveal
  a between-sample violation accepted by the optimisation grid. This is a useful
  numerical audit, not a formal continuous-time certificate (0 for
  `:SINGLE_PHASE`).
- `filter_resonance_hz` — undamped scalar LCL resonance estimate (Hz), or `NaN`
  when the filter is reduced, matrix-valued, or lacks two positive inductive arms.
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
    p_cap_loss::Float64
    p_dc::Float64
    p_filter_loss::Float64
    v_int_mag::Vector{Float64}
    v_filter_mag::Vector{Float64}
    i_mag::Vector{Float64}
    i_grid_mag::Vector{Float64}
    i_filter_shunt_mag::Vector{Float64}
    i_neutral::Float64
    i_zero::Float64
    i_positive::Float64
    i_negative::Float64
    ripple::Float64
    dv2::Float64
    dv_mid::Float64
    v_mid_mean::Float64
    q_mid_balance::Float64
    i_cap::Float64
    i_cap_thermal::Float64
    i_cap_upper::Float64
    i_cap_lower::Float64
    i_cap_thermal_upper::Float64
    i_cap_thermal_lower::Float64
    i_cap_switching::Float64
    i_cap_switching_reserved::Float64
    pwm_reserve_margin::Float64
    dv_switching_rms::Float64
    dv_switching_pp::Float64
    pwm_modulation_margin::Float64
    pwm_iterations::Int
    switching_margin::Float64
    filter_resonance_hz::Float64
    bus::Dict{String,Any}
    solve::SolveStatus
end

solve_status(result::InverterResult) = result.solve

solve_diagnostics(result::InverterResult) =
    (topology=result.topology, p_loss=result.p_loss,
     p_cap_loss=result.p_cap_loss,
     p_filter_loss=result.p_filter_loss, ripple=result.ripple,
     neutral_current=result.i_neutral, zero_sequence=result.i_zero,
     negative_sequence=result.i_negative, dc_ripple=result.dv2,
     midpoint_ripple=result.dv_mid, midpoint_mean=result.v_mid_mean,
     cap_current=result.i_cap, cap_thermal_current=result.i_cap_thermal,
     cap_switching_current=result.i_cap_switching,
     cap_switching_reserved=result.i_cap_switching_reserved,
     pwm_reserve_margin=result.pwm_reserve_margin,
     dc_switching_ripple_rms=result.dv_switching_rms,
     dc_switching_ripple_pp=result.dv_switching_pp,
     pwm_modulation_margin=result.pwm_modulation_margin,
     pwm_iterations=result.pwm_iterations,
     switching_margin=result.switching_margin,
     filter_resonance_hz=result.filter_resonance_hz)

# Handles the hook publishes for post-solve reporting. Power/current/voltage
# expressions are in MODEL units; `sb`/`vb`/`ib` convert them back to SI.
struct _InvHandles
    p_poc; q_poc; p_conv; q_conv; p_loss; p_cap_loss; p_dc; p_filter_loss
    vrint::Vector{Any}; viint::Vector{Any}
    cri::Vector{Any}; cii::Vector{Any}
    vrmid::Vector{Any}; vimid::Vector{Any}
    gri::Vector{Any}; gii::Vector{Any}
    shre::Vector{Any}; shim::Vector{Any}
    ripple_re; ripple_im
    in_re; in_im          # SI neutral-current components (or nothing)
    seq                    # model-unit ((I0re,I0im), (I1...), (I2...)) or nothing
    dre; dim              # SI 2ω ripple-phasor components (or nothing)
    nre; nim              # SI split-link midpoint RMS phasor (or nothing)
    nmean; qbal           # SI mean midpoint shift (V) and balancing charge (C)
    icap_sq; icap_thermal_sq
    icap_upper_sq; icap_lower_sq
    icap_thermal_upper_sq; icap_thermal_lower_sq
    sb::Float64; vb::Float64; ib::Float64
end

_start_or(v, default) = (s = JuMP.start_value(v); s === nothing ? default : s)

# Stamp the inverter's internal-node model into `ctx`. Returns `_InvHandles`.
function _stamp_inverter!(ctx, inv::AdvancedInverter)
    m = _opf_model(ctx)
    vr, vi = _opf_voltage_maps(ctx)
    bus = inv.bus; phases = inv.phase_terminals; neutral = inv.neutral
    nph = length(phases)

    # Per-unit base factors (1.0 in SI mode), looked up per bus like the engine.
    bs = _opf_bases(ctx)
    sb = bs === nothing ? 1.0 : bs.s_base
    vb = bs === nothing ? 1.0 : get(bs.v_base, bus, 1.0)
    ib = bs === nothing ? 1.0 : get(bs.i_base, bus, 1.0)
    zb = bs === nothing ? 1.0 : get(bs.z_base, bus, 1.0)

    Rf, Xf = _series_filter_matrices(inv, nph, zb; side=:converter)
    Rg, Xg = _series_filter_matrices(inv, nph, zb; side=:grid)
    has_lcl = _lcl_active(inv)
    gmid_si, bmid_si = _mid_filter_admittance(inv)
    gmid, bmid = gmid_si * zb, bmid_si * zb # S → per-unit admittance

    is_3ph_topo = inv.topology in _THREE_PHASE_TOPOLOGIES

    vrint = Vector{Any}(undef, nph); viint = Vector{Any}(undef, nph)
    cri   = Vector{Any}(undef, nph); cii   = Vector{Any}(undef, nph)
    vrmid = Vector{Any}(undef, nph); vimid = Vector{Any}(undef, nph)
    gri   = Vector{Any}(undef, nph); gii   = Vector{Any}(undef, nph)
    shre  = Vector{Any}(undef, nph); shim  = Vector{Any}(undef, nph)
    dvr   = Vector{Any}(undef, nph); dvi   = Vector{Any}(undef, nph)

    # Allocate all phase currents before writing KVL. A non-zero neutral
    # impedance creates the common term Z_n ΣI_x in every phase-to-neutral
    # equation, so each phase equation depends on all currents.
    for (k, ph) in enumerate(phases)
        if neutral === nothing
            dvr[k] = vr[(bus, ph)]; dvi[k] = vi[(bus, ph)]
            vseedr = _start_or(vr[(bus, ph)], 0.0)
            vseedi = _start_or(vi[(bus, ph)], 0.0)
        else
            dvr[k] = JuMP.@expression(m, vr[(bus, ph)] - vr[(bus, neutral)])
            dvi[k] = JuMP.@expression(m, vi[(bus, ph)] - vi[(bus, neutral)])
            vseedr = _start_or(vr[(bus, ph)], 0.0) -
                     _start_or(vr[(bus, neutral)], 0.0)
            vseedi = _start_or(vi[(bus, ph)], 0.0) -
                     _start_or(vi[(bus, neutral)], 0.0)
        end
        vrint[k] = JuMP.@variable(m, base_name = "vrint_$(inv.id)_$(ph)")
        viint[k] = JuMP.@variable(m, base_name = "viint_$(inv.id)_$(ph)")
        cri[k]   = JuMP.@variable(m, base_name = "cri_$(inv.id)_$(ph)")
        cii[k]   = JuMP.@variable(m, base_name = "cii_$(inv.id)_$(ph)")
        JuMP.set_start_value(vrint[k], vseedr)
        JuMP.set_start_value(viint[k], has_lcl ? vseedi : 0.0)
        if has_lcl
            vrmid[k] = JuMP.@variable(m, base_name = "vrmid_$(inv.id)_$(ph)")
            vimid[k] = JuMP.@variable(m, base_name = "vimid_$(inv.id)_$(ph)")
            gri[k] = JuMP.@variable(m, base_name = "gri_$(inv.id)_$(ph)")
            gii[k] = JuMP.@variable(m, base_name = "gii_$(inv.id)_$(ph)")
            JuMP.set_start_value(vrmid[k], vseedr)
            JuMP.set_start_value(vimid[k], vseedi)
        else
            vrmid[k] = dvr[k]; vimid[k] = dvi[k]
            gri[k] = cri[k]; gii[k] = cii[k]
        end
    end
    for k in 1:nph
        if has_lcl
            shre[k] = JuMP.@expression(m, gmid*vrmid[k] - bmid*vimid[k])
            shim[k] = JuMP.@expression(m, bmid*vrmid[k] + gmid*vimid[k])
        else
            shre[k] = 0.0; shim[k] = 0.0
        end
    end
    sumcr = JuMP.@expression(m, sum(cri[k] for k in 1:nph))
    sumci = JuMP.@expression(m, sum(cii[k] for k in 1:nph))
    filter_dr, filter_di = _filter_voltage_drops(
        Rf, Xf, cri, cii, nph, neutral !== nothing)
    grid_dr, grid_di = _filter_voltage_drops(
        Rg, Xg, gri, gii, nph, neutral !== nothing)

    Pconv = zero(JuMP.QuadExpr); Qconv = zero(JuMP.QuadExpr)
    Ppoc  = zero(JuMP.QuadExpr); Qpoc  = zero(JuMP.QuadExpr)
    isq_sum = zero(JuMP.QuadExpr)     # Σ |I|²  (for c_loss)
    imag_sum = zero(JuMP.AffExpr)     # Σ |I|   (for a_loss)

    for (k, ph) in enumerate(phases)
        # Reduced mode: Vint--Zf--POC. Explicit LCL mode inserts the midpoint,
        # grid-side primitive Zg, and the damped capacitor branch. KCL separates
        # converter current from grid current at the midpoint.
        if has_lcl
            JuMP.@constraint(m, vrint[k] - vrmid[k] == filter_dr[k])
            JuMP.@constraint(m, viint[k] - vimid[k] == filter_di[k])
            JuMP.@constraint(m, vrmid[k] - dvr[k] == grid_dr[k])
            JuMP.@constraint(m, vimid[k] - dvi[k] == grid_di[k])
            JuMP.@constraint(m, cri[k] == gri[k] + shre[k])
            JuMP.@constraint(m, cii[k] == gii[k] + shim[k])
        else
            JuMP.@constraint(m, vrint[k] - dvr[k] == filter_dr[k])
            JuMP.@constraint(m, viint[k] - dvi[k] == filter_di[k])
        end

        # Converter-side power uses converter current; POC power uses grid-arm
        # current. They differ when the midpoint branch carries current.
        Pconv += JuMP.@expression(m, vrint[k]*cri[k] + viint[k]*cii[k])
        Qconv += JuMP.@expression(m, viint[k]*cri[k] - vrint[k]*cii[k])
        Ppoc  += JuMP.@expression(m, dvr[k]*gri[k] + dvi[k]*gii[k])
        Qpoc  += JuMP.@expression(m, dvi[k]*gri[k] - dvr[k]*gii[k])

        # Only the grid-side arm reaches the POC bus KCL.
        BMOPFTools.add_terminal_injection!(ctx, bus, ph, gri[k], gii[k])
        if neutral !== nothing
            BMOPFTools.add_terminal_injection!(ctx, bus, neutral, -gri[k], -gii[k])
        end

        # |I|² (for the quadratic loss term and the current cap) and the current
        # magnitude |I| via an auxiliary. The magnitude only feeds the LINEAR
        # loss term, and how tightly it is pinned matters only there:
        #   • a_loss ≠ 0 → the shifted smooth norm
        #     `sqrt(isq + ε²) - ε`, ε=1 mA. An epigraph can be inflated by an
        #     unrelated charging/price objective; the equality cannot. Smoothing
        #     removes the singular Jacobian of the exact norm at I=0.
        #   • a_loss == 0 → the original epigraph is retained as a benign
        #     regularising constraint; its slack cannot affect the zero coefficient.
        isq = JuMP.@expression(m, cri[k]^2 + cii[k]^2)
        isq_sum += isq
        if inv.i_max !== nothing
            JuMP.@constraint(m, isq <= (inv.i_max / ib)^2)
        end
        if inv.i_grid_max !== nothing
            JuMP.@constraint(m, gri[k]^2 + gii[k]^2 <= (inv.i_grid_max / ib)^2)
        end
        if inv.a_loss != 0.0
            imag_eps = _IMAG_EPS_SI / ib
            im = JuMP.@variable(m, base_name = "imag_$(inv.id)_$(ph)",
                                lower_bound = imag_eps)
            JuMP.set_start_value(im, imag_eps)
            JuMP.@constraint(m, im^2 == isq + imag_eps^2)
            imag_sum += im - imag_eps
        else
            im = JuMP.@variable(m, base_name = "imag_$(inv.id)_$(ph)", lower_bound = 0.0)
            JuMP.@constraint(m, im^2 >= isq)
            imag_sum += im
        end

        # Internal EMF magnitude bounds (skipped under grid-forming, which pins the
        # magnitude via v_gfm below). The v_int box applies to every topology; the
        # scalar modulation cap only to SINGLE_PHASE (three-phase topologies use the
        # sampled switching-polytope instead).
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

    # The neutral semiconductor leg of a four-leg bridge carries -ΣI_phase.
    # Include that fourth conducting path in the same fitted per-leg loss curve.
    # A split-link converter has no fourth leg; its neutral heating is accounted
    # for in the capacitor thermal allocation instead.
    if inv.topology == :FOUR_LEG
        isq_n = JuMP.@expression(m, sumcr^2 + sumci^2)
        isq_sum += isq_n
        if inv.a_loss != 0.0
            imag_eps = _IMAG_EPS_SI / ib
            imn = JuMP.@variable(m, base_name = "imag_$(inv.id)_neutral",
                                 lower_bound = imag_eps)
            JuMP.set_start_value(imn, imag_eps)
            JuMP.@constraint(m, imn^2 == isq_n + imag_eps^2)
            imag_sum += imn - imag_eps
        else
            imn = JuMP.@variable(m, base_name = "imag_$(inv.id)_neutral", lower_bound = 0.0)
            JuMP.@constraint(m, imn^2 >= isq_n)
            imag_sum += imn
        end
    end

    # Non-conjugated Σ V_int·I (converter side) — the 2ω oscillating-power
    # phasor. Its magnitude is a sinusoidal amplitude; capacitor thermal loading
    # below converts that amplitude to RMS with the explicit 1/sqrt(2) factor.
    rip_re, rip_im = _ripple_components(vrint, viint, cri, cii)

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
            BMOPFTools.add_terminal_injection!(ctx, bus, ph, b * vip, -b * vrp)
            Qpoc += JuMP.@expression(m, b * (vrp^2 + vip^2))
            if neutral !== nothing
                BMOPFTools.add_terminal_injection!(ctx, bus, neutral, -b * vip, b * vrp)
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

    # ── Three-phase switching-polytope feasibility (time-sampled) ─────────────
    in_re = nothing; in_im = nothing; seq = nothing
    dre = nothing; dim = nothing; nre = nothing; nim = nothing
    nmean = nothing; qbal = nothing
    icap_sq = nothing; icap_thermal_sq = nothing
    icap_upper_sq = nothing; icap_lower_sq = nothing
    icap_thermal_upper_sq = nothing; icap_thermal_lower_sq = nothing
    P_cap_loss = 0.0
    if is_3ph_topo
        w = 2pi * inv.f
        v_dc = inv.v_dc; c_dc = inv.c_dc
        # Converter output (pole) voltages and 2ω oscillating power, in SI.
        UreSI = [JuMP.@expression(m, vrint[k]*vb) for k in 1:3]
        UimSI = [JuMP.@expression(m, viint[k]*vb) for k in 1:3]
        S2re_SI = JuMP.@expression(m, rip_re * sb)
        S2im_SI = JuMP.@expression(m, rip_im * sb)

        # 2ω bus-ripple phasor D = j·S2/(2ω·C_eq·Vdc). Unequal split
        # half-banks use their exact series equivalent C_u C_l/(C_u+C_l).
        if inv.topology == :SPLIT_DC
            cu, cl = _split_capacitances(inv)
            csum = cu + cl
            ceq = cu * cl / csum
            denom = 2w * ceq * v_dc
        else
            cu = cl = csum = NaN
            denom = 2w * c_dc * v_dc
        end
        dre = JuMP.@variable(m, base_name = "dre_$(inv.id)")
        dim = JuMP.@variable(m, base_name = "dim_$(inv.id)")
        JuMP.@constraint(m, dre == -S2im_SI / denom)
        JuMP.@constraint(m, dim ==  S2re_SI / denom)
        inv.dv2_max !== nothing && JuMP.@constraint(m, dre^2 + dim^2 <= inv.dv2_max^2)

        # Split-link midpoint fundamental ripple plus mean displacement. Neutral
        # current divides in proportion to half-bank capacitance; optional
        # balancing charge shifts the mean away from the equal-charge value.
        if inv.topology == :SPLIT_DC
            nre = JuMP.@expression(m,  sumci*ib / (w * csum))
            nim = JuMP.@expression(m, -sumcr*ib / (w * csum))
            if inv.dv_mid_max !== nothing
                JuMP.@constraint(m, nre^2 + nim^2 <= inv.dv_mid_max^2)
            end
            natural_mean = v_dc * (cu - cl) / (2csum)
            if inv.q_mid_balance_max === nothing
                nmean = natural_mean
            else
                qmax = inv.q_mid_balance_max
                qbal = JuMP.@variable(m, base_name = "qbal_$(inv.id)",
                    lower_bound = -qmax, upper_bound = qmax)
                nmean = JuMP.@variable(m, base_name = "nmean_$(inv.id)")
                JuMP.set_start_value(qbal, 0.0)
                JuMP.set_start_value(nmean, natural_mean)
                JuMP.@constraint(m, nmean == natural_mean + 2qbal/csum)
                if inv.v_mid_mean_max !== nothing
                    JuMP.@constraint(m, nmean^2 <= inv.v_mid_mean_max^2)
                end
            end
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
                        sg*(sqrt(2)*((UreSI[k]+nre)*c1 - (UimSI[k]+nim)*s1) +
                            nmean) <= railh)
                end
            end
        end

        # Neutral current: 3-wire carries none; 4-wire limited by In_max.
        if inv.topology == :THREE_LEG
            JuMP.@constraint(m, sum(cri[k] for k in 1:3) == 0)
            JuMP.@constraint(m, sum(cii[k] for k in 1:3) == 0)
        elseif inv.In_max !== nothing
            # Standalone neutral rating. Optional for :SPLIT_DC, where i_cap_max
            # can bound |I_n| instead (via the capacitor thermal budget below).
            JuMP.@constraint(m, sumcr^2 + sumci^2 <= (inv.In_max / ib)^2)
        end
        # Neutral current (SI), for reporting: I_n = −Σ I_phase.
        in_re = JuMP.@expression(m, -sum(cri[k] for k in 1:3) * ib)
        in_im = JuMP.@expression(m, -sum(cii[k] for k in 1:3) * ib)

        # Fortescue components expose the two physically different unbalance
        # channels: I0 sets neutral current, while I2 interacts with positive-
        # sequence voltage to produce the dominant 2ω DC-link power pulsation.
        seq = _sequence_components(cri, cii)
        for (name, comp, limit) in (("zero", seq[1], inv.i_zero_max),
                                    ("positive", seq[2], inv.i_positive_max),
                                    ("negative", seq[3], inv.i_negative_max))
            if limit !== nothing
                sre = JuMP.@variable(m, base_name = "i$(name)_re_$(inv.id)")
                sim = JuMP.@variable(m, base_name = "i$(name)_im_$(inv.id)")
                JuMP.@constraint(m, sre == comp[1])
                JuMP.@constraint(m, sim == comp[2])
                JuMP.@constraint(m, sre^2 + sim^2 <= (limit / ib)^2)
            end
        end

        # ── Capacitor RMS ripple current: SIMULTANEOUS thermal allocation ───────
        # The bank always carries the 2ω bus current, and in the split link the
        # fundamental neutral current as well. They sit at different frequencies,
        # so they combine in RMS and must be allocated together — charging the
        # whole bank rating to neutral current double-counts the capacitors.
        #
        #   I_2ω,rms = |S̃|/(√2·Vdc) = k·|D|,   k = denom/(√2·Vdc)
        #
        # Writing it through the ripple phasor D keeps this QUADRATIC in variables
        # that already exist, rather than quartic in the bilinear S̃.
        k2w = denom / (sqrt(2) * v_dc)
        i2w_sq = JuMP.@expression(m, k2w^2 * (dre^2 + dim^2))
        wn, w2, wsw = inv.cap_thermal_weights
        if inv.topology == :SPLIT_DC
            neutral_sq = JuMP.@expression(m, in_re^2 + in_im^2)
            au, al = cu/csum, cl/csum
            icap_upper_sq = JuMP.@expression(m,
                i2w_sq + au^2*neutral_sq + inv.i_sw^2)
            icap_lower_sq = JuMP.@expression(m,
                i2w_sq + al^2*neutral_sq + inv.i_sw^2)
            icap_thermal_upper_sq = JuMP.@expression(m,
                w2*i2w_sq + wn*au^2*neutral_sq + wsw*inv.i_sw^2)
            icap_thermal_lower_sq = JuMP.@expression(m,
                w2*i2w_sq + wn*al^2*neutral_sq + wsw*inv.i_sw^2)
            if inv.i_cap_max !== nothing
                JuMP.@constraint(m, icap_thermal_upper_sq <= inv.i_cap_max^2)
                if cu != cl
                    JuMP.@constraint(m, icap_thermal_lower_sq <= inv.i_cap_max^2)
                end
            end
            if inv.i_cap_upper_max !== nothing
                JuMP.@constraint(m, icap_thermal_upper_sq <= inv.i_cap_upper_max^2)
            end
            if inv.i_cap_lower_max !== nothing
                JuMP.@constraint(m, icap_thermal_lower_sq <= inv.i_cap_lower_max^2)
            end
            P_cap_loss = JuMP.@expression(m,
                (inv.esr_dc_upper*icap_thermal_upper_sq +
                 inv.esr_dc_lower*icap_thermal_lower_sq) / sb)
        else
            icap_sq = JuMP.@expression(m, i2w_sq + inv.i_sw^2)
            icap_thermal_sq = JuMP.@expression(m, w2*i2w_sq + wsw*inv.i_sw^2)
            if inv.i_cap_max !== nothing
                JuMP.@constraint(m, icap_thermal_sq <= inv.i_cap_max^2)
            end
            P_cap_loss = JuMP.@expression(m, inv.esr_dc*icap_thermal_sq / sb)
        end
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
    P_dc = JuMP.@expression(m, Pconv + P_loss + P_cap_loss)
    P_filter_loss = JuMP.@expression(m, Pconv - Ppoc)

    # Single-phase standalone 2ω ripple magnitude bound (VA → per-unit power).
    if inv.p_ripple_max !== nothing && !is_3ph_topo
        rr = JuMP.@variable(m, base_name = "rip_re_$(inv.id)")
        ri = JuMP.@variable(m, base_name = "rip_im_$(inv.id)")
        JuMP.@constraint(m, rr == rip_re)
        JuMP.@constraint(m, ri == rip_im)
        JuMP.@constraint(m, rr^2 + ri^2 <= (inv.p_ripple_max / sb)^2)
    end

    return _InvHandles(Ppoc, Qpoc, Pconv, Qconv, P_loss, P_cap_loss, P_dc, P_filter_loss,
                       vrint, viint, cri, cii,
                       vrmid, vimid, gri, gii, shre, shim, rip_re, rip_im,
                       in_re, in_im, seq, dre, dim, nre, nim, nmean, qbal,
                       icap_sq, icap_thermal_sq,
                       icap_upper_sq, icap_lower_sq,
                       icap_thermal_upper_sq, icap_thermal_lower_sq,
                       sb, vb, ib)
end

function stamp_device!(ctx, inverter::AdvancedInverter; kwargs...)
    inverter.pwm_strategy == :NONE || throw(ArgumentError(
        "pwm_strategy=:$(inverter.pwm_strategy) uses the conservative outer " *
        "reserve iteration provided by solve_advanced_inverter; direct stamping " *
        "requires pwm_strategy=:NONE and an explicit i_sw reserve"))
    return _stamp_inverter!(ctx, inverter)
end

link_device!(model, inverter::AdvancedInverter, handles, sb, grid::TimeGrid) = nothing

function extract_device(inverter::AdvancedInverter, h::_InvHandles,
                        status::SolveStatus)
    solved = status.publishable
    sb = h.sb; vb = h.vb; ib = h.ib
    scaled(e, scale) = solved ? JuMP.value(e) * scale : NaN
    mag(a, b, scale) = solved ? hypot(JuMP.value(a), JuMP.value(b)) * scale : NaN
    nph = length(inverter.phase_terminals)
    seqmag(k) = h.seq === nothing ? 0.0 : mag(h.seq[k][1], h.seq[k][2], ib)
    sqmag(e) = e === nothing ? 0.0 :
        (solved ? sqrt(max(JuMP.value(e), 0.0)) : NaN)
    icu = sqmag(h.icap_upper_sq); icl = sqmag(h.icap_lower_sq)
    ictu = sqmag(h.icap_thermal_upper_sq); ictl = sqmag(h.icap_thermal_lower_sq)
    icap = inverter.topology == :SPLIT_DC ? max(icu, icl) : sqmag(h.icap_sq)
    icap_thermal = inverter.topology == :SPLIT_DC ? max(ictu, ictl) :
        sqmag(h.icap_thermal_sq)
    switching_margin = if !(inverter.topology in _THREE_PHASE_TOPOLOGIES)
        0.0
    elseif !solved
        NaN
    else
        ure = [JuMP.value(h.vrint[k]) * vb for k in 1:3]
        uim = [JuMP.value(h.viint[k]) * vb for k in 1:3]
        dr = JuMP.value(h.dre); di = JuMP.value(h.dim)
        nr = h.nre === nothing ? nothing : JuMP.value(h.nre)
        ni = h.nim === nothing ? nothing : JuMP.value(h.nim)
        nm = h.nmean === nothing ? 0.0 : JuMP.value(h.nmean)
        _switching_margin(inverter, ure, uim, dr, di, nr, ni, nm)
    end
    return (
        p_poc=scaled(h.p_poc, sb), q_poc=scaled(h.q_poc, sb),
        p_conv=scaled(h.p_conv, sb), q_conv=scaled(h.q_conv, sb),
        p_loss=scaled(h.p_loss, sb), p_cap_loss=scaled(h.p_cap_loss, sb),
        p_dc=scaled(h.p_dc, sb),
        p_filter_loss=scaled(h.p_filter_loss, sb),
        v_int_mag=[mag(h.vrint[k], h.viint[k], vb) for k in 1:nph],
        v_filter_mag=[mag(h.vrmid[k], h.vimid[k], vb) for k in 1:nph],
        i_mag=[mag(h.cri[k], h.cii[k], ib) for k in 1:nph],
        i_grid_mag=[mag(h.gri[k], h.gii[k], ib) for k in 1:nph],
        i_filter_shunt_mag=[mag(h.shre[k], h.shim[k], ib) for k in 1:nph],
        i_neutral=h.in_re === nothing ? 0.0 :
            (solved ? hypot(JuMP.value(h.in_re), JuMP.value(h.in_im)) : NaN),
        i_zero=seqmag(1), i_positive=seqmag(2), i_negative=seqmag(3),
        ripple=mag(h.ripple_re, h.ripple_im, sb),
        dv2=h.dre === nothing ? 0.0 :
            (solved ? hypot(JuMP.value(h.dre), JuMP.value(h.dim)) : NaN),
        dv_mid=h.nre === nothing ? 0.0 :
            (solved ? hypot(JuMP.value(h.nre), JuMP.value(h.nim)) : NaN),
        v_mid_mean=h.nmean === nothing ? 0.0 : scaled(h.nmean, 1.0),
        q_mid_balance=h.qbal === nothing ? 0.0 : scaled(h.qbal, 1.0),
        i_cap=icap, i_cap_thermal=icap_thermal,
        i_cap_upper=icu, i_cap_lower=icl,
        i_cap_thermal_upper=ictu, i_cap_thermal_lower=ictl,
        switching_margin=switching_margin,
        filter_resonance_hz=_filter_resonance_hz(inverter),
    )
end

"""
    solve_advanced_inverter(net, inverter; objective=:max_export, kwargs...)
        -> InverterResult

Stamp `inverter` into `net` and solve, demonstrating the experimental internal-node
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
function _solve_advanced_inverter_once(net::Dict{String,Any}, inverter::AdvancedInverter;
                                       objective::Symbol=:max_export,
                                       p_set::Union{Float64,Nothing}=nothing,
                                       q_set::Union{Float64,Nothing}=nothing,
                                       per_unit::Bool=false,
                                       s_base::Float64=1e6,
                                       optimizer=Ipopt.Optimizer,
                                       verbose::Bool=false,
                                       solver_options=(),
                                       pwm_audit_inverter::AdvancedInverter=inverter,
                                       pwm_iterations::Int=0)
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
        q_set === nothing || JuMP.@constraint(_opf_model(ctx), h.q_poc == q_set / h.sb)
        if objective == :max_export
            JuMP.@objective(_opf_model(ctx), Max, h.p_poc)
        else
            JuMP.@constraint(_opf_model(ctx), h.p_poc == p_set / h.sb)
            JuMP.@objective(_opf_model(ctx), Min, h.p_loss + h.p_cap_loss)
        end
    end

    ctx = build_opf_model(net; per_unit=per_unit, s_base=s_base,
                          add_objective=false, model_hook! = hook!,
                          optimizer=optimizer, verbose=verbose)
    _set_solver_options!(_opf_model(ctx), solver_options)
    enforce_kcl!(ctx)
    JuMP.optimize!(_opf_model(ctx))

    outcome = _solve_outcome(_opf_model(ctx))
    status = string(outcome.termination_status)
    h = handles[]
    device_result = extract_device(inverter, h, SolveStatus(outcome))

    pwm = if pwm_audit_inverter.pwm_strategy != :NONE && _publishable(outcome)
        ure = [JuMP.value(h.vrint[k]) * h.vb for k in 1:3]
        uim = [JuMP.value(h.viint[k]) * h.vb for k in 1:3]
        ire = [JuMP.value(h.cri[k]) * h.ib for k in 1:3]
        iim = [JuMP.value(h.cii[k]) * h.ib for k in 1:3]
        dr = h.dre === nothing ? 0.0 : JuMP.value(h.dre)
        di = h.dim === nothing ? 0.0 : JuMP.value(h.dim)
        nr = h.nre === nothing ? 0.0 : JuMP.value(h.nre)
        ni = h.nim === nothing ? 0.0 : JuMP.value(h.nim)
        nm = h.nmean === nothing ? 0.0 : JuMP.value(h.nmean)
        _pwm_ripple_audit(pwm_audit_inverter, ure, uim, ire, iim,
                          dr, di, nr, ni, nm)
    elseif pwm_audit_inverter.pwm_strategy != :NONE
        (i_rms=NaN, dv_rms=NaN, dv_pp=NaN, modulation_margin=NaN)
    else
        (i_rms=0.0, dv_rms=0.0, dv_pp=0.0, modulation_margin=NaN)
    end
    pwm_target = pwm_audit_inverter.pwm_strategy == :NONE ? NaN :
        hypot(pwm_audit_inverter.i_sw, pwm.i_rms)
    pwm_reserve_margin = isnan(pwm_target) ? NaN : inverter.i_sw - pwm_target

    result = _extract_result(ctx, outcome)
    return InverterResult(status, inverter.topology,
                          device_result.p_poc, device_result.q_poc,
                          device_result.p_conv, device_result.q_conv,
                          device_result.p_loss, device_result.p_cap_loss,
                          device_result.p_dc,
                          device_result.p_filter_loss,
                          device_result.v_int_mag, device_result.v_filter_mag,
                          device_result.i_mag, device_result.i_grid_mag,
                          device_result.i_filter_shunt_mag,
                          device_result.i_neutral, device_result.i_zero,
                          device_result.i_positive, device_result.i_negative,
                          device_result.ripple, device_result.dv2,
                          device_result.dv_mid, device_result.v_mid_mean,
                          device_result.q_mid_balance,
                          device_result.i_cap, device_result.i_cap_thermal,
                          device_result.i_cap_upper, device_result.i_cap_lower,
                          device_result.i_cap_thermal_upper,
                          device_result.i_cap_thermal_lower,
                          pwm.i_rms, inverter.i_sw, pwm_reserve_margin,
                          pwm.dv_rms, pwm.dv_pp, pwm.modulation_margin,
                          pwm_iterations,
                          device_result.switching_margin,
                          device_result.filter_resonance_hz,
                          result["bus"], SolveStatus(outcome))
end
