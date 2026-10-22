# powerflow/helm.jl
#
# HELM — Holomorphic Embedding Load-flow Method — for unbalanced multiconductor
# (4-wire and beyond) distribution power flow, on the augmented nodal admittance
# matrix (`ybus_augmented`: ideal switches / ideal transformers as constraint
# rows). Nothing is Kron-reduced: every explicit neutral, PE or extra conductor
# keeps its own node, and the recursion below is generic over node indices, so
# the solver itself has no conductor-count assumption.
#
# The power-flow equations are embedded in a complex parameter `s`:
#
#   • voltage sources stay FIXED for all s — either as a Dirichlet boundary
#     (`source_reference = :earth`, eliminated from the unknowns) or as
#     bordered phase-to-source-neutral constraint rows
#     (`source_reference = :neutral`, the source neutral then floats with the
#     network's earthing);
#   • loads enter through one of two embeddings:
#       `:load_scaling`   — every load part is scaled by s; s = 0 is the
#                            energized NO-LOAD network and s is a load
#                            multiplier (the classic HELM germ);
#       `:impedance_germ` — constant-impedance parts sit in the matrix, and
#                            every other part is blended from an equivalent
#                            admittance y_g at s = 0 to its true law at s = 1:
#                            I(s) = (1 − s)·y_g·ΔV + s·I_law(ΔV). The germ is
#                            then a loaded linear network, which determines
#                            nodes held only by loads (e.g. an unearthed star
#                            point) that the no-load germ cannot.
#
# Each voltage (and coupling / source current) is a holomorphic function of s,
# expanded as a power series whose coefficients solve ONE constant linear system
# per order:
#
#   K_UU · x[n] = rhs(x[0..n−1])
#
# with the same LU factorization for the germ and every order.
#
# Load laws. A sub-load across (pos, neg) with ΔV = V_pos − V_neg draws
#
#   I = Σ_γ κ_γ · |ΔV|^γ / conj(ΔV)
#
# (γ = 0 constant power, γ = 1 constant current, γ = 2 constant impedance,
# other γ exponential; κ_γ = a_γ − j·b_γ from P = Σ a_γ|ΔV|^γ, Q = Σ b_γ|ΔV|^γ).
# The conj(·) is handled with the conjugate-reflection trick: the reflected
# series d̃[k] = conj(d[k]) is the expansion of conj(ΔV(conj s)), itself
# holomorphic, and agrees with conj(ΔV(s)) on the real s axis. Hence
# m(s) = ΔV(s)·d̃(s) is holomorphic and equals |ΔV|² on the real axis, so
# |ΔV|^γ = m^{γ/2} is holomorphic too (m[0] = |ΔV₀|² > 0). The reciprocal 1/d̃
# and the power m^{γ/2} both have incremental recursions (the latter J. C. P.
# Miller's formula), so the order-n RHS needs only coefficients 0..n−1.
# Constraint rows have zero RHS at every order above the germ — the
# ideal-element identities hold term-by-term and therefore exactly in the
# summed solution.
#
# The series is evaluated at s = 1 by Padé analytic continuation (Wynn epsilon,
# powerflow/pade.jl). Finite-order Padé spread, physical mismatch, and
# coefficient-tail growth are numerical diagnostics: they do not by themselves
# prove non-existence of a power-flow solution. A Domb–Sykes coefficient-ratio
# extrapolation is reported as a `singularity_estimate`; it needs validation
# against continuation or an analytically known system before being interpreted
# as a voltage-collapse loading multiplier (and it is only a load multiplier
# under `:load_scaling`).
#
# Scope: lines, shunts, capacitors, every transformer subtype the passive Ybus
# stamps, ideal switches/transformers as constraint rows, and loads of every
# BMOPF model (constant power / current / impedance, ZIP, exponential) in WYE,
# DELTA and SINGLE_PHASE connections. Generators and IBRs have no fixed
# power-flow law in the data model (their set-points are OPF variables or
# control-law outputs) and raise an `ArgumentError` rather than being silently
# dropped. Ideal transformers keep their no-load (magnetising) shunt and
# neutral-grounding branch: those are restored on top of `ybus_augmented`,
# which drops them together with the singular series stamp.

# The BMOPFTools compatibility imports used here are isolated in
# `src/upstream.jl`; see the contributing guide for their upstream replacement.

"""
    HelmResult

Result of [`helm_series`](@ref) — the HELM power-flow solution and its
diagnostics.

Fields:
- `V`           — `Dict{(bus,terminal),ComplexF64}` node-to-earth voltages at
                  `s = 1` (every declared terminal, including source/fixed,
                  aliased, and earth-referenced ones).
- `w`           — unscaled bordered-row solution entries; the PHYSICAL current
                  of `couplings[j]` is `couplings[j].scale * w[j]`.
- `couplings`   — the `IdealCoupling`s of the underlying
                  `ybus_augmented` matrix, in `w` order.
- `coeffs`      — series coefficients, one row per unknown (non-fixed nodes,
                  then couplings, then source constraint rows) and one column
                  per order (diagnostic).
- `converged`   — `true` iff the s = 1 nonlinear current mismatch is within
                  tolerance.
- `status`      — `:converged` | `:series_diverged` (the finite coefficient
                  tail is growing) | `:max_order_reached`.
- `residual`    — max |current mismatch| (A) over every non-Dirichlet node
                  (source-bus terminals included when they are unknowns) and
                  every constraint row, evaluated at the returned solution.
- `n_order`     — highest series order computed.
- `pade_spread` — absolute difference between the last two Wynn-epsilon Padé
                  estimates for each row of `coeffs`. Values have the units of
                  their corresponding row.
- `coefficient_tail_norms` — max-norm of the final coefficient orders used by
                  the status classifier.
- `coefficient_tail_ratios` — adjacent ratios of those tail norms. A tail with
                  every ratio above one is classified as `:series_diverged`.
- `singularity_estimate` — heuristic Domb–Sykes estimate of the nearest
                  coefficient-dominating singularity in the embedding
                  parameter. Under `embedding = :load_scaling` that parameter
                  is a load multiplier; it is never a certified loading margin.
                  `NaN` means the series is too short or featureless.
- `source_injection` — `Dict{(bus,terminal),ComplexF64}` current (A) injected
                  INTO the network by the voltage-source boundary at each node
                  it acts on: the slack phase currents, and, with
                  `source_reference = :earth`, the earthing current of the
                  fixed source neutral.
- `embedding`, `source_reference` — the embedding options used.

`result.load_margin` remains a compatibility alias for
`result.singularity_estimate`; new code should use the latter name.
"""
struct HelmResult <: AbstractSolveResult
    V::Dict{_Node,ComplexF64}
    w::Vector{ComplexF64}
    couplings::Vector{IdealCoupling}
    coeffs::Matrix{ComplexF64}
    converged::Bool
    status::Symbol
    residual::Float64
    n_order::Int
    pade_spread::Vector{Float64}
    coefficient_tail_norms::Vector{Float64}
    coefficient_tail_ratios::Vector{Float64}
    singularity_estimate::Float64
    source_injection::Dict{_Node,ComplexF64}
    embedding::Symbol
    source_reference::Symbol
end

function Base.getproperty(r::HelmResult, name::Symbol)
    name === :load_margin && return getfield(r, :singularity_estimate)
    getfield(r, name)
end

Base.propertynames(::HelmResult, private::Bool=false) =
    (fieldnames(HelmResult)..., :load_margin)

function solve_status(result::HelmResult)
    publishable = result.converged && result.status === :converged
    _result_solve_status(string(result.status), publishable;
        primal_status=publishable ? "FEASIBLE_POINT" : "NO_SOLUTION")
end

solve_diagnostics(result::HelmResult) =
    (residual=result.residual, order=result.n_order,
     pade_spread=result.pade_spread,
     coefficient_tail_norms=result.coefficient_tail_norms,
     coefficient_tail_ratios=result.coefficient_tail_ratios,
     singularity_estimate=result.singularity_estimate,
     load_margin=result.load_margin,
     embedding=result.embedding,
     source_reference=result.source_reference)

Base.show(io::IO, r::HelmResult) =
    print(io, "HelmResult($(length(r.V)) terminals, status=$(r.status), " *
              "order=$(r.n_order), residual=$(round(r.residual; sigdigits=3)) A, " *
              "singularity_estimate=$(round(r.singularity_estimate; sigdigits=4)))")

# One non-constant-impedance sub-load in the recursion: endpoint positions in
# the unknown vector (0 = fixed/earth endpoint, with `v0p`/`v0n` carrying the
# order-0 voltage), the current-law terms (κ, γ), and the growing series of ΔV
# (`d`), its reflection (`dt`), the reciprocal of the reflection (`e`), the
# squared modulus m = d·dt (`m`) and its powers m^{γ/2} (`p`, one per γ ≠ 0).
mutable struct _HelmLoad
    up::Int
    un::Int
    v0p::ComplexF64
    v0n::ComplexF64
    terms::Vector{Tuple{ComplexF64,Float64}}
    d::Vector{ComplexF64}
    dt::Vector{ComplexF64}
    e::Vector{ComplexF64}
    m::Vector{ComplexF64}
    p::Dict{Float64,Vector{ComplexF64}}
    name::String
end

# One voltage-source constraint row (`source_reference = :neutral`):
# Σ coeffs·V[nodes] = E, with the source current as its bordered unknown.
struct _HelmSourceRow
    nodes::Vector{Int}
    coeffs::Vector{Float64}
    E::ComplexF64
    terminal::_Node
end

# ── load laws ────────────────────────────────────────────────────────────────

# Current-law terms of a sub-load excluding its constant-impedance part:
# I_drawn = Σ κ·|ΔV|^γ / conj(ΔV). The P family contributes κ = a, the Q family
# κ = −j·b (conj(P + jQ) = P − jQ).
function _helm_current_terms(sl::_SubLoad)
    acc = Dict{Float64,ComplexF64}()
    add!(γ, κ) = iszero(κ) || (acc[γ] = get(acc, γ, 0.0im) + κ)
    for (t, rot) in ((sl.pt, 1.0 + 0.0im), (sl.qt, -1.0im))
        if t.nl === nothing
            add!(0.0, rot * t.cc)
            add!(1.0, rot * t.cs)
        else
            base, γ = t.nl
            add!(Float64(γ), rot * base / sl.Vnom^γ)
        end
    end
    Tuple{ComplexF64,Float64}[(κ, γ) for (γ, κ) in sort!(collect(acc); by=first)]
end

# Drawn current of the law terms at a given ΔV (convergence check / germ).
_helm_law_current(terms, dv::ComplexF64) =
    sum(κ * abs(dv)^γ for (κ, γ) in terms; init=0.0im) / conj(dv)

# ── unsupported physics ──────────────────────────────────────────────────────

# Device families with no fixed power-flow law in the BMOPF data model. HELM
# refuses them instead of solving a network that silently omits them.
function _helm_unsupported_devices(net::Dict{String,Any})
    out = String[]
    for kind in ("generator", "ibr")
        devs = get(net, kind, nothing)
        devs isa AbstractDict && !isempty(devs) || continue
        push!(out, "$(length(devs)) $(kind)(s) ($(join(sort!(collect(String.(keys(devs)))), ", ")))")
    end
    out
end

# ── boundary (fixed nodes / source rows) ─────────────────────────────────────

function _helm_fix!(fixed::Dict{Int,ComplexF64}, aug::AugYbusResult,
                    nd::_Node, v::ComplexF64, sid::String)
    gi = get(aug.index, nd, 0)
    if gi == 0
        abs(v) <= 1e-9 || throw(ArgumentError(
            "voltage source '$sid' fixes $(nd) to $(v) V, but that terminal " *
            "is perfectly grounded (0 V)"))
        return
    end
    if haskey(fixed, gi)
        isapprox(fixed[gi], v; atol=1e-9) || throw(ArgumentError(
            "conflicting fixed voltages at node $(nd) (aliased terminals): " *
            "$(fixed[gi]) vs $(v)"))
        return
    end
    fixed[gi] = v
    return
end

# Iterate the voltage sources as (sid, bus, configuration-checked phasor list).
function _helm_sources(net::Dict{String,Any})
    out = Tuple{String,String,Vector{Tuple{String,ComplexF64}}}[]
    vss = get(net, "voltage_source", Dict())
    for sid in sort!(collect(String.(keys(vss))))
        vs = vss[sid]
        bus = string(get(vs, "bus", ""))
        cfg = uppercase(string(get(vs, "configuration", "WYE")))
        cfg in ("WYE", "SINGLE_PHASE") || throw(ArgumentError(
            "voltage source '$sid': configuration '$cfg' is not supported by " *
            "HELM (WYE / SINGLE_PHASE phase-to-neutral references only)"))
        tm    = Vector{String}(string.(get(vs, "terminal_map", String[])))
        v_mag = Float64.(get(vs, "v_magnitude", Float64[]))
        v_ang = Float64.(get(vs, "v_angle",     Float64[]))
        ph = Tuple{String,ComplexF64}[]
        for (k, t) in enumerate(tm)
            (length(v_mag) >= k && length(v_ang) >= k) || continue
            push!(ph, (t, v_mag[k] * cis(v_ang[k])))
        end
        push!(out, (sid, bus, ph))
    end
    out
end

# `source_reference = :earth`: voltage-source phase terminals at their reference
# phasors, source-bus neutral-labelled terminals at 0 V (mirrors the OPF's
# `_add_source_constraints!` semantics).
function _helm_fixed_nodes(net::Dict{String,Any}, aug::AugYbusResult)
    fixed = Dict{Int,ComplexF64}()
    nlabels = _neutral_labels(net)
    buses = get(net, "bus", Dict())
    for (sid, bus, ph) in _helm_sources(net)
        for (t, v) in ph
            _helm_fix!(fixed, aug, (bus, t), v, sid)
        end
        busd = get(buses, bus, Dict{String,Any}())
        nt = _neutral_terminal(busd)
        for t in get(busd, "terminal_names", String[])
            tt = string(t)
            (tt == nt || tt in nlabels) || continue
            gi = get(aug.index, (bus, tt), 0)
            (gi == 0 || haskey(fixed, gi)) && continue
            fixed[gi] = 0.0 + 0.0im
        end
    end
    fixed
end

# `source_reference = :neutral`: each source phasor is imposed between the phase
# terminal and the source bus's neutral terminal (or earth, if the bus has no
# neutral) as a bordered constraint row. The source neutral is NOT tied to
# earth by the source: its potential follows the network's earthing (perfect
# grounds, grounding shunts, transformer neutral impedances).
function _helm_source_rows(net::Dict{String,Any}, aug::AugYbusResult)
    rows = _HelmSourceRow[]
    nlabels = _neutral_labels(net)
    buses = get(net, "bus", Dict())
    for (sid, bus, ph) in _helm_sources(net)
        busd = get(buses, bus, Dict{String,Any}())
        nt = _neutral_terminal(busd)
        for (t, E) in ph
            if t == nt || t in nlabels
                abs(E) <= 1e-9 || throw(ArgumentError(
                    "voltage source '$sid' assigns $(E) V to the neutral terminal " *
                    "'$t'; with source_reference = :neutral the source neutral is " *
                    "the reference of the other phasors"))
                continue
            end
            acc = Dict{Int,Float64}()
            gp = get(aug.index, (bus, t), 0)
            gp == 0 && throw(ArgumentError(
                "voltage source '$sid' drives ($bus, $t), which is perfectly grounded"))
            acc[gp] = 1.0
            if nt !== nothing && nt != t
                gn = get(aug.index, (bus, string(nt)), 0)
                gn != 0 && (acc[gn] = get(acc, gn, 0.0) - 1.0)
            end
            nds = sort!(collect(keys(acc)))
            filter!(i -> abs(acc[i]) > 1e-12, nds)
            isempty(nds) && throw(ArgumentError(
                "voltage source '$sid': ($bus, $t) is aliased to its own neutral"))
            push!(rows, _HelmSourceRow(nds, [acc[i] for i in nds], E, (bus, t)))
        end
    end
    rows
end

# ── restored shunt branches of ideal transformers ───────────────────────────

# `ybus_augmented` drops the whole Yprim of a zero-leakage two-winding
# transformer it promotes to a constraint row (non-unity ratio) or node alias
# (unity ratio) — including its magnetising shunt and its neutral-grounding
# impedance, which are finite and physically present. Re-stamp exactly that
# shunt-only part: the transformer's own Yprim with the series leakage zeroed.
function _helm_ideal_xfmr_shunts(net::Dict{String,Any}, aug::AugYbusResult,
                                 config, ideal_xfmrs::Symbol)
    n = length(aug.nodes)
    I = Int[]; J = Int[]; V = ComplexF64[]
    xfmrs = get(net, "transformer", Dict())
    xfmrs isa AbstractDict || return spzeros(ComplexF64, n, n)
    zmin = _ybus_xfmr_z_min(_domain_thresholds(config))
    for subtype in ("single_phase", "wye_delta", "delta_wye")
        sub = get(xfmrs, subtype, Dict())
        sub isa AbstractDict || continue
        for xid in sort!(collect(String.(keys(sub))))
            x = sub[xid]
            _xfmr_is_ideal(x, zmin) || continue
            unity = isapprox(_xfmr_turns_ratio(x), 1.0; atol=1e-9)
            # Non-unity ideal units under `:stamp` already carry the shunt-only
            # block from `ybus_passive`; restoring it again would double it.
            (unity || ideal_xfmrs === :constrain) || continue
            xs = copy(x)
            for f in ("r_series_from", "x_series_from", "r_series_to", "x_series_to")
                xs[f] = 0.0
            end
            nodes, Y = with_logger(NullLogger()) do
                transformer_yprim(xs, subtype)
            end
            gi = [get(aug.index, nd, 0) for nd in nodes]
            for a in eachindex(gi), b in eachindex(gi)
                (gi[a] == 0 || gi[b] == 0 || iszero(Y[a, b])) && continue
                push!(I, gi[a]); push!(J, gi[b]); push!(V, Y[a, b])
            end
        end
    end
    sparse(I, J, V, n, n)
end

# ── singular-system diagnosis ────────────────────────────────────────────────

# Best-effort explanation for a singular K_UU: BFS the structural pattern of the
# germ matrix and report node groups with no path to any voltage reference. A
# connected-but-singular case (e.g. an ungrounded delta common mode reachable
# only through voltage DIFFERENCES) gets the generic message.
function _helm_floating_error(K::SparseMatrixCSC, nodes::Vector{_Node},
                              seeds, embedding::Symbol)
    n = length(nodes)
    ntot = size(K, 1)
    adj = [Int[] for _ in 1:ntot]
    rows = rowvals(K)
    for j in 1:ntot, p in nzrange(K, j)
        i = rows[p]
        i == j && continue
        push!(adj[i], j); push!(adj[j], i)
    end
    seen = falses(ntot)
    for f in seeds
        seen[f] && continue
        stack = [f]; seen[f] = true
        while !isempty(stack)
            u = pop!(stack)
            for v in adj[u]
                seen[v] || (seen[v] = true; push!(stack, v))
            end
        end
    end
    hint = embedding === :load_scaling ?
        " If the affected nodes are held only by loads (e.g. an unearthed " *
        "star point), use `embedding = :impedance_germ`." : ""
    floating = [nodes[i] for i in 1:n if !seen[i]]
    if isempty(floating)
        error("HELM: the augmented system is singular although every node " *
              "reaches a voltage reference — typically a floating common mode " *
              "(e.g. an ungrounded delta winding or source neutral with no earth " *
              "path). Add a grounding element or reference for the affected " *
              "subnetwork." * hint)
    else
        preview = join(string.(floating[1:min(6, length(floating))]), ", ")
        error("HELM: $(length(floating)) node(s) have no path to any voltage " *
              "source or earth reference (islanded / floating): $preview" *
              (length(floating) > 6 ? ", …" : "") * "." * hint)
    end
end

# ── coefficient diagnostics (Domb–Sykes + classifier tail) ─────────────────

# Domb–Sykes: when a component is dominated by an algebraic singularity at λ,
# its coefficient ratios behave as
#   |c[k]| / |c[k−1]|  ≈  (1/λ)·(1 + b/k).
# A finite tail can be affected by other real/complex singularities and
# numerical noise, so this is a heuristic singularity estimate, not a collapse
# certificate or a guaranteed loading margin.
function _helm_singularity_estimate(coeffs::AbstractMatrix{ComplexF64}, n_nodes::Int)
    N = size(coeffs, 2) - 1          # highest order
    N < 8 && return NaN
    λmin = Inf
    for i in 1:n_nodes
        c = [abs(coeffs[i, j]) for j in 1:N+1]
        scale = max(maximum(c), 1e-300)
        c[end] <= 1e-14 * scale && continue     # converged tail: no signal
        ks = Int[]; rs = Float64[]
        for k in max(3, N - 9):N                # tail ratios r_k = |c_k|/|c_{k−1}|
            (c[k] > 0 && c[k+1] > 0) || continue
            push!(ks, k); push!(rs, c[k+1] / c[k])
        end
        length(ks) < 5 && continue
        xs = 1.0 ./ ks
        sx = sum(xs); sxx = sum(abs2, xs)
        sr = sum(rs); sxr = sum(xs .* rs)
        nfit = length(ks)
        den = nfit * sxx - sx^2
        den <= 0 && continue
        A = (sr * sxx - sx * sxr) / den         # intercept = 1/λ*
        A > 0 || continue
        λ = 1.0 / A
        0 < λ < λmin && (λmin = λ)
    end
    isfinite(λmin) ? λmin : NaN
end

function _helm_tail_diagnostics(coeffs::AbstractMatrix{ComplexF64})
    norms = [maximum(abs, view(coeffs, :, j); init=0.0)
             for j in axes(coeffs, 2)]
    first_order = max(firstindex(norms), lastindex(norms) - 6)
    tail_norms = Float64.(norms[first_order:end])
    tail_ratios = [tail_norms[j] / max(tail_norms[j-1], 1e-300)
                   for j in 2:length(tail_norms)]
    tail_norms, tail_ratios
end

# ── series recursion of one load ─────────────────────────────────────────────

# Extend `ld`'s series with order k of ΔV (taken from `xprev = x[k]`) and return
# the order-k coefficient of its law current Σ κ·m^{γ/2}·(1/d̃).
function _helm_load_step!(ld::_HelmLoad, xprev::AbstractVector{ComplexF64},
                          k::Int, Vscale::Float64)
    dp = ld.up != 0 ? xprev[ld.up] : (k == 0 ? ld.v0p : 0.0im)
    dn = ld.un != 0 ? xprev[ld.un] : (k == 0 ? ld.v0n : 0.0im)
    dk = dp - dn
    push!(ld.d, dk)
    push!(ld.dt, conj(dk))
    if k == 0
        abs(dk) <= 1e-9 * Vscale && error(
            "HELM: load '$(ld.name)' is connected across a node pair with zero " *
            "germ voltage — its S/ΔV* term is singular at the germ.")
        push!(ld.e, 1.0 / ld.dt[1])
    else
        acc = 0.0im
        for j in 1:k
            acc += ld.dt[j+1] * ld.e[k-j+1]
        end
        push!(ld.e, -ld.e[1] * acc)
    end
    if !isempty(ld.p)
        mk = 0.0im
        for j in 0:k
            mk += ld.d[j+1] * ld.dt[k-j+1]
        end
        push!(ld.m, k == 0 ? complex(real(mk)) : mk)   # m[0] = |ΔV₀|² exactly
        m0 = real(ld.m[1])
        for (γ, p) in ld.p
            α = γ / 2
            if k == 0
                push!(p, complex(m0^α))
            else
                # Miller: p[k] = Σ_{j=1..k} (α·j − (k−j))·m[j]·p[k−j] / (k·m[0])
                acc = 0.0im
                for j in 1:k
                    acc += (α * j - (k - j)) * ld.m[j+1] * p[k-j+1]
                end
                push!(p, acc / (k * m0))
            end
        end
    end
    c = 0.0im
    for (κ, γ) in ld.terms
        if γ == 0.0
            c += κ * ld.e[k+1]
        else
            p = ld.p[γ]
            acc = 0.0im
            for j in 0:k
                acc += p[j+1] * ld.e[k-j+1]
            end
            c += κ * acc
        end
    end
    c
end

# ── the solver ───────────────────────────────────────────────────────────────

"""
    helm_series(net; config=_DEFAULT_CONFIG, switches=:alias,
                ideal_xfmrs=:constrain, max_order=40, tol=1e-8,
                embedding=:load_scaling, source_reference=:earth) -> HelmResult

Solve the unbalanced multiconductor power flow of `net` with the Holomorphic
Embedding Load-flow Method on the augmented nodal admittance matrix. Every
explicit conductor (neutral, PE, …) keeps its own node.

Deterministic and non-iterative: no initial guess, one LU factorization for
the germ and every series order. A converged result is the operational branch
continuously connected to the germ. A non-converged result reports whether its
finite coefficient tail grew (`:series_diverged`) or merely exhausted the
requested order (`:max_order_reached`); neither is a proof that no power-flow
solution exists.

Keywords:
- `switches`, `ideal_xfmrs` — forwarded to `ybus_augmented`;
  `switches = :constrain` additionally yields every switch-conductor current.
- `max_order` — highest series order (default 40).
- `tol` — relative convergence tolerance on the s = 1 current mismatch,
  scaled by the largest load-current magnitude.
- `embedding` — `:load_scaling` (default; the germ is the energized no-load
  network and the embedding parameter is a load multiplier) or
  `:impedance_germ` (constant-impedance parts in the matrix, other load parts
  blended from an equivalent admittance at the germ — needed when a node is
  held only by loads, e.g. an unearthed star point).
- `source_reference` — `:earth` (default; source phase terminals fixed to
  their phasors and the source-bus neutral to 0 V, as in the OPF) or
  `:neutral` (phasors imposed phase-to-source-neutral; the source neutral
  floats with the network's earthing, so neutral-earthing impedances at the
  source are honoured).

Supported loads: constant power, constant current, constant impedance, ZIP
and exponential models in WYE, DELTA and SINGLE_PHASE connections. Networks
with generators or IBRs raise an `ArgumentError`: they have no fixed
power-flow law in the data model and are not silently dropped.
"""
function helm_series(net::Dict{String,Any}; config=_DEFAULT_CONFIG,
                     switches::Symbol=:alias, ideal_xfmrs::Symbol=:constrain,
                     max_order::Int=40, tol::Float64=1e-8,
                     embedding::Symbol=:load_scaling,
                     source_reference::Symbol=:earth)::HelmResult
    max_order >= 2 || throw(ArgumentError("max_order must be at least 2"))
    embedding in (:load_scaling, :impedance_germ) || throw(ArgumentError(
        "embedding must be :load_scaling or :impedance_germ, got $(repr(embedding))"))
    source_reference in (:earth, :neutral) || throw(ArgumentError(
        "source_reference must be :earth or :neutral, got $(repr(source_reference))"))

    unsupported = _helm_unsupported_devices(net)
    isempty(unsupported) || throw(ArgumentError(
        "HELM has no power-flow law for " * join(unsupported, "; ") *
        ": generator and IBR set-points are OPF variables or control-law " *
        "outputs in the BMOPF model. Remove them (or model fixed injections as " *
        "negative loads) to solve the remaining network with HELM."))

    aug = ybus_augmented(net; config, switches, ideal_xfmrs)
    n = length(aug.nodes)
    m = length(aug.couplings)

    # ── boundary ──────────────────────────────────────────────────────────────
    fixed = source_reference === :earth ? _helm_fixed_nodes(net, aug) :
                                          Dict{Int,ComplexF64}()
    srows = source_reference === :neutral ? _helm_source_rows(net, aug) :
                                            _HelmSourceRow[]
    r = length(srows)
    (isempty(fixed) && isempty(srows)) && throw(ArgumentError(
        "HELM needs at least one voltage source (no voltage reference found)"))

    # ── base matrix: augmented Ybus + restored transformer shunts + source rows
    Yxs = _helm_ideal_xfmr_shunts(net, aug, config, ideal_xfmrs)
    ntot = n + m + r
    y0 = m > 0 ? aug.couplings[1].scale :
         _helm_median_scale(abs.(Vector(diag(aug.K)))[1:n])
    Is = Int[]; Js = Int[]; Vs = ComplexF64[]
    for (j, sr) in enumerate(srows)
        row = n + m + j
        for (ni, cf) in zip(sr.nodes, sr.coeffs)
            push!(Is, row); push!(Js, ni); push!(Vs, y0 * cf)
            push!(Is, ni); push!(Js, row); push!(Vs, y0 * cf)
        end
    end
    Kbase = blockdiag(aug.K, spzeros(ComplexF64, r, r)) +
            blockdiag(Yxs, spzeros(ComplexF64, m + r, m + r)) +
            sparse(Is, Js, Vs, ntot, ntot)

    Fidx = sort!(collect(keys(fixed)))
    VF = ComplexF64[fixed[i] for i in Fidx]
    isfree = trues(n); isfree[Fidx] .= false
    U = findall(isfree)
    nU = length(U)
    posU = zeros(Int, n)
    for (p, i) in enumerate(U); posU[i] = p; end
    keep = vcat(U, collect(n+1:ntot))
    nkeep = nU + m + r
    Erhs = zeros(ComplexF64, nkeep)                  # source rows: y0·E at the germ
    for (j, sr) in enumerate(srows)
        Erhs[nU + m + j] = y0 * sr.E
    end
    Vsrc = maximum(abs, vcat(VF, [sr.E for sr in srows]); init=0.0)

    # ── loads: constant-Z admittances and law-current sub-loads ──────────────
    loads = get(net, "load", Dict())
    subs = Tuple{String,_SubLoad}[]
    for lid in sort!(collect(String.(keys(loads))))
        append!(subs, (lid, sl) for sl in _load_subloads(loads[lid], net))
    end
    endpoints(sl) = (get(aug.index, sl.pos, 0),
                     sl.neg === nothing ? 0 : get(aug.index, sl.neg, 0))

    Iz = Int[]; Jz = Int[]; Vz = ComplexF64[]
    for (_, sl) in subs
        gp, gn = endpoints(sl)
        _stamp_pair!(Iz, Jz, Vz, gp, gn, _subload_yz(sl))
    end
    YZ = sparse(Iz, Jz, Vz, n, n)
    padn(Y) = blockdiag(Y, spzeros(ComplexF64, m + r, m + r))

    law = [(lid, sl, _helm_current_terms(sl)) for (lid, sl) in subs]
    filter!(t -> !isempty(t[3]), law)

    # Germ matrix additions (Ygerm, into K) and the s-proportional admittance
    # moved to the RHS (Yrhs: order-n injection −Yrhs·V[n−1]).
    Ygerm, Yrhs = if embedding === :load_scaling
        spzeros(ComplexF64, n, n), YZ
    else
        Yg = _helm_germ_admittance(Kbase, YZ, law, endpoints, keep, Fidx, VF,
                                   Erhs, posU, n, Vsrc, padn,
                                   _neutral_labels(net))
        YZ + Yg, -Yg
    end

    Kg = Kbase + padn(Ygerm)
    Kk = Kg[keep, keep]
    b0 = Erhs .- Kg[keep, Fidx] * VF
    seeds = vcat(Fidx, collect(n+m+1:ntot))
    Flu = try
        lu(Kk)
    catch err
        err isa LinearAlgebra.SingularException &&
            _helm_floating_error(Kg, aug.nodes, seeds, embedding)
        rethrow()
    end

    # ── germ ─────────────────────────────────────────────────────────────────
    x0 = Flu \ b0
    all(isfinite, x0) || _helm_floating_error(Kg, aug.nodes, seeds, embedding)
    coeffs = Vector{ComplexF64}[x0]
    Vscale = max(1.0, maximum(abs, view(x0, 1:nU); init=0.0), Vsrc)

    Yrhs_uu = Yrhs[U, U]
    Yrhs_uf = Yrhs[U, Fidx]
    hasYrhs = nnz(Yrhs_uu) > 0 || nnz(Yrhs_uf) > 0

    hl = _HelmLoad[]
    for (lid, sl, terms) in law
        gp, gn = endpoints(sl)
        v0p = gp == 0 ? 0.0im : get(fixed, gp, 0.0im)
        v0n = gn == 0 ? 0.0im : get(fixed, gn, 0.0im)
        pw = Dict{Float64,Vector{ComplexF64}}(γ => ComplexF64[]
                                              for (_, γ) in terms if γ != 0.0)
        push!(hl, _HelmLoad(gp == 0 ? 0 : posU[gp], gn == 0 ? 0 : posU[gn],
                            v0p, v0n, terms, ComplexF64[], ComplexF64[],
                            ComplexF64[], ComplexF64[], pw, lid))
    end

    # ── order-n recursion: one back-substitution per order ───────────────────
    n_order = 0
    for ord in 1:max_order
        k = ord - 1                       # the ΔV order the RHS consumes
        xprev = coeffs[end]               # x[k]
        rhs = zeros(ComplexF64, nkeep)
        for ld in hl
            c = _helm_load_step!(ld, xprev, k, Vscale)   # drawn-current coefficient
            ld.up != 0 && (rhs[ld.up] -= c)
            ld.un != 0 && (rhs[ld.un] += c)
        end
        if hasYrhs
            rhs[1:nU] .-= Yrhs_uu * view(xprev, 1:nU)
            ord == 1 && (rhs[1:nU] .-= Yrhs_uf * VF)
        end
        xn = Flu \ rhs
        push!(coeffs, xn)
        n_order = ord
        # early stop: the tail is numerically negligible
        maximum(abs, xn; init=0.0) <= 1e-3 * tol * Vscale && break
    end

    Cmat = reduce(hcat, coeffs)           # nkeep × (n_order+1)

    # ── Padé evaluation at s = 1, per component ──────────────────────────────
    vals = Vector{ComplexF64}(undef, nkeep)
    pade_spread = Vector{Float64}(undef, nkeep)
    for i in 1:nkeep
        vals[i], pade_spread[i] = _wynn_epsilon(cumsum(Cmat[i, :]))
    end

    # ── assemble the full state ──────────────────────────────────────────────
    Vfull = zeros(ComplexF64, n)
    for (p, i) in enumerate(U); Vfull[i] = vals[p]; end
    for (i, v) in fixed;        Vfull[i] = v;       end
    w = vals[nU+1:nU+m]
    ws = vals[nU+m+1:end]
    z = vcat(Vfull, w, ws)

    Vdict = Dict{_Node,ComplexF64}()
    for (nd, gi) in aug.index
        Vdict[nd] = gi == 0 ? 0.0 + 0.0im : Vfull[gi]
    end

    # ── physical convergence check: nonlinear current mismatch at s = 1 ──────
    # Node rows: (Y + restored shunts)·V + coupling and source columns − load
    # injection; constraint rows: a·V − (E for source rows). Dirichlet-fixed
    # nodes are the slack — excluded (their mismatch is the source injection).
    Kz = Kbase * z
    mis = Vector{ComplexF64}(Kz[1:n])
    Iscale = 1.0
    for (_, sl) in subs
        gp, gn = endpoints(sl)
        dv = (gp == 0 ? 0.0im : Vfull[gp]) - (gn == 0 ? 0.0im : Vfull[gn])
        a = abs(dv)
        a > 0 || continue
        Iload = conj(_subload_S(sl, a)) / conj(dv)      # drawn into the load
        Iscale = max(Iscale, abs(Iload))
        gp != 0 && (mis[gp] += Iload)                   # injection = −draw
        gn != 0 && (mis[gn] -= Iload)
    end
    residual = 0.0
    for i in 1:n
        haskey(fixed, i) && continue
        residual = max(residual, abs(mis[i]))
    end
    for j in 1:m
        residual = max(residual, abs(Kz[n+j]))
    end
    for (j, sr) in enumerate(srows)
        residual = max(residual, abs(Kz[n+m+j] - y0 * sr.E))
    end
    converged = residual <= tol * Iscale

    # ── boundary currents ────────────────────────────────────────────────────
    source_injection = Dict{_Node,ComplexF64}()
    for i in Fidx
        source_injection[aug.nodes[i]] = mis[i]
    end
    for (j, sr) in enumerate(srows), (ni, cf) in zip(sr.nodes, sr.coeffs)
        nd = aug.nodes[ni]
        source_injection[nd] = get(source_injection, nd, 0.0im) - y0 * cf * ws[j]
    end

    # ── status taxonomy + auditable coefficient diagnostics ─────────────────
    coefficient_tail_norms, coefficient_tail_ratios =
        _helm_tail_diagnostics(Cmat)
    status = if converged
        :converged
    else
        (!isempty(coefficient_tail_ratios) &&
         all(>(1.0), coefficient_tail_ratios)) ? :series_diverged :
                                                 :max_order_reached
    end

    singularity_estimate = _helm_singularity_estimate(Cmat, nU)

    HelmResult(Vdict, w, aug.couplings, Cmat, converged, status,
               residual, n_order, pade_spread, coefficient_tail_norms,
               coefficient_tail_ratios, singularity_estimate, source_injection,
               embedding, source_reference)
end

# Conditioning scale of the source constraint rows when the network has no
# ideal couplings to borrow one from: the median positive diagonal magnitude,
# floored at 1 S (the `ybus_augmented` heuristic).
function _helm_median_scale(d::Vector{Float64})
    filter!(>(0.0), d)
    isempty(d) ? 1.0 : max(1.0, sort!(d)[cld(length(d), 2)])
end

# Equivalent germ admittance of every law sub-load, y_g = I_law(Vref)/Vref at a
# real reference voltage. Any non-zero y_g yields the same s = 1 solution; a
# good one keeps the series short. Vref is |ΔV| of the no-load (constant-Z
# loaded) network where that network is solvable and the pair is energized,
# else the sub-load's nominal voltage, else the largest source phasor (×√3
# for a pair of two non-earth, non-neutral terminals).
function _helm_germ_admittance(Kbase, YZ, law, endpoints, keep, Fidx, VF, Erhs,
                               posU, n, Vsrc, padn, nlabels)
    x_nl = nothing
    try
        K0 = Kbase + padn(YZ)
        x_nl = lu(K0[keep, keep]) \ (Erhs .- K0[keep, Fidx] * VF)
        all(isfinite, x_nl) || (x_nl = nothing)
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
    end
    Vfix = Dict{Int,ComplexF64}(zip(Fidx, VF))
    vnode(gi) = gi == 0 ? 0.0im :
                posU[gi] != 0 ? (x_nl === nothing ? NaN * im : x_nl[posU[gi]]) :
                Vfix[gi]
    Ig = Int[]; Jg = Int[]; Vg = ComplexF64[]
    for (_, sl, terms) in law
        gp, gn = endpoints(sl)
        vref = abs(vnode(gp) - vnode(gn))
        if !(isfinite(vref) && vref > 1e-6 * max(Vsrc, 1.0))
            vref = sl.Vnom != 1.0 ? sl.Vnom :
                   (gn != 0 && !(sl.neg[2] in nlabels) ? sqrt(3.0) : 1.0) * Vsrc
        end
        yg = _helm_law_current(terms, complex(vref)) / vref
        iszero(yg) && continue
        _stamp_pair!(Ig, Jg, Vg, gp, gn, yg)
    end
    sparse(Ig, Jg, Vg, n, n)
end

# ── result-dict wrapper ──────────────────────────────────────────────────────

"""
    solve_pf_helm(net; config=_DEFAULT_CONFIG, switches=:alias,
                  ideal_xfmrs=:constrain, max_order=40, tol=1e-8,
                  embedding=:load_scaling, source_reference=:earth)
        -> Dict{String,Any}

Solve the power flow with HELM ([`helm_series`](@ref)) and return the standard
result dictionary (SI units, same `"bus"` shape as the OPF results, compatible
with `write_result`/`read_result`).

Top-level keys:
- `"termination_status"` — `"HELM_CONVERGED"` |
  `"HELM_SERIES_DIVERGED"` (the computed coefficient tail is growing) |
  `"HELM_MAX_ORDER"` (series order exhausted before the tolerance was met).
  Neither non-converged status is a proof of power-flow non-existence.
- `"feasible"` — `true` iff converged.
- `"solve_time"` — wall-clock seconds.
- `"bus"` — `bus_id => terminal => {vr, vi, vm [V], va [rad]}`. NaN-filled
  when HELM does not converge (the standard infeasible convention).
- `"coupling"` — `id => conductor => {kind, ir, ii [A], im [A]}`: the PHYSICAL
  current through each ideal coupling (closed-switch conductor with
  `switches = :constrain`, ideal-transformer winding core).
- `"source_injection"` — `bus_id => terminal => {ir, ii, im [A]}`: current
  injected into the network by the voltage-source boundary.
- `"helm"` — diagnostics: `order`, `residual` (A), raw `pade_spread`,
  `coefficient_tail_norms`, `coefficient_tail_ratios`, the heuristic
  `singularity_estimate`, and the `embedding` / `source_reference` used.
  `load_margin` is retained as a compatibility alias.
"""
function solve_pf_helm(net::Dict{String,Any}; config=_DEFAULT_CONFIG,
                       switches::Symbol=:alias, ideal_xfmrs::Symbol=:constrain,
                       max_order::Int=40, tol::Float64=1e-8,
                       embedding::Symbol=:load_scaling,
                       source_reference::Symbol=:earth)::Dict{String,Any}
    t0 = time()
    hr = helm_series(net; config, switches, ideal_xfmrs, max_order, tol,
                     embedding, source_reference)
    dt = time() - t0

    ok = hr.converged
    num(x) = ok ? x : NaN                    # infeasible ⇒ NaN-filled numerics

    bus = Dict{String,Any}()
    for (nd, v) in hr.V
        b, t = nd
        bt = get!(bus, b, Dict{String,Any}())
        bt[t] = Dict{String,Any}(
            "vr" => num(real(v)), "vi" => num(imag(v)),
            "vm" => num(abs(v)),  "va" => num(angle(v)))
    end

    coupling = Dict{String,Any}()
    for (j, c) in enumerate(hr.couplings)
        iw = c.scale * hr.w[j]
        ct = get!(coupling, c.id, Dict{String,Any}())
        ct[string(c.conductor)] = Dict{String,Any}(
            "kind" => string(c.kind),
            "ir" => num(real(iw)), "ii" => num(imag(iw)), "im" => num(abs(iw)))
    end

    source_injection = Dict{String,Any}()
    for ((b, t), i) in hr.source_injection
        bt = get!(source_injection, b, Dict{String,Any}())
        bt[t] = Dict{String,Any}(
            "ir" => num(real(i)), "ii" => num(imag(i)), "im" => num(abs(i)))
    end

    status = hr.status === :converged            ? "HELM_CONVERGED" :
             hr.status === :series_diverged      ? "HELM_SERIES_DIVERGED" :
                                                   "HELM_MAX_ORDER"

    Dict{String,Any}(
        "termination_status" => status,
        "feasible"           => ok,
        "solve_time"         => dt,
        "bus"                => bus,
        "coupling"           => coupling,
        "source_injection"   => source_injection,
        "helm" => Dict{String,Any}(
            "order"            => hr.n_order,
            "residual"         => hr.residual,
            "pade_spread"      => copy(hr.pade_spread),
            "coefficient_tail_norms" => copy(hr.coefficient_tail_norms),
            "coefficient_tail_ratios" => copy(hr.coefficient_tail_ratios),
            "singularity_estimate" => hr.singularity_estimate,
            "load_margin"      => hr.singularity_estimate,
            "embedding"        => string(hr.embedding),
            "source_reference" => string(hr.source_reference),
            "series_tail_norm" => maximum(abs, view(hr.coeffs, :, size(hr.coeffs, 2));
                                          init=0.0),
        ),
    )
end
