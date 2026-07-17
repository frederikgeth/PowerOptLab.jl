# powerflow/helm.jl
#
# HELM — Holomorphic Embedding Load-flow Method — for 4-wire unbalanced
# multiphase distribution power flow, on the augmented nodal admittance matrix
# (`ybus_augmented`: ideal switches / ideal transformers as constraint rows).
#
# The power-flow equations are embedded in a complex parameter `s`:
#
#   • voltage sources stay FIXED for all s (Dirichlet boundary, eliminated
#     from the unknown vector);
#   • every load is scaled by s: constant-power draws s·S, constant-impedance
#     admittances become s·Y_Z.
#
# At s = 0 the network is the energized NO-LOAD state (the germ — genuinely
# unbalanced, floating neutrals at their equilibrium; one linear solve). At
# s = 1 it is the actual problem. Each voltage (and coupling current) is a
# holomorphic function of s, expanded as a power series whose coefficients
# solve ONE constant linear system per order:
#
#   K_UU · x[n] = rhs(x[0..n−1])
#
# with the same LU factorization for the germ and every order. The nonlinear
# constant-power term S*/ΔV* is handled with the classic conjugate-reflection
# trick: the reflected series d̃[k] = conj(d[k]) is itself holomorphic, and its
# reciprocal e = 1/d̃ has an incremental convolution recursion, so the order-n
# RHS needs only coefficients 0..n−1. Constraint rows (ideal couplings) have
# identically zero RHS at every order — the ideal-element identities hold
# term-by-term and therefore exactly in the summed solution.
#
# The series is evaluated at s = 1 by Padé analytic continuation (Wynn epsilon,
# powerflow/pade.jl). Finite-order Padé spread, physical mismatch, and
# coefficient-tail growth are numerical diagnostics: they do not by themselves
# prove non-existence of a power-flow solution. A Domb–Sykes coefficient-ratio
# extrapolation is reported as a `singularity_estimate`; it needs validation
# against continuation or an analytically known system before being interpreted
# as a voltage-collapse loading multiplier.
#
# v1 scope: constant-power and constant-impedance load parts (incl. their ZIP
# fractions). Constant-current fractions and non-integer exponential loads are
# NOT holomorphic in this embedding and raise a validation error — use the
# `ybus_linearized` fixed-point map or the OPF path for those, or wait for the
# outer-loop follow-up. Voltage-source configurations: WYE / SINGLE_PHASE.
# NOTE: a transformer promoted to an ideal coupling is modelled as PURELY
# ideal — its no-load (magnetising) shunt and neutral-grounding branch are
# dropped along with the singular series stamp.

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
- `coeffs`      — `(n_unknown + n_couplings) × (order+1)` series coefficients
                  (diagnostic; row order = non-fixed nodes then couplings).
- `converged`   — `true` iff the s = 1 nonlinear current mismatch is within
                  tolerance.
- `status`      — `:converged` | `:series_diverged` (the finite coefficient
                  tail is growing) | `:max_order_reached`.
- `residual`    — max |current mismatch| (A) over non-source nodes and
                  constraint rows, evaluated at the returned solution.
- `n_order`     — highest series order computed.
- `pade_spread` — absolute difference between the last two Wynn-epsilon Padé
                  estimates for each row of `coeffs` (node rows, then coupling
                  rows). Values have the units of their corresponding row.
- `coefficient_tail_norms` — max-norm of the final coefficient orders used by
                  the status classifier.
- `coefficient_tail_ratios` — adjacent ratios of those tail norms. A tail with
                  every ratio above one is classified as `:series_diverged`.
- `singularity_estimate` — heuristic Domb–Sykes estimate of the nearest
                  coefficient-dominating singularity in the loading parameter.
                  It is not a certified loading margin. `NaN` means the series
                  is too short or featureless to extrapolate.

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
     load_margin=result.load_margin)

Base.show(io::IO, r::HelmResult) =
    print(io, "HelmResult($(length(r.V)) terminals, status=$(r.status), " *
              "order=$(r.n_order), residual=$(round(r.residual; sigdigits=3)) A, " *
              "singularity_estimate=$(round(r.singularity_estimate; sigdigits=4)))")

# One constant-power sub-load in the recursion: endpoint positions in the
# unknown vector (0 = fixed/earth endpoint, with `v0p`/`v0n` carrying the
# order-0 voltage), the conjugated constant power, and the growing reflected
# ΔV series `dt` and its reciprocal series `e`.
mutable struct _HelmPQ
    up::Int
    un::Int
    v0p::ComplexF64
    v0n::ComplexF64
    Sc::ComplexF64                # conj(P + jQ) of the constant-power part
    dt::Vector{ComplexF64}        # d̃[k] = conj(ΔV[k])
    e::Vector{ComplexF64}         # series of 1/d̃
    name::String                  # load id (for error messages)
end

# ── validation ───────────────────────────────────────────────────────────────

# HELM v1 supports the holomorphic load parts only: constant-power (cc) and
# constant-impedance (cW). Constant-current (cs) and non-integer exponential
# (nl) involve |ΔV| — not holomorphic in the standard embedding.
_helm_subload_ok(sl::_SubLoad) =
    sl.pt.cs == 0.0 && sl.qt.cs == 0.0 && sl.pt.nl === nothing && sl.qt.nl === nothing

# ── boundary (fixed nodes) ───────────────────────────────────────────────────

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

# Collect the Dirichlet boundary: voltage-source phase terminals at their
# reference phasors, source-bus neutral-labelled terminals at 0 V (mirrors the
# OPF's `_add_source_constraints!` semantics). Returns (fixed, source_buses).
function _helm_fixed_nodes(net::Dict{String,Any}, aug::AugYbusResult)
    fixed = Dict{Int,ComplexF64}()
    srcbuses = Set{String}()
    nlabels = _neutral_labels(net)
    buses = get(net, "bus", Dict())
    vss = get(net, "voltage_source", Dict())
    for sid in sort!(collect(String.(keys(vss))))
        vs = vss[sid]
        bus = string(get(vs, "bus", ""))
        push!(srcbuses, bus)
        cfg = uppercase(string(get(vs, "configuration", "WYE")))
        cfg in ("WYE", "SINGLE_PHASE") || throw(ArgumentError(
            "voltage source '$sid': configuration '$cfg' is not supported by " *
            "HELM (WYE / SINGLE_PHASE phase-to-neutral references only)"))
        tm    = Vector{String}(string.(get(vs, "terminal_map", String[])))
        v_mag = Float64.(get(vs, "v_magnitude", Float64[]))
        v_ang = Float64.(get(vs, "v_angle",     Float64[]))
        for (k, t) in enumerate(tm)
            (length(v_mag) >= k && length(v_ang) >= k) || continue
            _helm_fix!(fixed, aug, (bus, t), v_mag[k] * cis(v_ang[k]), sid)
        end
        # Source-bus neutral(s) → 0 V (system ground), as in the OPF.
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
    (fixed, srcbuses)
end

# ── singular-system diagnosis ────────────────────────────────────────────────

# Best-effort explanation for a singular K_UU: BFS the structural pattern of K
# and report node groups with no path to any fixed node. A connected-but-
# singular case (e.g. an ungrounded delta common mode reachable only through
# voltage DIFFERENCES) gets the generic message.
function _helm_floating_error(aug::AugYbusResult, fixed::Dict{Int,ComplexF64})
    n = length(aug.nodes)
    ntot = size(aug.K, 1)
    adj = [Int[] for _ in 1:ntot]
    rows = rowvals(aug.K)
    for j in 1:ntot, p in nzrange(aug.K, j)
        i = rows[p]
        i == j && continue
        push!(adj[i], j); push!(adj[j], i)
    end
    seen = falses(ntot)
    for f in keys(fixed)
        seen[f] && continue
        stack = [f]; seen[f] = true
        while !isempty(stack)
            u = pop!(stack)
            for v in adj[u]
                seen[v] || (seen[v] = true; push!(stack, v))
            end
        end
    end
    floating = [aug.nodes[i] for i in 1:n if !seen[i]]
    if isempty(floating)
        error("HELM: the augmented system is singular although every node " *
              "reaches a voltage reference — typically a floating common mode " *
              "(e.g. an ungrounded delta winding with no earth path). Add a " *
              "grounding element or reference for the affected subnetwork.")
    else
        preview = join(string.(floating[1:min(6, length(floating))]), ", ")
        error("HELM: $(length(floating)) node(s) have no path to any voltage " *
              "source or earth reference (islanded / floating): $preview" *
              (length(floating) > 6 ? ", …" : ""))
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

# ── the solver ───────────────────────────────────────────────────────────────

"""
    helm_series(net; config=_DEFAULT_CONFIG, switches=:alias,
                ideal_xfmrs=:constrain, max_order=40, tol=1e-8) -> HelmResult

Solve the 4-wire multiphase power flow of `net` with the Holomorphic Embedding
Load-flow Method on the augmented nodal admittance matrix.

Deterministic and non-iterative: no initial guess, one LU factorization for
the germ and every series order. A converged result is the operational branch
continuously connected to the no-load state. A non-converged result reports
whether its finite coefficient tail grew (`:series_diverged`) or merely
exhausted the requested order (`:max_order_reached`); neither is a proof that
no power-flow solution exists.

Keywords:
- `switches`, `ideal_xfmrs` — forwarded to `ybus_augmented`;
  `switches = :constrain` additionally yields every switch-conductor current.
- `max_order` — highest series order (default 40).
- `tol` — relative convergence tolerance on the s = 1 current mismatch,
  scaled by the largest load-current magnitude.

Requirements/limitations (v1): at least one WYE / SINGLE_PHASE voltage source;
loads restricted to their constant-power + constant-impedance (ZIP Z/P) parts —
constant-current fractions and non-integer exponential models raise an
`ArgumentError` naming the offending loads.
"""
function helm_series(net::Dict{String,Any}; config=_DEFAULT_CONFIG,
                     switches::Symbol=:alias, ideal_xfmrs::Symbol=:constrain,
                     max_order::Int=40, tol::Float64=1e-8)::HelmResult
    max_order >= 2 || throw(ArgumentError("max_order must be at least 2"))

    aug = ybus_augmented(net; config, switches, ideal_xfmrs)
    n = length(aug.nodes)
    m = length(aug.couplings)

    # ── loads: validate + split into Y_Z stamps and constant-P convolutions ──
    loads = get(net, "load", Dict())
    bad = String[]
    subs = Tuple{String,_SubLoad}[]
    for lid in sort!(collect(String.(keys(loads))))
        sls = _load_subloads(loads[lid], net)
        any(sl -> !_helm_subload_ok(sl), sls) && push!(bad, lid)
        append!(subs, (lid, sl) for sl in sls)
    end
    isempty(bad) || throw(ArgumentError(
        "HELM v1 supports constant-power and constant-impedance load parts " *
        "only; load(s) $(join(bad, ", ")) have constant-current fractions or " *
        "non-integer exponential models (|ΔV| terms are not holomorphic in " *
        "the load-scaling embedding). Use the ybus_linearized fixed-point " *
        "map or the OPF power flow for these."))

    # ── boundary ──────────────────────────────────────────────────────────────
    fixed, srcbuses = _helm_fixed_nodes(net, aug)
    isempty(fixed) && throw(ArgumentError(
        "HELM needs at least one voltage source (no fixed-voltage boundary found)"))

    Fidx = sort!(collect(keys(fixed)))
    VF = ComplexF64[fixed[i] for i in Fidx]
    isfree = trues(n); isfree[Fidx] .= false
    U = findall(isfree)
    nU = length(U)
    posU = zeros(Int, n)
    for (p, i) in enumerate(U); posU[i] = p; end
    keep = vcat(U, collect(n+1:n+m))
    nkeep = nU + m

    K = aug.K
    Kk = K[keep, keep]
    b0 = -(K[keep, Fidx] * VF)

    Flu = try
        lu(Kk)
    catch err
        err isa LinearAlgebra.SingularException && _helm_floating_error(aug, fixed)
        rethrow()
    end

    # ── germ: energized no-load state ────────────────────────────────────────
    x0 = Flu \ b0
    all(isfinite, x0) || _helm_floating_error(aug, fixed)
    coeffs = Vector{ComplexF64}[x0]
    Vscale = max(1.0, maximum(abs, x0; init=0.0), maximum(abs, VF; init=0.0))

    # ── constant-Z: s·Y_Z on the RHS recursion ───────────────────────────────
    Iz = Int[]; Jz = Int[]; Vz = ComplexF64[]
    for (_, sl) in subs
        gp = get(aug.index, sl.pos, 0)
        gn = sl.neg === nothing ? 0 : get(aug.index, sl.neg, 0)
        _stamp_pair!(Iz, Jz, Vz, gp, gn, _subload_yz(sl))
    end
    YZ = sparse(Iz, Jz, Vz, n, n)
    YZ_uu = YZ[U, U]
    YZ_uf = YZ[U, Fidx]

    # ── constant-P sub-loads ─────────────────────────────────────────────────
    pqs = _HelmPQ[]
    for (lid, sl) in subs
        Sp = complex(sl.pt.cc, sl.qt.cc)
        iszero(Sp) && continue
        gp = get(aug.index, sl.pos, 0)
        gn = sl.neg === nothing ? 0 : get(aug.index, sl.neg, 0)
        up = gp == 0 ? 0 : posU[gp]
        un = gn == 0 ? 0 : posU[gn]
        v0p = gp == 0 ? 0.0im : (haskey(fixed, gp) ? fixed[gp] : 0.0im)
        v0n = gn == 0 ? 0.0im : (haskey(fixed, gn) ? fixed[gn] : 0.0im)
        push!(pqs, _HelmPQ(up, un, v0p, v0n, conj(Sp),
                           ComplexF64[], ComplexF64[], lid))
    end

    # ── order-n recursion: one back-substitution per order ───────────────────
    n_order = 0
    for ord in 1:max_order
        k = ord - 1                       # the ΔV order the RHS consumes
        xprev = coeffs[end]               # x[k]
        rhs = zeros(ComplexF64, nkeep)

        for pq in pqs
            # extend d̃ with order-k: unknowns from xprev, fixed only at k = 0
            dp = pq.up != 0 ? xprev[pq.up] : (k == 0 ? pq.v0p : 0.0im)
            dn = pq.un != 0 ? xprev[pq.un] : (k == 0 ? pq.v0n : 0.0im)
            push!(pq.dt, conj(dp - dn))
            if k == 0
                abs(pq.dt[1]) <= 1e-9 * Vscale && error(
                    "HELM: load '$(pq.name)' is connected across a de-energized " *
                    "node pair (zero no-load voltage) — its S/ΔV* term is " *
                    "singular at the germ.")
                push!(pq.e, 1.0 / pq.dt[1])
            else
                acc = 0.0 + 0.0im
                for j in 1:k
                    acc += pq.dt[j+1] * pq.e[k-j+1]
                end
                push!(pq.e, -pq.e[1] * acc)
            end
            # drawn current s·conj(S)·e(s): order-ord coefficient = Sc·e[k]
            c = pq.Sc * pq.e[k+1]
            pq.up != 0 && (rhs[pq.up] -= c)
            pq.un != 0 && (rhs[pq.un] += c)
        end

        # constant-Z: injection −s·Y_Z·V(s) → order-ord term −(Y_Z·V)[k]
        if nnz(YZ_uu) > 0 || nnz(YZ_uf) > 0
            xu_prev = view(xprev, 1:nU)
            rhs[1:nU] .-= YZ_uu * xu_prev
            ord == 1 && (rhs[1:nU] .-= YZ_uf * VF)
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

    # ── assemble the full voltage state ──────────────────────────────────────
    Vfull = zeros(ComplexF64, n)
    for (p, i) in enumerate(U); Vfull[i] = vals[p]; end
    for (i, v) in fixed;        Vfull[i] = v;       end
    w = vals[nU+1:end]

    Vdict = Dict{_Node,ComplexF64}()
    for (nd, gi) in aug.index
        Vdict[nd] = gi == 0 ? 0.0 + 0.0im : Vfull[gi]
    end

    # ── physical convergence check: nonlinear current mismatch at s = 1 ──────
    # Node rows: Y·V + (coupling columns)·w − load injection; constraint rows:
    # a·V (must vanish). Source-bus rows are the slack — excluded.
    mis = Vector{ComplexF64}(K[1:n, 1:n] * Vfull)
    for (j, c) in enumerate(aug.couplings)
        for (ni, cf) in zip(c.nodes, c.coeffs)
            mis[ni] += c.scale * cf * w[j]
        end
    end
    Iscale = 1.0
    for (_, sl) in subs
        gp = get(aug.index, sl.pos, 0)
        gn = sl.neg === nothing ? 0 : get(aug.index, sl.neg, 0)
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
        (haskey(fixed, i) || aug.nodes[i][1] in srcbuses) && continue
        residual = max(residual, abs(mis[i]))
    end
    for c in aug.couplings                              # constraint rows: a·V = 0
        acc = sum(cf * Vfull[ni] for (ni, cf) in zip(c.nodes, c.coeffs); init=0.0im)
        residual = max(residual, c.scale * abs(acc))
    end
    converged = residual <= tol * Iscale

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
               coefficient_tail_ratios, singularity_estimate)
end

# ── result-dict wrapper ──────────────────────────────────────────────────────

"""
    solve_pf_helm(net; config=_DEFAULT_CONFIG, switches=:alias,
                  ideal_xfmrs=:constrain, max_order=40, tol=1e-8)
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
- `"helm"` — diagnostics: `order`, `residual` (A), raw `pade_spread`,
  `coefficient_tail_norms`, `coefficient_tail_ratios`, and the heuristic
  `singularity_estimate`. `load_margin` is retained as a compatibility alias.
"""
function solve_pf_helm(net::Dict{String,Any}; config=_DEFAULT_CONFIG,
                       switches::Symbol=:alias, ideal_xfmrs::Symbol=:constrain,
                       max_order::Int=40, tol::Float64=1e-8)::Dict{String,Any}
    t0 = time()
    hr = helm_series(net; config, switches, ideal_xfmrs, max_order, tol)
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

    status = hr.status === :converged            ? "HELM_CONVERGED" :
             hr.status === :series_diverged      ? "HELM_SERIES_DIVERGED" :
                                                   "HELM_MAX_ORDER"

    Dict{String,Any}(
        "termination_status" => status,
        "feasible"           => ok,
        "solve_time"         => dt,
        "bus"                => bus,
        "coupling"           => coupling,
        "helm" => Dict{String,Any}(
            "order"            => hr.n_order,
            "residual"         => hr.residual,
            "pade_spread"      => copy(hr.pade_spread),
            "coefficient_tail_norms" => copy(hr.coefficient_tail_norms),
            "coefficient_tail_ratios" => copy(hr.coefficient_tail_ratios),
            "singularity_estimate" => hr.singularity_estimate,
            "load_margin"      => hr.singularity_estimate,
            "series_tail_norm" => maximum(abs, view(hr.coeffs, :, size(hr.coeffs, 2));
                                          init=0.0),
        ),
    )
end
