# powerflow/pade.jl
#
# Analytic continuation of a power series to its expansion boundary via Wynn's
# epsilon algorithm — the workhorse of the HELM power-flow solver.
#
# Given the partial sums s_j = Σ_{k≤j} c[k] of a series at the evaluation point,
# the epsilon table
#
#     ε^{(j)}_{-1} = 0,   ε^{(j)}_0 = s_j,
#     ε^{(j)}_{k+1} = ε^{(j+1)}_{k-1} + 1 / (ε^{(j+1)}_k − ε^{(j)}_k)
#
# has even columns ε^{(j)}_{2k} equal to the [j+k / k] Padé approximants of the
# series, evaluated at the point. When the finite diagonal sequence stabilizes,
# it can recover analytic continuations whose plain partial sums diverge. A
# large finite-order spread is useful numerical evidence of instability, but is
# not by itself a proof that the underlying function or power-flow solution
# does not exist.
#
# Self-contained: no dependencies beyond Base complex arithmetic.

"""
    _wynn_epsilon(psums) -> (value, spread)

Accelerate a sequence of partial sums with Wynn's epsilon algorithm.

Returns the highest-order even-column epsilon estimate (`value`) and `spread`,
the absolute difference between the last two even-column estimates — a
finite-order convergence indicator. `spread ≈ 0` means the available Padé
sequence has stabilized; a large spread means the continuation is numerically
unresolved at this order, not that non-existence has been proved.

Guards: a (near-)zero difference between adjacent entries means the sequence
has already converged — the algorithm returns that entry with `spread = 0`
rather than dividing by ~0 (the standard Wynn degeneracy guard). A single
partial sum returns `(psums[1], Inf)` (no acceleration information).
"""
function _wynn_epsilon(psums::AbstractVector{<:Complex})
    n = length(psums)
    n == 0 && throw(ArgumentError("empty partial-sum sequence"))
    n == 1 && return (ComplexF64(psums[1]), Inf)

    e_prev = zeros(ComplexF64, n + 1)      # ε_{k-1} column (ε_{-1} ≡ 0)
    e_curr = ComplexF64.(psums)            # ε_k column (k = 0)
    estimates = ComplexF64[e_curr[end]]    # even-column history

    k = 0
    while length(e_curr) > 1
        e_next = Vector{ComplexF64}(undef, length(e_curr) - 1)
        for i in 1:length(e_curr)-1
            d = e_curr[i+1] - e_curr[i]
            scale = abs(e_curr[i+1]) + abs(e_curr[i])
            if abs(d) <= eps(Float64) * scale + floatmin(Float64)
                # Degenerate step: the sequence has converged at this entry.
                return (e_curr[i+1], 0.0)
            end
            e_next[i] = e_prev[i+1] + 1.0 / d
        end
        e_prev = e_curr
        e_curr = e_next
        k += 1
        iseven(k) && push!(estimates, e_curr[end])
    end

    value = estimates[end]
    spread = length(estimates) >= 2 ? abs(estimates[end] - estimates[end-1]) : Inf
    (value, spread)
end

"""
    _pade_sum(coeffs) -> (value, spread)

Evaluate a power series `Σ coeffs[k+1]·sᵏ` at `s = 1` by analytic continuation:
build the partial sums and accelerate them with [`_wynn_epsilon`](@ref).

This can converge to the underlying function's value at `s = 1` even when
`s = 1` lies outside the series' radius of convergence, and returns the last
two-approximant `spread` so callers can inspect finite-order stability.
"""
function _pade_sum(coeffs::AbstractVector{<:Complex})
    isempty(coeffs) && throw(ArgumentError("empty coefficient sequence"))
    psums = cumsum(ComplexF64.(coeffs))
    _wynn_epsilon(psums)
end
