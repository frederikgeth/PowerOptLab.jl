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
# series, evaluated at the point. For a series with a genuine limit (Stahl's
# convergence-in-capacity for the HELM voltage series) the even-column diagonal
# converges to it — including for series whose plain partial sums DIVERGE, as
# long as the underlying function is analytic at the evaluation point. When no
# limit exists (HELM: no power-flow solution / voltage collapse) the even
# columns stagnate at a large spread instead, which is exactly the signal the
# solver's failure taxonomy uses.
#
# Self-contained: no dependencies beyond Base complex arithmetic.

"""
    _wynn_epsilon(psums) -> (value, spread)

Accelerate a sequence of partial sums with Wynn's epsilon algorithm.

Returns the highest-order even-column epsilon estimate (`value`) and `spread`,
the absolute difference between the last two even-column estimates — the
convergence indicator: `spread ≈ 0` means the (Padé) continuation has locked
in; a `spread` that stagnates at a large value across increasing series order
means the series has no limit at the evaluation point.

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

This converges to the underlying function's value at `s = 1` whenever that
value exists — even when `s = 1` lies outside the series' radius of convergence
(the HELM situation for heavily loaded feeders) — and reports a stagnating
`spread` when it does not (no power-flow solution).
"""
function _pade_sum(coeffs::AbstractVector{<:Complex})
    isempty(coeffs) && throw(ArgumentError("empty coefficient sequence"))
    psums = cumsum(ComplexF64.(coeffs))
    _wynn_epsilon(psums)
end
