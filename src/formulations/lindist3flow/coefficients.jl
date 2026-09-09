"""
    cross_voltage_coefficients(vbar_phi, vbar_psi)

Return the general fixed-angle affine closure for
`v_phi * conj(v_psi)` as `constant + coefficient_phi*w_phi +
coefficient_psi*w_psi`. Both reference phasors must be nonzero.
"""
function cross_voltage_coefficients(vbar_phi::Number, vbar_psi::Number)
    vp, vq = ComplexF64(vbar_phi), ComplexF64(vbar_psi)
    iszero(vp) && throw(ArgumentError("vbar_phi must be nonzero"))
    iszero(vq) && throw(ArgumentError("vbar_psi must be nonzero"))
    a = conj(vq) / (2conj(vp))
    b = vp / (2vq)
    c = vp * conj(vq) - a * abs2(vp) - b * abs2(vq)
    CrossVoltageCoefficients(c, a, b)
end

"""
    evaluate_cross_voltage(c, w_phi, w_psi)

Evaluate the closure of [`cross_voltage_coefficients`](@ref) at squared
magnitudes `w_phi`, `w_psi`. At the reference point it reproduces
`vbar_phi * conj(vbar_psi)` exactly; away from it, it is the first-order Taylor
expansion of `sqrt(w_phi*w_psi)*exp(im*dtheta)` with the angle difference held
at its reference value.
"""
evaluate_cross_voltage(c::CrossVoltageCoefficients, w_phi::Real, w_psi::Real) =
    c.constant + c.coefficient_phi * w_phi + c.coefficient_psi * w_psi

"""
    winding_voltage_coefficients(d, vbar)

Return affine coefficients for `|d*v|^2`, where `d` is a real incidence row
and `vbar` supplies the nonzero terminal reference phasors.
"""
function winding_voltage_coefficients(d::AbstractVector{<:Real},
                                      vbar::AbstractVector{<:Number})
    length(d) == length(vbar) || throw(DimensionMismatch("d and vbar must have equal length"))
    n = length(d)
    coefficients = zeros(Float64, n)
    constant = 0.0
    for phi in 1:n, psi in 1:n
        scale = Float64(d[phi]) * Float64(d[psi])
        iszero(scale) && continue
        c = cross_voltage_coefficients(vbar[phi], vbar[psi])
        constant += scale * real(c.constant)
        coefficients[phi] += scale * real(c.coefficient_phi)
        coefficients[psi] += scale * real(c.coefficient_psi)
    end
    AffineScalarCoefficients(constant, coefficients)
end

"""
    evaluate_affine(c, w)

Evaluate the real affine scalar `c.constant + c.coefficients' * w`, the form
produced by [`winding_voltage_coefficients`](@ref) for a squared winding
voltage. `w` is indexed by the terminals the incidence row spans.
"""
evaluate_affine(c::AffineScalarCoefficients, w::AbstractVector{<:Real}) =
    c.constant + dot(c.coefficients, w)

const _L3F_OPEN_DELTA_PAIRS = Dict(
    "ABBC" => ((1, 2), (2, 3)),
    "BCAC" => ((2, 3), (1, 3)),
    "CABA" => ((3, 1), (2, 1)),
)

"""
    regulator_gain_matrix(configuration, tap_ratio;
                          connection="ABBC", regulator_type="B")

Return the real gain matrix ``A`` in ``v_from = A*v_to`` for a fixed ideal
WYE, closed-delta, or open-delta step-voltage-regulator bank. The matrices
follow Bazrafshan, Gatsis, and Zhu (PSCC 2018, Table I). `tap_ratio` contains
three ratios for WYE/closed delta and two for open delta. Following the BMOPF
regulator contract, ANSI type A uses the declared ratio directly and type B
uses its reciprocal as the effective ratio.
"""
function regulator_gain_matrix(configuration, tap_ratio;
                               connection="ABBC", regulator_type="B")
    cfg = uppercase(replace(String(configuration), '-' => '_'))
    ratio = tap_ratio isa AbstractVector ? Float64.(tap_ratio) : [Float64(tap_ratio)]
    rt = uppercase(String(regulator_type))
    rt in ("A", "B") || throw(ArgumentError("regulator_type must be A or B"))
    all(x -> isfinite(x) && x > 0.0, ratio) ||
        throw(ArgumentError("tap ratios must be positive and finite"))
    effective = rt == "A" ? ratio : inv.(ratio)

    if cfg == "WYE"
        length(effective) == 3 || throw(DimensionMismatch(
            "WYE regulator bank requires three tap ratios"))
        return Matrix(Diagonal(effective))
    elseif cfg in ("CLOSED_DELTA", "CLOSEDDELTA")
        length(effective) == 3 || throw(DimensionMismatch(
            "closed-delta regulator bank requires three tap ratios"))
        rab, rbc, rca = effective
        return [rab 1-rab 0.0; 0.0 rbc 1-rbc; 1-rca 0.0 rca]
    elseif cfg in ("OPEN_DELTA", "OPENDELTA")
        length(effective) == 2 || throw(DimensionMismatch(
            "open-delta regulator bank requires two tap ratios"))
        pairs = get(_L3F_OPEN_DELTA_PAIRS, uppercase(String(connection)), nothing)
        pairs === nothing && throw(ArgumentError(
            "open-delta connection must be ABBC, BCAC, or CABA"))
        shared = only(intersect(collect(pairs[1]), collect(pairs[2])))
        A = Matrix{Float64}(I, 3, 3)
        for (r, pair) in zip(effective, pairs)
            other = pair[1] == shared ? pair[2] : pair[1]
            A[other, :] .= 0.0
            A[other, other] = r
            A[other, shared] = 1-r
        end
        return A
    end
    throw(ArgumentError("unsupported regulator configuration '$configuration'"))
end

"""
    connection_power_map(D, vbar)

Construct `H = diag(vbar) * D' * diag(D*vbar)^(-1)`, mapping absorbed
physical-channel powers to terminal powers. A zero reference winding voltage
is rejected.
"""
function connection_power_map(D::AbstractMatrix{<:Real},
                              vbar::AbstractVector{<:Number})
    size(D, 2) == length(vbar) ||
        throw(DimensionMismatch("D columns must match vbar length"))
    vb = ComplexF64.(vbar)
    ubar = ComplexF64.(D * vb)
    any(iszero, ubar) && throw(ArgumentError("connection has a zero reference winding voltage"))
    H = Diagonal(vb) * transpose(Float64.(D)) * Diagonal(inv.(ubar))
    matrix = Matrix{ComplexF64}(H)
    ConnectionPowerMap(matrix, real.(matrix), imag.(matrix), ubar)
end

"""
    line_drop_coefficients(Z, vbar_from)

Construct the direct complex-number LinDist3Flow coefficients
`M=2real(conj(Z).*gamma)` and `N=-2imag(conj(Z).*gamma)`, with
`gamma[phi,psi]=vbar[phi]/vbar[psi]`.
"""
function line_drop_coefficients(Z::AbstractMatrix{<:Number},
                                vbar_from::AbstractVector{<:Number})
    n = length(vbar_from)
    size(Z) == (n, n) || throw(DimensionMismatch("Z must be square and match vbar_from"))
    vb = ComplexF64.(vbar_from)
    any(iszero, vb) && throw(ArgumentError("line reference phasors must be nonzero"))
    M = zeros(Float64, n, n)
    N = zeros(Float64, n, n)
    for phi in 1:n, psi in 1:n
        value = conj(ComplexF64(Z[phi, psi])) * vb[phi] / vb[psi]
        M[phi, psi] = 2real(value)
        N[phi, psi] = -2imag(value)
    end
    LineDropCoefficients(M, N)
end
