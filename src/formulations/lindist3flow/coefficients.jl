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

evaluate_affine(c::AffineScalarCoefficients, w::AbstractVector{<:Real}) =
    c.constant + dot(c.coefficients, w)

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

"""
    regular_polygon_coefficients(K; inner=true)

Return equally spaced unit normals and the radius multiplier for a regular
polygon approximation of a circle. `inner=true` uses `cos(pi/K)`, so every
polygon point lies inside the requested circle.
"""
function regular_polygon_coefficients(K::Integer; inner::Bool=true)
    K >= 3 || throw(ArgumentError("K must be at least 3"))
    theta = 2pi .* (0:Int(K)-1) ./ Int(K)
    normals = hcat(cos.(theta), sin.(theta))
    RegularPolygonCoefficients(normals, inner ? cos(pi / Int(K)) : 1.0)
end
