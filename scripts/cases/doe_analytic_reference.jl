module DOEAnalyticReference

export resistive_voltage, negative_sequence, phase_voltages,
       equal_export_limit, equal_import_limit, polyhedral_box_headroom

"""
High-voltage branch of P = V(V - Vs)/R, with positive P denoting export.
Assumes a resistive uncoupled phase, zero Q, fixed source angle and grounded
neutral. This is an independent scalar power-flow reference, not a DOE solver.
"""
function resistive_voltage(p::Real; source=230.0, resistance=0.4)
    source > 0 && resistance > 0 || throw(ArgumentError("positive source and resistance required"))
    discriminant = source^2 + 4resistance * p
    discriminant >= 0 || throw(DomainError(p, "no real power-flow solution"))
    return (source + sqrt(discriminant)) / 2
end

phase_voltages(p; kwargs...) = resistive_voltage.(p; kwargs...) .* cis.([0, -2pi/3, 2pi/3])
negative_sequence(v) = abs((v[1] + cis(-2pi/3)*v[2] + cis(2pi/3)*v[3]) / 3)

# On this balanced, independent-phase fixture the voltage-magnitude box has
# maximum negative sequence ΔV/3, attained at any nontrivial binary corner.
# Convexity of the norm of a linear phasor map bounds every interior point.
function equal_export_limit(; source=230.0, resistance=0.4, vmax=260.0, vneg=1.0)
    v = min(vmax, source + 3vneg)
    return v * (v - source) / resistance
end

function equal_import_limit(; source=230.0, resistance=0.4, vmin=200.0, vneg=1.0)
    v = max(vmin, source - 3vneg)
    v >= source/2 || throw(ArgumentError("reference requires the high-voltage branch"))
    return v * (source - v) / resistance
end

"""Exact residual b - max(Ap : lower ≤ p ≤ upper), one entry per row."""
function polyhedral_box_headroom(A, b, lower, upper)
    all(lower .<= upper) || throw(ArgumentError("box bounds reversed"))
    return b - max.(A, 0) * upper - min.(A, 0) * lower
end

end
