# Analytic reference cases for DOE development

<!-- doe-executable -->

All Julia blocks run in order from the repository root and are executed by CI.
Run `julia --project=. scripts/run_doe_tutorials.jl` to reproduce them.

These references test physical quantities independently of the nonlinear solver.
They distinguish a mathematical statement about one simple feeder from the
generic finite-point evidence returned by the DOE API.

## A resistive phase with a known solution

Assume a fixed 230 V source, a 0.4 Ω phase conductor, no mutual impedance,
zero reactive injection, and neutrals grounded at both ends. There are no loads
or controllable support devices. Positive power denotes export. The receiving
voltage stays aligned with its source phase on the high-voltage branch:

```math
P=VI=\frac{V(V-V_s)}{R},\qquad
V(P)=\frac{V_s+\sqrt{V_s^2+4RP}}{2}.
```

The 200 V minimum excludes the other algebraic branch. For import, substitute
negative `P`; a real solution requires `P ≥ -Vs²/(4R)`. This is an exact scalar
AC solution under these assumptions, not a linear voltage approximation.

```julia
using PowerOptLab
include("scripts/cases/doe_range_benchmark.jl")
include("scripts/cases/doe_analytic_reference.jl")
using .DOEAnalyticReference

@assert resistive_voltage(19500) == 260
@assert resistive_voltage(-1702.5) == 227
case = doe_benchmark_case()
bound = solve_operating_envelope(case.nets, case.connection_points)
@assert isapprox(only(bound.total_capacity), 58500; atol=0.1)
```

## Prove the dropout violation without an infeasible OPF

Let `a = exp(j2π/3)`. For phase-to-neutral phasors the negative-sequence voltage
is `V₂ = (Vₐ + a²Vᵦ + aV𝚌)/3`. With magnitudes `(230,260,260)` at the source
phase angles, its magnitude is 10 V. The declared limit is 1 V.

```julia
dropout_voltage = phase_voltages([0, 19500, 19500])
@assert isapprox(negative_sequence(dropout_voltage), 10; atol=1e-12)
```

This is a physical violating witness with a 9 V excess. An Ipopt
`LOCALLY_INFEASIBLE` status alone would only be a candidate requiring further
investigation; `ITERATION_LIMIT` would be unresolved numerical evidence.

## A case-specific continuous-box certificate

For equal independent export limits, every phase magnitude lies in
`[230, 230+ΔV]`. Negative-sequence magnitude is a convex function of the three
magnitudes: a norm of a linear map. Its maximum over this magnitude box is
attained at a vertex and equals `ΔV/3`. Monotonicity of `V(P)` maps the entire
power box into that magnitude box. Thus `ΔV ≤ 3 V` is both necessary and
sufficient for the 1 V negative-sequence limit on this fixture.

The exact equal export limit is `233×3/0.4 = 1747.5 W` per phase, or 5242.5 W
in total. The corresponding import reference is `227×3/0.4 = 1702.5 W` per
phase. Both satisfy the phase-voltage limits.

```julia
corners = solve_operating_envelope(case.nets, case.connection_points;
    security=:corners, control_policy=IssuePlusLocalLaws())
@assert isapprox(only(corners.total_capacity), 3equal_export_limit(); atol=0.15)
@assert !corners.diagnostics[1]["global_certificate"]

points = [[0.0, 1.0, 1.0], [0.25, 0.6, 0.9], [1.0, 1.0, 1.0]]
check = verify_operating_envelope(case.nets, case.connection_points, corners;
    utilizations=points)
@assert all(check.feasible)
for (point, context) in zip(points, check.context_results[1])
    expected = phase_voltages(point .* equal_export_limit())
    actual = [complex(context.snapshot["bus"]["b1"][string(i)]["vr"],
                      context.snapshot["bus"]["b1"][string(i)]["vi"]) for i in 1:3]
    @assert maximum(abs.(actual - expected)) < 1e-3
end
```

The certificate belongs to this derivation and its assumptions. It does not
change the generic solver's `global_certificate=false`. Floating neutrals,
mutual impedance, nonzero reactance, voltage-dependent demand, and nonlinear
inverter laws invalidate this particular proof and require other references.

## An affine reference for general box geometry

For `Ap ≤ b` and `ℓ ≤ p ≤ u`, exact containment is equivalent to
`A⁺u + A⁻ℓ ≤ b`, where `A⁺ = max(A,0)` and `A⁻ = min(A,0)` elementwise.
Each row maximizes independently over the box. This is useful for testing sign
conventions, import/export geometry, and algorithms using linear sensitivities.
It certifies the stated affine model only; AC approximation error remains a
separate question.

```julia
A = [1.0 2.0; -1.0 2.0; 1.0 -1.0]
b = [10.0, 5.0, 4.0]
lower, upper = [-1.0, -2.0], [2.0, 1.5]
headroom = polyhedral_box_headroom(A, b, lower, upper)
@assert headroom ≈ [5, 1, 0]
@assert any(polyhedral_box_headroom(A, b, lower, 2upper) .< 0)
```

As a research exercise, add neutral impedance without grounding the receiving
neutral, then quantify where the independent-phase reference ceases to agree.
Do not retain its containment claim after changing its physical assumptions.
