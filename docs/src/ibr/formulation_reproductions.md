# Reproducing forward operation with optimisation

The [formulation comparison table](../components/inverter_operating_constraints.md)
separates hardware feasibility from controller decisions. These cases test when
optimisation can reproduce the forward algebraic model. Run them with
`julia --project=. scripts/formulations/ibr_reproductions.jl`.

## Independent feeder equilibria

Three balanced feeders have source voltages of 215, 230 and 259 V and 0.05 Ω
phase resistance. The inverter has nonzero filter impedance and fixed losses.
Both formulations use the same smooth average-voltage curves, source budget,
sensing reference and balanced converter-current allocation. Hardware limits
are inactive; the 259 V case has an active volt-watt ceiling. Maximising export
then reproduces the forward dispatch to the declared approximation tolerance.

```@example ibr_reproductions
using PowerOptLab
include(joinpath(pkgdir(PowerOptLab), "scripts", "formulations", "ibr_reproductions.jl"))
rows = ibr_formulation_reproductions()
for row in rows.feeder
    @assert row.forward.solve.publishable && row.constrained.solve.publishable
    @assert row.p_error_W < .2 && row.q_error_var < .2
    @assert row.current_error_A < .001 && row.voltage_error_V < .001
end
[(source_V=r.voltage, delta_P_W=r.p_error_W, delta_Q_var=r.q_error_var,
  delta_I_A=r.current_error_A, delta_V_V=r.voltage_error_V) for r in rows.feeder]
```

These are different solved feeder equilibria, not a replay at one model's
voltage. Small differences remain because the forward path regularises the
power-to-current conversion (here `power_voltage_floor=1e-3` V) and uses smooth
conservative selectors/backoff. The optimisation retains the curve equalities
and inequalities directly. This test is conditional, not a theorem for
arbitrary unbalanced networks, current targets or binding hardware constraints.

## Saturated proportional current limiting

A lossless 6 kVA three-leg plant has ``I_{max}=1.1(6000)/(3\cdot230)`` A.
At fixed 230 and 207 V POC voltage, its raw request is 6000 W + j2640 var.
The optimisation **adds** the controller's ray restriction
``I_{phase}=t I_{phase,raw}``, ``0\le t\le1``, then maximises ``t`` under the
physical capability inequalities. This additional relationship reproduces the
proportional limiter; it is not an implication of a current inequality.

An independent lossless oracle is
``t=\min(1,S_{max}/|S_{raw}|,3VI_{max}/|S_{raw}|)``. The optimiser agrees with
this scale to ``10^{-7}``; the forward model agrees within 0.4 W and 0.2 var,
including its conservative smoothing/backoff.

```@example ibr_reproductions
for row in rows.proportional
    @assert row.forward.solve.publishable && row.optimisation.status.publishable
    @assert abs(row.optimisation.scale-row.analytic_scale) < 1e-7
    @assert row.p_error_W < .4 && row.q_error_var < .2
end
[(POC_V=r.voltage, scale=r.optimisation.scale, analytic=r.analytic_scale,
  delta_P_W=r.p_error_W, delta_Q_var=r.q_error_var) for r in rows.proportional]
```

By contrast, `VoltVarWattConstraints` retains Q=2640 var and curtails more P
under current saturation. The [207 V comparison](@ref inverter-saturation-comparison)
and its regression intentionally preserve that difference. Maximising P alone
does not reproduce a proportional limiter.

## Sequence-priority allocation at a fixed measurement

The final cases have ``U_1=230`` V, ``U_2=10e^{j0.2}`` V, no zero sequence,
and a 6 A conductor limit. A purely reactive negative-sequence admittance
requests approximately 3 A, alongside a 3000 W positive-sequence request.
For both `priority=:power` and `:balance`, the reserved sequence alone fits
inside the conductor rating. Optimisation holds that sequence fixed and
maximises the other sequence's scale. It agrees with both `evaluate_exact`
and an independent phase-by-phase quadratic disk intersection.

```@example ibr_reproductions
for row in rows.sequence_priority
    @assert row.optimisation.status.publishable
    @assert abs(row.optimisation.scale-row.analytic_scale) < 1e-7
    @assert row.current_error_A < 1e-6
end
[(priority=r.priority, scale=r.optimisation.scale,
  delta_I_A=r.current_error_A) for r in rows.sequence_priority]
```

This is a fixed-measurement allocation reproduction. It does not establish
equivalence of independently solved sequence-control network equilibria, the
reserved-sequence backoff when that sequence alone is infeasible, or carrier
reserve closure. The test suite separately checks full feeder reproduction,
physical margins, failed-solve masking and the intentional strict-Q mismatch.
