# Reproducing forward operation with optimisation

Start with the [P–Q tutorial](../tutorials/inverter_constraints_and_control.md)
for the meaning of a feasible region, an operating relationship and an objective.

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

## Unbalanced feeders and LCL current targets

The extended cases solve both models on a feeder with ``U_1=230`` V and
``U_2=8e^{j0.2}`` V at the source and 0.05 Ω phase resistance. Their POC voltages
are solved separately. The LCL cases use nonzero series impedances, a 30 µF
midpoint capacitor and damping resistance. Both models have the same physical
plant, a 3000 W DC source budget and a constant 600 var command. The reduced
filter comparison removes the midpoint branch.

| Case | Expected comparison | What must match |
|:--|:--|:--|
| Unbalanced reduced filter | Agreement | Balanced converter current and total POC Q |
| LCL, grid-current target | Agreement | Balanced grid current and total POC Q |
| LCL, converter-current target | Agreement after matching the Q location | Balanced converter current and Q formed from POC voltage and converter current |
| LCL, converter target versus strict POC Q | Deliberate difference | The two models impose different Q relationships; current-allocation choice alone does not align them |

The forward command is formed using its sensing voltage and its selected target
current:

```math
Q_{command}=\Im\sum_k U_{POC,k}I_{target,k}^*.
```

Total POC exchange instead uses grid-arm current. In these examples there is
no separate POC shunt, so
``Q_{POC}=\Im\sum_k U_{POC,k}I_{grid,k}^*``. The LCL midpoint capacitor makes
``I_{converter}\ne I_{grid}``. Consequently, a 600 var converter-target command
can coexist with approximately 2103 var of POC injection. This is not a
reactive-tracking error at the declared current target.

For the converter-target reproduction, the example stamps the general free-Q
policy, adds the command's Q relationship at the matching location and maximises
POC export. The converter's *internal-terminal* Q is a third quantity because
its voltage includes filter drops. The tests independently verify all three
power definitions, feeder KVL and the capacitor relation
``I_{converter}-I_{grid}=Y_{mid}U_{mid}``.

```@example ibr_reproductions
network = ibr_network_reproductions()
for row in network.pairs
    @assert row.forward.solve.publishable && row.constrained.solve.publishable
    if row.expected == :agreement
        @assert row.p_error_W < .2 && row.q_error_var < .2
        @assert row.current_error_A < .001 && row.voltage_error_V < .001
    else
        @assert row.q_error_var > 1000.
    end
end
[(case=r.id, expected=r.expected, delta_P_W=r.p_error_W,
  delta_Q_var=r.q_error_var) for r in network.pairs]
```

The matched cases agree within 0.2 W/var, 1 mA and 1 mV. The unmatched Q
locations differ by about 1.5 kvar despite both solves being feasible. Here
maximum export also reproduces the source-limited active dispatch; that
selection is conditional on the same source budget and inactive other limits.

## Saturation and infeasibility are different outcomes

The remaining cases deliberately exercise three different situations:

| Case | Forward model | Optimisation model | Independent explanation |
|:--|:--|:--|:--|
| Unbalanced 207 V feeder with a binding conductor limit | Var-priority apparent-power allocation followed by current backoff reduces P and Q | Strict POC Q=2640 var, maximum P retains Q and curtails more P | Different operating relationships on the same physical plant |
| Requested Q=7200 var on a lossless 6000 VA plant | Limiter backs off the request | Strict Q service is infeasible | ``|S|\ge|Q|=7200>6000`` VA |
| Grounded LCL shunt, three-leg bridge, nonzero POC ``U_0`` and balanced grid current | A corresponding exact zero-sum grid-current command has the same circuit obstruction | The tested balanced-grid constraints are infeasible | Capacitor zero-sequence current has no path consistent with both imposed current sums |

```@example ibr_reproductions
saturated = network.saturated
@assert saturated.forward.solve.publishable && saturated.constrained.solve.publishable
@assert saturated.constrained.plant.q_poc-saturated.forward.plant.q_poc > 15.
@assert saturated.forward.plant.p_poc-saturated.constrained.plant.p_poc > 5.
@assert !network.impossible_service.constrained.solve.publishable
@assert !network.zero_sequence.balanced.solve.publishable
@assert network.zero_sequence.released.solve.publishable
(forward=(P_W=saturated.forward.plant.p_poc,Q_var=saturated.forward.plant.q_poc),
 strict_q=(P_W=saturated.constrained.plant.p_poc,Q_var=saturated.constrained.plant.q_poc),
 impossible_Q=network.impossible_service.constrained.termination_status,
 incompatible_balance=network.zero_sequence.balanced.termination_status)
```

The last case fixes POC phasors to isolate the circuit obstruction, adding a
3 V zero-sequence voltage to the unbalanced source. Balanced grid current has
``I_{0,grid}=0``. Identical grid-arm impedances then give
``U_{0,mid}=U_{0,POC}=3`` V. The grounded capacitor requires
``I_{0,shunt}=Y_{mid}U_{0,mid}\ne0``, while three-leg topology requires
``I_{0,converter}=0``. Midpoint KCL would require
``0=0+I_{0,shunt}``, a contradiction. This is an analytical obstruction,
not a conclusion inferred only from a local solver status.

Switching the allocation to balanced **converter** current restores a path:
the grid arm carries ``I_{0,grid}=-Y_{mid}U_0/(1+Y_{mid}Z_{grid})``. The paired
solve and test verify that complex current. This changes the operating
assumption while retaining the physical plant. It does not make every
infeasible service request feasible.

These cases complement the fixed-measurement sequence-priority example above.
They do not establish arbitrary unbalanced sequence-controller equilibrium
equivalence, dynamic stability, or exact-versus-smooth network fixed-point
equivalence. The quantities checked here are separately solved smooth-model
operating points and independently evaluated circuit/capability identities.
