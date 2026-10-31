# Four-leg per-phase control counterfactual

This study changes the three-wire bridge into an explicitly connected four-wire
bridge. It asks what independent phase reactive requests do when a physical
neutral return is available. It is a separate research-control scenario, not a
claim that this output policy meets Australian three-phase connection rules.

The phase rating (6 kVA total, 10 A per phase), 750 V / 300 µF DC link, 1 A
capacitor RMS budget, loss coefficients and phase filter are unchanged. The
additional neutral has an explicit 10 A RMS rating and the same R/X as one phase
filter arm. These are declared topology assumptions; no hardware sizing is
optimized. Neutral conduction adds its existing physical loss terms.

Each phase evaluates the Australia A VV and VW knots using its own PN voltage.
Q bases and rated P bases are one third of the total request. The physical
allocator maximizes a common service fraction for the three Q requests, fixes
that fraction, then maximizes total active export with equal POC phase active
powers. Equal phase powers generally require unequal currents and a nonzero
neutral return. One shared DC source constrains the whole inverter. Phase and
neutral RMS bounds, total POC VA, DC ripple/capacitor current and continuous
modulation rails remain active. No per-phase VA rating is inferred.

Eleven predeclared cases cover 0/50/100% PV at 1.00/1.03 pu source voltage,
0/1/2% aligned source VUF, balanced versus 60/25/15% load shares, and two
disconnected cases. Energized zero-PV cases retain reactive service and loss
imports. Disconnection opens the entire AC model; DC standby is outside scope.
Every case uses the common fixed numerical profile without hidden retries.

```sh
julia --project=. scripts/australian_pv/test_four_leg.jl
julia --project=. scripts/australian_pv/four_leg.jl /tmp/au_four_leg
julia --project=. scripts/australian_pv/audit_four_leg.jl /tmp/au_four_leg
julia --project=. scripts/australian_pv/test_four_leg_records.jl /tmp/au_four_leg
python scripts/australian_pv/plot_four_leg.py /tmp/au_four_leg
```

The driver saves all planned-case statuses, iteration diagnostics, physical
checks, device metrics and phasors. A failed solve remains explicit and cannot
issue an accepted command. The separate saved-point audit reconstructs SI
physical and phase-policy checks from those phasors and reruns fixed-PQ network
injections. The latter uses the same network engine; small-feeder package tests
separately check Ohm's law, phase/neutral KCL and full LCL power losses.

These are quasi-static local nonconvex references. They establish neither
global service optimality nor controller dynamics, response time, switching
peaks, thermal lifetime or Australian connection conformity. The 253 V customer
criterion in the coordinated benchmark is not imposed by this local VV+VW
policy. A converged four-leg point can still exceed it.
