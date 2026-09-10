# Equations, physical boundaries and numerical formulations

## Energy is an account with a defined boundary

Let ``p^c_t,p^d_t\ge0`` be **AC-side** charging and discharging power in W.
Positive network injection means export:

```math
P_t=p^d_t-p^c_t,\qquad
E_{t+1}=E_t+\Delta t_t(\eta_c p^c_t-p^d_t/\eta_d).
```

``E`` is Wh and ``\Delta t`` is hours. The fixed efficiencies encompass the
chosen AC-to-stored-energy boundary, including whatever conversion and battery
losses were included in their calibration. Adding a converter loss model on top
without changing that boundary double-counts losses. The PE account is not
coulomb counting and does not certify electrochemical limits.

The initial energy is fixed. All subsequent states respect the energy bounds.
At departure interval ``d``, ``E_{d+1}\ge E^{req}``. Optional `energy_final`
fixes ``E_{T+1}``; for `StorageDevice` it overrides the cyclic convention.

With fixed powers and efficiencies, energy is affine within each interval.
Endpoint bounds therefore bound this PE energy throughout the interval. This
statement does not transfer automatically to temperature, voltage, or a
state-dependent nonlinear loss model.

## Electrical port

For each phase-to-neutral voltage ``V_{\phi n}`` and injected current ``I_\phi``,

```math
S_\phi=V_{\phi n}I_\phi^*,\quad
P=\sum_\phi\Re S_\phi,\quad Q=\sum_\phi\Im S_\phi.
```

The current is stamped into phase KCL and its negative into the neutral. Optional
limits are imposed as smooth polynomial inequalities:

```math
P^2+Q^2\le\bar S^2,\qquad
(\Re I_\phi)^2+(\Im I_\phi)^2\le\bar I^2,\qquad
\left|\sum_\phi I_\phi\right|^2\le\bar I_n^2.
```

The implementation normalizes these three squared-norm inequalities by their
own equipment ratings, giving a right-hand side of one. Otherwise a small
per-unit rating squared can be comparable to the solver's absolute feasibility
tolerance and permit several watts of excess on a kW-scale charger.

The apparent-power circle alone cannot constrain phase currents whose complex
powers cancel. The default equal-power policy imposes ``P_\phi=P_1`` and
``Q_\phi=Q_1``. `:independent` deliberately models independent phase dispatch;
it is appropriate only for equipment supporting that control freedom.

When disconnected, **each real and imaginary phase current is fixed to zero**.
Zero aggregate P and Q alone do not express an open connection. An EVSE standby
load or empty-station reactive service would be another physical object.

This is a fundamental-frequency steady-state port. It does not resolve current
harmonics, switching transients, contactors, or grid-forming behavior.

## Three operating-mode formulations

Let ``a=p^c/\bar p^c`` and ``b=p^d/\bar p^d`` for positive ratings. A zero rating
fixes its power to zero and removes the need for a complementarity pair.

| `operation` | Rows | Meaning |
|:--|:--|:--|
| `:independent` | ``a,b\ge0`` and power bounds | Outer relaxation; artificial simultaneous operation may be optimal |
| `:relaxed` (default) | additionally ``ab\le\tau`` | Smooth NLP with finite mode leakage; ``\tau=10^{-8}`` by default |
| `:complementarity` | ``0\le a\perp b\ge0`` | Exact mathematical mode exclusion, solved to numerical MPCC tolerances |

The implementation creates explicit normalized variables for the exact pair,
linked to the network power variables. This is required by the tested
MathOptComplements bridge and keeps the mode scale independent of the network
VA base. For the relaxed path, the normalized product is imposed directly.

Because ``\min(a,b)\le\sqrt{\tau}``, a small product can still permit appreciable
power when ratings are large. Numerical feasibility residuals add to the chosen
``\tau``. At unit efficiency, simultaneous operation can be degenerate; at
negative prices, inefficient simultaneous operation can be strictly attractive.
[The mode tutorial](modes.md) derives a counterexample.

## CCOpt configuration

The optional environment preserves the library's immutable BMOPFTools revision:

```sh
julia --project=. scripts/instantiate_pinned.jl
julia scripts/formulations/setup.jl /tmp/pol-ev-optional
julia --project=/tmp/pol-ev-optional scripts/ev/ccopt.jl
```

The destination for `setup.jl` must be new or empty. CCOpt remains optional;
production dependencies are not replaced by the optional environment.

```julia
using PowerOptLab, JuMP, CCOpt, MathOptComplements, NLPModelsJuMP
result = solve_multiperiod_opf(nets, [ev];
    optimizer=CCOpt.Optimizer,
    configure! = MathOptComplements.Bridges.add_all_bridges,
    solver_options=(tol=1e-8,
        relaxation_update=CCOpt.ProportionalRelaxationUpdate(sigma_min=1e-12),
        max_iter=1000,))
```

Here `ev.operation == :complementarity`. `configure!` runs before optimizer
attachment. The package calls `optimize!` once; CCOpt owns its relaxation,
penalty/homotopy and active-set algorithms. See the
[CCOpt implementation](https://github.com/madsuite-org/CCOpt.jl) and
[2026 preprint](https://arxiv.org/abs/2604.18726).

The tested CCOpt configuration supplies both an NLP stopping tolerance and a
smaller floor on the chosen relaxation update (`relaxation_update.sigma_min`). Tightening only the NLP tolerance can
leave the solver at an acceptable-tolerance status while its relaxed MPCC
subproblem uses a larger floor. The optional regression reports physical
residuals in addition to status; these settings are not universal tuning rules.
The [zero-dispatch investigation](ccopt_investigation.md) also varies the
backend's bound relaxation. Its effect must not be confused with the
continuation floor or stopping tolerance.

An exact product equality passed to an ordinary NLP solver has MPCC degeneracy;
calling it smooth does not restore standard constraint qualifications. Distinct
MPCC stationarity notions matter; a successful MOI termination is not a
stationarity certificate. See
[Nurkanović, Pozharskiy and Diehl](https://arxiv.org/abs/2312.11022).

## Optional charging acceptance

`EV.charge_acceptance` is a nonnegative, nonincreasing `PWLFunction` mapping
PE energy fraction ``z=E/E^{max}`` to an AC power ceiling in W. Its declared
units must be `:unitless` and `:W`. It must cover the usable energy domain.
It is an envelope, not a forced charging trajectory or pilot response law.

For each connected interval the compiler enforces

```math
p^c_t\le f(z_t),\qquad p^c_t\le f(z_{t+1}).
```

Because the PE state is affine within an interval and the canonical cap is
nonincreasing, both endpoints bound its minimum along that path. This is a
conservative **piecewise-constant dispatch** restriction. It does not claim to
integrate a continually changing maximum charging rate exactly.

For concave PWL ``f=\min_k(a_kz+b_k)``, `:auto` lowers these upper bounds to
ordinary affine inequalities. No minimum operator, smoothing, or integer
variable is required. Nonconcave envelopes require an explicit smooth
formulation or `ComplementarityGraph`; integer and hull encodings are rejected
by the EV API. Smooth bounds use conservative error correction by default.
See [the approximation tutorial](approximations.md).

## Results and independent checks

PE dispatch includes `energy_wh`, the compatible `soc`, AC powers, per-phase
P/Q, real/imaginary/magnitude phase currents, neutral magnitude, mode policy,
and normalized relaxation tolerance. Sessions also expose vehicle and outlet
ids and occupancy. Snapshot `custom_injection` includes these device injections.

`diagnostics` evaluates energy balance (Wh), energy and power bounds,
departure shortfall (Wh), apparent-power excess (VA), current excess (A),
disconnected current (A), phase-sharing errors (W/var), simultaneous power (W),
and normalized mode product/minimum. Sessions with an acceptance curve also
report violation of the **original** endpoint envelope in W.

A missing optional rating produces `nothing` for its violation diagnostic.
Failed solves return NaN trajectories; an inapplicable diagnostic may still be
zero or `nothing`. Check solve status before aggregating diagnostics. These
checks are independent arithmetic on the candidate, not validation against a
measured vehicle or an AC time-domain replay.
