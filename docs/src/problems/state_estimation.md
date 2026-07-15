# State estimation

> **Kind:** Problem specification · **Maturity:** promotion candidate · **Direction:** inverse · **Temporal:** single-snapshot

[`solve_state_estimation`](@ref) is a *different problem specification* over the
same network physics: given noisy measurements of an energised network, find the
bus voltage state that best fits them in a weighted-least-squares (WLS) sense. It
is the maximum-likelihood estimate under independent, zero-mean, Gaussian
measurement errors with a diagonal covariance (weight ``1/\sigma_i^2``).

It reuses the BMOPFTools device model but with three changes, all expressed
through the public seams:

1. **No operational bounds.** The estimation network is a *physics-only* model —
   buses, lines, transformers, shunts, and a voltage source. It carries no loads,
   generators, IBRs, or operational limits.
2. **Free injections instead of fixed loads.** A `model_hook!` adds a free
   injection current at each measured bus (so KCL closes with the voltages free to
   fit the data).
3. **A residual objective instead of generation cost.** The hook sets the WLS
   objective ``\sum_i (z_i - h_i(\text{state}))^2 / \sigma_i^2``.

## The estimation network is a contract

An operational net biases the estimate: a retained load adds its fixed injection
*on top of* the estimated one, and a retained voltage/thermal limit constrains
the state to a feasible region that has nothing to do with the measurements. On
the worked example below, retaining the loads moves the estimate by ~20 V and
inflates the objective from ≈0 to ≈355; a 950 V bus limit pins the estimate to
950 V — both while still reporting `LOCALLY_SOLVED`.

`solve_state_estimation` therefore **rejects** a net that carries injecting
devices (`load`, `generator`, `ibr`, `dc_source`) or operational limits
(`v_min`/`v_max`/…, line `i_max`/`s_max`/`s_rating`). Pass
`allow_operational=true` to override with a warning when you deliberately want to
estimate against such a model.

## Injection coverage is explicit

Absence of telemetry is **not** evidence of zero injection. Every non-source
phase terminal must be *either* a measured injection *or* explicitly declared
zero-injection:

- A measured injection needs **both** `:pinj` and `:qinj` (an injection is a
  complex quantity; a lone `P` or `Q` is ill-posed). That pair gives the bus a
  free injection current.
- Everything else must be listed in `zero_injection` (a bus id, expanded over its
  phase terminals, or a `(bus, terminal)` pair) — the classic zero-injection
  pseudo-measurement.

An un-declared, un-measured bus is an **error**, not a silent zero injection.

## Measurements and reference convention

Each [`Measurement`](@ref) is an SI scalar with a standard deviation `sigma`
(WLS weight ``1/\sigma^2``), validated at construction (finite value, `sigma>0`,
supported `kind`):

- `:vmag` — voltage magnitude across `(bus, terminal)` → `reference`, in volts.
- `:pinj` — active power injected into the network at that terminal pair, watts.
- `:qinj` — reactive power injection, vars.

All three share **one** reference terminal per measurement (`reference`,
defaulting to the solve's `neutral`): a smart-meter reading is phase-to-neutral
for voltage *and* power. Pass `reference=nothing` to reference a terminal to
ground, or an explicit terminal name for a bespoke return path. This resolves the
earlier inconsistency where `:vmag` was terminal-to-ground while `:pinj`/`:qinj`
were phase-to-neutral — masked only because the test feeders perfectly ground
every neutral.

## Observability — convergence is not uniqueness

The WLS fit is **nonconvex**: the squared `P`/`Q` residuals are quartic and the
voltage-magnitude equality is nonconvex, so Ipopt returns a *local* stationary
point (and low-voltage solutions exist). A converged solve does **not** prove the
state is unique.

The result carries an `observability` diagnostic: the Jacobian of the
measurement + zero-injection equations with respect to the rectangular node
voltages is formed at the returned point (reusing BMOPFTools' `ybus_passive` for
`I = Y·V`), and `observable = rank == n_states`. It reports `redundancy` (surplus
equations), the smallest singular value, and the condition number, and
`solve_state_estimation` warns on rank deficiency. This is a **local numerical**
identifiability check (grounded-neutral networks; it does not model floating
neutral displacement), not a global uniqueness proof.

!!! warning "Local optimum"
    Like the IVQ battery solve, this is a
    nonconvex problem solved to a local stationary point. For a promotion-grade
    estimator, add multistart or a Gauss–Newton/normal-equations method, and
    check `observability.observable` and `primal_status` before trusting a result.

## Residuals

`residuals[i].standardized = residual/σ` is the **σ-normalised raw residual** — a
scale-free residual, *not* the classical leverage-adjusted normalised residual
``r^N_i = r_i/\sqrt{S_{ii}}`` (with ``S`` the residual covariance) used for
bad-data identification. A χ² goodness-of-fit test and largest-normalised-residual
bad-data processing are not yet implemented.

## Worked example

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf

# A physics-only net: buses, lines, source; no loads or operational limits.
net = parse_bmopf("""
{"bus":{
    "src": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus1":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus2":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],
     "v_magnitude":[1000.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.5}},
 "line":{
    "l1":{"bus_from":"src","bus_to":"bus1","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
    "l2":{"bus_from":"bus1","bus_to":"bus2","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
"""; from_string=true)

# Both load buses are measured (a :pinj+:qinj pair each), so both get a free
# injection; no bus is left to a silent zero-injection assumption.
meas = [
    Measurement(kind=:vmag, bus="bus1", value=979.5,     sigma=2.0),
    Measurement(kind=:pinj, bus="bus1", value=-20_000.0, sigma=400.0),
    Measurement(kind=:qinj, bus="bus1", value=0.0,       sigma=400.0),
    Measurement(kind=:vmag, bus="bus2", value=969.2,     sigma=2.0),
    Measurement(kind=:pinj, bus="bus2", value=-20_000.0, sigma=400.0),
    Measurement(kind=:qinj, bus="bus2", value=0.0,       sigma=400.0),
]

se = solve_state_estimation(net, meas)

se.primal_status            # "FEASIBLE_POINT" — trust the estimate only then
se.bus["bus1"]["1"]["vm"]   # estimated |V| at bus1 (V); NaN if not FEASIBLE_POINT
se.residuals                # per-measurement measured/estimated/residual/standardized
se.observability            # (observable, n_states, rank, redundancy, min_singular, cond)
se.objective                # optimal weighted-residual sum
```

Here the six measurements exceed the four voltage components, and — because they
are placed so the state is observable (`se.observability.observable == true`,
`redundancy == 2`) — the fused estimate filters measurement noise. Redundancy
alone is not observability: *placement* and Jacobian rank decide whether the
state is determined.

## Not yet supported

Branch-flow / branch-current / phasor (PMU) / source measurements; bad-data
detection and identification (χ², largest-`rᴺ`); floating/displaced-neutral
observability; multistart / global search.

## Literature

- Schweppe & Wildes, *Power system static-state estimation, Part I* (1970) — the
  foundational exact WLS model.
- Baran & Kelley, *State estimation for real-time monitoring of distribution
  systems* (1994) — node-voltage DSSE with power/voltage/current measurements.
- Monticelli & Garcia, *Reliable bad data processing for real-time state
  estimation* (1983) — classical normalised-residual bad-data processing.
- Dehghanpour et al., *A survey on state estimation techniques and challenges in
  smart distribution systems* (2019) — pseudo-measurements, observability,
  topology, and meter placement in modern DSSE.

See the API reference for [`Measurement`](@ref), [`solve_state_estimation`](@ref),
and [`StateEstimationResult`](@ref).
