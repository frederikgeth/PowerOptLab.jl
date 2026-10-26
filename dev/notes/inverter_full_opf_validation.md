# Full inverter OPF validation plan

Recorded before implementation of the larger validation example.

Existing capability: `ConstrainedInverter` can be stamped into a shared native
network OPF with explicit requests and an external objective. Standalone solving
is one-device convenience; it is not the limit of the physical formulation.
Keep current topology, smooth-curve, fixed-DC and PWM limitations explicit.

1. Build a reproducible seven-bus radial four-wire network with three inverters,
   unequal constant-PQ phase loads, floating downstream neutrals and full mutual
   line impedances. Keep its independent physical data alongside the native net.
2. Compare separately built forward-control and maximum-export OPF models under
   matched POC Q/current relationships, source availability and inactive hardware
   saturation. Include nonlinear voltage feedback, filter/loss physics and
   multiple inverter policies. Demonstrate agreement rather than assume it.
3. Solve a genuinely coordinated OPF: free Q, unequal active-export weights,
   reactive-use penalty, phase-to-neutral voltage bounds and a feeder conductor
   rating. Verify at least one nontrivial network constraint binds and that the
   result differs from the independent local-controller solution.
4. Replay both controlled and optimised P/Q dispatch with an independently coded
   four-wire backward/forward sweep. Check complex voltages and branch currents,
   KCL, KVL, source/load/inverter/loss balance and inverter capability margins.
   Replay validates equilibrium consistency; it does not prove OPF global
   optimality, uniqueness, stability or controller equivalence of the free-Q OPF.
5. Exercise power-base changes, failure masking and an overloaded/infeasible
   service variant where useful. Repair only demonstrated integration gaps.
6. Update full-OPF documentation with executable composition and validation
   evidence, including explicit voltage references and avoiding duplicate native
   inverter ports. Run affected regressions and a complete documentation build.

Integration audit: add a public `solve_status(::JuMP.Model)` adapter for safe
external-model extraction, and make `inverter_handles` accept the handles
returned by stamping a `ControlledDevice`, as it already does for constrained
devices. Preserve existing result/status and raw-plant-handle behaviour.

## Completed and verified

- Added public `solve_status(::JuMP.Model)` and the controlled-wrapper handle
  adapter, with success, infeasible and unstarted status coverage. Updated the
  existing reproduction/generalised-generator examples to use the public status
  method rather than a private adapter.
- Added `scripts/formulations/inverter_network_opf.jl` with a seven-bus,
  four-wire network, full mutual impedances, eighteen unequal phase loads and
  three inverters. Five separately built models cover matched forward controls,
  matching maximum export, unmanaged high availability, coordinated free-Q
  dispatch and a feasible fixed-PQ benchmark. The matching model has 206 variables.
- Matched controls agree to maximum 0.1167 W, 0.0362 var, 0.000177 V,
  0.000166 A converter current and 0.000409 A branch current. Source budgets
  bind in the constrained comparison; hardware saturation is inactive there.
- Coordinated dispatch binds both 244 V phase-to-neutral and 12 A feeder-head
  limits, as well as converter S bounds. Unmanaged controls at the same high
  availability reach about 247.94 V / 13.29 A. The coordinated objective improves
  on a feasible 2000 W / zero-Q-per-inverter benchmark.
- Native `vpn_min`/`vpn_max` and independently stamped voltage inequalities
  agree, as do 10 kVA and 1 MVA coordinate bases. An independent SI radial sweep
  reproduces all five solutions' complex voltages/currents near floating-point
  precision. KCL, KVL and source/load/inverter/line power balance are checked.
- Added 294 checks for public status integration and the full case, including
  physical losses, capability margins, failed-solve masking, replay convergence,
  coordinate invariance and native/hook formulation agreement. The broader
  affected suite passed 3,974 checks with one optional OpenDSS skip.
- Full Documenter build passed. The public-API composition recipe and all
  comparison/replay assertions execute in the new guide, and the rendered
  numerical results and local links were checked. Existing documentation size
  warnings remain.

The readiness table in `docs/src/ibr/full_network_opf.md` separates working
staged composition from missing convenience wrappers and unsupported operating
laws. The independent oracle is explicitly scoped to the radial, no-shunt,
balanced-current fixture; no global-optimality, uniqueness or stability claim
is inferred from these agreement checks.
