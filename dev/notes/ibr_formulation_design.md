# IBR simulation and optimisation formulation design

Initial design, recorded before implementation on `codex/ibr-constraint-control`.

## Purpose and boundaries

Share the existing `AdvancedInverter` physics between forward prescribed-current
control and optimisation of feasible operation. Hardware limits remain hard.
Operating relationships, current allocation, DC source bounds, numerical
encoding and the selection objective are separate modelling choices.

Preserve existing `SequenceController`, limiter and `VoltVarWattConstraints`
defaults. A missing volt-var curve in the latter continues to mean Q = 0.
Introduce a general operating policy whose default explicitly leaves Q free.
Do not silently interpret an optimisation result as an autonomous controller.

## API and representation

- Introduce an abstract operating-policy interface. Retain `VoltVarWattConstraints`
  as the voltage-law implementation, and add `InverterOperatingConstraints` for
  free Q, fixed Q and signed fixed PF, with optional P/Q bounds. It may reuse a
  positive policy for watt ceilings but rejects an ambiguous second Q law.
- Represent DC source limits with `DCSourceBounds`. The legacy default retains
  the request's upper availability bound; `p_min=0` explicitly represents a
  non-charging source. Bounds apply to DC power, so PV can still import AC
  power for inverter losses. Inconsistent requested service is infeasible;
  malformed/nonfinite/reversed configured bounds are rejected.
- Add current allocation to `ConstrainedInverter`: `:balanced_converter`
  (existing default), `:balanced_grid`, or `:free`. Enforce the two sequence
  equalities at the selected location without duplicating plant topology rows.
- Support THREE_LEG, FOUR_LEG and SPLIT_DC with the existing plant's physical
  neutral and capacitor requirements. A single-phase voltage-policy interface,
  GFM operating laws and automatic PWM reserve closure remain outside this slice.
- Preserve direct stamping without replacing an external objective. Retain
  max-export and fixed-dispatch feasibility convenience solves. Report selected
  policy/allocation/source and physical margins; use `nothing` for an absent Q
  target, never a fabricated current-limiter scale.
- Finish opt-in rating normalisation for capability rows used by the supported
  topologies, while retaining the old numerical default for forward controllers.

## Simulation reproduction contract

Reproduction is conditional on matching the mathematical problem:

1. Unsaturated forward curves and maximum-export constraints can agree when
   source availability, sensing, power/current locations and current allocation
   match. Test both fixed POC voltages and independently solved feeder equilibria.
   Account explicitly for the forward model's P-to-I voltage floor and smoothing.
2. A proportional limiter is reproducible by maximising a common scale on a
   prescribed current direction under the same capability inequalities. This
   ray restriction is part of the reproduced controller, not part of the general
   optimisation policy.
3. Sequence priority is reproducible by fixing the admissible reserved sequence
   and maximising the other sequence's scale. This checks the allocation at a
   fixed measurement; it does not prove equivalence of arbitrary network equilibria.
4. A strict volt-var equality and proportional current backoff intentionally
   disagree under current saturation. Keep a regression for this distinction.
5. Raw exact PWL, encoded smooth relations, network equilibrium, hardware
   feasibility and solver status have separate residuals. Do not call a replay
   at the same voltage an independent network fixed-point validation.

## Test matrix

- Validate finite/order/sign/domain constraints, PF endpoints and invalid zero
  PF, unsupported topology/controller combinations, missing neutral sensing,
  request contradictions and unpublishable output masking.
- Analytic PQ capability, source direction, positive/negative/zero dispatch,
  rated/available watt ceilings, losses, zero availability and fixed PF signs.
- Balanced and free sequence allocation; converter versus grid balance under
  an LCL filter; three-leg zero sum; four-leg neutral and split-bank constraints.
- Physical current, S, voltage, ripple and capacitor margins at different power
  bases. Use independent expected values, not only stamp/extract identities.
- Deterministic fixed-measurement limiter reproductions, unsaturated feeder
  reproductions, and deliberate non-equivalence at current saturation.
- Existing advanced-inverter and controller regressions, API ambiguity/export
  checks, documentation examples and a full documentation build.

## Implementation sequence

1. Operating policies, source and allocation API; adapt stamping/extraction.
2. Topology validation and rating-scaled physical constraints.
3. Edge-case, capability and allocation tests; repair observed failures.
4. Independent optimisation reproductions and executable comparative examples.
5. Concept-by-concept documentation table and updated verification/API references.
6. Run regressions, inspect results and commit the completed slice.

Deferred: new curve graph backends in complete inverter models, continuous-time
switching constraint generation, PWM closure orchestration, variable DC rails,
GFM droop, new fleet/DOE scheduling APIs and continuous hardware design. Record
these as planned in the comparison table rather than implying implementation.

## Implemented and verified

Completed the six steps above. The implementation also normalises operating
power relationships and fixed-dispatch rows to the inverter rating: an
incompatible charging request exposed false solver success with poorly scaled
rows. A regression now rejects that case at 6 kVA, 1 MVA and 100 MVA bases.

Validation on Julia 1.12.6 / Ipopt:

- 3,303 checks passed across control lowering, advanced inverter, generalised
  generator, forward control, operating constraints, formulation variants,
  reproductions, controller numerics/reviews, fleet studies, experiment and
  sizing runners, closed-loop evidence and battery suites. One optional
  OpenDSSDirect comparison was skipped because that test dependency was absent.
- Export-definition and method-ambiguity checks passed.
- `scripts/formulations/ibr_reproductions.jl` ran independently. Three feeder
  comparisons differed by less than 0.1 W/var; proportional limiting by less
  than 0.32 W / 0.14 var; both fixed-measurement sequence-priority comparisons
  agreed within 1 µA. Assertions use the documented, slightly wider tolerances.
- Full Documenter build, executable examples and rendered comparison/reproduction
  links passed. Existing API-page/search-index size warnings remain.

The conceptual comparison table is in
`docs/src/components/inverter_operating_constraints.md`; the executable
explanation is in `docs/src/ibr/formulation_reproductions.md`. Deferred variants
remain explicitly identified rather than being inferred from plant support.
