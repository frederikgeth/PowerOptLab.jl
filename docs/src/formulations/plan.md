# Configurable function formulations: development plan

## Purpose

Deliver tools that enable researchers to investigate nonsmooth inverter controls
and operating-envelope formulations. PowerOptLab stages adaptable implementations,
worked examples and transparent diagnostics. This branch does not prescribe a
research protocol, declare a winning formulation, or require an engine-grade proof
and regression programme before an experimental idea can be used.

The same canonical function should support alternative representations, physical
error accounting and independent evaluation. Researchers choose the cases,
optimizers, options, tolerances, parameter sweeps and interpretation of outcomes.
External libraries may own continuation/homotopy. PowerOptLab will not add its own
continuation or automatic retry algorithm in this work.

## Delivery sequence

### A. Mathematical building blocks and physical budgets

- Retain the bounded PWL definitions and all existing representations.
- Add physical-error-to-width selection and explicit flat-extension evaluation.
- Add purpose-specific shifted-lower and upper magnitude primitives and their
  signed error bounds; keep exact positive-domain roots distinct.
- Supply simple composition and local equilibrium-sensitivity tools with explicit
  assumptions. Avoid automatic claims of whole-network safety from primitive error.
- Decouple complementarity normalization from the network coordinate base. Let
  callers choose the normalization and inspect physical residuals.
- Provide documented public extension methods so new smoothing families can be
  staged without editing a central list of built-in types.

### B. Configurable experiment infrastructure

- Cases are callbacks that construct models; methods supply representations,
  optimizers, options and configuration hooks. User-owned configuration records
  define sweeps, with no fixed feeder, source-voltage list or error target.
- Record build/solve times, raw solver outcomes, candidate diagnostics, errors and
  unsupported combinations. Preserve non-success candidates as diagnostics.
- Acceptance is an optional caller callback, separate from raw solver status.
- Export result data and caller-selected source fingerprints, with explicit
  version/schema metadata. Do not serialize solver/model objects as evidence.
- Keep previous fixed examples as reproducible illustrations of one configuration.

### C. Controller integration and example adapters

- Let existing PiecewiseLinearLaw users select softplus or local C2 independently
  for each curve. Preserve legacy constructors and default behavior.
- Keep numeric smooth evaluation faithful to the stamped formulation, including
  flat tails beyond the first/last control breakpoint.
- Expose adapters for real controlled-inverter studies. A generic case callback
  also supports feeder/DOE experiments without a new rigid campaign API.
- Keep complementarity and exact/hull graph integration explicit: substituting
  them for a smooth controller changes problem class and available guarantees.

### D. Researcher documentation

- An executable getting-started path from canonical curve to custom experiment.
- A mathematical guide to smoothing, signed error, composition, norm direction,
  conditioning, complementarity and graph/hull semantics.
- A controller tutorial, including how to vary curves and solver settings and
  diagnose non-publishable candidates.
- An extension tutorial showing a researcher adding a smoothing family and a
  custom case/metric without modifying the engine.
- Ground claims in primary sources. Separate analytic guarantees, local estimates,
  example observations and unassessed properties. Explain assumptions pedagogically.

## Verification proportionate to staging

Use focused checks for new contracts, model construction, representative solves,
status preservation and executable examples. Existing package regressions remain
useful for compatibility. Do not require every research solver/case combination to
converge or encode one campaign's numerical outcomes as universal acceptance
criteria. Exhaustive derivative/error proofs, platform-wide reliability campaigns,
performance guarantees and final API stabilization belong to a later proposal for
BMOPFTools' engine.

## Boundaries

A complete hybrid dynamical-system simulator, automatic whole-network safety
certification, a new MPCC/global solver, and general multivariate spline fitting
are outside this delivery. Composition, stationarity and equilibrium estimates
must state what they actually assess. Reuse external solver and graph libraries.

## Additional increment: shared primitives in advanced IBR controls

The PWL integration is complete. This refactor moves the existing algebraic
min/max/positive-part/clip operations and purpose-specific controller norms onto
shared primitives. Preserve established physical widths, upper-versus-lower norm
direction, exact hardware inequalities, expression sharing and default control
semantics. Keep exact implicit roots where a positive-domain restriction or
smoothing would change the feasible set. Reuse the shared magnitude primitive for
the advanced inverter's linear-current loss term where applicable.

Add focused derivative, scaling and same-state numeric/stamped comparisons; run
the existing inverter accuracy regressions. Document precisely which operations
are shared and which remain domain-specific. Commit working increments, then
request the user's review once these changes and the current A–D scope are ready.

Implementation is now delivered: algebraic hinge/selector/clip builders and their
signed contracts; shared controller extrema, capability norms and scale
composition; shifted unbalance and advanced-inverter loss norms backed by
BMOPFTools. Existing physical widths, exact roots, two-component norm provenance
and hardware constraints are preserved. The controller tutorial includes an
operation-by-operation integration map and explains the remaining domain logic.
Sixty focused selector checks, the existing controller start/base regressions and
the advanced inverter tests pass on Julia 1.10.11 and 1.12.6. All tutorials execute
in the documentation build. The final Julia 1.12.6 complete suite passes 5,238
functional assertions plus package-quality checks. This completes the additional
IBR increment and is the review point for the current branch scope. Engine
migration and broader research campaigns remain separate work.

## Existing foundation

The initial commits provide PWL softplus/C2/complementarity/exact-graph/hull
representations, an analytic scalar feeder, source-fingerprinted example evidence,
and optional backend CI. They pass 5,129 functional checks and 224 optional checks
on Julia 1.10 and 1.12. The initial CCOpt snapshot is an example observation
(0/20 strict successes), not a library-wide claim or an acceptance policy.

## Progress

- A delivered: public hinge extension methods, physical budget selection,
  magnitude/root constructions, signed affine error accounting, local sensitivity
  and independent complementarity normalization.
- B delivered: callback cases/methods/configurations, raw outcomes and candidate
  audits, caller-owned assessments, errors and versioned TOML export. Optional
  backend diagnostic-getter failures do not discard primal evidence; certificates
  are distinguished from candidate operating points.
- C delivered: selectable `PiecewiseLinearLaw` families, faithful flat-tail numeric
  evaluation/stamping, scalar analytic and physical inverter case adapters.
- D delivered: executable experiment, physical-error, controller and extension
  tutorials, with a dedicated API reference and primary scientific sources.

Validation of this increment: the full Julia 1.12.6 suite passes 5,176 functional
assertions plus package-quality checks. After the final diagnostic/certificate
handling changes, all 49 focused checks pass on Julia 1.10.11 and 1.12.6. The 1,135
primitive/model checks also pass on Julia 1.10.11. The expanded optional suite
passes 231 checks on Julia 1.12.6. Documentation builds with all examples executed.
These are compatibility and focused integration checks, not a new solver
reliability campaign or an engine-migration gate.


## Semantic control intent and compilation guide

The follow-up scope separates physical volt-var/watt intent from numerical
encoding. Deliver immutable `VoltVarWattIntent`, independent per-curve
`VoltVarWattEncoding`, lowering into existing smooth policies, and explicitly
bounded curve ports for smooth/exact/hull/complementarity experiments. Preserve
flat-tail semantics when restricting a port's voltage domain. Prepare immutable
hinge coefficients once and reuse identical model-local curve blocks.

Provide an executable explanation of numeric replay, staged and plain JuMP
nonlinear operators, solver bridges and external CCOpt homotopy. Expose physical
budgets, coordinate and complementarity scales, graph methods, solver options,
cache behavior and supported full-controller versus scalar-curve boundaries.
Demonstrate all encodings and Ipopt/MadNLP/HiGHS/CCOpt paths, plus structural reuse
counts. Update the physical inverter example to use semantic intent. Validate
focused semantics/domain/reuse checks, existing formulation/controller checks,
all documentation examples and optional solver integration. Whole-controller MIP
or MPCC compilation and general expression-graph optimization remain future
extensions, rather than being implied by scalar graph support.
