# Function formulation layer: implementation plan

## Objective and boundary

Stage a small reusable mathematical layer in PowerOptLab before proposing stable
pieces for BMOPFTools. Preserve canonical function semantics while selecting a
smooth surrogate, an exact graph encoding, or an explicitly labelled outer
relaxation. External libraries may own homotopy/continuation; PowerOptLab will not
implement a continuation or retry algorithm for this addition.

The first executable slice covers bounded, continuous, scalar PWL curves (including
hinges, deadbands, and clamps), fixed smoothing parameters, a shared scalar
feeder/controller equilibrium, and physical-unit error accounting. Existing
inverter defaults are unchanged. Arbitrary discontinuous steps, hysteresis,
multivariate surfaces, automatic error propagation through a network, and changes
to the inverter's purpose-specific norm formulas are subsequent work.

## Commit sequence and acceptance

1. **Canonical curves and contracts.** Validated immutable data; physical domain
   and unit labels; exact numeric oracle; softplus and compact C2 approximations;
   analytic derivative and signed error bounds. Tests at knots, endpoints, narrow
   bands, and rescaled coordinates. Reuse BMOPFTools softplus evaluation.
2. **JuMP representations.** Model builders with explicit coordinate scaling;
   reuse BMOPFTools' expression builder on staged contexts; analytic scalar
   operators for standalone models; native MOI complementarity; bounded vertex
   hull; optional PiecewiseLinearOpt exact graph extension. Returned handles carry
   semantics and allow independent exact-graph discrepancy assessment. Check
   derivative integration, constant/affine cases, and hull-versus-graph witnesses.
3. **Executable comparison and external solvers.** A small voltage/controller
   equilibrium with analytic segment reference; compare at matched output-error
   budgets over starts and physical bases. Optional isolated CCOpt / MIP environment
   with a reproducible setup and CI job. Record solver status, graph discrepancy,
   electrical residual, approximation allowance, and relaxation semantics
   independently; never promote a relaxed candidate to a physical certificate.
4. **Documentation and validation.** Executable tutorial, API documentation,
   assumptions, scientific sources, tests on supported Julia versions, full package
   regressions and documentation build. Record remaining research questions.

## Mathematical contracts

- All canonical breakpoints, values, widths and absolute error budgets use the
  declared physical units. Coordinate scales convert physical values to solver
  variables and preserve the underlying problem.
- Softplus and compact C2 smoothing operate on the same clamped PWL hinge
  expansion. Signed coefficient sums bound signed function error; zero coefficients
  are removed. Smoothing is a surrogate graph, not a physical safety certificate.
- A continuous vertex hull is the exact convex hull of the bounded PWL graph in
  (input, output) coordinates, but composing hulls with other equations generally
  gives only an outer relaxation of the full model.
- Complementarity and segment selection encode the canonical graph mathematically;
  numerical feasibility, complementarity and solver optimality remain separate.
- Finite domains, unit labels, smoothness, error direction, and representation kind
  are inspectable. Unsupported combinations fail explicitly.

## Reuse and references

BMOPFTools supplies stable telescoping softplus. PiecewiseLinearOpt supplies exact
PWL graph formulations. JuMP/MOI supplies complementarity constraints, and CCOpt
with MathOptComplements supplies their external solution machinery. Closed-form
local polynomial patches and vertex hull equations are small enough to maintain
here without a general spline or global-optimization dependency.

- Chen and Mangasarian (1996): https://doi.org/10.1007/BF00249052
- Chen (2012): https://doi.org/10.1007/s10107-012-0569-0
- Huchette and Vielma: https://arxiv.org/abs/1708.00050
- Nurkanovic, Pozharskiy and Diehl (2024): https://arxiv.org/abs/2312.11022
- CCOpt (2026 preprint): https://arxiv.org/abs/2604.18726

## Progress

- Branch created from merged main f65cb32; canonical primitives and all five
  representations implemented, with optional PiecewiseLinearOpt delegation.
- Shared electrical reference, executable tutorial, isolated external-backend
  setup, source-fingerprinted snapshot, and dedicated optional CI job implemented.
- Julia 1.12: 5,129 functional assertions pass; documentation examples build.
  Julia 1.10 validation is in progress.
- Exact graph: 20/20 strict successes and physical reference checks. CCOpt 0.1.0:
  0/20 strict successes, 8/20 microunit canonical-equation checks. The latter is
  explicitly a characterization backend; its numerical limitations remain work
  rather than being hidden by the optional integration test result.
- Subsequent research: physical complementarity normalization and stationarity,
  purpose-specific norm/root contracts, and realistic inverter/fleet integration.
  No migration to BMOPFTools is proposed until those interfaces are established.
