# Research program

The research objective is to make distribution-network operating decisions with
uncertain measurements and network models. This requires methods for identifying
plausible models, selecting informative measurements and checking the resulting
operating limits. The current package provides several of these methods
separately; their combination remains a research objective.

The proposed sequence is:

```text
telemetry and metadata
        ↓
model forensics → informative, safe interventions
        ↓                    ↓
plausible model ensemble ← new evidence
        ↓
verified operating decisions
```

Several network models may fit the same measurements. The useful questions are
which parameters remain unidentified, whether the alternative models change an
operating decision, and which additional measurements would distinguish them.

## Current methods and limitations

The following methods are implemented. Their assumptions determine how the
results can be used in a combined estimation and operating study.

| Study | Available methods | Main limitations |
|---|---|---|
| Constrained state estimation | Neutral-explicit state, exact equations versus stochastic residuals, tangent-space observability, selected local covariance, constraint multipliers | Network identity and topology are fixed at compile time; there is no topology-hypothesis search |
| Operational parameter estimation | Shared line lengths and transformer taps over multiple snapshots | The line prototype is single-terminal per end; it does not jointly infer full mutually coupled impedances, phase maps, switch states, or grounding models |
| Inverse Carson reconstruction | Enumeration of overhead-construction candidates, compatible-candidate retention, local rank and profile intervals, materialized primitive four-wire models | It uses sequence linecode data, not joint operational telemetry, and does not combine candidates with the state estimator |
| Operating envelopes | Shared allocation across explicit forecast/model scenarios; bound-point, corner, and custom-utilization checks | Scenario models must be supplied by the caller; corner checks are not a global robust-feasibility certificate |
| Differentiated controller experiment | Local DiffOpt sensitivity of smooth lower-level PV/controller responses | Derivatives are conditional on the returned local solution and stable active set; discrete modes and alternative solution branches are not covered |
| Post-OPF operability | Native static snapshot checks, scaling provenance, directional sensitivities, natural/pseudo-arclength continuation, local fold refinement and sufficient fixed-point certificates | Generator/IBR/controller residual extensions and general global branch discovery are outside this scope |
| HELM power flow | No-load-connected branch with physical residual and series diagnostics | It does not discover every power-flow branch or certify non-existence when its finite series fails |

A compatible model fits the specified data; local observability concerns a
particular operating point; corner feasibility concerns the tested combinations.
The individual guides define these conditions and the numerical checks used.

## Priority 1: joint network-model identification

Joint identification would combine the existing estimation methods to compare
alternative network models. The uncertain model may include
switch and neutral continuity, phase labels, transformer terminal maps and taps,
line construction or mutually coupled impedance, source impedance, grounding,
and meter sign, scale, or reference errors.

The comparison should report:

- ranked discrete hypotheses and compatible continuous parameter regions;
- locally unidentified directions and observational equivalence classes;
- groups that are different physically but equivalent for a declared decision;
- model discrepancies that affect the intended operating decision; and
- candidate evidence that could distinguish material alternatives.

This requires comparing discrete configurations as well as fitting continuous
parameters over the same operating points. The four-wire equations and
measurement assumptions must be consistent across the comparisons. A local
optimum alone does not establish unique identification.

## Priority 2: differentiable experiment design

For a set of plausible network models, the next question is which measurements
or operating-point changes best distinguish them: meter placement or synchronization, inverter active/reactive
probing, a tap or capacitor action, or a suitable operating interval.

Candidate objectives include increasing the smallest identifiable singular
value, reducing a selected parameter or decision covariance, and separating
discrete hypotheses. Every intervention must also satisfy the same voltage,
current, converter, and disturbance limits used for operation.

The existing methods supply residual Jacobians, local observability measures
and sensitivities of selected operating responses. A complete experiment-design
solver would also need to enforce equipment limits and check the validity of
those sensitivities near changes in operating mode or solution branch.

## Priority 3: counterexample-guided decision verification

The main limitation of the current operating-envelope calculations is that
feasibility at the simultaneous bound, or even at every box corner, does not
prove feasibility throughout a nonconvex AC utilization box.

The current assumptions and numerical checks are documented in the
[DOE problem reference](problems/operating_envelope.md). The
detailed implementation audit and benchmark rationale remain in the [DOE
quantification scientific audit](problems/doe_quantification_review.md).
The accompanying [DOE literature evidence
register](problems/doe_literature_evidence.md) separates established results,
finite case evidence, reanalysis, and proposed experiments.

The next step is an outer allocation and inner counterexample search over both
actual customer utilization and the plausible model ensemble. Each discovered
violation becomes a new scenario for the allocation problem. Results should
distinguish at least:

- a falsification search that found a violating point;
- a search-stable result for which no violation was found under a declared
  multistart and stopping policy; and
- a certificate backed by a valid global relaxation bound.

A robust-feasibility certificate requires a valid bound covering the full
utilisation and uncertainty set. A search that finds no violation establishes
only the result of that search.

For a fixed polyhedral feasible region, non-enumerative box-containment
machinery already exists in the literature. The open target here is narrower:
use the nonlinear four-wire AC equations as the reference and test sensitivity-based
screening without treating random utilization samples as a certificate, and
eventually compare it with a valid linear robust inner-set baseline. A
four-wire assessment of sensitivity filtering remains further research.

## Foundational track: solution validity and differentiability

The inference and intervention program depends on local nonlinear solves. The
[static operability checker](problems/operability.md) already provides local
Jacobian diagnostics, finite-difference sensitivity checks, natural-parameter and
pseudo-arclength continuation, local fold refinement, and sufficient fixed-point
certificates within its native static scope. General multiple-solution discovery,
controller/active-set transition handling, and KKT or continuation validation of
DiffOpt sensitivities remain research extensions.

A reported sensitivity should identify the solution branch and operating mode,
its local conditioning, and any nearby change in active constraints. Alternative
local solutions may produce different sensitivities and operating decisions.

The [Post-OPF voltage operability reference](problems/operability.md) and
[study tutorial](tutorials/post_opf_operability_study.md) document the implemented
processor. They separate endpoint feasibility,
local Jacobian evidence, phase/sequence sensitivity, energized-base-connected
reachability, and path-specific margin rather than compressing them into one
stability index.

## Further research topics

Possible extensions include:

- **Multiconductor relaxations and global bounds:** develop SDP, SOCP, or QC
  relaxations for network identification and operating-envelope bounds;
  compare them against the neutral-explicit IVR model.
- **Hybrid controller equilibria:** represent discrete taps, switched
  capacitors, deadbands, hysteresis, and interacting inverter controls when the
  question is reachability, multiplicity, or cycling—not merely optimal dispatch.
- **Fundamental-plus-harmonic optimization:** add frequency-indexed primitives
  when harmonic voltage, neutral heating, resonance, or emission uncertainty
  changes a hosting-capacity or operating-envelope decision.

Moving-horizon estimation, decision-focused forecasting, safe learned local
controllers, decomposition, and protection-aware envelopes are also plausible
extensions when they close a demonstrated gap in this loop.

## Research scope

The emphasis is on distribution-network estimation, controlled experiments and
operating decisions with uncertain models. Additional device models, learning
methods or solution algorithms should address a specific question in that
programme. Their treatment of topology, grounding, uncertainty and model error
should be stated alongside the results.
