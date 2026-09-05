# DOE literature: evidence and interpretation

> **Review date:** 5 September 2026 · **Scope:** literature used to define the
> claims, examples, and research boundaries of the PowerOptLab DOE framework.

Dynamic operating envelope (DOE) does not name one mathematical object across
the literature. It can mean a time-varying import/export limit, a feasible
operating region, a decoupled range intended to tolerate independent customer
use, a market-shaped allocation, or a pointwise request-and-approval process.
PowerOptLab therefore records the *claim* separately from the DOE label.

The particularly strong object used in robust range studies is a **decoupled
box DOE**. Given a utilization set ``\mathcal U``, uncertainty set ``\Xi``, and
declared control policy, every represented realization must admit a network
state satisfying the selected operational model. An ordinary bound-point DOE
does not make that full-range claim.

## Evidence classes

This review uses four evidence classes:

| Label | Meaning |
|---|---|
| Primary result | A claim established or numerically demonstrated in the cited research paper |
| Review evidence | A synthesis or field result reported through a review; the underlying trial report has not necessarily been checked |
| Reanalysis | A calculation or interpretation made from published values |
| Research hypothesis | A testable proposal that is not established by the cited work |

Numbers from a finite numerical experiment remain evidence about that case.
They are not promoted to universal properties, calibrated probabilities, or
mathematical certificates.

## Method taxonomy

| Object | Participant freedom | Typical evidence | PowerOptLab interpretation |
|---|---|---|---|
| Bound-point allocation | One simultaneous operating point | Local OPF solution and replay | `security=:bound_point`; not a range guarantee |
| Finitely tested range | Independent use at enumerated or sampled points | Corner, scenario, or adaptive-search results | Tested-point evidence; may falsify but cannot certify the nonlinear box |
| Certified decoupled box | Independent use throughout a Cartesian product | Valid inner-set or containment certificate | Future baseline for fixed polyhedral models; unavailable for the present nonlinear AC model |
| Coupled operating region | Aggregator coordinates a point in a polytope, ellipse, or sampled region | Region construction plus disaggregation policy | More expressive, but changes communication, privacy, and control responsibilities |
| Market-shaped envelope | Limits depend on bids, products, or customer value | Market/network co-optimization | Capacity or volume alone is not the objective |
| Request–approve–curtail | A requested setpoint is accepted or modified | Feasible pointwise AC dispatch | Stronger coordination with less autonomous range freedom |

Control recourse is orthogonal to this taxonomy. Any row may assume one fixed
issue-time setting, a scenario-dependent setting, an implemented local law, or
pointwise perfect recourse. [`PerfectRecourse`](@ref) remains supported so that
published anticipative formulations can be replicated without presenting them
as operationally non-anticipative.

## Core evidence register

| Evidence | Supported conclusion | Important boundary |
|---|---|---|
| [Lankeshwara et al. (2025)](https://doi.org/10.1016/j.rser.2025.115696) | Reviews DOE models, allocation objectives, uncertainty approaches, implementation architectures, and field activity. | A review supports taxonomy and synthesis; trial figures cited through it are secondary evidence. |
| [Barzegar et al. (2026)](https://doi.org/10.1016/j.apenergy.2026.128172) | Provides a current perspective on DOE advances, projects, regulation, and unresolved implementation questions. | Reported project outcomes should be traced to project reports before being treated as primary field evidence. |
| [Moring and Mathieu (2023)](https://doi.org/10.1109/SmartGridComm57358.2023.10333939) | Demonstrates that an SOCP relaxation can materially over- or underestimate operating envelopes; AC replay of the 56-bus relaxed result is unsafe. | One set of cases establishes a fidelity failure mode, not that every SOCP DOE is inexact. |
| [Liu and Braslavsky (2022)](https://doi.org/10.1109/ACCESS.2022.3203062) | Demonstrates that moving customers inside equal limits can worsen cross-phase voltage and violate the network in an unbalanced feeder. | It is a published counterexample to general monotonicity, not proof that all bound-point envelopes fail. |
| [Liu and Braslavsky (2024)](https://doi.org/10.1109/TPWRS.2023.3308104) | Defines a decoupled hyperrectangle inside a fixed linearized feasible region and gives a Motzkin-based polyhedral containment construction. | The certificate applies to the declared polyhedral model. AC-model error remains separate; local solution of the final bilinear expansion affects maximality, not containment of an MTT-feasible returned box. |
| [Liu, Braslavsky, and Mahdavi (2023)](https://doi.org/10.35833/MPCE.2023.000653) | Develops norm-bounded robust linear feasible regions under impedance and load uncertainty. | The large capacity sensitivity to impedance is one feeder result under unmatched uncertainty radii and norms. |
| [Liu and Ma (2023)](https://doi.org/10.1109/PESGM52003.2023.10253037) | Constructs an impedance-robust feasible region and evaluates it against random impedance realizations. | Reported 0%, 62%, and 100% figures are finite sample-containment rates, not probabilities or certificates. |
| [Liu and Braslavsky (2024)](https://doi.org/10.17775/CSEEJPES.2024.03220) | Uses operating-point sensitivities to reduce the vertices represented in a nonlinear unbalanced formulation and reports strong computational scaling. | Thirty thousand random utilization tests are useful evidence but do not prove that the filtered vertices dominate the full box. |
| [Attarha et al. (2024)](https://doi.org/10.1016/j.epsr.2024.110639) | Shows that bid-aware shaped envelopes can produce greater market benefit than capacity-oriented DOE baselines in the studied trial data. | This is a counterexample to envelope size being a sufficient value proxy, not evidence of general anticorrelation; the formulation also uses richer information and soft limits. |
| [Badmus and Pandey (2024)](https://doi.org/10.1109/CDC56724.2024.10886763) | Demonstrates a request–approve–curtail alternative using a three-phase nonlinear AC formulation and different curtailment norms. | Pointwise local nonlinear feasibility is not a global guarantee, and its fairness quantities are not directly comparable with envelope-width fairness. |

## Quantitative anchors and their permitted use

The following values are useful motivating examples when their case boundaries
remain attached:

- In the 56-bus case of Moring and Mathieu, the equal nonlinear envelope is
  161.5 kW per node and the SOCP result is 627.9 kW; AC replay produces severe
  overvoltage. The 3.89 ratio is a case result, not a generic correction factor.
- In one Australian-network case, a 10% impedance uncertainty set reduces the
  reported optimized-reactive-power total from 112.3 kW to 43.1 kW, while a
  20% load uncertainty set gives 109.0 kW. This motivates a matched uncertainty
  experiment; the different physical uncertainties and set geometries prevent
  a universal dominance conclusion.
- Liu and Ma report that 62 of 100 randomly generated Australian-network
  feasible regions contain their computed robust region, compared with zero of
  100 for the initial region. This reveals both the importance of parameter
  uncertainty and the local-algorithm limitation.
- The shaped-envelope case reports greater benefit for a smaller bid-aware
  envelope than for its best capacity-oriented DOE comparator. The defensible
  conclusion is that capacity is not sufficient to predict delivered value.

## Open questions, separated from established results

The most useful next research questions are:

1. Can a nonlinear or four-wire DOE receive a checkable box-containment
   certificate without assuming a unique power-flow solution?
2. Under what network, grounding, and operating-point conditions does a
   sensitivity-filtered vertex set dominate all omitted vertices?
3. How should parameter uncertainty and participant-utilization coverage be
   composed without conflating a robust feasible region with a robust box?
4. How much apparent capacity is attributable to anticipative network-control
   recourse, and how much survives under issue-time or causal local control?
5. Under matched experimental designs, which uncertainty sources most affect
   safety and useful capacity?
6. Which object is required over time: a sequence of independent limits, a
   finite-horizon feasible trajectory set, a viability kernel, or a controlled
   invariant set?
7. Which stakeholder outcome—access, delivered energy, curtailment, value, or
   reliability—is the subject of a fairness claim, and over what horizon?

The four-wire sensitivity-filtering experiment, intertemporal envelopes, and a
polyhedral robust baseline are separate research tracks. They are not claims of
the present implementation.
