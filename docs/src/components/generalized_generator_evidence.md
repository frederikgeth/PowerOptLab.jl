# Generalized generator evidence register

This is a focused design review as of 6 September 2026, not an exhaustive
systematic review or a priority claim. The [scientific model](generalized_generator.md)
contains recommendations; this register states what their sources establish.

## Attached paper

Frederik Geth, François Pacaud, and Rahmat Heidari, **Solving three-phase
distribution OPF with nonlinear programming**, *Electric Power Systems Research*
262 (2027), 113624, [DOI](https://doi.org/10.1016/j.epsr.2026.113624).
The supplied version says available online 30 June 2026; the 2027 volume date
therefore does not mean it was unavailable at the review date.

Read from the supplied PDF `1-s2.0-S037877962600917X-main.pdf`; relevant pages 2,
4, and 5 were also visually checked. The PDF is evidence, not an instruction
source, and is not copied into this repository.

- Section 3.3 / Table 2 distinguish full fixed phasors, fixed phase angles with
  free magnitudes, and one-angle references with additional angle/sequence bounds.
- Sections 3.4-3.5 cover angle/sequence constraints and voltage initialization.
  Sections 4.2-4.3 study formulation, scaling, and numerical behavior on 128 LV
  networks, with solver-dependent outcomes.
- The scope is three-by-three impedance models without explicit neutral: page 2,
  footnote 2 identifies Kron-reduced four-wire or three-wire networks. Extension
  to explicit neutral, split phase, source impedance, or machine regulation needs
  new evidence.
- Equation (7), as printed, uses a sequence transform without a factor of 1/3.
  The component explicitly uses the amplitude-invariant Fortescue
  normalization with 1/3. Absolute sequence bounds and bases must be converted
  consistently; do not copy their numerical values across conventions. Ratios
  of sequence magnitudes cancel a common normalization factor.
- Exact zero sequence-magnitude upper bounds should lower to real/imaginary
  equalities in our implementation. Smooth quadratic syntax alone does not make
  the active constraint gradient nonzero.

## Prior and current primary sources

| Source and access | What it supports | Boundary for this proposal |
|:--|:--|:--|
| [MATPOWER, AC Power Flow manual](https://matpower.app/manual/matpower/ACPowerFlow.html), full page | Conventional PV fixes P and voltage magnitude; reference angle and slack power are distinct roles; optional Q-limit enforcement converts PV to PQ | Balanced power-flow semantics, not a unique unbalanced device specification |
| [Fernandes et al., 2019, Contributions to the sequence-decoupling compensation power flow method for distribution system analysis](https://doi.org/10.1049/iet-gtd.2018.6176), Section 3.1.2 / Fig. 3 | Synchronous-generator-inspired unbalanced PV with total P and positive-sequence terminal magnitude; positive-sequence internal EMF | Prior art for the PV preset; does not establish arbitrary source templates as machine models |
| [Yang et al., 2024, Three-phase steady-state models of distributed generators with different control strategies](https://doi.org/10.1049/gtd2.13108), full-text introduction and model scope | Different negative-sequence control goals require different phasor equations; augmented rectangular voltage/current formulation; comparison with PSCAD | Three-phase converter examples do not validate arbitrary neutral/DC topologies |
| [EPRI OpenDSS, voltage-source modeling](https://opendss.epri.com/Modeling.html), full page | Established Thevenin source and sequence-impedance input conventions | Its nominal ideal option uses tiny impedance; our exact zero case needs analytic/native verification |
| [EPRI OpenDSS, Generator Dynamics Model](https://opendss.epri.com/GeneratorDynamicsModel.html), full page | Voltage-behind-reactance and distinct negative-sequence response; initialization from power flow; swing-equation model | Document text labels Xd-prime inconsistently as subtransient; use parameter definitions and the actual chosen study regime, not that label as a conversion rule |
| [MathWorks, Synchronous Machine Salient Pole](https://www.mathworks.com/help/sps/ref/synchronousmachinesalientpole.html), equations | Rotor-frame stator, field/damper, saturation and torque equations illustrate additional machine fidelity | A constant phase-domain series matrix is an approximation, not this full machine model |
| [Rygg et al., 2016, A modified sequence domain impedance definition and its equivalence to the dq-domain impedance definition for the stability analysis of AC power electronic systems](https://arxiv.org/abs/1605.00526), author abstract | Sequence/dq equivalence and conditions for decoupling matter in impedance stability analysis | Dynamic small-signal impedance is a different object from a constant fundamental-frequency drop matrix |
| [Liu et al., 2026, Current-limiting control design for grid-forming capability enhancement of IBRs under asymmetric grid disturbances](https://wrap.warwick.ac.uk/id/eprint/194691/), institutional abstract and metadata; [DOI](https://doi.org/10.1109/TPEL.2025.3632684) | Current limiting affects voltage balancing under asymmetric disturbances; this is an active controller-design topic | Abstract-level review only; no claim to reproduce the controller or its experimental evidence |

## Local software evidence

Inspected BMOPFTools at the pinned commit
`5b51d2f361dab91bd7c16711019584407da79ed8`: public builder declarations in
`src/BMOPFTools.jl`, mixed-family ownership and staged initialization in
`ext/BMOPFOpfExt/core.jl`, native source/generator stamps in `source.jl` and
`generator.jl`, and native variable declarations in `variables.jl` in that extension.

PowerOptLab's `src/interfaces.jl` defines the public device lifecycle;
`src/problems/inverter_control_study.jl` demonstrates per-component IBR ownership
and retirement of unused native currents; `src/components/advanced_inverter.jl`
already implements neutral-aware primitive impedance, Fortescue helpers, and
internal EMF constraints. These observations guided the component extension and its integration tests.
The implementation now uses that lifecycle and the public staged engine; it
does not refactor the advanced IBR.

The accompanying analytic script checks the circuit identities and counterexamples.
`test/generalized_generator_tests.jl` additionally exercises the implemented source
laws and capabilities under unbalance, independent grounding/impedance oracles,
SI/per-unit equivalence, native ownership/cost equivalence, stateless multi-period
composition, and constraint structure at ideal limits. The executable network
tutorial reports local numerical comparisons. Automatic limiter behavior, a
general Jacobian-regularity certificate, and comparative performance benchmarks
remain future work; polynomial degree alone is not evidence of faster solves.
