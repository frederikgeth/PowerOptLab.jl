# References and related software

These sources support the distinctions used throughout the guides. A theorem
about a primitive, or a solver benchmark on another problem class, does not by
itself establish accuracy or scalability for a composed AC network model.

## Mathematical foundations

- **Convex sets, epigraphs and hypographs.** S. Boyd and L. Vandenberghe,
  *Convex Optimization* (Cambridge University Press, 2004), Chapters 2–3.
  [Author-hosted book](https://web.stanford.edu/~boyd/cvxbook/index.html).
  This is the basis for the exact supporting-line formulations in
  [bounds and relations](bounds_and_relations.md). A convex hypograph requires a
  concave function on the relevant domain; a graph equality is more restrictive.
- **Smooth complementarity functions.** C. Chen and O. L. Mangasarian,
  “A class of smoothing functions for nonlinear and mixed complementarity
  problems,” *Computational Optimization and Applications* 5, 97–138 (1996).
  [DOI](https://doi.org/10.1007/BF00249052).
  Integrated-density constructions motivate smooth positive-part functions.
  The specific compact quartic patch and its constants are derived in
  [the mathematical guide](theory.md).
- **Nonconvex smoothing theory.** X. Chen, “Smoothing methods for nonsmooth,
  nonconvex minimization,” *Mathematical Programming* 134, 71–99 (2012).
  [DOI](https://doi.org/10.1007/s10107-012-0569-0).
  Function-value approximation, gradient consistency and limiting stationarity
  are different properties; decreasing a width alone is not a convergence proof.
- **Mixed-integer PWL formulations.** J. Huchette and J. P. Vielma,
  “Nonconvex piecewise linear functions: Advanced formulations and simple
  modeling tools,” *Operations Research* 71(5), 1835–1856 (2023).
  [Author preprint](https://arxiv.org/abs/1708.00050),
  [DOI](https://doi.org/10.1287/opre.2019.1973).
  Consult this for exact graph encodings, formulation size and relaxation strength.
- **Complementarity stationarity and algorithms.** A. Nurkanović,
  A. Pozharskiy and M. Diehl, “Solving mathematical programs with complementarity
  constraints arising in nonsmooth optimal control” (2024 revision).
  [Preprint](https://arxiv.org/abs/2312.11022).
  This explains the constraint-qualification difficulties and distinct
  stationarity concepts relevant to MPCCs. Its optimal-control benchmarks do not
  determine performance on distribution-system models.
- **External MPCC solver.** A. Pozharskiy, F. Pacaud, M. Diehl and
  A. Nurkanović, “CCOpt: an Open-Source Solver for Large-Scale Mathematical
  Programs with Complementarity Constraints” (2026).
  [Preprint](https://arxiv.org/abs/2604.18726).
  CCOpt provides specialized relaxation, penalty and active-set methods built on
  MadNLP. PowerOptLab delegates these algorithms to the external solver.

## Modeling and solver libraries

| Resource | Role in these examples |
|:--|:--|
| [BMOPFTools.jl](https://github.com/frederikgeth/BMOPFTools.jl) | Distribution-system model construction and stable telescoping softplus evaluation |
| [JuMP nonlinear modeling](https://jump.dev/JuMP.jl/stable/manual/nonlinear/) | Expression tracing, scalar operators, derivatives and nonlinear model construction |
| [MathOptInterface solution semantics](https://jump.dev/MathOptInterface.jl/stable/manual/solutions/) | Termination, primal/dual result status and certificates |
| [PiecewiseLinearOpt.jl](https://github.com/jump-dev/PiecewiseLinearOpt.jl) | Exact bounded PWL graph encodings; formulation selection remains explicit |
| [CCOpt.jl](https://github.com/madsuite-org/CCOpt.jl) | MPCC optimizer and documented JuMP/MathOptComplements bridge setup |
| [Ipopt](https://coin-or.github.io/Ipopt/) and [MadNLP](https://madsuite.org/MadNLP.jl/stable/) | Smooth nonlinear optimization; local outcomes on nonconvex models |
| [HiGHS](https://highs.dev/) | Linear and mixed-integer linear examples |

PowerOptLab keeps canonical curve data, physical units, occurrence-specific
relations and candidate diagnostics above these libraries. A compatible scalar
encoding does not establish solver support for every other constraint in a model.
See [compilation](compilation.md) for the pathways and
[solver choice](tolerances_and_solvers.md) for whole-model limitations.

## Power-system tolerance semantics

EPRI's OpenDSS documentation describes
[checking the need for control action](https://opendss.epri.com/Checkingtheneedforcontrolaction.html),
[properties of functions that control](https://opendss.epri.com/PropertiesoffunctionsthatCONTROL.html),
[properties of functions that limit](https://opendss.epri.com/PropertiesoffunctionsthatLIMITac.html),
and [common properties](https://opendss.epri.com/Commonproperties.html).
The [DSS-Extensions InvControl reference](https://dss-extensions.org/dss-format/InvControl.html)
is a useful implementation-oriented companion; identify the implementation and
version when reproducing a study.

Read each tolerance with its comparison quantity, normalization base and time or
iteration index. A stopping test for successive voltages is not a static voltage
uncertainty band; a tolerance on a power limit is not necessarily a tolerance on
dispatched power. The [worked tolerance guide](tolerances_and_solvers.md) derives
explicit residual accounting and distinguishes these tests from an intentionally
modeled band around a control curve.
