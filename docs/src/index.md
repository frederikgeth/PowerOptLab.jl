# PowerOptLab.jl

PowerOptLab is a research laboratory for **four-wire distribution-network
decisions when the network state or model is uncertain**.

## Start here

Follow [Your first PowerOptLab study](getting_started.md) to install from a
checkout, build a complete network, schedule a battery, and check status and
physical residuals. The guide uses synthetic data and the standard solver, and
its Julia examples execute during documentation builds.

For an existing checkout, install the pinned dependencies from its root:

```sh
julia --project=. scripts/instantiate_pinned.jl
```

Then choose a study topic below.

## Study topics

The chapters bring models, worked studies, API references and evidence together
around the question you want to answer.

| Study topic | Question | Start here |
|:--|:--|:--|
| State and model inference | What states, shared parameters or line constructions are supported by the available evidence? | [Choose an inference study](estimation/studies.md) |
| Device modeling and scheduling | Which generator, storage or EV assumptions matter, and how are stored states linked over time? | [Choose device models](components/index.md) |
| Inverter control and hardware | How do local laws, converter topology and ratings change network and hardware outcomes? | [Inverter studies](ibr/index.md) |
| Dynamic operating envelopes | What capacity can be allocated and what does the declared verification establish? | [DOE in 15 minutes](tutorials/doe_getting_started.md) |
| Power flow and operability | What operating point, branch, static margin or model reduction is supported? | [Choose a power-flow or operability study](algorithms/index.md) |
| Numerical formulations and solvers | How do encodings, bounds and tolerances affect a declared mathematical model? | [Formulation studies](formulations/index.md) |

Use the [API index](api.md) to find a callable interface or the
[symbol index](symbol_index.md) to look up a name. The
[capability status catalogue](concepts.md#Capability-status) retains the maturity
and direction labels for the component and problem prototypes.

## Research context

The organizing questions are:

1. What states and models are compatible with telemetry and metadata?
2. What safe measurement or intervention would distinguish material
   alternatives?
3. What operating decision is justified by the remaining alternatives?

The package builds on the
[BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl) reference
current–voltage OPF engine and reuses its neutral-explicit device physics,
per-unit handling, and result extraction through public extension seams.

The complete loop is not implemented yet. Current capabilities provide its
foundations as separate, explicit research prototypes; the documentation does
not treat proposed joint model-forensics, active-probing, or general
uncertainty-set certification APIs as existing functionality.
[Concepts](concepts.md) explains the shared modeling
contracts; the [Research program](research_program.md) distinguishes current
capabilities from the intended integrated workflow.

## Extending the package

See [Concepts](concepts.md) for how each capability layers over the BMOPFTools
staged API, and [Contributing](contributing.md) for how to add your own and the
promotion path back to the BMOPF spec. The topic guides link equations, worked
examples, API references and validation evidence.
