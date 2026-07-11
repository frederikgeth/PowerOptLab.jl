"""
    PowerOptLab

A staging ground for component models and problem specifications built on top of
the BMOPFTools reference OPF engine, using its public extension seams
(`model_hook!` / `solution_hook!` and the staged `build_opf_model` /
`enforce_kcl!` / `generation_cost` / `extract_result` API) rather than forking
the engine.

The BMOPFTools engine deliberately stays small and correct — it ships one
four-wire rectangular current–voltage OPF and refuses to be a "model zoo".
PowerOptLab is where the zoo lives: experimental devices (storage, EVs) and
alternative problem specifications (multi-period co-optimisation, weighted
least-squares state estimation) that reuse the engine's device physics,
per-unit handling, and result extraction unchanged. Anything here that matures
into accepted practice can later be folded back into the spec and the engine.

## Modules of capability

- **Devices** ([`StorageDevice`](@ref), [`EVDevice`](@ref)) — battery/EV inverter
  ports stamped into a solve as current injections with an inter-temporal
  state-of-charge state.
- **Multi-period OPF** ([`solve_multiperiod_opf`](@ref)) — several network
  snapshots co-optimised in one model with storage/EV state linking each step to
  the next.
- **State estimation** ([`solve_state_estimation`](@ref)) — weighted
  least-squares estimation of the network state from noisy measurements, reusing
  the same device physics as a pure measurement-fitting problem.
- **Dynamic operating envelopes** ([`solve_operating_envelope`](@ref)) —
  per-connection-point DER export limits that respect the network's voltage and
  thermal constraints, recomputed per interval.
- **Advanced inverter** ([`AdvancedInverter`](@ref)) — a prototype internal-AC-node
  IBR with an output filter, internal-EMF/DC-modulation bounds, grid-forming
  operation, converter losses, and double-frequency ripple limits.

Everything is SI at the interface; per-unit conditioning inside the solve is
handled via the engine's `ctx.bases`.
"""
module PowerOptLab

using BMOPFTools
using JuMP
using Ipopt

include("devices.jl")
include("multiperiod.jl")
include("state_estimation.jl")
include("operating_envelope.jl")
include("advanced_inverter.jl")

# Devices
export StorageDevice, EVDevice

# Multi-period OPF
export solve_multiperiod_opf, MultiperiodResult

# State estimation
export Measurement, solve_state_estimation, StateEstimationResult

# Dynamic operating envelopes
export ConnectionPoint, solve_operating_envelope, OperatingEnvelopeResult

# Advanced inverter (prototype internal-node IBR)
export AdvancedInverter, solve_advanced_inverter, InverterResult

end # module PowerOptLab
