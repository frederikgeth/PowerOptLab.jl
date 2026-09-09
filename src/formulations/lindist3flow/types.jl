"""
    L3FOptions(; kwargs...)

Options for the staged LinDist3Flow BMOPF formulation. It deliberately supports
only affine physics with linear and native second-order-cone bounds. The
canonical network remains a fixed-angle, lossless LinDist3Flow approximation;
it is not an exact AC model. Unsupported
physics is reported by [`check_l3f_applicability`](@ref), never silently
dropped. BMOPF input and results are always SI; `per_unit=true` (the default)
uses BMOPFTools' public classic-base preparation for the independent model's
working coordinates, with system power base `s_base`.

| Field | Values | Meaning |
|:------|:-------|:--------|
| `current_limit_policy` | `:voltage_product` | Ampacity is the native rotated-SOC bound `p²+q² ≤ w Iᵐᵃˣ²`, using the live squared terminal or winding voltage. |
| `validate_nonlinear` | `Bool` (`true`) | Replay the optimized dispatch through BMOPFTools' nonlinear power flow and report the voltage difference. |
| `reference_policy` | `:auto`, `:explicit`, `:source_propagated` | Which linearization point to use. `:auto` prefers a supplied `reference` and otherwise propagates the source phasors; `:explicit` requires a supplied `reference`; `:source_propagated` always uses the propagated flat profile and ignores a supplied `reference`. |
| `kron_reduce` | `Bool` (`true`) | Kron-reduce an explicit-neutral input on a copy. When `false`, an explicit neutral is an error. |
| `require_neutral_provenance` | `Bool` (`false`) | Require recorded `_meta["kron_reduction"]` provenance, or an explicit reference, before building. |
| `unsupported` | `:reject` (default), `:lower`, `:approximate` | How to treat data outside the supported vocabulary. `:reject` reports it and refuses. `:lower` applies canonical L3F-preserving rewrites (switches, capacitors, line shunts, and justified single-phase transformer elements), reported as `L.L3F.*` at severity `:info`. `:approximate` additionally applies experimental lossy projections (constant-current and exponential load laws, adjustable taps, or delta gauge choices), reported as `A.L3F.*` at severity `:warning`; unsupported voltage bounds remain errors. |
| `objective` | `:cost`, `:feasibility`, `:source_import` | Linear per-channel energy cost, a zero objective, or total source active injection. |
| `per_unit` | `Bool` (`true`) | Optimization coordinates only. Input and results are SI either way. |
| `s_base` | `Real` (`1e6`) | System VA base for the per-unit working copy. |
"""
struct L3FOptions
    current_limit_policy::Symbol
    validate_nonlinear::Bool
    reference_policy::Symbol
    kron_reduce::Bool
    require_neutral_provenance::Bool
    unsupported::Symbol
    objective::Symbol
    per_unit::Bool
    s_base::Float64
end

function L3FOptions(;
        current_limit_policy::Symbol=:voltage_product,
        validate_nonlinear::Bool=true,
        reference_policy::Symbol=:auto,
        kron_reduce::Bool=true,
        require_neutral_provenance::Bool=false,
        unsupported::Symbol=:reject,
        objective::Symbol=:cost,
        per_unit::Bool=true,
        s_base::Real=1e6)
    current_limit_policy == :voltage_product ||
        throw(ArgumentError("only current_limit_policy=:voltage_product is defined"))
    reference_policy in (:auto, :explicit, :source_propagated) ||
        throw(ArgumentError("unknown reference_policy"))
    unsupported in (:reject, :lower, :approximate) ||
        throw(ArgumentError("unsupported must be :reject, :lower, or :approximate"))
    objective in (:cost, :feasibility, :source_import) ||
        throw(ArgumentError("objective must be :cost, :feasibility, or :source_import"))
    isfinite(s_base) && s_base > 0 ||
        throw(ArgumentError("s_base must be finite and > 0"))
    L3FOptions(current_limit_policy, validate_nonlinear, reference_policy, kron_reduce,
        require_neutral_provenance, unsupported, objective, per_unit, Float64(s_base))
end

"""Copy an option set while overriding selected fields by keyword."""
function _l3f_with_options(options::L3FOptions; kwargs...)
    names = fieldnames(L3FOptions)
    values = NamedTuple{names}(Tuple(getfield(options, name) for name in names))
    L3FOptions(; merge(values, (; kwargs...))...)
end

"""A stable, structured applicability diagnostic emitted by the L3F compiler."""
struct L3FFinding
    code::String
    severity::Symbol
    component::Symbol
    id::Union{Nothing,String}
    message::String
    evidence::Dict{String,Any}
end

L3FFinding(code, severity, component, id, message) =
    L3FFinding(String(code), Symbol(severity), Symbol(component),
        isnothing(id) ? nothing : String(id), String(message), Dict{String,Any}())

"""
    L3FApplicabilityReport

Result of [`check_l3f_applicability`](@ref). `status` is `:applicable` or
`:inapplicable`; warnings and `:info` findings remain available even when a
build is permitted. `lowered` records whether
`L3FOptions(unsupported=:lower|:approximate)` rewrote anything, in which case
the `L.L3F.*` and `A.L3F.*` findings say exactly what.
"""
struct L3FApplicabilityReport
    status::Symbol
    findings::Vector{L3FFinding}
    roots::Vector{String}
    islands::Vector{Vector{String}}
    kron_reduced::Bool
    lowered::Bool
end

L3FApplicabilityReport(status, findings, roots, islands, kron_reduced) =
    L3FApplicabilityReport(status, findings, roots, islands, kron_reduced, false)

"""`true` when the report carries no findings at all, warnings included. This is
strictly stronger than [`is_l3f_applicable`](@ref), which tolerates warnings."""
Base.isempty(r::L3FApplicabilityReport) = isempty(r.findings)

"""
    is_l3f_applicable(report) -> Bool

Whether [`build_l3f_opf`](@ref) will accept the network. `false` exactly when
the report contains at least one `:error` finding; warnings never block a build.
"""
is_l3f_applicable(r::L3FApplicabilityReport) = r.status == :applicable

"""Raised by `build_l3f_opf` when the applicability report contains an error."""
struct L3FInapplicableError <: Exception
    report::L3FApplicabilityReport
end

function Base.showerror(io::IO, e::L3FInapplicableError)
    errors = filter(f -> f.severity == :error, e.report.findings)
    print(io, "L3F-BMOPF is inapplicable")
    for finding in errors
        print(io, "\n  ", finding.code, ": ", finding.message)
    end
end

"""Immutable complex reference phasors used to construct L3F coefficients."""
struct L3FReferenceState
    voltage::Dict{Tuple{String,String},ComplexF64}
    provenance::Symbol
    source_hash::String
    nonlinear_status::Symbol
end

"""Affine coefficients `c + a*w_phi + b*w_psi` for one cross-voltage entry."""
struct CrossVoltageCoefficients
    constant::ComplexF64
    coefficient_phi::ComplexF64
    coefficient_psi::ComplexF64
end

"""Affine real scalar `constant + coefficients' * w`."""
struct AffineScalarCoefficients
    constant::Float64
    coefficients::Vector{Float64}
end

"""Reference-phasor channel-to-terminal complex-power allocation map."""
struct ConnectionPowerMap
    matrix::Matrix{ComplexF64}
    real_part::Matrix{Float64}
    imag_part::Matrix{Float64}
    reference_winding_voltage::Vector{ComplexF64}
end

"""Real LinDist3Flow voltage-drop coefficient matrices."""
struct LineDropCoefficients
    active::Matrix{Float64}
    reactive::Matrix{Float64}
end


"""
    L3FOrientedLine

One two-port device oriented away from its island's source, as published in
`L3FBuild.topology`. `parent_map`/`child_map` are the aligned conductor maps in
that orientation, and `reversed` records whether it is opposite to the input's
`bus_from`/`bus_to`.
"""
struct L3FOrientedLine
    id::String
    family::Symbol
    subtype::String
    parent::String
    child::String
    parent_map::Vector{String}
    child_map::Vector{String}
    reversed::Bool
end

"""
    L3FBuild

Staged LinDist3Flow LP/SOCP containing semantic variable/constraint maps, SI
and working-coordinate BMOPF snapshots/reference states, scaling bases, and
applicability evidence.
"""
struct L3FBuild
    model::JuMP.Model
    variables::Dict{Symbol,Any}
    constraints::Dict{Symbol,Any}
    reference::L3FReferenceState
    working_reference::L3FReferenceState
    applicability::L3FApplicabilityReport
    options::L3FOptions
    network::Dict{String,Any}
    working_network::Dict{String,Any}
    topology::Vector{L3FOrientedLine}
    bases::Any
end

"""
    L3FResult

Result of [`solve_l3f_opf`](@ref). Numerical fields are published only for an
optimal feasible solve under PowerOptLab's standard result contract.
"""
struct L3FResult <: AbstractSolveResult
    buses::Dict{String,Any}
    lines::Dict{String,Any}
    transformers::Dict{String,Any}
    generators::Dict{String,Any}
    sources::Dict{String,Any}
    objective::Float64
    formulation::Dict{String,Any}
    reference::L3FReferenceState
    applicability::L3FApplicabilityReport
    validation::Dict{String,Any}
    network::Dict{String,Any}
    solve::SolveStatus
end

solve_status(r::L3FResult) = r.solve
solve_diagnostics(r::L3FResult) = (
    objective=r.objective,
    formulation=r.formulation,
    applicability=r.applicability.status,
    validation=get(r.validation, "status", "not_run"),
)
