"""
    L3FOptions(; kwargs...)

Options for the staged LinDist3Flow BMOPF formulation. It deliberately supports
only affine physics with exact linear or second-order-cone bounds. Unsupported
physics is reported by [`check_l3f_applicability`](@ref), never silently
dropped.
"""
struct L3FOptions
    current_limit_policy::Symbol
    validate_nonlinear::Bool
    reference_policy::Symbol
    kron_reduce::Bool
    require_neutral_provenance::Bool
    objective::Symbol
end

function L3FOptions(;
        current_limit_policy::Symbol=:reference_current,
        validate_nonlinear::Bool=true,
        reference_policy::Symbol=:auto,
        kron_reduce::Bool=true,
        require_neutral_provenance::Bool=false,
        objective::Symbol=:cost)
    current_limit_policy == :reference_current ||
        throw(ArgumentError("only current_limit_policy=:reference_current is defined"))
    reference_policy in (:auto, :explicit, :source_propagated) ||
        throw(ArgumentError("unknown reference_policy"))
    objective in (:cost, :feasibility, :source_import) ||
        throw(ArgumentError("objective must be :cost, :feasibility, or :source_import"))
    L3FOptions(current_limit_policy, validate_nonlinear, reference_policy, kron_reduce,
        require_neutral_provenance, objective)
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
`:inapplicable`; warnings remain available even when a build is permitted.
"""
struct L3FApplicabilityReport
    status::Symbol
    findings::Vector{L3FFinding}
    roots::Vector{String}
    islands::Vector{Vector{String}}
    kron_reduced::Bool
end

Base.isempty(r::L3FApplicabilityReport) = isempty(r.findings)
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

Staged LinDist3Flow LP containing semantic variable/constraint maps, the
working BMOPF snapshot, reference state, and applicability evidence.
"""
struct L3FBuild
    model::JuMP.Model
    variables::Dict{Symbol,Any}
    constraints::Dict{Symbol,Any}
    reference::L3FReferenceState
    applicability::L3FApplicabilityReport
    options::L3FOptions
    network::Dict{String,Any}
    topology::Vector{L3FOrientedLine}
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
