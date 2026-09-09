"""
    L3FOptions(; kwargs...)

Options for the staged LinDist3Flow BMOPF prototype. Version 0.1 deliberately
supports the minimal radial, constant-power, grounded-wye slice. Unsupported
physics is reported by [`check_l3f_applicability`](@ref), never silently
dropped.
"""
struct L3FOptions
    polygon_sides::Int
    current_limit_policy::Symbol
    current_load_policy::Symbol
    exponential_load_policy::Symbol
    droop_policy::Symbol
    include_fixed_losses::Bool
    validate_nonlinear::Bool
    reference_policy::Symbol
    trust_region_voltage::Union{Nothing,Float64}
    trust_region_control::Union{Nothing,Float64}
    kron_reduce::Bool
    require_neutral_provenance::Bool
    objective::Symbol
end

function L3FOptions(;
        polygon_sides::Integer=24,
        current_limit_policy::Symbol=:reference_current,
        current_load_policy::Symbol=:local_tangent,
        exponential_load_policy::Symbol=:local_tangent,
        droop_policy::Symbol=:reject,
        include_fixed_losses::Bool=false,
        validate_nonlinear::Bool=true,
        reference_policy::Symbol=:auto,
        trust_region_voltage=nothing,
        trust_region_control=nothing,
        kron_reduce::Bool=true,
        require_neutral_provenance::Bool=false,
        objective::Symbol=:cost)
    polygon_sides >= 3 || throw(ArgumentError("polygon_sides must be at least 3"))
    current_limit_policy == :reference_current ||
        throw(ArgumentError("only current_limit_policy=:reference_current is defined"))
    current_load_policy in (:local_tangent, :band_chord, :local_taylor) ||
        throw(ArgumentError("unknown current-load approximation policy"))
    exponential_load_policy in (:local_tangent, :band_chord, :local_taylor) ||
        throw(ArgumentError("unknown exponential-load approximation policy"))
    droop_policy == :reject ||
        throw(ArgumentError("the prototype only supports droop_policy=:reject"))
    reference_policy in (:auto, :explicit, :source_propagated) ||
        throw(ArgumentError("unknown reference_policy"))
    objective in (:cost, :feasibility, :source_import) ||
        throw(ArgumentError("objective must be :cost, :feasibility, or :source_import"))
    tv = isnothing(trust_region_voltage) ? nothing : Float64(trust_region_voltage)
    tc = isnothing(trust_region_control) ? nothing : Float64(trust_region_control)
    tv === nothing || tv > 0 || throw(ArgumentError("trust_region_voltage must be positive"))
    tc === nothing || tc > 0 || throw(ArgumentError("trust_region_control must be positive"))
    L3FOptions(Int(polygon_sides), current_limit_policy, current_load_policy,
        exponential_load_policy, droop_policy, include_fixed_losses,
        validate_nonlinear, reference_policy, tv, tc, kron_reduce,
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

"""Normals and inradius multiplier for a regular inner polygon."""
struct RegularPolygonCoefficients
    normals::Matrix{Float64}
    radius_scale::Float64
end

struct L3FOrientedLine
    id::String
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
