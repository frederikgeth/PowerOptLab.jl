"""Common supertype for network devices contributed by PowerOptLab."""
abstract type AbstractDevice end

"""Common supertype for scalar telemetry consumed by estimation problems."""
abstract type AbstractMeasurement end

"""
Common supertype for structured solver results. Use [`solve_status`](@ref) and
[`solve_diagnostics`](@ref) instead of interpreting result-specific strings.
"""
abstract type AbstractSolveResult end

"""
    TimeGrid(durations_h)
    TimeGrid(periods, duration_h=1.0)

Validated period durations in hours. `durations_h[t]` weights both the period's
rate objective and its state transition, so nonuniform horizons do not silently
assume one-hour snapshots.
"""
struct TimeGrid
    durations_h::Vector{Float64}
    function TimeGrid(durations_h::AbstractVector{<:Real})
        isempty(durations_h) && throw(ArgumentError("TimeGrid needs at least one period"))
        values = Float64.(durations_h)
        all(x -> isfinite(x) && x > 0, values) || throw(ArgumentError(
            "TimeGrid durations must be finite and > 0 hours"))
        new(values)
    end
end

function TimeGrid(periods::Integer, duration_h::Real=1.0)
    periods >= 1 || throw(ArgumentError("TimeGrid periods must be >= 1"))
    TimeGrid(fill(duration_h, periods))
end

Base.length(grid::TimeGrid) = length(grid.durations_h)
Base.getindex(grid::TimeGrid, period::Integer) = grid.durations_h[period]
Base.iterate(grid::TimeGrid, state...) = iterate(grid.durations_h, state...)

function _resolve_time_grid(periods::Integer, dt_h::Real,
                            time_grid::Union{Nothing,TimeGrid})
    grid = time_grid === nothing ? TimeGrid(periods, dt_h) : time_grid
    length(grid) == periods || throw(ArgumentError(
        "TimeGrid has $(length(grid)) periods, expected $periods"))
    grid
end

"""A JuMP model and the BMOPFTools contexts built into that shared model."""
struct MultiContext{M,C}
    model::M
    contexts::Vector{C}
end

Base.length(multi::MultiContext) = length(multi.contexts)
Base.getindex(multi::MultiContext, period::Integer) = multi.contexts[period]

"""
    build_multi_context(nets; hook_factory, kwargs...) -> MultiContext

Build every network snapshot into one JuMP model through BMOPFTools' staged
API. `hook_factory(t)` returns the `model_hook!` for snapshot `t`. The builder
centralizes the shared optimizer, solver options, time index, and common unit
settings; callers add linking constraints/objectives and then enforce KCL. Pass
an existing `model` when shared variables/operators must be created first.
"""
function build_multi_context(nets::AbstractVector;
                             model=nothing,
                             hook_factory::Function = _ -> (_ -> nothing),
                             per_unit::Bool=true,
                             s_base::Real=1e6,
                             optimizer=Ipopt.Optimizer,
                             verbose::Bool=false,
                             solver_options=(),
                             context_options::NamedTuple=NamedTuple())
    periods = length(nets)
    periods >= 1 || throw(ArgumentError("need at least one snapshot"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError(
        "s_base must be finite and > 0"))
    model = model === nothing ? JuMP.Model(optimizer) : model
    verbose || JuMP.set_silent(model)
    _set_solver_options!(model, solver_options)
    contexts = [build_opf_model(nets[t]; model=model, t_index=t,
                    per_unit=per_unit, s_base=Float64(s_base),
                    add_objective=false, model_hook! = hook_factory(t),
                    context_options...)
                for t in 1:periods]
    MultiContext(model, contexts)
end

"""Normalized, user-facing status returned by [`solve_status`](@ref)."""
struct SolveStatus
    termination_status::String
    primal_status::String
    has_primal::Bool
    feasible::Bool
    optimal::Bool
    publishable::Bool
end

SolveStatus(outcome::SolveOutcome) = SolveStatus(
    string(outcome.termination_status), string(outcome.primal_status),
    outcome.has_primal, outcome.feasible, outcome.optimal, _publishable(outcome))

function _result_solve_status(termination_status::AbstractString,
                              publishable::Bool;
                              primal_status::AbstractString =
                                  publishable ? "FEASIBLE_POINT" : "UNKNOWN")
    SolveStatus(String(termination_status), String(primal_status), publishable,
                publishable, publishable, publishable)
end

"""Return the normalized [`SolveStatus`](@ref) for a structured solve result."""
function solve_status end

"""Return result-specific numerical diagnostics as a named tuple."""
solve_diagnostics(::AbstractSolveResult) = NamedTuple()

"""Return the stable identifier used to key a device's handles and results."""
device_id(device::AbstractDevice) = getfield(device, :id)

"""Validate a device against a horizon before adding variables or constraints."""
function validate_device end

"""Stamp one device into one BMOPFTools model context and return its handles."""
function stamp_device! end

"""Link one device's per-period handles across a [`TimeGrid`](@ref)."""
function link_device! end

"""Extract one device's published numerical result from its model handles."""
function extract_device end

"""Return the symbolic quantity represented by a scalar measurement."""
measurement_kind(m::AbstractMeasurement) = getfield(m, :kind)

"""Return a scalar measurement's SI value."""
measurement_value(m::AbstractMeasurement) = getfield(m, :value)

"""Return a scalar measurement's positive SI standard deviation."""
measurement_sigma(m::AbstractMeasurement) = getfield(m, :sigma)

function _validate_measurement_scalar(kind, value::Real, sigma::Real)
    isfinite(value) || throw(ArgumentError(
        "measurement $kind value must be finite (got $value)"))
    isfinite(sigma) && sigma > 0 || throw(ArgumentError(
        "measurement $kind sigma must be finite and > 0 (got $sigma)"))
    nothing
end
