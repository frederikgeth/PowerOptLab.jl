"""
    SolveOutcome

Normalized MathOptInterface solve status used by PowerOptLab entry points.
`has_primal` records whether the optimizer returned a candidate point,
`feasible` requires a fully feasible primal status, and `optimal` additionally
requires `OPTIMAL` or `LOCALLY_SOLVED` termination. Only `optimal` outcomes are
published as final numerical results by the research optimization wrappers.

`acceptable` identifies MOI's relaxed `ALMOST_OPTIMAL` /
`ALMOST_LOCALLY_SOLVED` outcomes. It is exposed for diagnostics but is not
silently promoted to an optimal result.
"""
struct SolveOutcome
    termination_status::JuMP.MOI.TerminationStatusCode
    primal_status::JuMP.MOI.ResultStatusCode
    result_count::Int
    has_primal::Bool
    feasible::Bool
    optimal::Bool
    acceptable::Bool
end

function _solve_outcome(model::JuMP.Model)
    termination = JuMP.termination_status(model)
    primal = JuMP.primal_status(model)
    count = JuMP.result_count(model)
    has_primal = count >= 1 && primal != JuMP.MOI.NO_SOLUTION
    feasible = has_primal && primal == JuMP.MOI.FEASIBLE_POINT
    optimal = feasible && termination in (JuMP.MOI.OPTIMAL, JuMP.MOI.LOCALLY_SOLVED)
    acceptable = has_primal &&
        termination in (JuMP.MOI.ALMOST_OPTIMAL, JuMP.MOI.ALMOST_LOCALLY_SOLVED) &&
        primal in (JuMP.MOI.FEASIBLE_POINT, JuMP.MOI.NEARLY_FEASIBLE_POINT)
    SolveOutcome(termination, primal, count, has_primal, feasible, optimal, acceptable)
end

_publishable(outcome::SolveOutcome) = outcome.optimal
_value_or_nan(outcome::SolveOutcome, value) =
    _publishable(outcome) ? JuMP.value(value) : NaN

_mask_unpublished(value::Bool) = value
_mask_unpublished(value::Number) = NaN
_mask_unpublished(value::AbstractDict) =
    Dict(key => _mask_unpublished(item) for (key, item) in value)
_mask_unpublished(value::AbstractVector) = [_mask_unpublished(item) for item in value]
_mask_unpublished(value) = value

function _extract_result(ctx, outcome::SolveOutcome)
    result = extract_result(ctx)
    _publishable(outcome) ? result : _mask_unpublished(result)
end

function _set_solver_options!(model::JuMP.Model, solver_options)
    options = solver_options isa NamedTuple ? pairs(solver_options) : solver_options
    for (name, value) in options
        JuMP.set_attribute(model, string(name), value)
    end
    return model
end
