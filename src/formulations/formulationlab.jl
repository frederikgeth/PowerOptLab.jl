# Compatibility for downstream studies; electrical formulations live in FormulationLab.
import FormulationLab
import FormulationLab: L3FOptions, L3FFinding, L3FApplicabilityReport, L3FInapplicableError, L3FReferenceState, CrossVoltageCoefficients, AffineScalarCoefficients, ConnectionPowerMap, LineDropCoefficients, L3FBuild, L3FResult, is_l3f_applicable, cross_voltage_coefficients, evaluate_cross_voltage, winding_voltage_coefficients, evaluate_affine, connection_power_map, line_drop_coefficients, regulator_gain_matrix, check_l3f_applicability, build_l3f_opf, l3f_model_class

function _formulationlab_replay(net; optimizer=nothing, solver_options=(), kwargs...)
    attributes = solver_options isa NamedTuple ? pairs(solver_options) : solver_options
    BMOPFTools.solve_pf(_kr_copy(net); optimizer=optimizer === nothing ? Ipopt.Optimizer : optimizer,
        solver_options=attributes, kwargs...)
end

function solve_l3f_opf(args...; powerflow=_formulationlab_replay, kwargs...)
    FormulationLab.solve_l3f_opf(args...; powerflow, kwargs...)
end

function l3f_reference_from_powerflow(args...; powerflow=_formulationlab_replay, kwargs...)
    FormulationLab.l3f_reference_from_powerflow(args...; powerflow, kwargs...)
end

function validate_l3f_solution(args...; powerflow=_formulationlab_replay, kwargs...)
    FormulationLab.validate_l3f_solution(args...; powerflow, kwargs...)
end

function solve_status(r::FormulationLab.L3FResult)
    s = FormulationLab.solve_status(r)
    SolveStatus(s.termination_status, s.primal_status, s.has_primal,
                s.feasible, s.optimal, s.publishable)
end
solve_diagnostics(r::FormulationLab.L3FResult) = FormulationLab.solve_diagnostics(r)
