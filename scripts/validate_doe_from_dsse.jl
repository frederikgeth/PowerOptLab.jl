#!/usr/bin/env julia

# Reproducible DSSE → DOE validation runner.
#
# Usage:
#   julia --project=. scripts/validate_doe_from_dsse.jl path/to/case_builder.jl
#
# The builder must define `doe_validation_case()`, returning a NamedTuple with:
#   physics_net       passive DSSE network (no injecting devices / limits)
#   operational_net   same snapshot with known DSSE P/Q loads and DER/STATCOMs
#   measurements      Vector{Measurement} for the DSSE solve
#   connection_points Vector{ConnectionPoint}; each must bind a single-phase IBR
#
# Required adapter: materialize_estimate(estimate, operational_template).
# Optional fields: `with_statcom_net`, `truth_net`, `doe_keywords`.
# `truth_net` defaults to `operational_net`.  The independent AC check fixes each
# bound IBR's active-power setpoint to the issued DOE and calls BMOPFTools.solve_pf.

using PowerOptLab
using JuMP
using BMOPFTools: solve_pf, augment_case, opf_model

function _max_voltage_error(estimate, truth)
    Set(keys(estimate)) == Set(keys(truth)) || error("voltage comparison bus sets differ")
    errors = Float64[]
    for (bus, terminals) in estimate
        Set(keys(terminals)) == Set(keys(truth[bus])) || error("voltage terminal sets differ at $bus")
        for (terminal, value) in terminals
            reference = truth[bus][terminal]
            all(isfinite(get(entry, field, NaN)) for entry in (value, reference)
                for field in ("vr", "vi")) || error("missing/nonfinite voltage at $bus/$terminal")
            push!(errors, hypot(value["vr"] - reference["vr"], value["vi"] - reference["vi"]))
        end
    end
    isempty(errors) && error("no voltages to compare")
    return maximum(errors)
end

function _fixed_dispatch_net(net, cps, allocation, doe_snapshot;
                             direction, control_audit)
    direction in (:export, :import) || error("unknown DOE direction")
    fixed = deepcopy(net)
    # The pinned solve_pf strips s_max before stamping devices: VAR_MAX/S_MAX
    # droop bases become zero. Fixing p_max would also change a P_MAX law.
    # Do not silently replace a prescribed controller with a different equation.
    profiles = get(fixed, "control_profile", Dict())
    for (id, inv) in get(fixed, "ibr", Dict())
        profile = get(profiles, get(inv, "control_profile", ""), Dict())
        any(haskey(profile, law) for law in ("volt_var", "volt_watt")) &&
            throw(ArgumentError("independent PF cannot preserve rating-normalized " *
                "Volt-VAr/Volt-Watt laws for IBR '$id' with the pinned engine; " *
                "use verify_operating_envelope for faithful same-formulation replay"))
    end
    for cp in cps
        cp.ibr_id === nothing && error("independent PF validation requires ibr_id for '$(cp.id)'")
        inv = get(get(fixed, "ibr", Dict{String,Any}()), cp.ibr_id, nothing)
        inv isa Dict || error("IBR '$(cp.ibr_id)' is absent from operational_net")
        uppercase(String(get(inv, "topology", ""))) == "SINGLE_PHASE" || error(
            "independent PF validation currently supports single-phase IBRs; '$(cp.id)' is not")
        pmin, pmax = get(inv, "p_min", nothing), get(inv, "p_max", nothing)
        pmin isa AbstractVector && pmax isa AbstractVector && length(pmin) == length(pmax) || error(
            "IBR '$(cp.ibr_id)' must expose matching p_min/p_max vectors")
        setpoint = (direction == :export ? 1 : -1) * allocation[cp.id] / length(pmin)
        inv["p_min"] = fill(setpoint, length(pmin))
        inv["p_max"] = fill(setpoint, length(pmax))
    end
    # Reproduce the DOE's chosen operating point for other controllable assets
    # (for example a STATCOM). Customer IBR Q remains governed by its mandatory
    # control law rather than being converted into a dispatchable Q decision.
    for item in control_audit
        item["native_classification"] in (:free, :mixed) || continue
        id, quantity = item["id"], item["quantity"]
        item["component"] == :ibr && quantity in (:active_power, :reactive_power) ||
            error("independent PF adapter cannot freeze $(item["component"])/$id/$quantity")
        item["native_classification"] == :free || error("mixed control classification cannot be replayed")
        inv = fixed["ibr"][id]
        phase = item["position"]
        field = quantity == :active_power ? "p" : "q"
        value = doe_snapshot["ibr"][id][string(phase)][field * "g"]
        isfinite(value) || error("nonfinite control $id/$quantity/$phase")
        for bound in ("_min", "_max")
            values_ = copy(inv[field * bound])
            values_[phase] = value
            inv[field * bound] = values_
        end
    end
    fixed, _ = augment_case(fixed)
    return fixed
end

function _require_pf_witness(result)
    get(result, "termination_status", "UNKNOWN") in ("LOCALLY_SOLVED", "OPTIMAL") ||
        error("independent PF did not converge: $(get(result, "termination_status", "UNKNOWN"))")
    get(result, "doe_pf_primal_status", "NO_SOLUTION") == "FEASIBLE_POINT" ||
        error("independent PF did not return a feasible primal witness")
    return result
end

function _solve_checked_pf(net)
    result = solve_pf(net; per_unit=false, solution_hook! = (ctx, result) ->
        (result["doe_pf_primal_status"] = string(JuMP.primal_status(opf_model(ctx)))))
    return _require_pf_witness(result)
end

function _validate_snapshot(label, net, cps; doe_keywords=NamedTuple(), voltage_tolerance_V=1e-3)
    doe = solve_operating_envelope(net, cps; doe_keywords...)
    settings = doe.diagnostics[1]["formulation_settings"]
    settings.volt_var_watt_eps == 2e-3 && settings.context_hook! === nothing ||
        error("independent PF API cannot reproduce custom smoothing or context hooks")
    check = verify_operating_envelope(net, cps, doe; utilizations=:bound_point)
    # An optimal envelope binds a network constraint, so re-solving at the
    # issued capacity sits exactly on that boundary: the normalized margin here
    # is around -5e-9 and the joint solve can fail to converge on one platform
    # while succeeding on another. A joint solve that does not converge on the
    # boundary is not evidence of infeasibility, so only a candidate violation
    # is an error. `:unresolved` is reported, not silently accepted.
    outcomes = [diagnostics["verification_outcome"] for diagnostics in check.diagnostics]
    any(outcome -> outcome == :candidate_violation, outcomes) && error(
        "$label DOE verification returned a candidate violation at the issued capacity")
    pf_net = _fixed_dispatch_net(net, cps,
        Dict(cp.id => doe.envelope[cp.id][1] for cp in cps), doe.snapshots[1];
        direction=doe.direction, control_audit=doe.diagnostics[1]["control_audit"])
    pf = _solve_checked_pf(pf_net)
    difference = _max_voltage_error(doe.snapshots[1]["bus"], pf["bus"])
    difference <= voltage_tolerance_V || error("$label DOE/PF voltages disagree by $difference V")
    return (doe=doe, verification=check, pf=pf,
            verification_outcomes=outcomes,
            max_doe_pf_voltage_difference_V=difference)
end

function run_doe_validation(case)
    required = (:physics_net, :operational_net, :measurements, :connection_points, :materialize_estimate)
    all(name -> hasproperty(case, name), required) || error("case is missing one of $(required)")
    truth = _solve_checked_pf(get(case, :truth_net, case.operational_net))
    estimate = solve_state_estimation(case.physics_net, case.measurements)
    estimate.primal_status == "FEASIBLE_POINT" || error("DSSE did not return a feasible estimate")
    dsse_error = _max_voltage_error(estimate.bus, truth["bus"])
    kwargs = get(case, :doe_keywords, NamedTuple())
    operational = case.materialize_estimate(estimate, deepcopy(case.operational_net))
    base = _validate_snapshot("base", operational, case.connection_points; doe_keywords=kwargs)

    println("DSSE-to-DOE validation")
    println("  DSSE maximum voltage error [V]: ", dsse_error)
    println("  base DOE total capacity [W]:    ", base.doe.total_capacity[1])
    println("  base DOE/PF max |ΔV| [V]:       ", base.max_doe_pf_voltage_difference_V)
    println("  base fixed-capacity verdict:     ", base.verification_outcomes)

    statcom = nothing
    if hasproperty(case, :with_statcom_net)
        operational_statcom = case.materialize_estimate(estimate, deepcopy(case.with_statcom_net))
        statcom = _validate_snapshot("STATCOM", operational_statcom, case.connection_points;
            doe_keywords=kwargs)
        println("  STATCOM DOE total capacity [W]: ", statcom.doe.total_capacity[1])
        println("  STATCOM capacity gain [W]:      ",
                statcom.doe.total_capacity[1] - base.doe.total_capacity[1])
        println("  STATCOM DOE/PF max |ΔV| [V]:    ", statcom.max_doe_pf_voltage_difference_V)
        println("  STATCOM fixed-capacity verdict:  ", statcom.verification_outcomes)
    end
    return (estimate=estimate, base=base, statcom=statcom, max_dsse_voltage_error_V=dsse_error)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: julia --project=. scripts/validate_doe_from_dsse.jl path/to/case_builder.jl")
    include(abspath(only(ARGS)))
    Base.invokelatest(run_doe_validation, Base.invokelatest(doe_validation_case))
end
