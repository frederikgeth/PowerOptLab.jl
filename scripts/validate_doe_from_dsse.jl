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
# Optional fields: `with_statcom_net`, `truth_net`, `doe_keywords`.
# `truth_net` defaults to `operational_net`.  The independent AC check fixes each
# bound IBR's active-power setpoint to the issued DOE and calls BMOPFTools.solve_pf.

using PowerOptLab
using BMOPFTools: solve_pf, augment_case

length(ARGS) == 1 || error("usage: julia --project=. scripts/validate_doe_from_dsse.jl path/to/case_builder.jl")
include(abspath(only(ARGS)))
isdefined(Main, :doe_validation_case) || error("case builder must define doe_validation_case()")

case = doe_validation_case()
required = (:physics_net, :operational_net, :measurements, :connection_points)
all(name -> hasproperty(case, name), required) || error("case is missing one of $(required)")

function _max_voltage_error(estimate, truth)
    errors = Float64[]
    for (bus, terminals) in estimate
        haskey(truth, bus) || continue
        for (terminal, value) in terminals
            haskey(truth[bus], terminal) || continue
            vm = get(value, "vm", NaN)
            tm = get(truth[bus][terminal], "vm", NaN)
            isfinite(vm) && isfinite(tm) && push!(errors, abs(vm - tm))
        end
    end
    return isempty(errors) ? NaN : maximum(errors)
end

function _fixed_dispatch_net(net, cps, allocation, doe_snapshot)
    fixed = deepcopy(net)
    customer_ibrs = Set{String}()
    for cp in cps
        cp.ibr_id === nothing && error("independent PF validation requires ibr_id for '$(cp.id)'")
        inv = get(get(fixed, "ibr", Dict{String,Any}()), cp.ibr_id, nothing)
        inv isa Dict || error("IBR '$(cp.ibr_id)' is absent from operational_net")
        uppercase(String(get(inv, "topology", ""))) == "SINGLE_PHASE" || error(
            "independent PF validation currently supports single-phase IBRs; '$(cp.id)' is not")
        pmin, pmax = get(inv, "p_min", nothing), get(inv, "p_max", nothing)
        pmin isa AbstractVector && pmax isa AbstractVector && length(pmin) == length(pmax) || error(
            "IBR '$(cp.ibr_id)' must expose matching p_min/p_max vectors")
        setpoint = allocation[cp.id] / length(pmin)
        inv["p_min"] = fill(setpoint, length(pmin))
        inv["p_max"] = fill(setpoint, length(pmax))
        push!(customer_ibrs, cp.ibr_id)
    end
    # Reproduce the DOE's chosen operating point for other controllable assets
    # (for example a STATCOM). Customer IBR Q remains governed by its mandatory
    # control law rather than being converted into a dispatchable Q decision.
    for (id, inv) in get(fixed, "ibr", Dict{String,Any}())
        id in customer_ibrs && continue
        qmin, qmax = get(inv, "q_min", nothing), get(inv, "q_max", nothing)
        result = get(get(doe_snapshot, "ibr", Dict{String,Any}()), id, nothing)
        qmin isa AbstractVector && qmax isa AbstractVector && result isa Dict || continue
        length(qmin) == length(qmax) || continue
        qset = Float64[]
        for phase in 1:length(qmin)
            entry = get(result, string(phase), nothing)
            entry isa Dict && haskey(entry, "qg") || empty!(qset); break
            push!(qset, Float64(entry["qg"]))
        end
        length(qset) == length(qmin) || continue
        inv["q_min"] = qset
        inv["q_max"] = qset
    end
    fixed, _ = augment_case(fixed)
    return fixed
end

function _validate_snapshot(label, net, cps; doe_keywords=NamedTuple())
    doe = solve_operating_envelope(net, cps; doe_keywords...)
    check = verify_operating_envelope(net, cps, doe; utilizations=:bound_point)
    all(check.feasible) || error("$label DOE failed its nonlinear fixed-capacity verification")
    pf_net = _fixed_dispatch_net(net, cps,
        Dict(cp.id => doe.envelope[cp.id][1] for cp in cps), doe.snapshots[1])
    pf = solve_pf(pf_net; per_unit=false)
    difference = _max_voltage_error(doe.snapshots[1]["bus"], pf["bus"])
    return (doe=doe, verification=check, pf=pf,
            max_doe_pf_voltage_difference_V=difference)
end

truth = solve_pf(get(case, :truth_net, case.operational_net); per_unit=false)
estimate = solve_state_estimation(case.physics_net, case.measurements)
estimate.primal_status == "FEASIBLE_POINT" || error("DSSE did not return a feasible estimate")
dsse_error = _max_voltage_error(estimate.bus, truth["bus"])
kwargs = get(case, :doe_keywords, NamedTuple())
base = _validate_snapshot("base", case.operational_net, case.connection_points; doe_keywords=kwargs)

println("DSSE-to-DOE validation")
println("  DSSE maximum voltage error [V]: ", dsse_error)
println("  base DOE total capacity [W]:    ", base.doe.total_capacity[1])
println("  base DOE/PF max |ΔV| [V]:       ", base.max_doe_pf_voltage_difference_V)
println("  base fixed-capacity verified:    ", all(base.verification.feasible))

if hasproperty(case, :with_statcom_net)
    statcom = _validate_snapshot("STATCOM", case.with_statcom_net, case.connection_points;
        doe_keywords=kwargs)
    println("  STATCOM DOE total capacity [W]: ", statcom.doe.total_capacity[1])
    println("  STATCOM capacity gain [W]:      ",
            statcom.doe.total_capacity[1] - base.doe.total_capacity[1])
    println("  STATCOM DOE/PF max |ΔV| [V]:    ", statcom.max_doe_pf_voltage_difference_V)
end
