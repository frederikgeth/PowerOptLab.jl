# Offline allocation reference at a fixed electrical boundary. This deliberately
# does not replace the autonomous VV+VW controller or claim global optimality.

"""
    InverterCurrentProjectionResult

Offline current allocation at fixed POC voltage and midpoint balancing charge.
`request_feasibility` is `:feasible`, `:violated` or `:unresolved`, based on an
independent reconstruction of the requested point. `plant_feasibility` is
`:feasible` when either that point or the accepted allocation is a witness;
otherwise it is `:unresolved`. A failed local solve never proves plant
infeasibility. These statements apply only to the declared boundary/source model.

`allocation` retains the raw solve status and physical report. `distance_A` is
the Euclidean distance across the three complex arm currents (NaN unless the
allocation is publishable). `requested_state` and `requested_physical` are
`nothing` if the filter transfer is singular or numerically indeterminate.
"""
struct InverterCurrentProjectionResult <: AbstractSolveResult
    target::Symbol
    requested_current::NTuple{3,ComplexF64}
    requested_state::Union{Nothing,NamedTuple}
    requested_physical::Union{Nothing,InverterFeasibilityReport}
    request_feasibility::Symbol
    plant_feasibility::Symbol
    distance_A::Float64
    allocation::ConstrainedInverterResult
end

solve_status(r::InverterCurrentProjectionResult) = solve_status(r.allocation)
physical_feasibility(r::InverterCurrentProjectionResult) = physical_feasibility(r.allocation)
solve_diagnostics(r::InverterCurrentProjectionResult) = (
    target=r.target, request_feasibility=r.request_feasibility,
    plant_feasibility=r.plant_feasibility, distance_A=r.distance_A,
    requested_physical=r.requested_physical, allocation=solve_diagnostics(r.allocation))

function _projection_phasors(values, label)
    length(values)==3 || throw(DimensionMismatch("$label must contain three RMS phasors"))
    result=ComplexF64.(collect(values))
    all(isfinite,result) || throw(ArgumentError("$label must be finite"))
    result
end

function _projection_phase_impedance(inv,side)
    R,X=_series_filter_matrices(inv,3,1.;side)
    Z=complex.(R,X)
    inv.neutral===nothing && return Z
    B=[Matrix{Float64}(I,3,3); -ones(1,3)]
    B'*Z*B
end

# Recover the requested full filter state without optimization auxiliaries.
# At a converter-arm target, (I + Y Zg) Vm = Vpoc + Zg Ic. A singular
# transfer need not make the command impossible; it prevents this direct audit.
function _projection_requested_state(inv,v,current,target,charge)
    Zc=_projection_phase_impedance(inv,:converter)
    Zg=_projection_phase_impedance(inv,:grid)
    Y=complex(_mid_filter_admittance(inv)...)
    if target===:grid
        ig=current
        vm=v+Zg*ig
        ic=ig+Y*vm
    else
        A=Matrix{ComplexF64}(I,3,3)+Y*Zg
        singular=svdvals(A)
        minimum(singular)>128eps(Float64)*max(1.,maximum(singular)) || return nothing
        ic=current
        vm=A\(v+Zg*ic)
        ig=ic-Y*vm
    end
    u=vm+Zc*ic
    all(isfinite,[u;vm;ic;ig]) || return nothing
    (;converter_voltage=u,converter_current=ic,grid_voltage=v,grid_current=ig,
        filter_voltage=vm,midpoint_charge=charge)
end

function _projection_request_feasibility(report)
    report===nothing && return :unresolved
    report.passed && return :feasible
    # A conservative lower enclosure below zero is not a violating witness.
    # The rail upper bound is a sampled headroom plus roundoff padding.
    violated=any(report.checks) do c
        !c.passed && isfinite(c.value) &&
            (c.name!="continuous_rail_headroom" || report.rail.upper_V < -c.tolerance)
    end
    violated ? :violated : :unresolved
end

# The plant audit alone does not assert that a network solver respected the
# imposed boundary. Certify those conditions separately in SI as well.
function _projection_boundary_report(report,voltage,fixed_voltage,charge,fixed_charge)
    checks=copy(report.checks)
    tolerance=InverterFeasibilityTolerances()
    for k in 1:3
        push!(checks,_inverter_limit_check("projection_fixed_voltage_$k",:V,
            abs(voltage[k]-fixed_voltage[k]),0.,0.,tolerance))
    end
    push!(checks,_inverter_limit_check("projection_fixed_midpoint_charge",:C,
        abs(charge-fixed_charge),0.,0.,tolerance))
    InverterFeasibilityReport(report.checked,report.passed && all(c.passed for c in checks),
        checks,report.rail)
end

function _projection_network(inv,v)
    terminals=copy(inv.phase_terminals)
    ground=String[]
    if inv.neutral!==nothing
        push!(terminals,inv.neutral)
        push!(ground,inv.neutral)
    end
    data=Dict("bus"=>Dict(inv.bus=>Dict("terminal_names"=>terminals,
        "perfectly_grounded_terminals"=>ground)),
        "voltage_source"=>Dict("projection_boundary"=>Dict("bus"=>inv.bus,
            "terminal_map"=>inv.phase_terminals,"v_magnitude"=>abs.(v),"v_angle"=>angle.(v))))
    parse_bmopf(JSON3.write(data);from_string=true)
end

"""
    project_inverter_current(inv, request; grid_voltage, requested_current,
        target=:converter, dc_source=DCSourceBounds(p_min=0.), midpoint_charge=0.,
        s_base=1e5, optimizer=Ipopt.Optimizer, verbose=false, solver_options=())

Find a local least-squares current allocation within all declared inverter limits
at **fixed** POC RMS voltage phasors. `target` selects converter-arm or grid-arm
current; grid-arm current excludes the separate POC shunt. Phase order follows
`inv.phase_terminals`. Voltage phasors refer to the inverter neutral (or earth
when absent). The metric is the sum of squared complex-current errors, with a
common numerical normalization that does not change its minimizer.

The reference includes filter/neutral coupling, losses, source bounds, sequence
and thermal limits, and adaptive full-cycle rail cuts. Default DC bounds model a
non-charging PV source; availability remains an upper bound. Split-link balancing
charge is fixed to the explicit `midpoint_charge` in **both** the requested-point
audit and the allocation. This is a conditional capability test, not optimization
over that actuator. Non-split topologies require zero charge.

Supports the same three-phase, non-grid-forming plants and centered-PWM scope as
[`ConstrainedInverter`](@ref). The requested point is reconstructed and audited
independently. A violated request together with a feasible allocation demonstrates
a command/capability conflict under these conditions. Failure to find an allocation
is inconclusive. This nonconvex offline reference does not certify a globally
nearest point, implement reactive-priority VV+VW, or validate closed-loop dynamics.
"""
function project_inverter_current(inv::AdvancedInverter,request::InverterControlRequest;
        grid_voltage,requested_current,target::Symbol=:converter,
        dc_source::DCSourceBounds=DCSourceBounds(p_min=0.),midpoint_charge::Real=0.,
        s_base::Real=1e5,optimizer=Ipopt.Optimizer,verbose::Bool=false,solver_options=())
    target in (:converter,:grid) || throw(ArgumentError("target must be :converter or :grid"))
    isfinite(s_base) && s_base>0 || throw(ArgumentError("s_base must be finite and positive"))
    isfinite(midpoint_charge) || throw(ArgumentError("midpoint_charge must be finite"))
    charge=Float64(midpoint_charge)
    if inv.topology!==:SPLIT_DC
        iszero(charge) || throw(ArgumentError("midpoint_charge requires SPLIT_DC"))
    elseif inv.q_mid_balance_max===nothing
        iszero(charge) || throw(ArgumentError("nonzero midpoint_charge requires balancing authority"))
    else
        abs(charge)<=inv.q_mid_balance_max || throw(ArgumentError("midpoint_charge exceeds balancing authority"))
    end
    device=ConstrainedInverter(inv,InverterOperatingConstraints();
        current_allocation=:free,dc_source)
    validate_device(device)
    v=_projection_phasors(grid_voltage,"grid_voltage")
    desired=_projection_phasors(requested_current,"requested_current")
    net=_projection_network(inv,v)
    validate_device(device,(net,);periods=1)
    state=_projection_requested_state(inv,v,desired,target,charge)
    report=state===nothing ? nothing : inverter_physical_feasibility(inv;state...,
        dc_min=dc_source.p_min,dc_max=_available_dc(device,request),
        pwm_thermal_headroom=inv.pwm_strategy===:NONE ? nothing : device.pwm_thermal_headroom)
    handles=Ref{_ConstrainedInverterHandles}()
    # A common scale preserves the SI metric and is independent of the OPF base.
    current_scale=inv.i_max===nothing ? inv.s_max/(3max(1.,maximum(abs,v))) : inv.i_max
    hook! = ctx -> begin
        handles[]=stamp_device!(ctx,device;request)
        h=handles[].plant; model=_opf_model(ctx)
        re,im=target===:converter ? (h.cri,h.cii) : (h.gri,h.gii)
        JuMP.@objective(model,Min,sum(
            ((re[k]*h.ib-real(desired[k]))/current_scale)^2+
            ((im[k]*h.ib-imag(desired[k]))/current_scale)^2 for k in 1:3))
        # Fix the actuator to the same value used by the direct request audit.
        h.qbal===nothing || JuMP.fix(h.qbal,charge;force=true)
    end
    ctx=build_opf_model(net;per_unit=true,s_base=Float64(s_base),add_objective=false,
        model_hook! = hook!,optimizer,verbose)
    _set_inverter_solver_options!(_opf_model(ctx),solver_options)
    _optimize_inverter_rails!(ctx,((inv,handles[].plant),))
    outcome=_solve_outcome(_opf_model(ctx)); status=SolveStatus(outcome)
    extracted=extract_device(device,handles[],status)
    network=_extract_result(ctx,outcome)
    h=handles[].plant
    physical=extracted.plant.physical_feasibility
    if status.has_primal
        actual_voltage=complex.(JuMP.value.(h.vrpoc),JuMP.value.(h.vipoc)).*h.vb
        actual_charge=h.qbal===nothing ? 0. : JuMP.value(h.qbal)
        physical=_projection_boundary_report(physical,actual_voltage,v,actual_charge,charge)
    end
    plant=merge(extracted.plant,(physical_feasibility=physical,))
    operation=merge(extracted.operation,(physical_feasible=physical.passed,
        physical_failed_checks=[c.name for c in physical.checks if !c.passed]))
    accepted=_physical_solve_status(status,physical)
    allocation=ConstrainedInverterResult(string(outcome.termination_status),plant,
        operation,extracted.converter_terminal,extracted.grid_phase_current,
        network["bus"],accepted)
    actual=target===:converter ? allocation.converter_terminal.phase_current : allocation.grid_phase_current
    distance=accepted.publishable ? sqrt(sum(abs2,collect(actual)-desired)) : NaN
    request_status=_projection_request_feasibility(report)
    plant_status=request_status===:feasible || accepted.publishable ? :feasible : :unresolved
    InverterCurrentProjectionResult(target,Tuple(desired),state,report,request_status,
        plant_status,distance,allocation)
end
