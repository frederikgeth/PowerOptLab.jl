# Fixed-measurement control reference. Does not run a feeder/chronological study.
# julia --project=. scripts/formulations/inverter_service_allocation.jl [output.json]
using PowerOptLab, SHA
import JSON3
module AustralianServiceReference
    include(joinpath(@__DIR__,"..","australian_pv","run.jl"))
end

function service_allocation_reference(output::Union{Nothing,AbstractString}=nothing)
    root=normpath(joinpath(@__DIR__,"..",".."))
    options=("tol"=>1e-9,"max_iter"=>1000,"mu_strategy"=>"adaptive","nlp_scaling_method"=>"none")
    a=cis(2pi/3);v=253 .* [1,a^2,a]
    synthetic=AdvancedInverter(id="synthetic",bus="poc",phase_terminals=["a","b","c"],neutral=nothing,
        topology=:THREE_LEG,s_max=6000.,i_max=10.,v_dc=620.,c_dc=300e-6,m_max=.96,
        r_filter=.13225,x_filter=1.3225,p_loss_fixed=20.,c_loss=.01,dv2_max=7.5,i_cap_max=1.)
    request=InverterControlRequest(p_available=5000.,p_rated=4800.,q_scale=6000.)
    cases=NamedTuple[]
    phasors(x)=[(real=real(z),imag=imag(z)) for z in x]
    finite(x)=isfinite(x) ? x : nothing
    function record(label,r;expected=true,measurement=v,grid=v)
        stage_rows=map(r.stages) do (name,s)
            status=solve_status(s)
            (name=name,termination=status.termination_status,primal=status.primal_status,
                optimal=status.optimal,accepted=status.publishable,
                raw_solver_status=s.operation.raw_solver_status,
                p_poc_W=finite(s.plant.p_poc),q_poc_var=finite(s.plant.q_poc),p_dc_W=finite(s.plant.p_dc),
                physical_checked=physical_feasibility(s).checked,physical_passed=physical_feasibility(s).passed,
                physical_checks=[(name=c.name,unit=c.unit,value=finite(c.value),
                    lower=finite(c.lower),upper=finite(c.upper),tolerance=c.tolerance,passed=c.passed)
                    for c in physical_feasibility(s).checks])
        end
        diagnostic=r.plant_diagnostic
        push!(cases,(label=label,expected_publishable=expected,accepted=solve_status(r).publishable,
            reactive_mode=r.reactive_mode,reactive_request_var=r.reactive_request_var,
            negative_scale=finite(r.negative_scale),plant_feasibility=r.plant_feasibility,
            policy_feasibility=r.policy_feasibility,sensed_voltage_V=phasors(measurement),
            grid_voltage_V=phasors(grid),
            current_command_A=r.current_command===nothing ? nothing : phasors(r.current_command),
            plant_diagnostic=diagnostic===nothing ? nothing :
                (accepted=solve_status(diagnostic).publishable,plant_feasibility=diagnostic.plant_feasibility,
                 p_poc_W=finite(diagnostic.allocation.plant.p_poc),q_poc_var=finite(diagnostic.allocation.plant.q_poc)),
            stages=stage_rows))
    end
    for (label,q,mode,expected) in (("synthetic_preserve_Q",-6000.,:within_request,true),
            ("synthetic_strict_ray_unresolved",-3600.,:within_request,false),
            ("synthetic_explicit_protective_override",-3600.,:protective_override,true))
        r=allocate_inverter_services(synthetic,request;grid_voltage=v,reactive_request=q,
            reactive_mode=mode,solver_options=options)
        record(label,r;expected)
    end
    aus=AustralianServiceReference
    plant=aus.plant("australian_reference","poc";settings=aus.TERMINAL_CONTROL_SETTINGS)
    # Deliberate local input: 1% negative sequence at 60 degrees, zero U0.
    unbalanced=v.+2.53cis(pi/3).*[1,a,a^2]
    for mode in aus.TERMINAL_MODES
        controlled=ControlledDevice(plant,aus.controller(mode;settings=aus.TERMINAL_CONTROL_SETTINGS,
            design=(kappa=0.,phi2=0.)))
        request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
        r=evaluate_voltage_aware(controlled,InverterControlMeasurement(unbalanced),request;
            grid_voltage=unbalanced,solver_options=options)
        record("australian_fixed_boundary_"*mode,r;measurement=unbalanced,grid=unbalanced)
    end
    controlled=ControlledDevice(plant,aus.controller("mean_vv_vw";
        settings=aus.TERMINAL_CONTROL_SETTINGS,design=(kappa=0.,phi2=0.)))
    night=InverterControlRequest(p_available=0.,p_rated=6000.,q_scale=6000.)
    r=evaluate_voltage_aware(controlled,InverterControlMeasurement(unbalanced),night;
        grid_voltage=unbalanced,solver_options=options)
    record("australian_fixed_boundary_energized_zero_PV",r;measurement=unbalanced,grid=unbalanced)
    hashes=aus.study_source_sha256()
    hashes[relpath(@__FILE__,root)]=bytes2hex(sha256(read(@__FILE__)))
    artifact=(scope="Solver-based local Q-P-I2 control reference at fixed measurements. No feeder equilibrium, response-time, global-optimality or EMT claim.",
        julia=string(VERSION),source_sha256=hashes,solver_options=Dict(options),s_base_VA=1e5,
        Project_sha256=bytes2hex(sha256(read(joinpath(root,"Project.toml")))),
        Manifest_sha256=bytes2hex(sha256(read(joinpath(root,"Manifest.toml")))),
        australian_hardware=Dict(string(k)=>getfield(plant,k) for k in fieldnames(AdvancedInverter)),
        synthetic_hardware=Dict(string(k)=>getfield(synthetic,k) for k in fieldnames(AdvancedInverter)),cases=cases)
    if output!==nothing
        mkpath(dirname(abspath(output)))
        open(output,"w") do io;JSON3.pretty(io,artifact);println(io);end
    end
    for c in cases
        println(c.label,": accepted=",c.accepted,", expected=",c.expected_publishable,
            ", plant=",c.plant_feasibility,", policy=",c.policy_feasibility)
    end
    all(c.accepted==c.expected_publishable for c in cases) || error("Unexpected allocation result; inspect saved stage diagnostics")
    artifact
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)<=1 || error("Usage: inverter_service_allocation.jl [output.json]")
    service_allocation_reference(isempty(ARGS) ? nothing : only(ARGS))
end
