# Reproduce the offline voltage-allocation counterexample with one fixed plant.
# julia --project=. scripts/formulations/inverter_current_projection.jl [output.json]
using PowerOptLab, Ipopt, SHA
import JSON3

function current_projection_reference(output::Union{Nothing,AbstractString}=nothing)
    root=normpath(joinpath(@__DIR__,"..",".."))
    balanced(x)=x.*[1,cis(-2pi/3),cis(2pi/3)]
    # Synthetic rail regression fixture, not Australian hardware calibration.
    inv=AdvancedInverter(id="projection",bus="poc",phase_terminals=["a","b","c"],
        neutral="n",topology=:THREE_LEG,s_max=6000.,i_max=10.,v_dc=620.,
        c_dc=300e-6,m_max=.96,r_filter=.13225,x_filter=1.3225,
        p_loss_fixed=20.,c_loss=.01,dv2_max=7.5,i_cap_max=1.)
    options=("tol"=>1e-9,"max_iter"=>1000,"mu_strategy"=>"adaptive","nlp_scaling_method"=>"none")
    inputs=((label="nominal_request",voltage=230.,available=5000.,requested_p=4800.),
        (label="rail_conflict",voltage=253.,available=5000.,requested_p=4800.),
        (label="energized_zero_pv",voltage=230.,available=0.,requested_p=0.))
    phasors(x)=[(;real=real(z),imag=imag(z)) for z in x]
    rows=map(inputs) do input
        v=balanced(input.voltage)
        command=balanced(input.requested_p/(3input.voltage))
        request=InverterControlRequest(p_available=input.available,p_rated=6000.,q_scale=6000.)
        r=project_inverter_current(inv,request;grid_voltage=v,requested_current=command,
            target=:grid,solver_options=options)
        accepted=solve_status(r).publishable
        p=r.allocation.plant
        (;input...,grid_voltage_V=phasors(v),requested_grid_current_A=phasors(command),
            request_feasibility=r.request_feasibility,plant_feasibility=r.plant_feasibility,
            requested_rail_upper_V=r.requested_physical.rail.upper_V,
            requested_failed_checks=[c.name for c in r.requested_physical.checks if !c.passed],
            termination=r.allocation.termination_status,accepted,
            distance_A=accepted ? r.distance_A : nothing,
            p_poc_W=accepted ? p.p_poc : nothing,q_poc_var=accepted ? p.q_poc : nothing,
            p_dc_W=accepted ? p.p_dc : nothing,
            allocated_grid_current_A=accepted ? phasors(r.allocation.grid_phase_current) : nothing,
            rail_lower_V=accepted ? physical_feasibility(r).rail.lower_V : nothing,
            physical_checks=[(;name=c.name,unit=c.unit,value=isfinite(c.value) ? c.value : nothing,
                margin=isfinite(c.margin) ? c.margin : nothing,passed=c.passed)
                for c in physical_feasibility(r).checks])
    end
    paths=String[]
    for (dir,_,files) in walkdir(joinpath(root,"src")),file in files
        endswith(file,".jl") && push!(paths,relpath(joinpath(dir,file),root))
    end
    push!(paths,relpath(@__FILE__,root))
    hashes=Dict(path=>bytes2hex(sha256(read(joinpath(root,path)))) for path in sort(paths))
    artifact=(scope="Local offline current projection at fixed voltage; averaged model, PWM disabled. Synthetic fixed hardware; no global-nearest or dynamic claim.",
        julia=string(VERSION),Project_sha256=bytes2hex(sha256(read(joinpath(root,"Project.toml")))),
        Manifest_sha256=bytes2hex(sha256(read(joinpath(root,"Manifest.toml")))),source_sha256=hashes,
        hardware=Dict(string(f)=>getfield(inv,f) for f in fieldnames(AdvancedInverter)),
        target="grid arm excluding POC shunt",dc_source=(p_min_W=0.,),s_base_VA=1e5,
        solver_options=Dict(options),cases=rows)
    if output!==nothing
        mkpath(dirname(abspath(output)))
        open(output,"w") do io
            JSON3.pretty(io,artifact)
            println(io)
        end
    end
    for row in rows
        println(row.label,": request=",row.request_feasibility,", plant=",row.plant_feasibility,
            ", accepted=",row.accepted,", P/Q=",(row.p_poc_W,row.q_poc_var))
    end
    # Save diagnostics before rejecting an incomplete reference campaign.
    all(r.accepted for r in rows) || error("Reference includes an unaccepted allocation; inspect saved diagnostics")
    artifact
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)<=1 || error("Usage: inverter_current_projection.jl [output.json]")
    current_projection_reference(isempty(ARGS) ? nothing : only(ARGS))
end
