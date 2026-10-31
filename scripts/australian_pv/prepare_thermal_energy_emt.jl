include("prepare_emt.jl")
module MixedThermalReference
    include("thermal_energy.jl")
end

function prepare_thermal_energy_emt(output,mixed_output=joinpath(ROOT,"results","thermal_energy"))
    repo=normpath(joinpath(ROOT,"..",".."))
    mixed_env=JSON3.read(read(joinpath(mixed_output,"environment.json"),String),Dict{String,Any})
    mixed_audit=JSON3.read(read(joinpath(mixed_output,"audit.json"),String),Dict{String,Any})
    mixed_audit["all_planned_accepted"] || error("audited mixed equilibria required")
    for (path,hash) in mixed_env["source_sha256"]
        bytes2hex(sha256(read(joinpath(repo,path))))==hash || error("stale mixed reference: $path")
    end
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    v=230 .* ComplexF64[1,cis(-2pi/3),cis(2pi/3)]
    cases=Dict{String,Any}[]
    base=plant("emt","poc";settings=TERMINAL_CONTROL_SETTINGS)
    fields=NamedTuple{fieldnames(AdvancedInverter)}(Tuple(getfield(base,k) for k in fieldnames(AdvancedInverter)))
    for (id,pwm,energy,four,q,fraction,negative) in (
            ("pwm_balanced",true,false,false,0.,1.,0im),
            ("energy_balanced",true,true,false,0.,1.,0im),
            ("energy_four_leg",true,true,true,-1800.,1.,0im),
            ("nonlinear_boundary",false,true,false,-1800.,.7,2.0+0im))
        inv=AdvancedInverter(;merge(fields,(topology=four ? :FOUR_LEG : :THREE_LEG,
            neutral=four ? "n" : nothing,In_max=four ? 10. : nothing,
            r_filter_neutral=four ? filter_r() : 0.,x_filter_neutral=four ? filter_x() : 0.,
            pwm_strategy=pwm ? :CENTERED : :NONE,f_sw=pwm ? 10000. : nothing,
            dc_energy_model=energy ? :nonlinear : :linearized))...)
        voltage=id=="nonlinear_boundary" ? v.+4.6.*[1,cis(2pi/3),cis(-2pi/3)] : v
        r=four ? allocate_phase_inverter_services(inv,request;grid_voltage=voltage,reactive_request=fill(q/3,3)) :
            allocate_inverter_services(inv,request;grid_voltage=voltage,reactive_request=q,
                active_fraction=fraction,negative_sequence_request=negative)
        solve_status(r).publishable || error("failed local EMT allocation: $id")
        record=emt_reference(inv,voltage,collect(r.allocation.grid_phase_current),
            collect(r.allocation.converter_terminal.phase_voltage),r.allocation.plant.p_dc;
            id,scope="fixed-sizing allocated reference; PWM=$pwm nonlinear-energy=$energy")
        record["expect_switched_hardware_accepted"]=pwm
        record["expect_averaged_hardware_accepted"]=true
        push!(cases,record)
    end
    phasors(v)=ComplexF64[complex(z[1],z[2]) for z in v]
    for caseid in ("energy_s1","energy_s2")
        case=only(c for c in MixedThermalReference.mixed_thermal_cases() if c.case_id==caseid)
        net=voltage_aware_network(parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json")),case)
        mixed=MixedThermalReference.mixed_thermal_fleet(net,case)
        point=JSON3.read(read(joinpath(mixed_output,"points",caseid*".json"),String),Dict{String,Any})
        truth=mixed_audit["cases"][caseid]["devices"]
        selected=Set(argmax(Dict(id=>r["capacitor_carrier_audit_current_A"] for (id,r) in truth if r["controller"]==label))
            for label in MixedThermalReference.MIXED_CONTROLLERS)
        for id in sort!(collect(keys(truth)))
            controlled=mixed.spec.devices[id];inv=controlled.device
            record=point["devices"][id];ts=point["bus"][inv.bus]
            voltage=[phasor(ts[p])-(inv.neutral===nothing ? 0im : phasor(ts[inv.neutral])) for p in inv.phase_terminals]
            exported=emt_reference(inv,voltage,phasors(record["grid_current"]),phasors(record["converter_voltage"]),
                truth[id]["p_dc_W"];id="mixed_$(caseid)_$id",scope="coupled mixed $caseid controller=$(mixed.labels[id]); customer $id")
            exported["controller"]=mixed.labels[id]
            exported["p_available_W"]=mixed.spec.requests[id].p_available
            exported["expect_switched_hardware_accepted"]=true
            exported["expect_averaged_hardware_accepted"]=true
            exported["save_waveforms"]=caseid=="energy_s2" && id in selected
            push!(cases,exported)
        end
    end
    calibration=JSON3.read(read(joinpath(@__DIR__,"thermal_calibration_negative.json"),String),Dict{String,Any})
    append!(cases,calibration["cases"])
    source=study_source_sha256()
    for file in ("four_leg.jl","voltage_aware.jl","thermal_energy.jl","prepare_emt.jl","prepare_thermal_energy_emt.jl")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    source[relpath(joinpath(@__DIR__,"thermal_calibration_negative.json"),repo)]=
        bytes2hex(sha256(read(joinpath(@__DIR__,"thermal_calibration_negative.json"))))
    mkpath(dirname(output))
    writejson(output,Dict("cases"=>cases,"source_sha256"=>source,"hardware"=>HARDWARE,
        "git_revision"=>strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String)),
        "scope"=>"fixed POC L-filter circuit references, including coupled mixed-controller equilibria; no service-loop/feeder dynamics",
        "dc_source"=>"constant input power at DC link, open carrier-frequency source; source current may ripple",
        "converter_loss"=>"fixed mean fitted RMS loss, no instantaneous semiconductor or switching-loss model",
        "timing"=>Dict("switching_Hz"=>10000.,"cycles"=>10,"discard_cycles"=>6,"steps_per_carrier"=>[32,64,128]),
        "validation_bands"=>Dict("fundamental_relative"=>.01,"fundamental_absolute_A"=>.01,
            "dc_2omega_relative"=>.02,"dc_2omega_absolute_V"=>.05,
            "mesh_current_A"=>.01,"mesh_dc_2omega_V"=>.01,"energy_relative"=>1e-6)))
    cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("usage: prepare_thermal_energy_emt.jl input.json audited_mixed_directory")
    prepare_thermal_energy_emt(ARGS[1],ARGS[2])
end
