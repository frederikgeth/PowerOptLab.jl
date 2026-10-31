include("engineering_records.jl")
module EngineeringBenchmarkEMTReferences
    include("prepare_emt.jl")
end

function prepare_engineering_benchmark_emt(folder,output)
    ispath(output) && error("preserve existing EMT inputs")
    read_record(path)=JSON3.read(read(path,String),Dict{String,Any})
    env=read_record(joinpath(folder,"environment.json"));audit=read_record(joinpath(folder,"audit.json"))
    audit["all_planned_accepted"] || error("complete matched benchmark audit required")
    repo=normpath(joinpath(ROOT,"..",".."))
    for (file,digest) in env["source_sha256"]
        bytes2hex(sha256(read(joinpath(repo,file))))==digest || error("stale benchmark source: $file")
    end
    bytes2hex(sha256(read(joinpath(folder,"environment.json"))))==audit["environment_sha256"] || error("benchmark environment changed")
    bytes2hex(sha256(read(joinpath(@__DIR__,"audit_engineering_benchmark.jl"))))==audit["audit_code_sha256"] || error("benchmark audit source changed")
    for (file,digest) in audit["record_sha256"]
        bytes2hex(sha256(read(joinpath(folder,file))))==digest || error("benchmark point changed: $file")
    end
    cases=Dict{String,Any}[]
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    phasors(x)=ComplexF64[complex(z...) for z in x]
    for condition in engineering_conditions()
        net=voltage_aware_network(base,condition);mixed=engineering_fleet(net,condition)
        for mode in ("local","coordinated")
            key=condition.case_id*"_"*mode
            point=read_record(joinpath(folder,"points",key*".json"))
            metrics=engineering_record_metrics(point,mixed;mode)
            metrics.accepted || error("reference reconstruction rejected: $key")
            for label in MIXED_CONTROLLERS
                candidates=sort!([r for r in metrics.rows if r.local_controller==label];
                    by=r->(-r.capacitor_thermal_A,r.device_id))
                selected=first(candidates);id=selected.device_id
                inv=mixed.spec.devices[id].device;request=mixed.spec.requests[id]
                record=point["devices"][id];ts=point["bus"][inv.bus]
                v=ComplexF64[phasor(ts[p])-(inv.neutral===nothing ? 0im : phasor(ts[inv.neutral])) for p in inv.phase_terminals]
                i=phasors(record["grid_current"]);u=phasors(record["converter_voltage"])
                actual=inverter_physical_feasibility(inv;converter_voltage=u,converter_current=i,
                    grid_voltage=v,grid_current=i,dc_min=0.,dc_max=request.p_available,
                    pwm_thermal_headroom=PowerOptLab._PWM_ALLOCATION_HEADROOM)
                actual.passed || error("actual source-availability reference rejected")
                reference=EngineeringBenchmarkEMTReferences.emt_reference(inv,v,i,u,selected.p_dc_W;
                    id=key*"_"*label*"_"*id,
                    scope="matched Australian fixed-POC snapshot; highest capacitor model stress in this controller group; released VV/VW commands in coordinated mode")
                # The older generic helper uses a 5000 W fixture source. This
                # campaign validates and preserves each actual site availability.
                reference["p_available_W"]=request.p_available
                reference["reference_physical_check_count"]=length(actual.checks)
                reference["parent_case_id"]=condition.case_id
                reference["parent_mode"]=mode;reference["local_controller"]=label
                reference["device_id"]=id;reference["parent_capacitor_model_A"]=selected.capacitor_thermal_A
                reference["parent_point_sha256"]=audit["record_sha256"][joinpath("points",key*".json")]
                reference["expect_averaged_hardware_accepted"]=true
                reference["expect_switched_hardware_accepted"]=true
                push!(cases,reference)
            end
        end
    end
    source=study_source_sha256()
    for file in ("thermal_energy.jl","four_leg.jl","voltage_aware.jl","switching_engineering.jl",
            "coordinated_engineering.jl","engineering_records.jl","prepare_emt.jl","emt_circuit.jl",
            "prepare_engineering_benchmark_emt.jl")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    mkpath(dirname(output))
    write_engineering_json(output,Dict("cases"=>cases,"source_sha256"=>source,"hardware"=>HARDWARE,
        "operating_limits"=>ENGINEERING_OPERATING_LIMITS,"sizing_changed"=>false,
        "benchmark_environment_sha256"=>audit["environment_sha256"],
        "benchmark_audit_sha256"=>bytes2hex(sha256(read(joinpath(folder,"audit.json")))),
        "benchmark_record_sha256"=>audit["record_sha256"],
        "selection_rule"=>"one highest model capacitor-current site per controller label, case and mode; exact ties resolved by device ID",
        "julia_version"=>string(VERSION),"git_revision"=>strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String)),
        "git_status"=>read(Cmd(["git","-C",repo,"status","--porcelain"]),String),
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Project.toml")))),
        "Manifest_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Manifest.toml")))),
        "timing"=>Dict("switching_Hz"=>10000.,"cycles"=>10,"discard_cycles"=>6,"steps_per_carrier"=>[32,64,128]),
        "validation_bands"=>Dict("fundamental_relative"=>.01,"fundamental_absolute_A"=>.01,
            "dc_2omega_relative"=>.02,"dc_2omega_absolute_V"=>.05,"mesh_current_A"=>.01,
            "mesh_dc_2omega_V"=>.01,"mesh_peak_current_A"=>.01,"mesh_dc_extrema_V"=>.01,"energy_relative"=>1e-6),
        "scope"=>"fixed POC L-filter snapshots, constant-power DC inputs and fixed mean converter loss; no feeder/controller loops, periodic-steady-state or all-time DC drift certificate"))
    println("Prepared ",length(cases)," matched-case EMT references")
    cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("usage: prepare_engineering_benchmark_emt.jl benchmark_folder output.json")
    prepare_engineering_benchmark_emt(ARGS[1],ARGS[2])
end
