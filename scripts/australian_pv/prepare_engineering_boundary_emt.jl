include("switching_engineering.jl")
module EngineeringEMTReferences
    include("prepare_emt.jl")
end

function engineering_boundary_references(output)
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    v=230 .* ComplexF64[1,cis(-2pi/3),cis(2pi/3)]
    cases=Dict{String,Any}[]
    for (id,four,guard) in (("combined",false,true),("lf_only_boundary",false,false),("four_leg_combined",true,true))
        fixture=Dict{String,Any}("ibr"=>Dict{String,Any}("pv"=>Dict{String,Any}("bus"=>"poc")))
        base=four ? first(values(four_leg_fleet(fixture,1.).devices)).device :
            plant("pv","poc";settings=TERMINAL_CONTROL_SETTINGS)
        inv=engineering_plant(base;limits=ENGINEERING_LEGACY_LIMITS,dc_min=747.,dc_max=753.,dc_guard=guard)
        r=four ? allocate_phase_inverter_services(inv,request;grid_voltage=v,reactive_request=[-600.,0.,0.],active_fraction=fill(.1,3),solver_options=ENGINEERING_OPTIONS) :
            allocate_inverter_services(inv,request;grid_voltage=v,reactive_request=-300.,active_fraction=.1,
                negative_sequence_request=1.0+0im,solver_options=ENGINEERING_OPTIONS)
        solve_status(r).publishable || error("boundary allocation failed: $id")
        p=r.allocation
        record=EngineeringEMTReferences.emt_reference(inv,v,collect(p.grid_phase_current),
            collect(p.converter_terminal.phase_voltage),p.plant.p_dc;id,
            scope="fixed-sizing LF/carrier voltage boundary; shared 14 A declared operating peak ceiling")
        record["expect_averaged_hardware_accepted"]=true
        record["expect_switched_hardware_accepted"]=guard
        record["dc_guard_enabled"]=guard
        record["allocation_stages"]=[Dict("stage"=>String(n),"termination"=>solve_status(s).termination_status,
            "raw_status"=>s.operation.raw_solver_status,"accepted"=>solve_status(s).publishable) for (n,s) in r.stages]
        record["allocation_request"]=Dict("p_available_W"=>5000.,"p_rated_W"=>6000.,"active_fraction"=>.1,
            "reactive_request_var"=>four ? [-600.,0.,0.] : -300.,"negative_sequence_request_A"=>four ? nothing : [1.,0.])
        push!(cases,record)
    end
    repo=normpath(joinpath(ROOT,"..",".."));source=study_source_sha256()
    for file in ("thermal_energy.jl","four_leg.jl","voltage_aware.jl","switching_engineering.jl",
            "prepare_engineering_boundary_emt.jl","prepare_emt.jl","emt_circuit.jl")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    mkpath(dirname(output))
    writejson(output,Dict("cases"=>cases,"source_sha256"=>source,"hardware"=>HARDWARE,
        "operating_limits"=>merge(ENGINEERING_LEGACY_LIMITS,(dc_min_V=747.,dc_max_V=753.)),"fixture_limits"=>"legacy",
        "solver_options"=>ENGINEERING_OPTIONS,"julia_version"=>string(VERSION),"sizing_changed"=>false,
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Project.toml")))),
        "Manifest_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Manifest.toml")))),
        "git_revision"=>strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String)),
        "git_status"=>read(Cmd(["git","-C",repo,"status","--porcelain"]),String),
        "timing"=>Dict("switching_Hz"=>10000.,"cycles"=>10,"discard_cycles"=>6,"steps_per_carrier"=>[32,64,128]),
        "validation_bands"=>Dict("fundamental_relative"=>.01,"fundamental_absolute_A"=>.01,
            "dc_2omega_relative"=>.02,"dc_2omega_absolute_V"=>.05,"mesh_current_A"=>.01,
            "mesh_dc_2omega_V"=>.01,"mesh_peak_current_A"=>.01,"mesh_dc_extrema_V"=>.01,"energy_relative"=>1e-6),
        "scope"=>"fixed POC L-filter circuits; no dynamic feeder/service loop or all-time drift certificate"))
    cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==1 || error("usage: prepare_engineering_boundary_emt.jl output.json")
    engineering_boundary_references(ARGS[1])
end
