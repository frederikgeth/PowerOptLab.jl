# Retain an initialization ablation with the same final hardware/constraints.
# Its raw outcome is evidence, never an accepted benchmark replacement.
include("run_engineering_benchmark.jl")

function run_engineering_numerical_control(benchmark,output;python="python3")
    isdir(output) && !isempty(readdir(output)) && error("retain previous control runs")
    read_record(file)=JSON3.read(read(file,String),Dict{String,Any})
    audit=read_record(joinpath(benchmark,"audit.json"))
    audit["all_planned_accepted"] || error("accepted parent benchmark required")
    condition=only(c for c in engineering_conditions() if c.case_id=="high_source")
    parent=joinpath(benchmark,"points","high_source_local.json")
    bytes2hex(sha256(read(parent)))==audit["record_sha256"][joinpath("points","high_source_local.json")] || error("parent point changed")
    point=read_record(parent)
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    net=voltage_aware_network(base,condition);mixed=engineering_fleet(net,condition)
    local_p=audit["cases"]["high_source_local"]["p_poc_W"]
    repo=normpath(joinpath(ROOT,"..",".."));source=study_source_sha256()
    for file in ("run_engineering_numerical_control.jl","run_engineering_benchmark.jl","engineering_records.jl",
            "coordinated_engineering.jl","switching_engineering.jl","thermal_energy.jl","four_leg.jl","voltage_aware.jl","network_si_audit.py")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    revision=strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String))
    status=read(Cmd(["git","-C",repo,"status","--porcelain"]),String)
    mkpath(output)
    write_engineering_json(joinpath(output,"environment.json"),Dict("source_sha256"=>source,
        "git_revision"=>revision,"git_status"=>status,"julia_version"=>string(VERSION),
        "parent_audit_sha256"=>bytes2hex(sha256(read(joinpath(benchmark,"audit.json")))),
        "parent_point_sha256"=>bytes2hex(sha256(read(parent))),"case"=>condition,"s_base_VA"=>VA_BASE,
        "solver_options"=>ENGINEERING_COORDINATED_OPTIONS,"initialize_neutral"=>false,
        "hardware_sha256"=>point["model_hardware_sha256"],
        "scope"=>"numerical initialization ablation; identical full-capability objective and physical/network limits; raw failures retained, never promoted by a successful SI check"))
    result=solve_engineering_coordinated(net,mixed;local_point=point,minimum_p_W=local_p,
        s_base=VA_BASE,initialize_neutral=false)
    candidate=engineering_point(result,condition;coordinated=true)
    candidate["model_hardware_sha256"]=Dict(id=>engineering_hardware_sha256(c.device) for (id,c) in mixed.spec.devices)
    write_engineering_json(joinpath(output,"point.json"),candidate)
    network_ok=engineering_network_audit(joinpath(output,"point.json"),joinpath(output,"network_audit.json");python)
    write_engineering_json(joinpath(output,"status.json"),Dict("raw"=>result.raw,"solve"=>result.solve,
        "raw_solver_status"=>result.raw_solver_status,"model_candidate_accepted"=>result.model_candidate_accepted,
        "network_accepted"=>network_ok,"attempts"=>result.attempts,"replay"=>result.replay,
        "temporary_neutral_constraints_remaining"=>result.temporary_neutral_constraints_remaining,
        "point_sha256"=>bytes2hex(sha256(read(joinpath(output,"point.json"))))))
    println("Initialization ablation raw=",result.raw.termination_status," accepted=",result.model_candidate_accepted)
    result.raw
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS) in (2,3) || error("usage: run_engineering_numerical_control.jl benchmark output [python]")
    run_engineering_numerical_control(ARGS[1],ARGS[2];python=length(ARGS)==3 ? ARGS[3] : "python3")
end
