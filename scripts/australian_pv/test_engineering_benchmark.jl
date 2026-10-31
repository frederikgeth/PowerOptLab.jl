using Test
include("audit_engineering_benchmark.jl")

const BENCHMARK_TEST_INPUT=get(ENV,"ENGINEERING_BENCHMARK_INPUT",joinpath(ROOT,"results","engineering_benchmark"))
const BENCHMARK_TEST_PYTHON=get(ENV,"ENGINEERING_BENCHMARK_PYTHON","python3")

function mutate_engineering_record!(folder,key,mutate!)
    read_record(path)=JSON3.read(read(path,String),Dict{String,Any})
    path=joinpath(folder,"points",key*".json");point=read_record(path)
    mutate!(point);write_engineering_json(path,point)
    digest=bytes2hex(sha256(read(path)))
    details=read_record(joinpath(folder,key*"_status.json"));details["point_sha256"]=digest
    write_engineering_json(joinpath(folder,key*"_status.json"),details)
    summary=read_record(joinpath(folder,"summary.json"));summary["cases"][key]["point_sha256"]=digest
    write_engineering_json(joinpath(folder,"summary.json"),summary)
end

function rejected_engineering_fixture(mutate!)
    mktempdir() do temp
        folder=joinpath(temp,"campaign");cp(BENCHMARK_TEST_INPUT,folder;force=true)
        mutate!(folder)
        @test_throws ArgumentError audit_engineering_benchmark(folder;python=BENCHMARK_TEST_PYTHON)
        report=JSON3.read(read(joinpath(folder,"audit.json"),String),Dict{String,Any})
        @test !report["all_planned_accepted"]
    end
end

@testset "Matched engineering campaign evidence and corruption gates" begin
    report=mktempdir() do temp
        folder=joinpath(temp,"campaign");cp(BENCHMARK_TEST_INPUT,folder)
        audit_engineering_benchmark(folder;python=BENCHMARK_TEST_PYTHON)
    end
    @test report["all_planned_accepted"]
    @test report["accepted"]==report["planned"]==8
    @test report["network_conditions"]==8*826
    @test report["physical_checks"]>3500
    @test length(report["record_sha256"])==8
    for (label,mutate!) in (
        "current corruption"=>(folder->mutate_engineering_record!(folder,"nominal_local",
            p->(p["devices"][first(sort!(collect(keys(p["devices"]))))]["grid_current"][1][1]+=.02))),
        "planning relaxation"=>(folder->mutate_engineering_record!(folder,"nominal_local",
            p->(p["planning_voltage_max_V"]=400.))),
        "missing inverter"=>(folder->mutate_engineering_record!(folder,"nominal_local",
            p->delete!(p["devices"],first(keys(p["devices"]))))),
        "hardware mismatch"=>(folder->mutate_engineering_record!(folder,"nominal_local",
            p->(p["model_hardware_sha256"][first(keys(p["devices"]))]="0"^64))),
        "raw iteration limit promoted"=>(folder->mutate_engineering_record!(folder,"nominal_local",
            p->(p["solver_status"]["termination_status"]="ITERATION_LIMIT"))),
        "wrong operating condition"=>(folder->mutate_engineering_record!(folder,"nominal_local",
            p->(p["case"]["source_pu"]=1.03))),
        "failed local stage"=>(folder->mutate_engineering_record!(folder,"nominal_local",
            p->(first(values(p["allocation_stages"]))[1]["solve"]["publishable"]=false))))
        @testset "$label" begin
            rejected_engineering_fixture(mutate!)
        end
    end
    for (label,key,value) in (("temporary row survived","temporary_neutral_constraints_remaining",1),
            ("final objective changed","weight",.1))
        @testset "$label" begin
            rejected_engineering_fixture() do folder
                path=joinpath(folder,"nominal_coordinated_status.json")
                record=JSON3.read(read(path,String),Dict{String,Any})
                key=="weight" ? (record["attempts"][end][key]=value) : (record[key]=value)
                write_engineering_json(path,record)
            end
        end
    end
    @testset "stale source" begin
        rejected_engineering_fixture() do folder
            path=joinpath(folder,"environment.json")
            record=JSON3.read(read(path,String),Dict{String,Any})
            record["source_sha256"][first(keys(record["source_sha256"]))]="0"^64
            write_engineering_json(path,record)
        end
    end
    @testset "recorded dependency environment corrupted" begin
        rejected_engineering_fixture() do folder
            open(joinpath(folder,"dependency_manifest.toml"),"a") do io
                write(io,"\n# corrupted fixture\n")
            end
        end
    end
end
