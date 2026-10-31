# Saved-point tests keep publication false after any failed or incomplete audit.
include("audit_coordinated.jl")

@testset "Coordinated saved-record acceptance" begin
    output=isempty(ARGS) ? joinpath(ROOT,"results","coordinated") : ARGS[1]
    function fixture(f)
        mktempdir() do temp
            path=joinpath(temp,"records");cp(output,path)
            f(path)
        end
    end
    fixture() do path
        r=audit_coordinated(path)
        @test r["all_planned_accepted"] && r["accepted_saved_points"]==24
        @test r["physical_checks_recomputed"]>0 && r["assertions"]>1000
    end
    fixture() do path
        env=JSON3.read(read(joinpath(path,"environment.json"),String),Dict{String,Any})
        id=first(env["planned_cases"])["case_id"]
        file=joinpath(path,"points",id*".json")
        point=JSON3.read(read(file,String),Dict{String,Any})
        first(values(point["devices"]))["converter_current"][1][1]+=100.
        writejson(file,point)
        @test_throws ArgumentError audit_coordinated(path)
        report=JSON3.read(read(joinpath(path,"audit.json"),String))
        @test !report.all_planned_accepted && report.reason=="AUDIT_FAILED"
    end
    fixture() do path
        file=joinpath(path,"summary.json")
        summary=JSON3.read(read(file,String),Dict{String,Any})
        id=first(keys(summary["cases"]))
        delete!(summary["cases"],id)
        summary["completed"]-=1;summary["accepted"]-=1;summary["all_planned_accepted"]=false
        writejson(file,summary)
        r=audit_coordinated(path)
        @test !r["all_planned_accepted"] && r["accepted_saved_points"]==23
        @test r["cases"][id]["reason"]=="MISSING"
    end
    fixture() do path
        file=joinpath(path,"environment.json")
        env=JSON3.read(read(file,String),Dict{String,Any})
        env["source_sha256"][first(keys(env["source_sha256"]))]="invalid"
        writejson(file,env)
        @test_throws ArgumentError audit_coordinated(path)
        @test !JSON3.read(read(joinpath(path,"audit.json"),String)).all_planned_accepted
    end
end
