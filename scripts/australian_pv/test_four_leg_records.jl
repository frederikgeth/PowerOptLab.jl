include("audit_four_leg.jl")

@testset "Four-leg saved evidence rejects corruption" begin
    source=isempty(ARGS) ? joinpath(ROOT,"results","four_leg") : ARGS[1]
    mktempdir() do tmp
        copy=joinpath(tmp,"evidence");cp(source,copy)
        path=joinpath(copy,"points","matched_s1_p1.json")
        original=read(path,String);point=JSON3.read(original,Dict{String,Any})
        first(values(point["devices"]))["converter_current"][1][1]+=20.
        writejson(path,point)
        @test_throws ArgumentError audit_four_leg(copy)
        @test !JSON3.read(read(joinpath(copy,"audit.json"),String))["all_planned_accepted"]
        write(path,original)
        point=JSON3.read(original,Dict{String,Any});delete!(point["devices"],first(keys(point["devices"])))
        writejson(path,point)
        @test_throws ArgumentError audit_four_leg(copy)
        write(path,original)
        envpath=joinpath(copy,"environment.json")
        original_env=read(envpath,String);env=JSON3.read(original_env,Dict{String,Any})
        env["source_sha256"][first(keys(env["source_sha256"]))]="stale"
        writejson(envpath,env)
        @test_throws ArgumentError audit_four_leg(copy)
        write(envpath,original_env)
        summarypath=joinpath(copy,"summary.json");summary=JSON3.read(read(summarypath,String),Dict{String,Any})
        summary["cases"]["matched_s1_p1"]["p_poc_W"]+=10.
        writejson(summarypath,summary)
        @test_throws ArgumentError audit_four_leg(copy)
    end
end
