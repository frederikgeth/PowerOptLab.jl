include("audit_thermal_energy.jl")
@testset "Mixed thermal/energy evidence corruption gates" begin
    source=isempty(ARGS) ? joinpath(ROOT,"results","thermal_energy") : ARGS[1]
    mktempdir() do tmp
        folder=joinpath(tmp,"records");cp(source,folder)
        path=joinpath(folder,"points","energy_s2.json");original=read(path,String)
        point=JSON3.read(original,Dict{String,Any})
        first(values(point["devices"]))["converter_current"][1][1]+=20.
        writejson(path,point)
        @test_throws ArgumentError audit_mixed_thermal(folder)
        @test !JSON3.read(read(joinpath(folder,"audit.json"),String))["all_planned_accepted"]
        write(path,original)
        point=JSON3.read(original,Dict{String,Any});delete!(point["devices"],first(keys(point["devices"])))
        writejson(path,point)
        @test_throws ArgumentError audit_mixed_thermal(folder)
        write(path,original)
        path=joinpath(folder,"environment.json");original=read(path,String)
        env=JSON3.read(original,Dict{String,Any});env["source_sha256"][first(keys(env["source_sha256"]))]="stale"
        writejson(path,env)
        @test_throws ArgumentError audit_mixed_thermal(folder)
        write(path,original)
        path=joinpath(folder,"summary.json");summary=JSON3.read(read(path,String),Dict{String,Any})
        summary["cases"]["energy_s2"]["p_poc_W"]+=1.
        writejson(path,summary)
        @test_throws ArgumentError audit_mixed_thermal(folder)
    end
end
