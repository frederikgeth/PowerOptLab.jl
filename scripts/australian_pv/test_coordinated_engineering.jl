using Test
include("coordinated_engineering.jl")

function engineering_solver_fixture(available)
    fixture=Dict{String,Any}("ibr"=>Dict{String,Any}("pv"=>Dict{String,Any}("bus"=>"poc")))
    original=four_leg_fleet(fixture,1.)
    controlled=ControlledDevice(engineering_plant(original.devices["pv"].device),original.devices["pv"].controller)
    a=cis(2pi/3)
    v=250 .* ComplexF64[1,a^2,a] .+ 2.5cis(.3) .* ComplexF64[1,a,a^2]
    net=PowerOptLab._projection_network(controlled.device,v)
    net["load"]=Dict{String,Any}()
    net["ibr"]=Dict{String,Any}("pv"=>Dict{String,Any}("bus"=>"poc",
        "topology"=>"FOUR_LEG","terminal_map"=>["a","b","c","n"]))
    spec=ControlledInverterFleetSpec(Dict("pv"=>controlled),
        Dict("pv"=>InverterControlRequest(p_available=available,p_rated=6000.,q_scale=6000.)))
    net,(spec=spec,modes=Dict("pv"=>:energized),labels=Dict("pv"=>"phase_vv_vw"))
end

@testset "Coordinated engineering: complete neutral release and coordinate bases" begin
    for available in (1250.,5000.)
        reference=nothing
        for base in (1e4,1e5)
            net,mixed=engineering_solver_fixture(available)
            r=solve_engineering_coordinated(net,mixed;s_base=base)
            @test r.raw.publishable
            @test r.model_candidate_accepted
            @test r.solve.publishable==r.model_candidate_accepted
            @test r.replay.accepted
            @test r.temporary_neutral_constraints_remaining==0
            @test first(r.attempts).stage===:neutral_initialization
            @test last(r.attempts).stage===:full_capability
            @test last(r.attempts).weight==0.
            @test all(a.raw.publishable for a in r.attempts)
            @test r.devices["pv"].plant.physical_feasibility.passed
            @test r.devices["pv"].plant.p_dc<=available+1e-4
            @test r.devices["pv"].plant.p_poc>1000.
            @test JuMP.objective_sense(PowerOptLab._opf_model(r.model_context))===JuMP.MOI.MAX_SENSE
            reference===nothing || @test r.devices["pv"].plant.p_poc≈reference atol=.005
            reference=r.devices["pv"].plant.p_poc
            bad=SolveStatus("ITERATION_LIMIT","INFEASIBLE_POINT",true,false,false,false)
            saved=engineering_point(merge(r,(raw=bad,)),(case_id="rejected_fixture",);coordinated=true)
            @test !saved["solver_status"].publishable
            @test all(isfinite(t["vr"])&&isfinite(t["vi"]) for ts in values(saved["bus"]) for t in values(ts))
            @test all(all(isfinite,pair) for d in values(saved["devices"]) for pair in d["grid_current"])
        end
    end
end

@testset "Coordinated engineering: objective and rejected-record gates" begin
    net,mixed=engineering_solver_fixture(1250.)
    @test_throws ErrorException solve_engineering_coordinated(net,mixed;continuation_weights=(1.,))
    @test_throws ErrorException solve_engineering_coordinated(net,mixed;objective_scale=-1.)
    @test_throws ErrorException solve_engineering_coordinated(net,mixed;solver_options=(obj_scaling_factor=1.,))
    check=InverterLimitCheck("rejected",:A,NaN,-Inf,Inf,1e-5,-Inf,false)
    report=InverterFeasibilityReport(true,false,[check],(lower_V=-Inf,upper_V=Inf))
    safe=JSON3.read(JSON3.write(_engineering_json_safe(report)))
    @test safe.passed===false
    @test safe.checks[1].value===nothing
    @test safe.checks[1].lower===nothing
    @test safe.rail.lower_V===nothing
    inv=mixed.spec.devices["pv"].device;a=cis(2pi/3)
    voltage=250 .* ComplexF64[1,a^2,a] .+ 2.5cis(.3) .* ComplexF64[1,a,a^2]
    seed=engineering_night_seed(inv,voltage)
    current=ComplexF64[complex(z...) for z in seed["grid_current"]]
    power=voltage.*conj.(current)
    @test maximum(abs.(real.(power).+inv.p_loss_fixed/3))<1e-10
    @test abs(sum(current))<1e-12
end
