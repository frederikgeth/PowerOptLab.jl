# Additional checks for the deliberately modified customer connections.
include("run.jl")
using Test, JuMP
function replay(net, result)
    fixed=deepcopy(net); empty!(fixed["ibr"])
    for (id,device) in result.devices
        bus=result.spec.devices[id].device.bus
        ts=result.network["bus"][bus]; vn=phasor(ts["n"])
        current=device.grid_phase_current
        @test abs(sum(current)) < 1e-6 # three-wire converter cannot inject I0
        for (i,phase) in enumerate(["a","b","c"])
            s=(phasor(ts[phase])-vn)*conj(current[i])
            fixed["load"]["replay_"*id*phase]=Dict("bus"=>bus,
                "terminal_map"=>[phase,"n"],"configuration"=>"SINGLE_PHASE",
                "model"=>"constant_power","v_nom"=>[230.],
                "p_nom"=>[-real(s)],"q_nom"=>[-imag(s)])
        end
    end
    r=solve_pf(fixed;per_unit=true,solver_options=OPTIONS)
    @test r["termination_status"]=="LOCALLY_SOLVED"
    delta=maximum(abs(phasor(t)-phasor(r["bus"][b][p]))
        for (b,ts) in result.network["bus"] for (p,t) in ts)
    @test delta < 0.01
    delta
end
function audit()
    provenance=JSON3.read(read(joinpath(ROOT,"source","provenance.json"),String))
    for f in provenance.files
        @test bytes2hex(sha256(read(joinpath(ROOT,"source",f.copy)))) == f.sha256
    end
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    n=study_network(base,1.,1.)
    for source in values(n["voltage_source"])
        v=source["v_magnitude"] .* cis.(source["v_angle"])
        a=cis(2pi/3)
        u1=(v[1]+a*v[2]+a^2*v[3])/3
        u2=(v[1]+a^2*v[2]+a*v[3])/3
        @test abs(u2)/abs(u1) ≈ 0.01
        @test abs(sum(v[1:3])) < 1e-8
    end
    @test length(n["load"])==36
    @test length(n["ibr"])==12
    @test sum(sum(ld["p_nom"]) for ld in values(n["load"]))≈18000.
    @test length(base["load"])==12 # input was not mutated
    write_bmopf(n,joinpath(ROOT,"networks","LV3_three_phase_pv.bmopf.json"))
    empty!(n["ibr"])
    write_bmopf(n,joinpath(ROOT,"networks","LV3_three_phase_no_pv.bmopf.json"))
    pf,report=validate_pf(parse_bmopf(joinpath(ROOT,"networks","LV3_three_phase_no_pv.bmopf.json")))
    report["customer_count"]=12; report["phase_load_records"]=36
    for mode in ["unity","worst_vv_vw","sequence_droop"]
        net=study_network(base,1.,1.03)
        result=solve_controlled_inverter_fleet(net,fleet(net,mode,1.);
            s_base=100000.,solver_options=OPTIONS)
        @test solve_status(result).publishable
        report[mode*"_fixed_PQ_replay_max_difference_V"]=replay(net,result)
    end
    writejson(joinpath(ROOT,"results","modified_validation.json"),report)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    @testset "Australian three-phase study" begin
        audit()
    end
end
