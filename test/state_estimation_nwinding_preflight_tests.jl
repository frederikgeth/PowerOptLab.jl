using Test, PowerOptLab

@testset "SE n-winding preflight guards" begin
    net = Dict{String,Any}(
        "bus"=>Dict("h"=>Dict("terminal_names"=>["a","n"],"perfectly_grounded_terminals"=>["n"]),
                    "l"=>Dict("terminal_names"=>["a","n"],"perfectly_grounded_terminals"=>["n"])),
        "voltage_source"=>Dict("s"=>Dict("bus"=>"h","terminal_map"=>["a"],"v_magnitude"=>[2400.],"v_angle"=>[0.])),
        "transformer"=>Dict("n_winding"=>Dict("t"=>Dict{String,Any}(
            "windings"=>[Dict{String,Any}("bus"=>"h","terminal_map"=>["a","n"],"v_nom"=>2400.,"r_winding"=>1.),
                         Dict{String,Any}("bus"=>"l","terminal_map"=>["a","n"],"v_nom"=>1200.,"r_winding"=>.25)],
            "x_sc"=>Dict("1_2"=>2.)))))
    @test state_estimator_preflight(net).supported
    for mutate in (
        t -> (t["windings"][1]["v_nom"] = 0.),
        t -> (t["windings"][2]["configuration"] = "INVALID"),
        t -> (t["windings"][2]["terminal_map"] = ["n"]),
        t -> (t["tap"] = 1.03),
        t -> push!(t["windings"],Dict{String,Any}()))
        bad=deepcopy(net); mutate(bad["transformer"]["n_winding"]["t"])
        @test !state_estimator_preflight(bad).supported
        @test_throws SEUnsupportedNetwork compile_state_estimator(bad)
    end
    net["transformer"]["n_winding"]["t"]["windings"][1]["i_max"]=100.
    @test any(o->o.field=="i_max",state_estimator_preflight(net).omitted)
end
