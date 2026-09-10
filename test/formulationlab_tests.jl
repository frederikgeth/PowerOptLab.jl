using Test, PowerOptLab, Clarabel
import FormulationLab

@testset "FormulationLab migration bridge" begin
    @test PowerOptLab.L3FOptions === FormulationLab.L3FOptions
    @test PowerOptLab.build_l3f_opf === FormulationLab.build_l3f_opf
    net = Dict{String,Any}(
        "bus"=>Dict("b"=>Dict("terminal_names"=>["a"])),
        "voltage_source"=>Dict("s"=>Dict("bus"=>"b","terminal_map"=>["a"],
            "v_magnitude"=>[230.0],"v_angle"=>[0.0])),
        "load"=>Dict("l"=>Dict("bus"=>"b","terminal_map"=>["a"],
            "configuration"=>"SINGLE_PHASE","model"=>"constant_power",
            "p_nom"=>[100.0],"q_nom"=>[20.0])))
    r=PowerOptLab.solve_l3f_opf(net,Clarabel.Optimizer;solver_options=(verbose=false,))
    @test PowerOptLab.solve_status(r).optimal
    @test PowerOptLab.solve_diagnostics(r).objective == r.objective
    @test r.sources["s"]["pg"] ≈ [100.0]
    reference=PowerOptLab.l3f_reference_from_powerflow(net;solver_options=(print_level=0,))
    @test reference.provenance == :power_flow
    replay=PowerOptLab.validate_l3f_solution(r)
    @test replay["status"] == "replayed"
end
