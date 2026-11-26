# The explicit-neutral Kron reduction moved to BMOPFTools (kron_reduce_neutrals,
# tested there). PowerOptLab re-exports it and keeps deprecated aliases.
import BMOPFTools

@testset "Kron reduction: BMOPFTools re-export and deprecated aliases" begin
    @test PowerOptLab.kron_reduce_neutrals === BMOPFTools.kron_reduce_neutrals
    net = Dict{String,Any}(
        "bus" => Dict("b" => Dict{String,Any}("terminal_names" => ["a", "n"],
                                              "perfectly_grounded_terminals" => ["n"])),
        "voltage_source" => Dict("s" => Dict{String,Any}("bus" => "b", "terminal_map" => ["a", "n"],
                                                       "v_magnitude" => [230.0], "v_angle" => [0.0])))
    expected = kron_reduce_neutrals(net)
    @test expected["bus"]["b"]["terminal_names"] == ["a"]
    for alias in (kron_reduce_bmopf, reduce_bmopf_neutrals)
        out = @test_deprecated alias(net)
        @test out["bus"]["b"] == expected["bus"]["b"]
        @test out["voltage_source"] == expected["voltage_source"]
    end
end
