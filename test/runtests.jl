using Test
using PowerOptLab
using JuMP, Ipopt

include("fixtures.jl")

@testset "PowerOptLab" begin
    include("multiperiod_tests.jl")
    include("ev_tests.jl")
    include("state_estimation_tests.jl")
    include("operating_envelope_tests.jl")
end
