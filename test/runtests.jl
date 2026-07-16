using Test
using TOML
using Aqua
using PowerOptLab
using JuMP, Ipopt

# JuMP + Ipopt are hard dependencies here, so they are always available. The
# HELM OpenDSS-parity tests optionally use OpenDSSDirect and skip themselves
# when it is not installed (mirrors the guard in BMOPFTools' own suite).
const _HAS_JUMP_IPOPT = true
const _HAS_ODS = !isnothing(Base.identify_package("OpenDSSDirect"))
if _HAS_ODS
    @eval using OpenDSSDirect
end

include("fixtures.jl")

@testset "Package quality" begin
    Aqua.test_all(PowerOptLab)
end

@testset "PowerOptLab" begin
    include("multiperiod_tests.jl")
    include("ev_tests.jl")
    include("state_estimation_tests.jl")
    include("constrained_state_estimation_tests.jl")
    include("parameter_estimation_tests.jl")
    include("inverse_carson_tests.jl")
    include("inverse_carson_benchmark_tests.jl")
    include("operating_envelope_tests.jl")
    include("advanced_inverter_tests.jl")
    include("battery_tests.jl")
    include("helm_tests.jl")
end
