using Test
using TOML
using Aqua
using PowerOptLab
using JuMP, Ipopt
import PiecewiseLinearOpt

# JuMP + Ipopt are hard dependencies here, so they are always available. The
# HELM OpenDSS-parity tests optionally use OpenDSSDirect and skip themselves
# when it is not installed (mirrors the guard in BMOPFTools' own suite).
const _HAS_JUMP_IPOPT = true

include("fixtures.jl")

@testset "Package quality" begin
    # On Julia 1.10, Aqua's persistent-task probe develops the package in a fresh
    # environment that cannot recover the source of the unregistered BMOPFTools
    # dependency. Keep that probe on newer Julia releases and retain every other
    # Aqua check on all supported versions.
    Aqua.test_all(PowerOptLab; persistent_tasks = VERSION >= v"1.11")
end

# OpenDSSDirect is a declared test dependency. Load it unconditionally for the
# oracle; unlike PowerOptLab's package runtime, the test target owns this dep.
using OpenDSSDirect
const _HAS_ODS = true

@testset "PowerOptLab" begin
    include("function_formulation_tests.jl")
    include("function_comparison_tests.jl")
    include("formulation_research_tests.jl")
    include("control_lowering_tests.jl")
    include("relation_lowering_tests.jl")
    include("formulation_observation_tests.jl")
    include("formulation_review_tests.jl")
    include("selector_primitive_tests.jl")
    include("multiperiod_tests.jl")
    include("ev_tests.jl")
    include("evse_tests.jl")
    include("ev_workplace_tests.jl")
    include("state_estimation_tests.jl")
    include("constrained_state_estimation_tests.jl")
    include("state_estimation_network_tests.jl")
    include("state_estimation_nwinding_preflight_tests.jl")
    include("parameter_estimation_tests.jl")
    include("inverse_carson_tests.jl")
    include("inverse_carson_benchmark_tests.jl")
    include("operating_envelope_tests.jl")
    include("doe_evidence_tests.jl")
    include("doe_cleanup_tests.jl")
    include("bilevel_tests.jl")
    include("advanced_inverter_tests.jl")
    include("generalized_generator_tests.jl")
    include("generator_data_tests.jl")
    include("inverter_control_tests.jl")
    include("inverter_control_numerics_tests.jl")
    include("inverter_control_study_tests.jl")
    include("inverter_control_experiment_tests.jl")
    include("inverter_control_sizing_tests.jl")
    include("closed_loop_evidence_tests.jl")
    include("operability_tests.jl")
    include("battery_tests.jl")
    include("helm_tests.jl")
    include("kron_reduction_tests.jl")
    include("formulationlab_tests.jl")
end
