include(joinpath(@__DIR__,"..","scripts","formulations","comparison.jl"))
@testset "Analytic electrical PWL comparison" begin
    r = PWLComparison.representations()
    @test PWLComparison.reference_equilibrium().current ≈ 40/3
    for vs in (210.,220.,230.,250.,260.), (si,so) in ((1.,1.),(230.,20.)),
        start in (230.,260.), rep in (r.softplus,r.local_c2)
        row = PWLComparison.run_case(rep;source_voltage=vs,
            voltage_scale=si,current_scale=so,start_voltage=start)
        @test row["strict_solver_success"]
        @test row["candidate_available"]
        if row["candidate_available"]
            @test row["electrical_residual_V"] <= 1e-6
            @test row["domain_violation_V"] <= 1e-6
            @test row["current_limit_violation_A"] <= 1e-6
            @test abs(row["surrogate_equation_error_A"]) <= 1e-6
            @test row["error_lower_A"]-1e-6 <= row["exact_graph_error_A"] <= row["error_upper_A"]+1e-6
            # For this monotone nonincreasing scalar law and R>=0 the equilibrium
            # current error is bounded by the controller's uniform output error.
            @test abs(row["reference_current_error_A"]) <= PWLComparison.ERROR_BUDGET_A+1e-6
        end
    end
    row = PWLComparison.run_case(r.hull)
    @test row["semantics"] == "outer_relaxation"
    @test !row["canonical_equations_satisfied"]
    @test row["current_A"] ≈ 16. atol=1e-6
    @test row["exact_graph_error_A"] ≈ 8. atol=1e-6
end
