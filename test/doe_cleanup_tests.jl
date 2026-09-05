@testset "DOE scrambled Halton and bounded dropout generation" begin
    inverse = PowerOptLab._doe_radical_inverse
    # Complement every binary digit, including the infinite trailing zeros:
    # phi_complement(n) = 1 - phi(n). The truncated implementation collided.
    expected = [1//2, 3//4, 1//4, 7//8, 3//8, 5//8, 1//8, 15//16]
    @test [inverse(n, 2, [1, 0]) for n in 1:8] ≈ expected
    @test length(unique(inverse(n, 2, [1, 0]) for n in 1:1024)) == 1024
    @test all(inverse(n, 3, [2, 1, 0]) ≈ 1 - inverse(n, 3) for n in 1:100)
    @test all(inverse(n, 5, collect(0:4)) == inverse(n, 5) for n in 1:100)
    diagnostics = PowerOptLab._doe_search_point_diagnostics
    @test !diagnostics(20, 2; dropout_depth=0,
        halton_scramble_seed=7)["halton_seed_stratified"]
    @test !diagnostics(2, 0; dropout_depth=0,
        halton_scramble_seed=7)["halton_seed_stratified"]
    @test diagnostics(1, 2; dropout_depth=0, sequence_offset=2,
        halton_scramble_seed=nothing)["halton_occupied_first_digit_strata"] == [2]
    @test diagnostics(1, 2; dropout_depth=0, sequence_offset=0,
        halton_scramble_seed=7)["halton_occupied_first_digit_strata"] == [1]
    @test PowerOptLab._doe_dropout_count(64, 64) == big(2)^64 - 1
    @test_throws ArgumentError PowerOptLab._doe_search_points(64, 0, 0;
        include_zero=false, include_bound=false, include_corners=false,
        max_exact_corners=10, dropout_depth=64, max_dropout_points=10)
end
