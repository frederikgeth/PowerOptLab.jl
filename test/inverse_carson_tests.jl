const IC_RAC = 4.5e-4
const IC_GMR = 4.0e-3
const IC_RAD = 5.0e-3

function _ic_test_candidate(id, geometry)
    n = geometry in (:horizontal_3, :triangle_3) ? 3 : 4
    bounds = if geometry in (:horizontal_3, :triangle_3)
        ([0.3, 6.0, 20.0], [0.8, 12.0, 90.0])
    elseif geometry == :horizontal_4
        ([0.15, 0.2, 6.0, 20.0], [0.5, 0.8, 12.0, 90.0])
    else
        ([0.3, 7.0, 0.3, 20.0], [0.8, 12.0, 2.5, 90.0])
    end
    OverheadCarsonCandidate(id=id, geometry=geometry,
        r_ac_ref=fill(IC_RAC, n), gmr=fill(IC_GMR, n),
        radius=fill(IC_RAD, n), lower=bounds[1], upper=bounds[2])
end

function _ic_observation_from(c, parameters; shunt=false, sigma_scale=1.0)
    seed = SequenceLineObservation(z0=1 + im, z1=1 + im,
        b0=shunt ? 1.0 : nothing, b1=shunt ? 1.0 : nothing,
        frequency=50.0, z_units=:ohm_per_m, b_units=:siemens_per_m,
        sigma=ones(shunt ? 6 : 4))
    u = (parameters .- c.lower) ./ (c.upper .- c.lower)
    y = PowerOptLab._ic_predict(c, seed, u)
    sigma = shunt ? [fill(1e-7 * sigma_scale, 4); fill(1e-10 * sigma_scale, 2)] :
            fill(1e-7 * sigma_scale, 4)
    SequenceLineObservation(z0=y[1] + im * y[2], z1=y[3] + im * y[4],
        b0=shunt ? y[5] : nothing, b1=shunt ? y[6] : nothing,
        frequency=50.0, z_units=:ohm_per_m, b_units=:siemens_per_m,
        sigma=sigma)
end

@testset "Inverse Carson" begin
    @testset "shunt-informed round trip and BMOPF materialization" begin
        c = _ic_test_candidate("horizontal", :horizontal_3)
        truth = [0.55, 9.15, 75.0]
        obs = _ic_observation_from(c, truth; shunt=true)
        result = solve_inverse_carson(obs, [c]; starts=6,
                                      solver_options=(tol=1e-9,))
        fit = only(result.fits)

        @test result.compatible_candidates == [c.id]
        @test fit.compatible
        @test fit.parameters ≈ truth atol=1e-4
        @test fit.jacobian_rank == 3
        @test fit.objective < 1e-8
        @test fit.Z_primitive == transpose(fit.Z_primitive)
        @test fit.C_primitive == transpose(fit.C_primitive)
        @test size(fit.Z_sequence) == (3, 3)
        @test size(something(fit.B_sequence)) == (3, 3)
        @test maximum(abs(fit.Z_sequence[i, j]) for i in 1:3 for j in 1:3
                      if i != j) > 0

        materialized = materialize_inverse_carson(fit, c)
        PowerOptLab.BMOPFTools.compile_linecode(materialized, c.id)
        lc = materialized["linecode"][c.id]
        Z = PowerOptLab.BMOPFTools._pattern_keys_to_matrix(lc, "R_series_") .+
            im .* PowerOptLab.BMOPFTools._pattern_keys_to_matrix(lc, "X_series_")
        @test Z ≈ fit.Z_primitive rtol=1e-12
    end

    @testset "series-only ambiguity and conductor-count rejection" begin
        horizontal = _ic_test_candidate("horizontal", :horizontal_3)
        triangle = _ic_test_candidate("triangle", :triangle_3)
        four_wire = _ic_test_candidate("four_wire", :horizontal_4)
        obs = _ic_observation_from(horizontal, [0.55, 9.15, 75.0])
        result = solve_inverse_carson(obs, [horizontal, triangle, four_wire];
                                      starts=6)

        @test Set(result.compatible_candidates) == Set(["horizontal", "triangle"])
        @test !result.fits[findfirst(f -> f.candidate_id == "four_wire",
                                    result.fits)].compatible
        @test all(f.jacobian_rank == 2 for f in result.fits if f.compatible)
        @test any(w -> contains(w, "ambiguous"), result.warnings)
        @test any(w -> contains(w, "rank-deficient"), result.warnings)
    end

    @testset "smooth derivative and no-match outcome" begin
        c = _ic_test_candidate("horizontal", :horizontal_3)
        obs = _ic_observation_from(c, [0.55, 9.15, 75.0])
        u = [0.4, 0.6, 0.7]
        f(v) = PowerOptLab._ic_objective(c, obs, v)
        gad = PowerOptLab.ForwardDiff.gradient(f, u)
        h = 1e-6
        gfd = [(f(u + h * (1:3 .== i)) - f(u - h * (1:3 .== i))) / (2h)
               for i in 1:3]
        @test gad ≈ gfd rtol=2e-5 atol=1e-5
        @test all(isfinite, gad)

        impossible = SequenceLineObservation(
            z0=0.2 * real(obs.z1) + im * imag(obs.z0), z1=obs.z1,
            frequency=obs.frequency, z_units=:ohm_per_m,
            sigma=fill(1e-7, 4))
        result = solve_inverse_carson(impossible, [c]; starts=4)
        @test isempty(result.compatible_candidates)
        @test any(w -> contains(w, "no candidate"), result.warnings)
    end

    @testset "input contracts" begin
        @test_throws ArgumentError SequenceLineObservation(
            z0=1 + im, z1=1 + im, b0=1.0, frequency=50.0)
        @test_throws DimensionMismatch SequenceLineObservation(
            z0=1 + im, z1=1 + im, frequency=50.0, sigma=ones(3))
        @test_throws ArgumentError OverheadCarsonCandidate(
            id="bad", geometry=:horizontal_3, r_ac_ref=fill(IC_RAC, 3),
            gmr=fill(IC_GMR, 3), radius=fill(IC_RAD, 3),
            lower=[0.8, 6.0, 20.0], upper=[0.3, 12.0, 90.0])
        @test_throws ArgumentError OverheadCarsonCandidate(
            id="bad-cap-radius", geometry=:horizontal_3,
            r_ac_ref=fill(IC_RAC, 3), gmr=fill(IC_GMR, 3),
            radius=fill(IC_RAD, 3), cap_radius=[IC_RAD, 0.0, IC_RAD],
            lower=[0.3, 6.0, 20.0], upper=[0.8, 12.0, 90.0])
        @test_throws ArgumentError OverheadCarsonCandidate(
            id="bad-temperature", geometry=:horizontal_3,
            r_ac_ref=fill(IC_RAC, 3), gmr=fill(IC_GMR, 3),
            radius=fill(IC_RAD, 3), alpha_20=0.01,
            lower=[0.3, 6.0, -100.0], upper=[0.8, 12.0, 90.0])
    end
end
