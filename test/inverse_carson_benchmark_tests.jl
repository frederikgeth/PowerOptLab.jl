const IC_BENCHMARK_DIR = joinpath(@__DIR__, "data", "inverse_carson")

function _ic_case_candidate(case; id=case["id"], conductor=nothing)
    construction = case["construction"]
    geometry = Symbol(construction["geometry"])
    n = construction["conductor_count"]
    parameters = Float64.(construction["parameters"])
    wire = conductor === nothing ? construction : conductor
    lower, upper = if geometry in (:horizontal_3, :triangle_3)
        ([0.2, 5.8, 0.0], [1.5, 21.5, 105.0])
    elseif geometry == :horizontal_4
        ([0.2, 0.2, 5.8, 0.0], [1.5, 1.5, 21.5, 105.0])
    else
        ([0.2, 5.8, 0.2, 0.0], [1.5, 21.5, 4.0, 105.0])
    end
    OverheadCarsonCandidate(
        id=id, geometry=geometry,
        r_ac_ref=fill(wire["r_ac_ref_ohm_per_m"], n),
        gmr=fill(wire["gmr_m"], n), radius=fill(wire["radius_m"], n),
        temperature_ref=wire["temperature_ref_degC"],
        alpha_20=wire["alpha_20_per_degC"],
        lower=lower, upper=upper, initial=parameters,
        angle=get(construction, "angle_deg", 21.67) * pi / 180)
end

function _ic_case_prediction(case, candidate)
    assumptions = case["assumptions"]
    seed = SequenceLineObservation(z0=1 + im, z1=1 + im,
        frequency=assumptions["frequency_hz"],
        earth_resistivity=assumptions["earth_resistivity_ohm_m"],
        z_units=:ohm_per_m, sigma=ones(4))
    parameters = Float64.(case["construction"]["parameters"])
    u = (parameters .- candidate.lower) ./ (candidate.upper .- candidate.lower)
    PowerOptLab._ic_predict(candidate, seed, u) .* 1000
end

@testset "Inverse Carson validation benchmarks" begin
    paper_data = TOML.parsefile(joinpath(IC_BENCHMARK_DIR,
                                         "paper_table_iv.toml"))
    paper_cases = paper_data["benchmark"]

    @testset "paper Table IV forward reproduction and provenance" begin
        @test length(paper_cases) == 5
        for case in paper_cases
            assumptions = case["assumptions"]
            uncertainty = case["uncertainty"]
            @test haskey(assumptions, "frequency_hz")
            @test haskey(assumptions, "earth_resistivity_ohm_m")
            @test haskey(assumptions, "transposition")
            @test haskey(assumptions, "neutral_treatment")
            @test haskey(case["observation"], "units")
            @test haskey(uncertainty, "measurement_uncertainty_status")

            candidate = _ic_case_candidate(case)
            predicted = _ic_case_prediction(case, candidate)
            reported = Float64.(case["observation"]["values"])
            # The paper prints sequence values to 4 decimals and some geometry
            # inputs to 2 decimals; two last-place units cover both roundings.
            @test predicted ≈ reported atol=2e-4 rtol=0
        end
    end

    @testset "rounded, noisy, wrong-catalog, and wrong-earth cases" begin
        case = only(filter(c -> contains(c["id"], "triangle_21_67"),
                           paper_cases))
        reported = Float64.(case["observation"]["values"])
        sigma = fill(2e-4, 4)
        correlation = [1.0 0.35 -0.15 0.0;
                       0.35 1.0 0.10 -0.10;
                      -0.15 0.10 1.0 0.25;
                       0.0 -0.10 0.25 1.0]
        covariance = sigma .* correlation .* transpose(sigma)
        noise = [0.35, -0.55, 0.25, -0.20] .* sigma
        noisy = reported .+ noise
        obs = SequenceLineObservation(
            z0=noisy[1] + im * noisy[2], z1=noisy[3] + im * noisy[4],
            frequency=case["assumptions"]["frequency_hz"],
            earth_resistivity=case["assumptions"]["earth_resistivity_ohm_m"],
            covariance=covariance)

        correct = _ic_case_candidate(case; id="mars")
        libra = Dict(
            "r_ac_ref_ohm_per_m" => 28.3e-9 / 49.48e-6,
            "gmr_m" => 2.18 * 1.5e-3,
            "radius_m" => 3.0 * 1.5e-3,
            "temperature_ref_degC" => 20.0,
            "alpha_20_per_degC" => 0.00403)
        wrong_catalog = _ic_case_candidate(case; id="libra", conductor=libra)
        result = solve_inverse_carson(obs, [correct, wrong_catalog]; starts=8)
        @test "mars" in result.compatible_candidates
        @test !("libra" in result.compatible_candidates)

        wrong_earth_obs = SequenceLineObservation(
            z0=reported[1] + im * reported[2],
            z1=reported[3] + im * reported[4], frequency=50.0,
            earth_resistivity=30.0, sigma=sigma)
        wrong_earth = solve_inverse_carson(wrong_earth_obs, [correct]; starts=6)
        @test isempty(wrong_earth.compatible_candidates)
    end

    @testset "shunt uncertainty widens height confidence" begin
        case = only(filter(c -> contains(c["id"], "horizontal_3wire"),
                           paper_cases))
        candidate = _ic_case_candidate(case; id="shunt-height")
        truth = Float64.(case["construction"]["parameters"])
        seed = SequenceLineObservation(z0=1 + im, z1=1 + im,
            b0=1.0, b1=1.0, frequency=50.0,
            z_units=:ohm_per_m, b_units=:siemens_per_m, sigma=ones(6))
        u = (truth .- candidate.lower) ./ (candidate.upper .- candidate.lower)
        y = PowerOptLab._ic_predict(candidate, seed, u)
        zsigma = fill(1e-7, 4)
        tight_sigma = [zsigma; 1e-10; 1e-10]
        loose_sigma = [zsigma; 1e-8; 1e-8]
        make_obs(s) = SequenceLineObservation(
            z0=y[1] + im * y[2], z1=y[3] + im * y[4],
            b0=y[5], b1=y[6], frequency=50.0,
            z_units=:ohm_per_m, b_units=:siemens_per_m,
            covariance=PowerOptLab.Diagonal(s .^ 2))
        tight = only(solve_inverse_carson(make_obs(tight_sigma), [candidate];
                                          starts=5).fits)
        loose = only(solve_inverse_carson(make_obs(loose_sigma), [candidate];
                                          starts=5).fits)
        tight_ci = something(tight.local_confidence_intervals)[2, :]
        loose_ci = something(loose.local_confidence_intervals)[2, :]
        @test loose_ci[2] - loose_ci[1] > tight_ci[2] - tight_ci[1]
    end

    @testset "independent OpenDSS artifact and model mismatch" begin
        dss_data = TOML.parsefile(joinpath(IC_BENCHMARK_DIR,
                                           "opendss_mars_triangle.toml"))
        case = only(dss_data["benchmark"])
        @test case["assumptions"]["earth_model"] == "deri"
        @test haskey(case["uncertainty"], "measurement_uncertainty_status")

        if _HAS_ODS
            OpenDSSDirect.Text.Command("Redirect " *
                joinpath(IC_BENCHMARK_DIR, case["provenance"]["input_file"]))
            OpenDSSDirect.LineGeometries.Name("mars_triangle")
            OpenDSSDirect.LineGeometries.RhoEarth(100.0)
            units = OpenDSSDirect.Lib.LineUnits_km
            Z = reshape(OpenDSSDirect.LineGeometries.Zmatrix(50.0, 1.0, units), 3, 3)
            C = reshape(OpenDSSDirect.LineGeometries.Cmatrix(50.0, 1.0, units), 3, 3)
            z_expected = ComplexF64[complex(v[1], v[2])
                                    for v in case["primitive"]["z_values"]]
            c_expected = Float64.(case["primitive"]["c_values"])
            @test vec(Z) ≈ z_expected atol=1e-10 rtol=0
            @test vec(C) ≈ c_expected atol=1e-10 rtol=0
            a = cis(2pi / 3)
            A = [1.0 1.0 1.0; 1.0 a^2 a; 1.0 a a^2]
            Z012 = inv(A) * Z * A
            B012 = inv(A) * (2pi * 50 .* C .* 1e-9) * A
            expected = Float64.(case["observation"]["values"])
            actual = [real(Z012[1, 1]), imag(Z012[1, 1]),
                      real(Z012[2, 2]), imag(Z012[2, 2]),
                      real(B012[1, 1]), real(B012[2, 2])]
            @test actual ≈ expected atol=1e-10 rtol=0

            construction = case["construction"]
            r75 = construction["r_ac_ref_ohm_per_m"] *
                  (1 + construction["alpha_20_per_degC"] * 55)
            x = [-1.1, 0.0, 1.1]
            y = [9.15, 9.15 + 1.1 * tan(21.67pi / 180), 9.15]
            bmopf = PowerOptLab.BMOPFTools.overhead_line_constants(
                fill(r75, 3), fill(construction["gmr_m"], 3),
                fill(construction["radius_m"], 3), x, y;
                frequency=50.0, earth_resistivity=100.0,
                earth_model="deri")
            uncertainty = case["uncertainty"]
            # DSS C-API has a small real-part convention difference, while
            # reactance and electrostatic capacitance agree much more closely.
            @test maximum(abs, real.(bmopf.Z .* 1000) .- real.(Z)) <=
                  uncertainty["deri_cross_implementation_real_z_tolerance"]
            @test maximum(abs, imag.(bmopf.Z .* 1000) .- imag.(Z)) <=
                  uncertainty["deri_cross_implementation_imag_z_tolerance"]
            @test maximum(abs, bmopf.C .* 1e12 .- C) <=
                  uncertainty["deri_cross_implementation_c_tolerance"]
        end

        values = Float64.(case["observation"]["values"])
        obs = SequenceLineObservation(
            z0=values[1] + im * values[2], z1=values[3] + im * values[4],
            frequency=50.0, earth_resistivity=100.0,
            sigma=fill(1e-4, 4))
        candidate = _ic_case_candidate(case; id="modified-carson-fit")
        mismatch = solve_inverse_carson(obs, [candidate]; starts=6)
        @test isempty(mismatch.compatible_candidates)
        @test maximum(abs, only(mismatch.fits).standardized_residual) > 3
    end
end
