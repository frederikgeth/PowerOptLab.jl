# Analytic companion to docs/src/tutorials/generalized_generator_models.md.
# Circuit identities only: this does not exercise a generalized-generator API.
using LinearAlgebra
using Printf
using Test

function analytic_examples()
    a = cis(2pi / 3)
    T = ComplexF64[1 1 1; 1 a a^2; 1 a^2 a] / 3
    A = inv(T)
    seq_z = ComplexF64[0.4 + 0.3im, 0.2 + 0.4im, 0.1 + 0.2im]
    Z = A * Diagonal(seq_z) * T
    e = 230 .* ComplexF64[1, a^2, a]
    i = A * ComplexF64[0, 10, 2]
    u = e - Z * i
    H = Hermitian((Z + Z') / 2)
    sint = sum(e .* conj.(i))
    spoc = sum(u .* conj.(i))

    @testset "Series primitive, sequences, power and rotation" begin
        @test T * e ≈ ComplexF64[0, 230, 0] atol=1e-12
        @test T * u ≈ ComplexF64[0, 228 - 4im, -0.2 - 0.4im] atol=1e-12
        @test sint - spoc ≈ dot(i, Z * i) atol=1e-10
        @test real(sint - spoc) ≈ real(dot(i, H * i)) atol=1e-10
        @test real(sint - spoc) ≈ 61.2 atol=1e-10
        @test sint ≈ 3sum((T * e) .* conj.(T * i)) atol=1e-10
        w = cis(0.37)
        @test w .* e - w .* u ≈ Z * (w .* i) atol=1e-12
        @test sum((w .* u) .* conj.(w .* i)) ≈ spoc atol=1e-10
        # The exact ideal primitive gives zero drop; allocation/variable-count
        # tests belong to the future component implementation.
        @test e - zeros(ComplexF64, 3, 3) * i == e
    end
    e_reg = real(seq_z[2] * 10) + sqrt(230^2 - imag(seq_z[2] * 10)^2)
    @test abs(e_reg - seq_z[2] * 10) ≈ 230
    @printf("Fixed EMF: |U1| = %.6f V; |U2| = %.6f V; loss = %.6f W\n",
        abs((T * u)[2]), abs((T * u)[3]), real(sint - spoc))
    @printf("Imposed-current terminal-voltage target: E1 = %.6f V\n", e_reg)

    @testset "Grounded source: distinct neutral and earth return" begin
        phase_i = ComplexF64[10, 0, 0]
        zn, zg = 0.2, 0.8
        vstar = -sum(phase_i) / (1 / zn + 1 / zg)
        in_ = vstar / zn
        ig = vstar / zg
        conductor_i = [phase_i; in_]
        @test vstar ≈ -1.6
        @test in_ ≈ -8
        @test ig ≈ -2
        @test sum(conductor_i) ≈ 2
        @test abs(sum(conductor_i) + ig) < 1e-12
        @test abs(3(T * phase_i)[1] + in_ + ig) < 1e-12
        v = ComplexF64[230, 230a^2, 230a, 0]
        emf = v[1:3] .- vstar # zero phase impedance, star-referenced EMFs
        source_power = sum(emf .* conj.(phase_i))
        terminal_power = sum(v .* conj.(conductor_i))
        ground_and_neutral_loss = zn * abs2(in_) + zg * abs2(ig)
        @test source_power - terminal_power ≈ ground_and_neutral_loss
        @test ground_and_neutral_loss ≈ 16
        @printf("Grounded source: In = %.1f A; Ig = %.1f A; four-wire sum = %.1f A; loss = %.1f W\n",
            real(in_), real(ig), real(sum(conductor_i)), ground_and_neutral_loss)

        # A nonzero PCC neutral makes the earth-return power correction visible.
        vn = 1.0
        shifted_star = (vn / zn - sum(phase_i)) / (1 / zn + 1 / zg)
        shifted_in = (shifted_star - vn) / zn
        shifted_ig = shifted_star / zg
        shifted_v = [v[1:3]; vn]
        shifted_j = [phase_i; shifted_in]
        phase_neutral_power = sum((v[1:3] .- vn) .* conj.(phase_i))
        @test sum(shifted_v .* conj.(shifted_j)) ≈ phase_neutral_power - vn * conj(shifted_ig)
        @test sum(shifted_v .* conj.(shifted_j)) ≈ 2291
    end

    @testset "Capability and topology counterexamples" begin
        seq_i = ComplexF64[0, 6, 6]
        phase_i = A * seq_i
        @test all(abs.(seq_i) .<= 8)
        @test phase_i ≈ ComplexF64[12, -6, -6] atol=1e-12
        @test maximum(abs.(phase_i)) > 10
        unequal = ComplexF64[230, 220a^2, 240a]
        @test abs((T * unequal)[1]) ≈ 10 / sqrt(3)
        @test abs((T * unequal)[3]) ≈ 10 / sqrt(3)
        @test abs(sum(unequal)) > 1 # cannot close a delta of ideal differences

        C = [1 0; 0 1; -1 -1]
        split_v = ComplexF64[120, -120, 0]
        split_i = ComplexF64[10, -6]
        split_u, split_j = C' * split_v, C * split_i
        @test split_j == ComplexF64[10, -6, -4]
        @test sum(split_j) == 0
        @test abs(split_v[1] - split_v[2]) == 240
        @test sum(split_u .* conj.(split_i)) == sum(split_v .* conj.(split_j)) == 1920

        # This phase/neutral reduction is valid only for the single return case.
        K = [Matrix{Float64}(I, 3, 3); -ones(1, 3)]
        Zc = Diagonal(ComplexF64[0.1, 0.2, 0.3, 0.4])
        @test K' * Zc * K ≈ Diagonal([0.1, 0.2, 0.3]) + 0.4ones(3, 3)
        @test sum(K * i) ≈ 0 atol=1e-12
    end
    println("All analytic checks passed; see tradeoffs.jl for network integration.")
end

analytic_examples()
