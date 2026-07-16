# Tests for the HELM power-flow solver stack.
#
# Layer 0: Padé / Wynn-epsilon analytic continuation on series with known
#          limits — including series that DIVERGE at the evaluation point but
#          have an analytic continuation there (the HELM situation).
# Later layers (added with the solver): analytic power-flow fixtures, OpenDSS
# parity, and numerical divergence diagnostics around known collapse points.
#
# HELM itself now lives in PowerOptLab (bespoke algorithm); it consumes the
# BMOPFTools engine through its public exports only.

using BMOPFTools: solve_pf, read_result, write_result, ybus_linearized, from_dss

@testset "helm: pade layer" begin

    _pade_sum = PowerOptLab._pade_sum
    _wynn_epsilon = PowerOptLab._wynn_epsilon

    @testset "geometric series inside the radius" begin
        # Σ zᵏ = 1/(1−z). Padé is exact for geometric series (rank-1 table),
        # so even slow-converging tails lock in immediately.
        for z in (0.5 + 0.0im, 0.9 + 0.0im, 0.6 + 0.3im)
            coeffs = [z^k for k in 0:12]
            v, spread = _pade_sum(coeffs)
            @test v ≈ 1 / (1 - z)  rtol=1e-10
            @test spread < 1e-8
        end
    end

    @testset "geometric series OUTSIDE the radius (analytic continuation)" begin
        # Σ 1.5ᵏ diverges at s=1, but the function 1/(1−z) continues to −2.
        coeffs = [(1.5 + 0.0im)^k for k in 0:10]
        v, spread = _pade_sum(coeffs)
        @test v ≈ -2.0 + 0.0im  rtol=1e-10
        @test spread < 1e-8
    end

    @testset "log(1+s) at s=1 (branch point at s=−1)" begin
        # c₀=0, c_k = (−1)^{k+1}/k → ln 2. Plain partial sums need ~10⁸ terms
        # for 1e-8; the epsilon table needs ~15.
        coeffs = ComplexF64[0.0; [(-1.0)^(k+1)/k for k in 1:15]]
        v, spread = _pade_sum(coeffs)
        @test v ≈ log(2.0)  atol=1e-9
        @test spread < 1e-6
    end

    @testset "sqrt-type branch point beyond s=1" begin
        # √(1−s/2) at s=1 = √(1/2); radius of convergence 2 > 1, but the
        # nearby branch point makes plain summation slow.
        cf = ComplexF64[1.0]
        for k in 1:20
            push!(cf, cf[end] * (0.5 - (k - 1)) / k * (-0.5))  # binomial(1/2,k)·(−1/2)ᵏ
        end
        v, spread = _pade_sum(cf)
        @test v ≈ sqrt(0.5)  atol=1e-10
        @test spread < 1e-8
    end

    @testset "no limit at the point: spread stagnates" begin
        # √(1−s) at s=1: the value exists (0) but the DERIVATIVE blows up;
        # harder: 1/√(1−s) has NO finite value at s=1. The epsilon table must
        # not pretend otherwise: spread stays large relative to convergence.
        cf = ComplexF64[1.0]
        for k in 1:24
            push!(cf, cf[end] * (k - 0.5) / k)   # coefficients of (1−s)^(−1/2)
        end
        _, spread = _pade_sum(cf)
        @test spread > 1e-4      # nothing like the ~1e-8 lock-in of true limits
    end

    @testset "degenerate/guard paths (pade)" begin
        # Already-converged sequence: zero differences short-circuit safely.
        v, spread = _wynn_epsilon(fill(3.0 + 1.0im, 6))
        @test v == 3.0 + 1.0im
        @test spread == 0.0
        # Single partial sum: no acceleration information.
        v, spread = _wynn_epsilon([2.0 + 0.0im])
        @test v == 2.0 + 0.0im && spread == Inf
        @test_throws ArgumentError _wynn_epsilon(ComplexF64[])
        @test_throws ArgumentError _pade_sum(ComplexF64[])
    end
end

# ── analytic solver layer ────────────────────────────────────────────────────

@testset "helm: analytic power flow" begin

    # 2-node single-phase 4-wire feeder with a closed-form solution:
    # source E on src.a (src.n perfectly grounded), line with phase resistance
    # R_a and neutral resistance R_n, constant-P load across (ld.a, ld.n).
    # Loop resistance R = R_a + R_n:
    #   ΔV² − E·ΔV + R·P = 0  ⇒  ΔV = (E + √(E² − 4RP))/2  (operational branch)
    # with collapse at P* = E²/(4R).
    function _two_node_net(P::Float64; model::String="constant_power",
                           E::Float64=230.0, Ra::Float64=0.5, Rn::Float64=0.5)
        net = Dict{String,Any}(
            "bus" => Dict{String,Any}(
                "src" => Dict{String,Any}("terminal_names" => ["a", "n"],
                                          "perfectly_grounded_terminals" => ["n"]),
                "ld"  => Dict{String,Any}("terminal_names" => ["a", "n"]),
            ),
            "voltage_source" => Dict{String,Any}(
                "vs" => Dict{String,Any}(
                    "bus" => "src", "terminal_map" => ["a"],
                    "configuration" => "WYE",
                    "v_magnitude" => [E], "v_angle" => [0.0]),
            ),
            "line" => Dict{String,Any}(
                "l1" => Dict{String,Any}(
                    "bus_from" => "src", "bus_to" => "ld",
                    "terminal_map_from" => ["a", "n"], "terminal_map_to" => ["a", "n"],
                    "R_series_1_1" => Ra, "X_series_1_1" => 0.0,
                    "R_series_2_2" => Rn, "X_series_2_2" => 0.0),
            ),
            "load" => Dict{String,Any}(
                "d1" => Dict{String,Any}(
                    "bus" => "ld", "terminal_map" => ["a", "n"],
                    "configuration" => "WYE", "model" => model,
                    "p_nom" => [P], "q_nom" => [0.0]),
            ),
        )
        model == "constant_power" || (net["load"]["d1"]["v_nom"] = [E])
        net
    end

    _dv_closed_form(P; E=230.0, R=1.0) = (E + sqrt(E^2 - 4R * P)) / 2

    # Residual oracle shared with ybus_linearized (the plan's cross-code-path
    # check): at the HELM solution, ‖Y_lin·V − i_comp(V)‖∞ ≈ 0 off-source.
    function _lin_residual(net, hr)
        lin = ybus_linearized(net; fold=:constant_z)
        n = length(lin.nodes)
        Vv = zeros(ComplexF64, n)
        for (nd, gi) in lin.index
            gi == 0 && continue
            Vv[gi] = hr.V[nd]
        end
        r = lin.Y * Vv .- lin.i_comp(Vv)
        srcb = Set(string(get(vs, "bus", "")) for (_, vs) in net["voltage_source"])
        maximum(abs(r[gi]) for (nd, gi) in lin.index if gi != 0 && !(nd[1] in srcb);
                init=0.0)
    end

    # Natural-parameter continuation for the nonlinear reference solve. Each
    # corrector starts from the last feasible voltage state, so this exercises a
    # connected loading path rather than a collection of independent flat-start
    # solves. The two-node fixture has an analytic nose, which supplies the
    # independent oracle for the continuation bracket.
    function _continuation_trace(lambdas, Pstar)
        previous = nothing
        trace = NamedTuple[]
        for lambda in lambdas
            hook! = previous === nothing ? nothing : ctx -> begin
                for ((bus, terminal), vr) in ctx.vars[:vr]
                    terminal_result = previous["bus"][bus][terminal]
                    JuMP.set_start_value(vr, terminal_result["vr"])
                    JuMP.set_start_value(ctx.vars[:vi][(bus, terminal)],
                                         terminal_result["vi"])
                end
            end
            result = solve_pf(_two_node_net(lambda * Pstar);
                per_unit=false, optimizer=Ipopt.Optimizer, model_hook! = hook!)
            push!(trace, (lambda=Float64(lambda), result=result))
            result["feasible"] && (previous = result)
        end
        trace
    end

    @testset "constant-P load: closed form + oracle residual" begin
        P = 2000.0
        hr = helm_series(_two_node_net(P))
        @test hr.status == :converged && hr.converged
        dv = hr.V[("ld", "a")] - hr.V[("ld", "n")]
        @test dv ≈ _dv_closed_form(P) rtol=1e-9
        @test abs(imag(dv)) < 1e-9
        # 4-wire detail: the neutral at the load rises above earth by the
        # return-path drop (+R_n·I, current flowing ld.n → src.n).
        Iload = P / real(dv)
        @test hr.V[("ld", "n")] ≈ 0.5 * Iload rtol=1e-9
        @test _lin_residual(_two_node_net(P), hr) < 1e-6
        # On this analytic saddle-node fixture, the Domb–Sykes singularity
        # estimate tracks P* = E²/4R = 13225 W ⇒ λ* ≈ 6.61.
        @test hr.singularity_estimate ≈ (230.0^2 / 4) / P rtol=0.1
        @test hr.load_margin == hr.singularity_estimate  # compatibility alias
        @test length(hr.pade_spread) == size(hr.coeffs, 1)
        @test all(x -> x >= 0, hr.pade_spread)
        @test length(hr.coefficient_tail_ratios) + 1 ==
              length(hr.coefficient_tail_norms)
    end

    @testset "constant-P at half the collapse loading: singularity estimate ≈ 2" begin
        Pstar = 230.0^2 / 4
        hr = helm_series(_two_node_net(Pstar / 2); max_order=60)
        @test hr.status == :converged
        @test hr.V[("ld", "a")] - hr.V[("ld", "n")] ≈ _dv_closed_form(Pstar / 2) rtol=1e-7
        @test hr.singularity_estimate ≈ 2.0 rtol=0.1
    end

    @testset "past analytic collapse: series-divergence diagnostic" begin
        Pstar = 230.0^2 / 4
        hr = helm_series(_two_node_net(1.2 * Pstar))
        @test hr.status == :series_diverged
        @test !hr.converged
        @test all(>(1.0), hr.coefficient_tail_ratios)
        @test 0.7 < hr.singularity_estimate < 1.0     # ≈ 1/1.2
    end

    @testset "natural-parameter continuation brackets analytic collapse" begin
        Pstar = 230.0^2 / 4
        lambdas = vcat(collect(0.50:0.02:0.98), [1.02, 1.06, 1.10])
        trace = _continuation_trace(lambdas, Pstar)
        feasible = [point.lambda for point in trace if point.result["feasible"]]
        infeasible = [point.lambda for point in trace if !point.result["feasible"]]
        @test maximum(feasible) ≈ 0.98
        @test minimum(infeasible) ≈ 1.02
        @test all(point.result["feasible"] for point in trace if point.lambda < 1)
        @test all(!point.result["feasible"] for point in trace if point.lambda > 1)

        # Prove the reference correctors were actually warm-started from their
        # predecessor, then check the whole feasible trace against the closed
        # form rather than only its endpoint.
        @test trace[2].result["initialisation"]["ld"]["a"]["vr_init"] ≈
              trace[1].result["bus"]["ld"]["a"]["vr"]
        for point in trace
            point.result["feasible"] || continue
            bus = point.result["bus"]["ld"]
            dv = complex(bus["a"]["vr"], bus["a"]["vi"]) -
                 complex(bus["n"]["vr"], bus["n"]["vi"])
            @test dv ≈ _dv_closed_form(point.lambda * Pstar) rtol=1e-6
        end

        below = helm_series(_two_node_net(0.98 * Pstar); max_order=80, tol=1e-6)
        above = helm_series(_two_node_net(1.02 * Pstar); max_order=80, tol=1e-6)
        @test below.status == :converged
        @test !above.converged
        @test above.status in (:series_diverged, :max_order_reached)
        @test below.singularity_estimate ≈ inv(0.98) rtol=0.1
        @test above.singularity_estimate ≈ inv(1.02) rtol=0.1
    end

    @testset "just below collapse still converges" begin
        Pstar = 230.0^2 / 4
        hr = helm_series(_two_node_net(0.9 * Pstar); max_order=60, tol=1e-6)
        @test hr.status == :converged
        @test hr.V[("ld", "a")] - hr.V[("ld", "n")] ≈ _dv_closed_form(0.9 * Pstar) rtol=1e-4
    end

    @testset "constant-Z load: linear, fast, oracle-exact" begin
        net = _two_node_net(2000.0; model="constant_impedance")
        hr = helm_series(net)
        @test hr.status == :converged
        @test _lin_residual(net, hr) < 1e-8
        # Hand check: y = P/Vnom² = 2000/230² S across (ld.a, ld.n).
        y = 2000.0 / 230.0^2
        dv = hr.V[("ld", "a")] - hr.V[("ld", "n")]
        @test dv ≈ 230.0 / (1 + y * 1.0) rtol=1e-9      # E·y_load/(y_load+1/R)… = E/(1+yR)
    end

    @testset "unbalanced 3-phase 4-wire: germ is not flat, neutral shifts" begin
        E = [240.0, 230.0, 220.0]
        net = Dict{String,Any}(
            "bus" => Dict{String,Any}(
                "src" => Dict{String,Any}("terminal_names" => ["a", "b", "c", "n"],
                                          "perfectly_grounded_terminals" => ["n"]),
                "ld"  => Dict{String,Any}("terminal_names" => ["a", "b", "c", "n"]),
            ),
            "voltage_source" => Dict{String,Any}(
                "vs" => Dict{String,Any}(
                    "bus" => "src", "terminal_map" => ["a", "b", "c"],
                    "configuration" => "WYE",
                    "v_magnitude" => E, "v_angle" => [0.0, -2π/3, 2π/3]),
            ),
            "line" => Dict{String,Any}(
                "l1" => Dict{String,Any}(
                    "bus_from" => "src", "bus_to" => "ld",
                    "terminal_map_from" => ["a", "b", "c", "n"],
                    "terminal_map_to"   => ["a", "b", "c", "n"],
                    Dict("R_series_$(k)_$(k)" => 0.4 for k in 1:4)...,
                    Dict("X_series_$(k)_$(k)" => 0.3 for k in 1:4)...),
            ),
            "load" => Dict{String,Any}(
                "d1" => Dict{String,Any}(
                    "bus" => "ld", "terminal_map" => ["a", "b", "c", "n"],
                    "configuration" => "WYE",
                    "p_nom" => [5000.0, 2000.0, 1000.0],
                    "q_nom" => [1000.0, 500.0, 0.0]),
            ),
        )
        hr = helm_series(net)
        @test hr.status == :converged
        # Germ (order-0 coefficients) carries the source unbalance exactly.
        @test hr.V[("src", "a")] ≈ 240.0 + 0im
        @test hr.V[("src", "b")] ≈ 230.0 * cis(-2π/3)
        # Unbalanced load ⇒ the load-bus neutral rises off earth.
        @test abs(hr.V[("ld", "n")]) > 1.0
        @test _lin_residual(net, hr) < 1e-5
    end

    @testset "switch interplay: :alias ≡ :constrain, w = feeder current" begin
        P = 2000.0
        net = _two_node_net(P)
        # Insert a closed switch src—mid, retarget the line mid—ld.
        net["bus"]["mid"] = Dict{String,Any}("terminal_names" => ["a", "n"])
        net["switch"] = Dict{String,Any}(
            "sw1" => Dict{String,Any}(
                "bus_from" => "src", "bus_to" => "mid", "status" => "closed",
                "terminal_map_from" => ["a", "n"], "terminal_map_to" => ["a", "n"]))
        net["line"]["l1"]["bus_from"] = "mid"

        ha = helm_series(net)                        # switches=:alias
        hc = helm_series(net; switches=:constrain)
        @test ha.status == :converged && hc.status == :converged
        for nd in (("ld", "a"), ("ld", "n"), ("mid", "a"))
            @test hc.V[nd] ≈ ha.V[nd] atol=1e-8
        end
        # Phase-conductor switch current = load current P/ΔV.
        dv = real(hc.V[("ld", "a")] - hc.V[("ld", "n")])
        iph = P / dv
        cph = hc.couplings[findfirst(c -> c.conductor == 1, hc.couplings)]
        wph = cph.scale * hc.w[findfirst(c -> c === cph, hc.couplings)]
        @test abs(abs(wph) - iph) < 1e-6 * iph
    end

    @testset "validation errors" begin
        # Constant-current fraction: not holomorphic in v1 → informative error.
        net = _two_node_net(2000.0; model="constant_current")
        @test_throws ArgumentError helm_series(net)

        # DELTA voltage source unsupported.
        net = _two_node_net(2000.0)
        net["voltage_source"]["vs"]["configuration"] = "DELTA"
        @test_throws ArgumentError helm_series(net)

        # No voltage source at all.
        net = _two_node_net(2000.0)
        delete!(net, "voltage_source")
        @test_throws ArgumentError helm_series(net)

        # Islanded bus: singular system, error names the floating node.
        net = _two_node_net(2000.0)
        net["bus"]["iso"] = Dict{String,Any}("terminal_names" => ["a"])
        err = try helm_series(net); nothing catch e; e end
        @test err isa ErrorException
        @test occursin("iso", sprint(showerror, err))
    end
end

# ── delta / line-to-line loads + result-dict wrapper ─────────────────────────

@testset "helm: delta and L-L loads, solve_pf_helm" begin

    # 3-phase 3-wire source bus with grounded neutral reference at the source.
    function _three_phase_net(; load::Dict{String,Any})
        Dict{String,Any}(
            "bus" => Dict{String,Any}(
                "src" => Dict{String,Any}("terminal_names" => ["a", "b", "c", "n"],
                                          "perfectly_grounded_terminals" => ["n"]),
                "ld"  => Dict{String,Any}("terminal_names" => ["a", "b", "c"]),
            ),
            "voltage_source" => Dict{String,Any}(
                "vs" => Dict{String,Any}(
                    "bus" => "src", "terminal_map" => ["a", "b", "c"],
                    "configuration" => "WYE",
                    "v_magnitude" => [230.0, 230.0, 230.0],
                    "v_angle" => [0.0, -2π/3, 2π/3]),
            ),
            "line" => Dict{String,Any}(
                "l1" => Dict{String,Any}(
                    "bus_from" => "src", "bus_to" => "ld",
                    "terminal_map_from" => ["a", "b", "c"],
                    "terminal_map_to"   => ["a", "b", "c"],
                    Dict("R_series_$(k)_$(k)" => 0.5 for k in 1:3)...,
                    Dict("X_series_$(k)_$(k)" => 0.25 for k in 1:3)...),
            ),
            "load" => Dict{String,Any}("d1" => load),
        )
    end

    _lin_residual3(net, hr) = begin
        lin = ybus_linearized(net; fold=:constant_z)
        Vv = zeros(ComplexF64, length(lin.nodes))
        for (nd, gi) in lin.index
            gi == 0 && continue
            Vv[gi] = hr.V[nd]
        end
        r = lin.Y * Vv .- lin.i_comp(Vv)
        maximum(abs(r[gi]) for (nd, gi) in lin.index if gi != 0 && nd[1] != "src";
                init=0.0)
    end

    @testset "SINGLE_PHASE line-to-line load: rotated closed form" begin
        # Load between phases a and b: loop R = 1.0, |E_ab| = 230·√3, and by
        # rotation invariance ΔV = e^{jθ}·dv with dv from the real quadratic.
        P = 5000.0
        net = _three_phase_net(load=Dict{String,Any}(
            "bus" => "ld", "terminal_map" => ["a", "b"],
            "configuration" => "SINGLE_PHASE", "p_nom" => [P], "q_nom" => [0.0]))
        # Make the a/b conductors purely resistive for the closed form.
        for k in 1:3
            net["line"]["l1"]["X_series_$(k)_$(k)"] = 0.0
        end
        hr = helm_series(net)
        @test hr.status == :converged
        Eab = 230.0 - 230.0 * cis(-2π/3)
        dv = hr.V[("ld", "a")] - hr.V[("ld", "b")]
        dv_cf = (abs(Eab) + sqrt(abs(Eab)^2 - 4 * 1.0 * P)) / 2
        @test abs(dv) ≈ dv_cf rtol=1e-9
        @test angle(dv) ≈ angle(Eab) atol=1e-9
        @test _lin_residual3(net, hr) < 1e-6
    end

    @testset "balanced DELTA load: oracle residual + symmetry" begin
        P = 4000.0
        net = _three_phase_net(load=Dict{String,Any}(
            "bus" => "ld", "terminal_map" => ["a", "b", "c"],
            "configuration" => "DELTA",
            "p_nom" => [P, P, P], "q_nom" => [500.0, 500.0, 500.0]))
        hr = helm_series(net)
        @test hr.status == :converged
        @test _lin_residual3(net, hr) < 1e-6
        # Balanced delta on a balanced source: phase voltage magnitudes equal.
        vms = [abs(hr.V[("ld", t)]) for t in ("a", "b", "c")]
        @test maximum(vms) - minimum(vms) < 1e-6 * maximum(vms)
    end

    @testset "unbalanced DELTA load: oracle residual" begin
        net = _three_phase_net(load=Dict{String,Any}(
            "bus" => "ld", "terminal_map" => ["a", "b", "c"],
            "configuration" => "DELTA",
            "p_nom" => [6000.0, 1500.0, 3000.0], "q_nom" => [1000.0, 0.0, -500.0]))
        hr = helm_series(net)
        @test hr.status == :converged
        @test _lin_residual3(net, hr) < 1e-6
    end

    @testset "solve_pf_helm result dict + write/read round-trip" begin
        P = 4000.0
        net = _three_phase_net(load=Dict{String,Any}(
            "bus" => "ld", "terminal_map" => ["a", "b", "c"],
            "configuration" => "DELTA", "p_nom" => [P, P, P],
            "q_nom" => [0.0, 0.0, 0.0]))
        res = solve_pf_helm(net)
        @test res["termination_status"] == "HELM_CONVERGED"
        @test res["feasible"] === true
        @test res["solve_time"] >= 0.0
        hr = helm_series(net)
        vb = res["bus"]["ld"]["a"]
        @test vb["vr"] ≈ real(hr.V[("ld", "a")])
        @test vb["vm"] ≈ abs(hr.V[("ld", "a")])
        @test vb["va"] ≈ angle(hr.V[("ld", "a")])
        @test res["bus"]["src"]["n"]["vm"] == 0.0      # grounded reference
        @test res["helm"]["order"] == hr.n_order
        @test isfinite(res["helm"]["residual"])

        # Round-trip through the standard result IO.
        path = joinpath(mktempdir(), "helm_result.json")
        write_result(res, path)
        back = read_result(path)
        @test back["termination_status"] == "HELM_CONVERGED"
        @test back["bus"]["ld"]["a"]["vm"] ≈ vb["vm"]
    end

    @testset "solve_pf_helm collapse: NaN-filled + margin < 1" begin
        Pstar = 230.0^2 / 4
        net = Dict{String,Any}(
            "bus" => Dict{String,Any}(
                "src" => Dict{String,Any}("terminal_names" => ["a", "n"],
                                          "perfectly_grounded_terminals" => ["n"]),
                "ld"  => Dict{String,Any}("terminal_names" => ["a", "n"]),
            ),
            "voltage_source" => Dict{String,Any}(
                "vs" => Dict{String,Any}(
                    "bus" => "src", "terminal_map" => ["a"],
                    "configuration" => "WYE",
                    "v_magnitude" => [230.0], "v_angle" => [0.0]),
            ),
            "line" => Dict{String,Any}(
                "l1" => Dict{String,Any}(
                    "bus_from" => "src", "bus_to" => "ld",
                    "terminal_map_from" => ["a", "n"], "terminal_map_to" => ["a", "n"],
                    "R_series_1_1" => 0.5, "X_series_1_1" => 0.0,
                    "R_series_2_2" => 0.5, "X_series_2_2" => 0.0),
            ),
            "load" => Dict{String,Any}(
                "d1" => Dict{String,Any}(
                    "bus" => "ld", "terminal_map" => ["a", "n"],
                    "configuration" => "WYE", "p_nom" => [1.3 * Pstar],
                    "q_nom" => [0.0]),
            ),
        )
        res = solve_pf_helm(net)
        @test res["termination_status"] == "HELM_SERIES_DIVERGED"
        @test res["feasible"] === false
        @test isnan(res["bus"]["ld"]["a"]["vm"])
        @test res["helm"]["singularity_estimate"] < 1.0     # ≈ 1/1.3
        @test res["helm"]["load_margin"] ==
              res["helm"]["singularity_estimate"]
        @test !isempty(res["helm"]["pade_spread"])
        @test all(>(1.0), res["helm"]["coefficient_tail_ratios"])
    end

    @testset "coupling currents in the result dict (switch :constrain)" begin
        P = 2000.0
        net = Dict{String,Any}(
            "bus" => Dict{String,Any}(
                "src" => Dict{String,Any}("terminal_names" => ["a", "n"],
                                          "perfectly_grounded_terminals" => ["n"]),
                "mid" => Dict{String,Any}("terminal_names" => ["a", "n"]),
                "ld"  => Dict{String,Any}("terminal_names" => ["a", "n"]),
            ),
            "voltage_source" => Dict{String,Any}(
                "vs" => Dict{String,Any}(
                    "bus" => "src", "terminal_map" => ["a"],
                    "configuration" => "WYE",
                    "v_magnitude" => [230.0], "v_angle" => [0.0]),
            ),
            "switch" => Dict{String,Any}(
                "sw1" => Dict{String,Any}(
                    "bus_from" => "src", "bus_to" => "mid", "status" => "closed",
                    "terminal_map_from" => ["a", "n"], "terminal_map_to" => ["a", "n"])),
            "line" => Dict{String,Any}(
                "l1" => Dict{String,Any}(
                    "bus_from" => "mid", "bus_to" => "ld",
                    "terminal_map_from" => ["a", "n"], "terminal_map_to" => ["a", "n"],
                    "R_series_1_1" => 0.5, "X_series_1_1" => 0.0,
                    "R_series_2_2" => 0.5, "X_series_2_2" => 0.0),
            ),
            "load" => Dict{String,Any}(
                "d1" => Dict{String,Any}(
                    "bus" => "ld", "terminal_map" => ["a", "n"],
                    "configuration" => "WYE", "p_nom" => [P], "q_nom" => [0.0]),
            ),
        )
        res = solve_pf_helm(net; switches=:constrain)
        @test res["termination_status"] == "HELM_CONVERGED"
        @test haskey(res["coupling"], "sw1")
        dv = 230.0 * 0 + (res["bus"]["ld"]["a"]["vr"] - res["bus"]["ld"]["n"]["vr"])
        iph = P / dv
        @test res["coupling"]["sw1"]["1"]["im"] ≈ iph rtol=1e-6
        @test res["coupling"]["sw1"]["1"]["kind"] == "switch"
    end
end

# ── OpenDSS parity + 3-way oracle (gated) ────────────────────────────────────

@testset "helm: OpenDSS parity" begin
    if !_HAS_ODS
        @test_skip "Requires OpenDSSDirect"
    else
        # OpenDSS node names use ".1/.2/.3/.4"; from_dss uses a/b/c/n.
        _TNH = Dict("a" => "1", "b" => "2", "c" => "3", "n" => "4",
                    "1" => "1", "2" => "2", "3" => "3", "4" => "4")

        # Solve `deck` in OpenDSS and with HELM; compare node-to-earth voltages
        # over all non-source terminals. Returns the max |ΔV| (V).
        function _helm_vs_ods(deck::String, srcbus::String; atol::Float64,
                              kwargs...)
            path = normpath(abspath(joinpath(@__DIR__, "data", "pf_comparison", deck)))
            OpenDSSDirect.dss("Clear")
            OpenDSSDirect.dss("Redirect \"$path\"")
            OpenDSSDirect.dss("Solve")
            names = lowercase.(OpenDSSDirect.Circuit.AllNodeNames())
            volts = OpenDSSDirect.Circuit.AllBusVolts()
            odsv = Dict(nm => v for (nm, v) in zip(names, volts))

            net = from_dss(path)
            hr = helm_series(net; kwargs...)
            @test hr.status == :converged

            maxerr = 0.0
            for (bid, b) in net["bus"], t in b["terminal_names"]
                bid == srcbus && continue
                k = "$(bid).$(get(_TNH, string(t), string(t)))"
                haskey(odsv, k) || continue
                maxerr = max(maxerr, abs(hr.V[(bid, string(t))] - odsv[k]))
            end
            @test maxerr < atol
            maxerr
        end

        # Constant-power decks (lines, explicit neutrals, delta loads, caps):
        # sub-0.1 V agreement. Excluded, with reasons:
        #   pf_3ph_line / pf_zip_3ph — the 4-wire-line import quirk (4×4
        #     linecode but the line terminal_map imports as a/b/c) leaves lb.n
        #     connected only through loads: structurally singular, and HELM's
        #     floating-node error correctly refuses it;
        #   pf_1ph_freeneutral — the load-bus neutral is held ONLY by the load,
        #     so the no-load germ is genuinely indeterminate (a fundamental
        #     HELM requirement: every node needs a linear path to a reference);
        #   pf_1ph_impedanceneutral — same 4-wire-line import quirk: the
        #     dropped neutral conductor forces ALL return current through the
        #     0.2 Ω grounding reactor (~2× the true neutral rise). The OPF
        #     comparisons pass on this deck only because they hand-build the
        #     net instead of importing it.
        _helm_vs_ods("pf_1ph_line.dss",   "src"; atol=0.1)
        _helm_vs_ods("pf_cap_wye.dss",    "src"; atol=0.1)
        _helm_vs_ods("pf_delta_load.dss", "src"; atol=0.1)
        # Perfect ground at lb.n hides most of the dropped-neutral quirk; the
        # ~0.2 V residual on phases b/c is the lost phase-neutral coupling
        # (import fidelity, not HELM).
        _helm_vs_ods("pf_1ph_perfectneutral.dss", "src"; atol=0.5)
        # Transformer decks: import fidelity dominates (same tolerance class as
        # the feasibility-OPF comparisons in powerflow_comparison_tests.jl).
        _helm_vs_ods("pf_yd_xfmr.dss", "hv"; atol=0.5)
        _helm_vs_ods("pf_dy_xfmr.dss", "hv"; atol=0.5)

        # ZIP decks carry constant-current fractions — HELM v1 must refuse them
        # loudly (documented limitation), not silently mis-solve.
        @testset "constant-I ZIP deck raises the v1 validation error" begin
            path = normpath(abspath(joinpath(@__DIR__, "data", "pf_comparison",
                                             "pf_zip_1ph.dss")))
            net = from_dss(path)
            @test_throws ArgumentError helm_series(net)
        end
    end
end

@testset "helm: 3-way oracle (HELM vs Ipopt solve_pf vs OpenDSS)" begin
    if !(_HAS_ODS && _HAS_JUMP_IPOPT)
        @test_skip "Requires OpenDSSDirect + JuMP/Ipopt"
    else
        path = normpath(abspath(joinpath(@__DIR__, "data", "pf_comparison",
                                         "pf_1ph_line.dss")))
        net = from_dss(path)

        hr = helm_series(net)
        @test hr.status == :converged
        res = solve_pf(net; optimizer=Ipopt.Optimizer)

        OpenDSSDirect.dss("Clear")
        OpenDSSDirect.dss("Redirect \"$path\"")
        OpenDSSDirect.dss("Solve")
        names = lowercase.(OpenDSSDirect.Circuit.AllNodeNames())
        volts = OpenDSSDirect.Circuit.AllBusVolts()
        odsv = Dict(nm => v for (nm, v) in zip(names, volts))
        _TNH = Dict("a" => "1", "b" => "2", "c" => "3", "n" => "4")

        for (bid, b) in net["bus"], t in b["terminal_names"]
            bid == "src" && continue
            tt = string(t)
            v_helm  = hr.V[(bid, tt)]
            tv      = res["bus"][bid][tt]
            v_ipopt = tv["vr"] + im * tv["vi"]
            @test abs(v_helm - v_ipopt) < 0.1          # HELM ≈ Ipopt PF
            k = "$(bid).$(get(_TNH, tt, tt))"
            haskey(odsv, k) && @test abs(v_helm - odsv[k]) < 0.1   # ≈ OpenDSS
        end
    end
end
