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

using BMOPFTools

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
                for (bus, bus_data) in BMOPFTools.opf_network(ctx)["bus"]
                    bid = String(bus)
                    for terminal in String.(bus_data["terminal_names"])
                        terminal_result = previous["bus"][bid][terminal]
                        vr = BMOPFTools.opf_object(ctx,
                            BMOPFTools.opf_bus_voltage_key(bid, terminal))
                        vi = BMOPFTools.opf_object(ctx,
                            BMOPFTools.opf_bus_voltage_key(bid, terminal;
                                                           component=:imag))
                        JuMP.set_start_value(vr, terminal_result["vr"])
                        JuMP.set_start_value(vi, terminal_result["vi"])
                    end
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
        @test solve_status(hr).publishable
        @test solve_diagnostics(hr).residual == hr.residual
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
                "bus_from" => "src", "bus_to" => "mid", "open_switch" => false,
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
        # Generators / IBRs have no fixed power-flow law: refused, not dropped.
        net = _two_node_net(2000.0)
        net["generator"] = Dict{String,Any}("g1" => Dict{String,Any}("bus" => "ld"))
        err = try helm_series(net); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("g1", sprint(showerror, err))
        net = _two_node_net(2000.0)
        net["ibr"] = Dict{String,Any}("pv1" => Dict{String,Any}("bus" => "ld"))
        @test_throws ArgumentError helm_series(net)

        # Unknown embedding options.
        @test_throws ArgumentError helm_series(_two_node_net(2000.0); embedding=:flat)
        @test_throws ArgumentError helm_series(_two_node_net(2000.0);
                                               source_reference=:ground)

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
                    "bus_from" => "src", "bus_to" => "mid", "open_switch" => false,
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

# ── load laws, source references, ideal-transformer earthing, n-wire ────────

@testset "helm: load laws, 4-wire references and n-wire networks" begin

    _line(f, t, tm, R, X) = Dict{String,Any}(
        "bus_from" => f, "bus_to" => t,
        "terminal_map_from" => tm, "terminal_map_to" => tm,
        Dict("R_series_$(k)_$(k)" => R for k in 1:length(tm))...,
        Dict("X_series_$(k)_$(k)" => X for k in 1:length(tm))...)
    _src3(E=230.0) = Dict{String,Any}(
        "bus" => "src", "terminal_map" => ["a", "b", "c"], "configuration" => "WYE",
        "v_magnitude" => [E, E, E], "v_angle" => [0.0, -2π/3, 2π/3])

    # Single-phase 2-node feeder, loop resistance R = 1 Ω (as in the analytic
    # layer), with an arbitrary load model.
    function _feeder(P; E=230.0, fields...)
        load = Dict{String,Any}("bus" => "ld", "terminal_map" => ["a", "n"],
                                "configuration" => "WYE", "p_nom" => [P],
                                "q_nom" => [0.0], "v_nom" => [E])
        for (k, v) in fields
            load[string(k)] = v
        end
        Dict{String,Any}(
            "bus" => Dict{String,Any}(
                "src" => Dict{String,Any}("terminal_names" => ["a", "n"],
                                          "perfectly_grounded_terminals" => ["n"]),
                "ld"  => Dict{String,Any}("terminal_names" => ["a", "n"])),
            "voltage_source" => Dict{String,Any}("vs" => Dict{String,Any}(
                "bus" => "src", "terminal_map" => ["a"], "configuration" => "WYE",
                "v_magnitude" => [E], "v_angle" => [0.0])),
            "line" => Dict{String,Any}("l1" => _line("src", "ld", ["a", "n"], 0.5, 0.0)),
            "load" => Dict{String,Any}("d1" => load))
    end
    _dv(hr) = hr.V[("ld", "a")] - hr.V[("ld", "n")]

    # Real scalar Newton on v² − E·v + R·P(v) = 0 (operational root), the
    # closed-form oracle for any |ΔV|-dependent load law on the 2-node feeder.
    function _feeder_root(Pv, dPv; E=230.0)
        v = E
        for _ in 1:60
            v -= (v^2 - E * v + Pv(v)) / (2v - E + dPv(v))
        end
        v
    end

    # ybus_linearized oracle, optionally skipping whole buses.
    function _lin_res(net, hr; skip=("src",))
        lin = ybus_linearized(net; fold=:constant_z)
        Vv = zeros(ComplexF64, length(lin.nodes))
        for (nd, gi) in lin.index
            gi == 0 && continue
            Vv[gi] = hr.V[nd]
        end
        r = lin.Y * Vv .- lin.i_comp(Vv)
        maximum(abs(r[gi]) for (nd, gi) in lin.index
                if gi != 0 && !(nd[1] in skip); init=0.0)
    end

    @testset "constant-current, exponential and full ZIP laws ($emb)" for emb in
            (:load_scaling, :impedance_germ)
        E = 230.0; P = 5000.0
        hr = helm_series(_feeder(P; model="constant_current"); embedding=emb)
        @test hr.status == :converged
        @test _dv(hr) ≈ E - P / E rtol=1e-11            # P·|v|/E ⇒ linear root

        hr = helm_series(_feeder(P; model="exponential", gamma_p=[0.7],
                                 gamma_q=[0.0]); embedding=emb)
        @test hr.status == :converged
        @test real(_dv(hr)) ≈ _feeder_root(v -> P * (v / E)^0.7,
                                           v -> 0.7P * (v / E)^(-0.3) / E) rtol=1e-11

        hr = helm_series(_feeder(P; model="zip", alpha_z=[0.3], alpha_i=[0.3],
                                 alpha_p=[0.4], beta_z=[0.0], beta_i=[0.0],
                                 beta_p=[1.0]); embedding=emb)
        @test hr.status == :converged
        @test real(_dv(hr)) ≈ _feeder_root(
            v -> P * (0.3 * (v / E)^2 + 0.3 * v / E + 0.4),
            v -> P * (0.6 * v / E^2 + 0.3 / E)) rtol=1e-11
    end

    @testset "unbalanced 4-wire ZIP/exponential loads: oracle residual" begin
        net = Dict{String,Any}(
            "bus" => Dict{String,Any}(
                "src" => Dict{String,Any}("terminal_names" => ["a", "b", "c", "n"],
                                          "perfectly_grounded_terminals" => ["n"]),
                "ld"  => Dict{String,Any}("terminal_names" => ["a", "b", "c", "n"])),
            "voltage_source" => Dict{String,Any}("vs" => _src3()),
            "line" => Dict{String,Any}(
                "l1" => _line("src", "ld", ["a", "b", "c", "n"], 0.4, 0.3)),
            "load" => Dict{String,Any}(
                "z" => Dict{String,Any}("bus" => "ld", "terminal_map" => ["a", "b", "c", "n"],
                    "configuration" => "WYE", "model" => "zip", "v_nom" => [230.0],
                    "p_nom" => [5000.0, 2000.0, 1000.0], "q_nom" => [1000.0, 500.0, 0.0],
                    "alpha_z" => [0.2], "alpha_i" => [0.5], "alpha_p" => [0.3],
                    "beta_z" => [0.6], "beta_i" => [0.1], "beta_p" => [0.3]),
                "x" => Dict{String,Any}("bus" => "ld", "terminal_map" => ["a", "b", "c"],
                    "configuration" => "DELTA", "model" => "exponential",
                    "v_nom" => [400.0], "p_nom" => [1500.0, 500.0, 900.0],
                    "q_nom" => [200.0, 0.0, 100.0],
                    "gamma_p" => [1.3], "gamma_q" => [2.5])))
        for emb in (:load_scaling, :impedance_germ)
            hr = helm_series(net; embedding=emb)
            @test hr.status == :converged
            @test _lin_res(net, hr) < 1e-8
        end
    end

    @testset "unearthed star point: needs the impedance germ" begin
        # 3-wire feeder; the wye load's star point is held only by the load.
        net = Dict{String,Any}(
            "bus" => Dict{String,Any}(
                "src" => Dict{String,Any}("terminal_names" => ["a", "b", "c", "n"],
                                          "perfectly_grounded_terminals" => ["n"]),
                "ld"  => Dict{String,Any}("terminal_names" => ["a", "b", "c", "n"])),
            "voltage_source" => Dict{String,Any}("vs" => _src3()),
            "line" => Dict{String,Any}("l1" => _line("src", "ld", ["a", "b", "c"], 0.3, 0.2)),
            "load" => Dict{String,Any}("d1" => Dict{String,Any}(
                "bus" => "ld", "terminal_map" => ["a", "b", "c", "n"],
                "configuration" => "WYE", "model" => "zip", "v_nom" => [230.0],
                "p_nom" => [4000.0, 3000.0, 2500.0], "q_nom" => [500.0, 0.0, 200.0],
                "alpha_z" => [0.7], "alpha_i" => [0.2], "alpha_p" => [0.1],
                "beta_z" => [1.0], "beta_i" => [0.0], "beta_p" => [0.0])))
        err = try helm_series(net); nothing catch e; e end
        @test err isa ErrorException
        @test occursin("impedance_germ", sprint(showerror, err))
        hr = helm_series(net; embedding=:impedance_germ)
        @test hr.status == :converged
        @test _lin_res(net, hr) < 1e-8
        @test abs(hr.V[("ld", "n")]) > 10.0             # genuine star-point shift
    end

    @testset "source_reference = :neutral honours source earthing impedance" begin
        # Unbalanced 4-wire feeder; source neutral earthed through 5 Ω, load
        # neutral through 10 Ω. Constant-Z loads make an independent linear
        # MNA solve (ybus_linearized + hand-built source rows) the oracle.
        net = Dict{String,Any}(
            "bus" => Dict{String,Any}(
                "src" => Dict{String,Any}("terminal_names" => ["a", "b", "c", "n"]),
                "ld"  => Dict{String,Any}("terminal_names" => ["a", "b", "c", "n"])),
            "voltage_source" => Dict{String,Any}("vs" => _src3()),
            "line" => Dict{String,Any}(
                "l1" => _line("src", "ld", ["a", "b", "c", "n"], 0.3, 0.2)),
            "shunt" => Dict{String,Any}(
                "gs" => Dict{String,Any}("bus" => "src", "terminal_map" => ["n"],
                                         "G_1_1" => 0.2, "B_1_1" => 0.0),
                "gl" => Dict{String,Any}("bus" => "ld", "terminal_map" => ["n"],
                                         "G_1_1" => 0.1, "B_1_1" => 0.0)),
            "load" => Dict{String,Any}("d1" => Dict{String,Any}(
                "bus" => "ld", "terminal_map" => ["a", "b", "c", "n"],
                "configuration" => "WYE", "model" => "constant_impedance",
                "v_nom" => [230.0], "p_nom" => [8000.0, 1000.0, 500.0],
                "q_nom" => [1000.0, 0.0, 0.0])))
        hr = helm_series(net; source_reference=:neutral)
        @test hr.status == :converged

        lin = ybus_linearized(net; fold=:constant_z)
        nn = length(lin.nodes)
        A = zeros(ComplexF64, 3, nn)
        Es = 230.0 .* cis.([0.0, -2π/3, 2π/3])
        for (k, t) in enumerate(("a", "b", "c"))
            A[k, lin.index[("src", t)]] = 1.0
            A[k, lin.index[("src", "n")]] = -1.0
        end
        M = [Matrix(lin.Y) transpose(A); A zeros(3, 3)]
        x = M \ vcat(zeros(ComplexF64, nn), Es)
        for (nd, gi) in lin.index
            gi == 0 && continue
            @test hr.V[nd] ≈ x[gi] atol=1e-8
        end
        @test abs(hr.V[("src", "n")]) > 1.0             # source neutral floats
        # Phase-to-neutral source voltage exact; source currents conserve.
        @test abs(hr.V[("src", "a")] - hr.V[("src", "n")]) ≈ 230.0 rtol=1e-12
        @test abs(sum(values(hr.source_injection))) < 1e-9
        # The default :earth reference still pins the source neutral at 0 V.
        @test helm_series(net).V[("src", "n")] == 0.0

        # With a perfectly grounded source neutral, both references coincide.
        g = deepcopy(net)
        g["bus"]["src"]["perfectly_grounded_terminals"] = ["n"]
        he = helm_series(g); hn = helm_series(g; source_reference=:neutral)
        @test hn.V[("ld", "n")] ≈ he.V[("ld", "n")] atol=1e-9
        @test hn.source_injection[("src", "a")] ≈ he.source_injection[("src", "a")] atol=1e-8
    end

    @testset "ideal transformer keeps its neutral earthing and magnetising shunt" begin
        # The LV side's only earth reference is the transformer's internal
        # neutral-grounding resistor (plus a second earth electrode at the load
        # so earth current actually flows). ybus_augmented drops the whole
        # Yprim of a promoted ideal unit; HELM restores the shunt-only part.
        function xnet(rs)
            Dict{String,Any}(
                "bus" => Dict{String,Any}(
                    "hv" => Dict{String,Any}("terminal_names" => ["a", "n"],
                                             "perfectly_grounded_terminals" => ["n"]),
                    "lv" => Dict{String,Any}("terminal_names" => ["a", "n"]),
                    "ld" => Dict{String,Any}("terminal_names" => ["a", "n"])),
                "voltage_source" => Dict{String,Any}("vs" => Dict{String,Any}(
                    "bus" => "hv", "terminal_map" => ["a"], "configuration" => "WYE",
                    "v_magnitude" => [6350.0], "v_angle" => [0.0])),
                "transformer" => Dict{String,Any}("single_phase" => Dict{String,Any}(
                    "t1" => Dict{String,Any}(
                        "bus_from" => "hv", "bus_to" => "lv",
                        "terminal_map_from" => ["a", "n"], "terminal_map_to" => ["a", "n"],
                        "v_nom_from" => 6350.0, "v_nom_to" => 230.0, "s_nom" => 50e3,
                        "r_series_from" => rs, "x_series_from" => rs,
                        "r_series_to" => 0.0, "x_series_to" => 0.0,
                        "r_neutral_to" => 2.0, "x_neutral_to" => 0.0,
                        "g_no_load" => 1e-4, "b_no_load" => -3e-4))),
                "line" => Dict{String,Any}("l1" => _line("lv", "ld", ["a", "n"], 0.2, 0.1)),
                "shunt" => Dict{String,Any}("e" => Dict{String,Any}(
                    "bus" => "ld", "terminal_map" => ["n"], "G_1_1" => 0.5, "B_1_1" => 0.0)),
                "load" => Dict{String,Any}("d1" => Dict{String,Any}(
                    "bus" => "ld", "terminal_map" => ["a", "n"], "configuration" => "WYE",
                    "p_nom" => [10e3], "q_nom" => [2e3])))
        end
        hi = helm_series(xnet(0.0))                     # promoted to a coupling
        @test hi.status == :converged
        @test length(hi.couplings) == 1
        # Earth KCL: current into earth through the neutral resistor returns
        # through the load electrode.
        @test abs(hi.V[("lv", "n")] / 2.0 + 0.5 * hi.V[("ld", "n")]) < 1e-9
        @test abs(hi.V[("ld", "n")]) > 0.1
        # Same unit with a tiny (stamped, non-ideal) leakage: same answer.
        hs = helm_series(xnet(1e-5))
        @test isempty(hs.couplings)
        for nd in (("lv", "a"), ("lv", "n"), ("ld", "a"), ("ld", "n"))
            @test hi.V[nd] ≈ hs.V[nd] atol=1e-4
        end
    end

    @testset "n-wire: 5-conductor feeder with separate N and PE" begin
        tm5 = ["a", "b", "c", "n", "pe"]
        R = [0.30 0.10 0.10 0.10 0.08; 0.10 0.30 0.10 0.10 0.08;
             0.10 0.10 0.30 0.10 0.08; 0.10 0.10 0.10 0.30 0.08;
             0.08 0.08 0.08 0.08 0.40]
        l5 = Dict{String,Any}("bus_from" => "src", "bus_to" => "ld",
                              "terminal_map_from" => tm5, "terminal_map_to" => tm5)
        for i in 1:5, j in 1:5
            l5["R_series_$(i)_$(j)"] = R[i, j]
            l5["X_series_$(i)_$(j)"] = 0.6 * R[i, j]
        end
        net = Dict{String,Any}(
            "bus" => Dict{String,Any}(
                "src" => Dict{String,Any}("terminal_names" => tm5,
                                          "perfectly_grounded_terminals" => ["n", "pe"]),
                "ld"  => Dict{String,Any}("terminal_names" => tm5)),
            "voltage_source" => Dict{String,Any}("vs" => _src3()),
            "line" => Dict{String,Any}("l1" => l5),
            "shunt" => Dict{String,Any}("pe_earth" => Dict{String,Any}(
                "bus" => "ld", "terminal_map" => ["pe"], "G_1_1" => 1.0, "B_1_1" => 0.0)),
            "load" => Dict{String,Any}(
                "d1" => Dict{String,Any}("bus" => "ld", "terminal_map" => ["a", "b", "c", "n"],
                    "configuration" => "WYE", "model" => "constant_current",
                    "v_nom" => [230.0], "p_nom" => [7000.0, 3000.0, 1000.0],
                    "q_nom" => [1000.0, 500.0, 0.0]),
                "leak" => Dict{String,Any}("bus" => "ld", "terminal_map" => ["c", "pe"],
                    "configuration" => "SINGLE_PHASE", "model" => "constant_impedance",
                    "v_nom" => [230.0], "p_nom" => [500.0], "q_nom" => [0.0])))
        hr = helm_series(net)
        @test hr.status == :converged
        @test _lin_res(net, hr) < 1e-8
        # Mutual coupling + the c–PE leakage lift PE off earth at the load
        # end, distinct from the neutral.
        @test abs(hr.V[("ld", "pe")]) > 0.1
        @test abs(hr.V[("ld", "n")] - hr.V[("ld", "pe")]) > 0.1
    end

    @testset "solve_pf_helm reports source injection and options" begin
        res = solve_pf_helm(_feeder(5000.0; model="constant_current");
                            embedding=:impedance_germ)
        @test res["termination_status"] == "HELM_CONVERGED"
        @test res["helm"]["embedding"] == "impedance_germ"
        @test res["helm"]["source_reference"] == "earth"
        dv = res["bus"]["ld"]["a"]["vr"] - res["bus"]["ld"]["n"]["vr"]
        @test res["source_injection"]["src"]["a"]["im"] ≈ 5000.0 * dv / 230.0 / dv rtol=1e-9
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

        # ZIP deck with constant-current fractions: the reflected-modulus
        # recursion carries |ΔV| holomorphically, so it is solved, not refused.
        _helm_vs_ods("pf_zip_1ph.dss", "src"; atol=1e-3)
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
