using Test, PowerOptLab, LinearAlgebra, JuMP, Ipopt
using BMOPFTools: parse_bmopf

@testset "SE regression: merged source and measurement aliases" begin
    function alias_estimation_net(; merged=false)
        net = parse_bmopf("""
        {"bus":{"src":{"terminal_names":["1"]},
                "b1":{"terminal_names":["1"]},"b2":{"terminal_names":["1"]}},
         "voltage_source":{"s":{"bus":"src","terminal_map":["1"],
             "v_magnitude":[230.0],"v_angle":[0.3]}},
         "linecode":{"lc":{"R_series_1_1":0.5}},
         "line":{"l1":{"bus_from":"src","bus_to":"b1","terminal_map_from":["1"],
             "terminal_map_to":["1"],"linecode":"lc","length":1.0},
                 "l2":{"bus_from":"b1","bus_to":"b2","terminal_map_from":["1"],
             "terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
        """; from_string=true)
        if merged
            # The source loses its representative name; the measured load alias
            # is not a representative either. Both must resolve through Ybus.index.
            for (alias, original) in (("aaa", "src"), ("zzz", "b2"))
                net["bus"][alias] = Dict{String,Any}("terminal_names" => ["1"])
                net["line"][alias] = Dict("bus_from" => original, "bus_to" => alias,
                    "terminal_map_from" => ["1"], "terminal_map_to" => ["1"],
                    "linecode" => "lc", "length" => 0.0)
            end
        end
        net
    end
    vs = 230cis(0.3)
    v1, v2 = vs - (3 + 2im), vs - (5 + 3im)
    powers = (v1 * conj((v1 - vs)/.5 + (v1 - v2)/.5),
              v2 * conj((v2 - v1)/.5))
    diagnostics = []
    for merged in (false, true)
        loadbus = merged ? "zzz" : "b2"
        ms = [Measurement(kind=:vmag, bus="src", reference=nothing, value=230., sigma=.5)]
        for (bus, voltage, power) in (("b1", v1, powers[1]), (loadbus, v2, powers[2]))
            append!(ms, [Measurement(kind=:vmag, bus=bus, reference=nothing, value=abs(voltage), sigma=.5),
                         Measurement(kind=:pinj, bus=bus, reference=nothing, value=real(power), sigma=10.),
                         Measurement(kind=:qinj, bus=bus, reference=nothing, value=imag(power), sigma=10.)])
        end
        # Supply only original/measurement names, so state extraction must also
        # find the voltage when the canonical representative has no entry.
        est = Dict(bus => Dict("1" => Dict("vr" => real(v), "vi" => imag(v)))
                   for (bus, v) in (("src", vs), ("b1", v1), (loadbus, v2)))
        d = PowerOptLab._observability(alias_estimation_net(; merged), ms, nothing,
                                       Set{Tuple{String,String}}(), est)
        push!(diagnostics, d)
        @test d.n_states == 4
        @test d.rank == 4
        @test d.observable
        @test d.tangent_dimension == 4
    end
    @test diagnostics[1].redundancy == diagnostics[2].redundancy
    @test diagnostics[1].min_singular ≈ diagnostics[2].min_singular rtol=1e-6
    @test diagnostics[1].cond ≈ diagnostics[2].cond rtol=1e-6

    grounded = alias_estimation_net(; merged=true)
    push!(grounded["bus"]["src"]["terminal_names"], "n")
    grounded["bus"]["src"]["perfectly_grounded_terminals"] = ["n"]
    source = grounded["voltage_source"]["s"]
    push!(source["terminal_map"], "n")
    push!(source["v_magnitude"], 0.); push!(source["v_angle"], 0.)
    d = PowerOptLab._observability(grounded, Measurement[], nothing,
                                  Set{Tuple{String,String}}(), Dict())
    @test d.n_states == d.rank == 4
    @test d.observable
end

function return_current_test_net(; phases=["1"], angles=[0.])
    terminals = vcat(phases, ["n"])
    lc = Dict{String,Any}("R_series_$(k)_$(k)" => .5 for k in eachindex(terminals))
    net = Dict{String,Any}(
        "bus" => Dict("src" => Dict("terminal_names" => terminals,
                                    "perfectly_grounded_terminals" => ["n"]),
                      "b" => Dict("terminal_names" => terminals)),
        "voltage_source" => Dict("s" => Dict("bus" => "src", "terminal_map" => phases,
            "v_magnitude" => fill(230., length(phases)), "v_angle" => angles)),
        "linecode" => Dict("lc" => lc),
        "line" => Dict("l" => Dict("bus_from" => "src", "bus_to" => "b",
            "terminal_map_from" => terminals, "terminal_map_to" => terminals,
            "linecode" => "lc", "length" => 1.)))
    parse_bmopf(String(PowerOptLab.JSON3.write(net)); from_string=true)
end

@testset "SE regression: injection starts respect complex reference voltage" begin
    for sb in (1., 1e6), vb in (1., 230.), angle in (.7, -2.1)
        # Test the production seeding routine on JuMP voltage objects, with a
        # displaced custom return. Check complex power, not individual formula terms.
        model = Model()
        vp = 230cis(angle) / vb
        vn = (3 - 2im) / vb
        @variable(model, pr, start=real(vp)); @variable(model, pi_, start=imag(vp))
        @variable(model, nr, start=real(vn)); @variable(model, ni, start=imag(vn))
        vr = Dict(("b", "phase") => pr, ("b", "return") => nr)
        vi = Dict(("b", "phase") => pi_, ("b", "return") => ni)
        power = -2000 - 400im
        current = PowerOptLab._se_injection_start(vr, vi, ("b", "phase"), "return",
                                                 real(power), imag(power), sb, 230/vb)
        @test sb * (vp - vn) * conj(current) ≈ power
        earth = PowerOptLab._se_injection_start(vr, vi, ("b", "phase"), nothing,
                                               real(power), imag(power), sb, 230/vb)
        @test sb * vp * conj(earth) ≈ power
    end
    for angle in (.7, -2.1), per_unit in (false, true)
        net = return_current_test_net(; angles=[angle])
        current = (-10 + 2im) * cis(angle)
        vp, vn = 230cis(angle) + .5current, -.5current
        power = (vp - vn) * conj(current)
        ms = [Measurement(kind=:vmag, bus="b", value=abs(vp-vn), sigma=.1),
              Measurement(kind=:pinj, bus="b", value=real(power), sigma=10.),
              Measurement(kind=:qinj, bus="b", value=imag(power), sigma=10.)]
        r = solve_state_estimation(net, ms; per_unit)
        @test r.primal_status == "FEASIBLE_POINT"
        @test r.objective < 1e-6
        @test complex(r.bus["b"]["1"]["vr"], r.bus["b"]["1"]["vi"]) ≈ vp atol=1e-4
        @test complex(r.bus["b"]["n"]["vr"], r.bus["b"]["n"]["vi"]) ≈ vn atol=1e-4
    end
end

function return_current_test_state(s, voltage)
    nf = length(s.free_state_map)
    x = zeros(2nf)
    for (node, k) in s.free_state_map
        x[k], x[nf+k] = reim(voltage[node])
    end
    x
end

@testset "SE regression: stochastic customer neutral return" begin
    net = return_current_test_net()
    # Independent Ohm's-law construction, including a displaced neutral.
    current = -10 + 2im
    vp, vn = 230 + .5current, -.5current
    power = (vp - vn) * conj(current)
    ms = [Measurement(kind=:pinj, bus="b", value=real(power), sigma=10.),
          Measurement(kind=:qinj, bus="b", value=imag(power), sigma=10.),
          Measurement(kind=:vmag, bus="b", value=abs(vp-vn) + .3, sigma=.1)]
    s = compile_state_estimator(net, ms)
    p = SEParameters(s, ms)
    x = return_current_test_state(s, Dict(("b", "1") => vp, ("b", "n") => vn))
    @test length(evaluate_state_estimator(s, p, x).constraints) == 2
    @test norm(evaluate_state_estimator(s, p, x).constraints) < 1e-10
    bad = copy(x); bad[s.free_state_map[("b", "n")]] += 1.
    @test norm(evaluate_state_estimator(s, p, bad).constraints) > 1.
    # Readings remain stochastic: changing their values cannot change physics.
    p.measurement_values[1] += 100.
    @test evaluate_state_estimator(s, p, bad).constraints ≈ evaluate_state_estimator(s, SEParameters(s, ms), bad).constraints
    p = SEParameters(s, ms)
    reference = solve_state_estimation(net, ms; per_unit=false)
    @test reference.primal_status == "FEASIBLE_POINT"
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator)
        r = solver(s, p, initial_state_estimator(s, p))
        @test r.status == :converged_unique
        @test norm(r.evaluation.constraints) < 1e-7
        ip = r.evaluation.current[s.node_index[("b", "1")]]
        inn = r.evaluation.current[s.node_index[("b", "n")]]
        @test ip + inn ≈ 0 atol=1e-7
        @test norm(r.evaluation.residual) > .01
        @test sum(abs2, r.evaluation.residual) ≈ reference.objective atol=1e-5 rtol=1e-5
        for terminal in ("1", "n")
            actual = r.evaluation.voltage[s.node_index[("b", terminal)]]
            expected = reference.bus["b"][terminal]
            @test actual ≈ complex(expected["vr"], expected["vi"]) atol=1e-4
        end
    end
end

@testset "SE regression: shared neutral and explicit earth returns" begin
    net = return_current_test_net(; phases=["1", "2"], angles=[.2, -2.])
    currents = [-10+2im, -4-3im]
    vn = -.5sum(currents)
    voltage = Dict(("b", "n") => vn)
    ms = Measurement[]
    for (k, (phase, angle)) in enumerate(zip(["1", "2"], [.2, -2.]))
        vp = 230cis(angle) + .5currents[k]
        voltage[("b", phase)] = vp
        power = (vp-vn) * conj(currents[k])
        for (kind, value) in ((:pinj, real(power)), (:qinj, imag(power)))
            push!(ms, Measurement(; kind, bus="b", terminal=phase, value, sigma=10.))
        end
    end
    s = compile_state_estimator(net, ms); p = SEParameters(s, ms)
    x = return_current_test_state(s, voltage)
    e = evaluate_state_estimator(s, p, x)
    @test length(e.constraints) == 2  # one sum, not one neutral equation per phase
    @test norm(e.constraints) < 1e-10
    bad = copy(x); bad[s.free_state_map[("b", "n")]] += .5
    @test norm(evaluate_state_estimator(s, p, bad).constraints) > .5
    C = constraint_jacobian(s, p, x)
    h = 1e-4
    basis = Matrix{Float64}(I, length(x), length(x))
    fd = hcat([(evaluate_state_estimator(s, p, x + h*basis[:, k]).constraints -
                evaluate_state_estimator(s, p, x - h*basis[:, k]).constraints)/(2h)
               for k in 1:length(x)]...)
    @test C ≈ fd atol=1e-7
    earth = [Measurement(kind=m.kind, bus=m.bus, terminal=m.terminal, reference=nothing,
                         value=m.value, sigma=m.sigma) for m in ms]
    se = compile_state_estimator(net, earth)
    @test isempty(evaluate_state_estimator(se, SEParameters(se, earth), x).constraints)

    grounded = deepcopy(net)
    grounded["bus"]["b"]["perfectly_grounded_terminals"] = ["n"]
    sg = compile_state_estimator(grounded, ms)
    @test size(constraint_jacobian(sg, SEParameters(sg, ms), initial_state_estimator(sg, SEParameters(sg, ms))), 1) == 0

    vm = [Measurement(kind=:vmag, bus="b", value=230., sigma=1.)]
    sv = compile_state_estimator(net, vm)
    @test isempty(evaluate_state_estimator(sv, SEParameters(sv, vm), x).constraints)

    # A measured phase-to-phase injection uses the same conservation law.
    line_to_line = [Measurement(kind=:pinj, bus="b", terminal="1", reference="2", value=-1000., sigma=10.)]
    sl = compile_state_estimator(net, line_to_line)
    pl = SEParameters(sl, line_to_line)
    vl = Dict(("b", "1") => 230cis(.2) + .5currents[1],
              ("b", "2") => 230cis(-2.) - .5currents[1], ("b", "n") => 0im)
    xl = return_current_test_state(sl, vl)
    @test length(evaluate_state_estimator(sl, pl, xl).constraints) == 2
    @test norm(evaluate_state_estimator(sl, pl, xl).constraints) < 1e-10

    # Repeated P/Q readings and original bus aliases must not duplicate or
    # split the shared-neutral component.
    aliased = deepcopy(net)
    aliased["bus"]["z"] = Dict{String,Any}("terminal_names" => ["1", "2", "n"])
    aliased["switch"] = Dict("sw" => Dict("bus_from" => "b", "bus_to" => "z",
        "terminal_map_from" => ["1", "2", "n"], "terminal_map_to" => ["1", "2", "n"], "status" => "closed"))
    aliases = vcat(ms, [Measurement(kind=m.kind, bus="z", terminal=m.terminal,
        value=m.value, sigma=m.sigma) for m in ms])
    sa = compile_state_estimator(aliased, aliases;formulation=:voltage)
    pa = SEParameters(sa, aliases)
    @test length(evaluate_state_estimator(sa, pa, x).constraints) == 2
    @test norm(evaluate_state_estimator(sa, pa, x).constraints) < 1e-10
end
