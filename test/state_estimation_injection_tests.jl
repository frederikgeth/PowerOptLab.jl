using Test, PowerOptLab, LinearAlgebra
using BMOPFTools: parse_bmopf

function merged_injection_net(; merged=true, grounded=true)
    bond = grounded ? ",\"perfectly_grounded_terminals\":[\"n\"]" : ""
    net = parse_bmopf("""
    {"bus":{"src":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
            "b":{"terminal_names":["1","n"]$bond}},
     "voltage_source":{"s":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.5,"R_series_2_2":0.5}},
     "line":{"l":{"bus_from":"src","bus_to":"b","terminal_map_from":["1","n"],
                    "terminal_map_to":["1","n"],"linecode":"lc","length":1.0}}}
    """; from_string=true)
    if merged
        net["bus"]["z"] = deepcopy(net["bus"]["b"])
        net["line"]["alias"] = Dict{String,Any}("bus_from"=>"b", "bus_to"=>"z",
            "terminal_map_from"=>["1","n"], "terminal_map_to"=>["1","n"], "linecode"=>"lc", "length"=>0.)
    end
    net
end

function injection_fd(f, x)
    h = 1e-4
    basis = Matrix{Float64}(I, length(x), length(x))
    hcat([(f(x+h*basis[:,k])-f(x-h*basis[:,k]))/(2h) for k in eachindex(x)]...)
end

@testset "DSSE merged-power profiling requires nonzero connection voltage" begin
    net = merged_injection_net()
    # Net power is zero but customer powers are not. At zero voltage finite
    # individual currents cannot reproduce these readings; a sum-only profile
    # must not publish an exact fit by implicitly using infinite currents.
    meters = [Measurement(kind=kind, bus=bus, reference=nothing, value=value, sigma=10.)
        for (bus, p, q) in (("b", -1000., -200.), ("z", 1000., 200.))
        for (kind, value) in ((:pinj, p), (:qinj, q))]
    s = compile_state_estimator(net, meters; formulation=:voltage); p = SEParameters(s, meters)
    @test_throws DomainError evaluate_state_estimator(s, p, [0., 0.])
    @test_throws DomainError residual_jacobian(s, p, [0., 0.])
    @test_throws DomainError evaluate_state_estimator(s, p, [p.voltage_min_model/2, 0.])
    @test all(isfinite, evaluate_state_estimator(s, p, [2p.voltage_min_model, 0.]).residual)
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator)
        r = solver(s, p, [0., 0.])
        @test r.status == :invalid_initial_domain
        @test !solve_status(r).publishable
    end
end

@testset "DSSE merged customers: sum powers and propagate uncertainty" begin
    net = merged_injection_net()
    # Independent two-customer total: Q = -460 Vi and
    # P = 2(Vr²+Vi²)-460 Vr on this 230 V, 0.5 ohm feeder.
    vi = 500/460
    vr = (230 + sqrt(230^2 - 4*(1250+vi^2)))/2
    ms = Measurement[]
    for (bus, power, sigma) in (("b", -1000-200im, 10.), ("z", -1500-300im, 20.))
        push!(ms, Measurement(kind=:pinj, bus=bus, reference=nothing, value=real(power), sigma=sigma),
                  Measurement(kind=:qinj, bus=bus, reference=nothing, value=imag(power), sigma=sigma))
    end
    s = compile_state_estimator(net, ms; formulation=:voltage); p = SEParameters(s, ms)
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator)
        r = solver(s, p, initial_state_estimator(s,p))
        @test r.status == :converged_unique
        @test r.evaluation.voltage[s.node_index[("b","1")]] ≈ complex(vr,vi) atol=1e-6
        @test norm(r.evaluation.residual) < 1e-6
        @test r.evaluation.predicted ≈ p.measurement_values atol=1e-5
    end
    aggregate = [Measurement(kind=:pinj, bus="b", reference=nothing, value=-2500., sigma=sqrt(500.)),
                 Measurement(kind=:qinj, bus="b", reference=nothing, value=-500., sigma=sqrt(500.))]
    legacy = solve_state_estimation(net, ms; per_unit=false)
    ag = solve_state_estimation(merged_injection_net(;merged=false), aggregate; per_unit=false)
    @test legacy.objective < 1e-5
    @test legacy.observability.redundancy == ag.observability.redundancy == 0
    @test legacy.observability.min_singular ≈ ag.observability.min_singular rtol=1e-6
    sa = compile_state_estimator(merged_injection_net(;merged=false), aggregate; formulation=:voltage)
    pa = SEParameters(sa, aggregate)
    x = [vr,vi]
    @test selected_state_covariance(s,p,x,[1,2]) ≈ selected_state_covariance(sa,pa,x,[1,2]) rtol=1e-6
    @test residual_jacobian(s,p,x) ≈ injection_fd(y->evaluate_state_estimator(s,p,y).residual,x) atol=1e-7

    # A noisy voltage reading and a repeated customer reading test profiling's
    # nonzero residuals, rather than only an interpolated clean solution.
    push!(ms, Measurement(kind=:vmag,bus="b",reference=nothing,value=abs(complex(vr,vi))+.2,sigma=.1),
              Measurement(kind=:pinj,bus="b",reference=nothing,value=-980.,sigma=10.))
    s = compile_state_estimator(net,ms; formulation=:voltage); p = SEParameters(s,ms)
    legacy = solve_state_estimation(net,ms;per_unit=false)
    for solver in (solve_compiled_state_estimator,solve_sparse_state_estimator,solve_ipopt_state_estimator)
        r=solver(s,p,initial_state_estimator(s,p))
        @test r.status == :converged_unique
        @test sum(abs2,r.evaluation.residual) ≈ legacy.objective atol=1e-5 rtol=1e-5
        v=legacy.bus["b"]["1"]
        @test r.evaluation.voltage[s.node_index[("b","1")]] ≈ complex(v["vr"],v["vi"]) atol=1e-5
        @test r.evaluation.predicted[1] ≈ r.evaluation.predicted[end]
        @test residual_jacobian(s,p,r.state) ≈ injection_fd(y->evaluate_state_estimator(s,p,y).residual,r.state) atol=1e-7
    end
    # Reordering customers and repeats changes neither the model nor the fit.
    sp=compile_state_estimator(net,reverse(ms); formulation=:voltage); pp=SEParameters(sp,reverse(ms))
    rp=solve_sparse_state_estimator(sp,pp,initial_state_estimator(sp,pp))
    @test sum(abs2,rp.evaluation.residual) ≈ legacy.objective atol=1e-5 rtol=1e-5

    # Mutable uncertainties must recompute the profile, not reuse stale weights.
    p.covariance_values[1] *= 2
    @test residual_jacobian(s,p,x) ≈ injection_fd(y->evaluate_state_estimator(s,p,y).residual,x) atol=1e-7
    incomplete = ms[[1,2,3]] # z's Q is unobserved; b's Q cannot identify total Q.
    si=compile_state_estimator(net,incomplete; formulation=:voltage); pi_=SEParameters(si,incomplete)
    @test norm(residual_jacobian(si,pi_,x)[2,:]) == 0
    @test evaluate_state_estimator(si,pi_,x).residual[2] == 0
end

@testset "DSSE references: reject conflicting returns independent of row order" begin
    net = merged_injection_net(;merged=false,grounded=false)
    ms = [Measurement(kind=:pinj,bus="b",reference="n",value=-2196.,sigma=10.),
          Measurement(kind=:qinj,bus="b",reference=nothing,value=-460.,sigma=10.)]
    for readings in (ms,reverse(ms))
        @test_throws ArgumentError solve_state_estimation(net,readings)
        @test !state_estimator_preflight(net,readings).supported
        @test_throws SEUnsupportedNetwork compile_state_estimator(net,readings; formulation=:voltage)
    end
    # Explicit earth and a named perfectly grounded neutral are equivalent.
    ground = merged_injection_net(;merged=false)
    equivalent = [ms[1],ms[2]]
    a=solve_state_estimation(ground,equivalent)
    b=solve_state_estimation(ground,reverse(equivalent))
    @test a.objective ≈ b.objective atol=1e-5
    @test a.bus["b"]["1"]["vm"] ≈ b.bus["b"]["1"]["vm"] atol=1e-6
    # Compiled preflight uses electrical equivalence too, not literal strings.
    for readings in (equivalent,reverse(equivalent))
        @test state_estimator_preflight(ground,readings).supported
    end

    # Different returns at a merged positive node cannot be disaggregated from
    # net Ybus current. Reject instead of giving each customer that same current.
    merged=merged_injection_net(;grounded=false)
    mixed=[Measurement(kind=:pinj,bus="b",reference="n",value=-1000.,sigma=10.),
           Measurement(kind=:qinj,bus="b",reference="n",value=-200.,sigma=10.),
           Measurement(kind=:pinj,bus="z",reference=nothing,value=-1500.,sigma=10.),
           Measurement(kind=:qinj,bus="z",reference=nothing,value=-300.,sigma=10.)]
    for readings in (mixed,reverse(mixed))
        @test state_estimator_preflight(merged,readings).supported
        @test solve_status(solve_state_estimation(merged,readings)).publishable
        augmented=compile_state_estimator(merged,readings)
        @test length(injection_state_indices(augmented))==2
    end
    overlap=[TerminalConnection(("b","1"),("b","n")),
             TerminalConnection(("b","n"),nothing)]
    @test state_estimator_preflight(net;injection_connections=overlap).supported
    @test_throws SEUnsupportedNetwork compile_state_estimator(net;injection_connections=overlap, formulation=:voltage)
    @test !state_estimator_preflight(net;injection_connections=[TerminalConnection(("b","1"),("b","absent"))]).supported
    @test !state_estimator_preflight(net;injection_connections=[TerminalConnection(("b","1"),("b","1"))]).supported
end

@testset "DSSE merged customers with an ungrounded return" begin
    net=merged_injection_net(;grounded=false)
    meters=[Measurement(kind=kind,bus=bus,reference="n",value=value,sigma=10.)
            for (bus,p,q) in (("b",-1000.,-200.),("z",-1500.,-300.))
            for (kind,value) in ((:pinj,p),(:qinj,q))]
    # Total phase+neutral series resistance is 1 ohm. Solve S=U*conj(U-230).
    ui=500/230
    ur=(230+sqrt(230^2-4*(2500+ui^2)))/2
    u=complex(ur,ui)
    phase=(230+u)/2; neutral=(230-u)/2
    s=compile_state_estimator(net,meters; formulation=:voltage); p=SEParameters(s,meters)
    legacy=solve_state_estimation(net,meters;per_unit=false)
    @test legacy.observability.observable
    @test legacy.observability.redundancy == 0
    for solver in (solve_compiled_state_estimator,solve_sparse_state_estimator,solve_ipopt_state_estimator)
        r=solver(s,p,initial_state_estimator(s,p))
        @test r.status == :converged_unique
        @test r.evaluation.voltage[s.node_index[("b","1")]] ≈ phase atol=1e-6
        @test r.evaluation.voltage[s.node_index[("b","n")]] ≈ neutral atol=1e-6
        @test norm(r.evaluation.constraints) < 1e-7
        @test norm(r.evaluation.residual) < 1e-6
        @test legacy.bus["b"]["1"]["vr"] ≈ real(phase) atol=1e-6
    end
end

@testset "DSSE connection declarations survive measurement dropout" begin
    net=merged_injection_net(;merged=false,grounded=false)
    connections=[TerminalConnection(("b","1"),("b","n"))]
    voltages=[Measurement(kind=:vr,bus="b",value=220.,sigma=.2),
              Measurement(kind=:vi,bus="b",value=2.,sigma=.2)]
    full=vcat(voltages,[Measurement(kind=:pinj,bus="b",value=-2196.,sigma=10.),
                        Measurement(kind=:qinj,bus="b",value=-460.,sigma=10.)])
    sf=compile_state_estimator(net,full;injection_connections=connections, formulation=:voltage)
    sv=compile_state_estimator(net,voltages;injection_connections=connections, formulation=:voltage)
    pf=SEParameters(sf,full); pv=SEParameters(sv,voltages)
    x=initial_state_estimator(sv,pv)
    @test constraint_jacobian(sf,pf,x) == constraint_jacobian(sv,pv,x)
    @test evaluate_state_estimator(sf,pf,x).constraints == evaluate_state_estimator(sv,pv,x).constraints
    @test observability_diagnostics(sv,pv,x).unobservable_dimension == 0
    for solver in (solve_compiled_state_estimator,solve_sparse_state_estimator,solve_ipopt_state_estimator)
        r=solver(sv,pv,x)
        @test r.status == :converged_unique
        @test r.evaluation.voltage[sv.node_index[("b","n")]] ≈ 5-im atol=1e-5
    end
    @test !state_estimator_preflight(net,full;injection_connections=TerminalConnection[]).supported
    # Declaring a second customer without readings leaves its power unknown.
    merged=merged_injection_net()
    cs=[TerminalConnection((bus,"1"),nothing) for bus in ("b","z")]
    only_b=[Measurement(kind=:pinj,bus="b",reference=nothing,value=-1000.,sigma=10.),
            Measurement(kind=:qinj,bus="b",reference=nothing,value=-200.,sigma=10.)]
    sm=compile_state_estimator(merged,only_b;injection_connections=cs, formulation=:voltage); pm=SEParameters(sm,only_b)
    @test norm(residual_jacobian(sm,pm,initial_state_estimator(sm,pm))) == 0
end

@testset "DSSE grounded neutral telemetry and diagnostic conventions" begin
    net=merged_injection_net(;merged=false,grounded=false)
    meters=[BranchMeasurement(kind=kind,line="l",side=:from,terminal="n",value=value,sigma=.1)
            for (kind,value) in ((:ire,-10.),(:iim,2.),(:imag,sqrt(104.)),(:pflow,0.),(:qflow,0.))]
    s=compile_state_estimator(net,meters; formulation=:voltage); p=SEParameters(s,meters)
    x=zeros(2length(s.free_state_map)); nf=length(s.free_state_map)
    for (terminal,v) in (("1",225+im),("n",5-im))
        k=s.free_state_map[("b",terminal)]; x[k]=real(v); x[nf+k]=imag(v)
    end
    @test norm(evaluate_state_estimator(s,p,x).residual) < 1e-9
    @test residual_jacobian(s,p,x) ≈ injection_fd(y->evaluate_state_estimator(s,p,y).residual,x) atol=1e-7
    device=ExactDeviceEquation(ConstantCurrentDevice([TerminalConnection(("b","1"),nothing)],ComplexF64[2]))
    ss=compile_state_estimator(merged_injection_net(;merged=false);exact_devices=[device], formulation=:voltage)
    ps=SEParameters(ss;exact_devices=[device])
    @test observability_diagnostics(ss,ps,[229.,0.]).condition_number == 1.
    series=solve_time_series_state_estimator(ss,[ps,deepcopy(ps),deepcopy(ps)],[230.,0.];max_iterations=1,initial_radius=1e-5)
    @test series.status == :time_series_stalled
    @test solve_diagnostics(series).requested_snapshots == 3
    @test series.completed_snapshots == 0
end
