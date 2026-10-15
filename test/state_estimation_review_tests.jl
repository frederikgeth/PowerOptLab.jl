using Test, PowerOptLab, LinearAlgebra, SparseArrays, Random, Logging
import BMOPFTools

const REVIEW_SE_SOLVERS=(
    solve_compiled_state_estimator,
    solve_sparse_state_estimator,
    solve_ipopt_state_estimator,
)

# Independent four-wire nodal Ohm's-law construction. No package admittance,
# device evaluator or power-flow solver generates the reference voltages/data.
function review_feeder(n; seed = 143)
    rng=MersenneTwister(seed)
    terms=["a", "b", "c", "n"]
    z=fill(0.03im, 4, 4)+Diagonal(fill(0.2+0.05im, 4))
    y=inv(0.05z)
    source=ComplexF64[230cis(a) for a in (0.0, -2pi/3, 2pi/3)]
    push!(source, 0.0)
    A=zeros(ComplexF64, 4n, 4n)
    rhs=zeros(ComplexF64, 4n)
    amps=ComplexF64[(3+4rand(rng)-im*(0.5+rand(rng)))*cis(-2pi*(k-1)/3) for b = 1:n, k = 1:3]
    bus=Dict{String,Any}("b0"=>Dict("terminal_names"=>terms, "perfectly_grounded_terminals"=>["n"]))
    lines=Dict{String,Any}()
    shunts=Dict{String,Any}()
    for b = 1:n
        rows=(4b-3):4b
        A[rows, rows]+=y
        if b==1
            rhs[rows]+=y*source
        else
            previous=(4b-7):(4b-4)
            A[previous, previous]+=y
            A[rows, previous]-=y
            A[previous, rows]-=y
        end
        rhs[rows]-=vcat(amps[b, :], -sum(amps[b, :]))
        if b%5==0
            A[4b, 4b]+=0.05
            shunts["g$b"]=Dict("bus"=>"b$b", "terminal_map"=>["n"], "G_1_1"=>0.05)
        end
        bus["b$b"]=Dict("terminal_names"=>terms)
        lines["l$b"]=Dict(
            "bus_from"=>"b$(b-1)",
            "bus_to"=>"b$b",
            "terminal_map_from"=>terms,
            "terminal_map_to"=>terms,
            "linecode"=>"lc",
            "length"=>0.05,
        )
    end
    code=Dict{String,Any}()
    for i = 1:4, j = 1:4
        code["R_series_$(i)_$(j)"]=real(z[i, j])
        code["X_series_$(i)_$(j)"]=imag(z[i, j])
    end
    net=Dict{String,Any}(
        "bus"=>bus,
        "line"=>lines,
        "linecode"=>Dict("lc"=>code),
        "shunt"=>shunts,
        "voltage_source"=>Dict(
            "s"=>Dict(
                "bus"=>"b0",
                "terminal_map"=>terms[1:3],
                "v_magnitude"=>fill(230.0, 3),
                "v_angle"=>[0.0, -2pi/3, 2pi/3],
            ),
        ),
    )
    volts=reshape(A\rhs, 4, n)
    meters=Measurement[]
    for b = 1:n, k = 1:3
        u=volts[k, b]-volts[4, b]
        power=u*conj(-amps[b, k])
        for (kind, value, sigma) in
            ((:pinj, real(power), 10.0), (:qinj, imag(power), 10.0), (:vmag, abs(u), 0.5))
            push!(
                meters,
                Measurement(;
                    kind,
                    bus = "b$b",
                    terminal = terms[k],
                    reference = "n",
                    value = value+sigma*randn(rng),
                    sigma,
                ),
            )
        end
    end
    (; net, meters, volts, amps)
end

@testset "Review: fully fixed voltages retain residual degrees of freedom" begin
    net=Dict{String,Any}(
        "bus"=>Dict("s"=>Dict("terminal_names"=>["p", "n"], "perfectly_grounded_terminals"=>["n"])),
        "voltage_source"=>Dict(
            "source"=>Dict(
                "bus"=>"s",
                "terminal_map"=>["p"],
                "v_magnitude"=>[230.0],
                "v_angle"=>[0.0],
            ),
        ),
    )
    meters=[
        Measurement(; kind, bus = "s", terminal = "p", reference = nothing, value, sigma = 1.0) for
        (kind, value) in ((:vr, 232.0), (:vi, 3.0))
    ]
    native=solve_state_estimation(net, meters; per_unit = false)
    s=compile_state_estimator(net, meters)
    p=SEParameters(s, meters)
    results=[native; [f(s, p, Float64[]) for f in REVIEW_SE_SOLVERS]]
    for r in results
        @test solve_status(r).publishable
        @test solve_diagnostics(r).objective ≈ 13.0
        @test solve_diagnostics(r).dof == 2
        @test solve_diagnostics(r).nuisance_dimension == 0
    end
    @test native.observability.redundancy == 2
end

@testset "Review: fixed voltages retain declared current states in the initial guess" begin
    net=Dict{String,Any}(
        "bus"=>Dict("s"=>Dict("terminal_names"=>["p", "n"], "perfectly_grounded_terminals"=>["n"])),
        "voltage_source"=>Dict(
            "source"=>Dict(
                "bus"=>"s",
                "terminal_map"=>["p"],
                "v_magnitude"=>[230.0],
                "v_angle"=>[0.0],
            ),
        ),
    )
    connections=[TerminalConnection(("s", "p"), ("s", "n"))]
    s=compile_state_estimator(net; injection_connections = connections, formulation = :current)
    p=SEParameters(s)
    initial=initial_state_estimator(s, p)
    @test isempty(s.free_state_map)
    @test length(initial)==2
    @test initial==zeros(2)
    @test injection_state_indices(s)[("s", "p")] == (real = 1, imag = 2)
    for solver in REVIEW_SE_SOLVERS
        r=solver(s, p, initial)
        @test r.status==:converged_underobserved
        @test r.tangent_dimension==2
        @test r.observable_dimension==0
        @test solve_diagnostics(r).dof==0
    end
    # A real current prior identifies those states without changing fixed V.
    set_state_prior!(p, StatePrior([1, 2], [3.0, -2.0], [0.1, 0.1]))
    for solver in REVIEW_SE_SOLVERS
        r=solver(s, p, initial)
        @test r.status==:converged_unique
        @test r.state ≈ [3.0, -2.0] atol=1e-7
        @test r.evaluation.voltage[s.node_index[("s", "p")]]==230+0im
        @test solve_diagnostics(r).dof==0
    end
end

@testset "Review: conflicting zero and unknown injection declarations fail closed" begin
    case=review_feeder(1)
    meters=filter(m->m.terminal=="a", case.meters)
    connection=TerminalConnection(("b1", "a"), ("b1", "n"))
    for formulation in (:auto, :voltage, :current), terminal in ("a", "n")
        options=(; formulation, zero_injection = [("b1", terminal)])
        report=state_estimator_preflight(case.net, meters; options...)
        @test !report.supported
        @test any(f->f.code==:conflicting_injection_declaration, report.findings)
        @test_throws SEUnsupportedNetwork compile_state_estimator(case.net, meters; options...)
        # The same conflict exists without P/Q telemetry when the path persists.
        @test !state_estimator_preflight(case.net; options..., injection_connections = [connection]).supported
    end
    # Electrical aliases cannot evade the declaration check.
    net=deepcopy(case.net)
    net["bus"]["alias"]=deepcopy(net["bus"]["b1"])
    net["switch"]=Dict(
        "bond"=>Dict(
            "bus_from"=>"alias",
            "bus_to"=>"b1",
            "terminal_map_from"=>["a"],
            "terminal_map_to"=>["a"],
            "status"=>"closed",
        ),
    )
    @test !state_estimator_preflight(net, meters; zero_injection = [("alias", "a")]).supported
    # A known device is not an unknown injection; its exact KCL remains valid.
    exact=[ExactDeviceEquation(ConstantCurrentDevice([connection], ComplexF64[2]))]
    @test state_estimator_preflight(
        case.net,
        meters;
        zero_injection = [("b1", "a")],
        exact_devices = exact,
    ).supported
    @test !state_estimator_preflight(
        case.net,
        meters;
        zero_injection = [("b1", "a")],
        exact_devices = exact,
        injection_connections = [connection],
    ).supported
    net_meter=NetInjectionMeasurement(
        kind = :pinj,
        bus = "b1",
        terminal = "a",
        reference = "n",
        value = 0.0,
        sigma = 10.0,
    )
    @test state_estimator_preflight(case.net, [net_meter]; zero_injection = [("b1", "a")]).supported
end

@testset "Review: loaded 5/20/40-bus four-wire backend agreement" begin
    for n in (5, 20, 40)
        @testset "$n downstream buses" begin
            case=review_feeder(n)
            s=compile_state_estimator(case.net, case.meters)
            p=SEParameters(s, case.meters)
            initial=initial_state_estimator(s, p)
            results=[f(s, p, initial; optimality_tolerance = 1e-6) for f in REVIEW_SE_SOLVERS]
            @test all(r->solve_status(r).publishable, results)
            @test all(r->r.constraint_rank==results[1].constraint_rank, results)
            native=solve_state_estimation(case.net, case.meters; per_unit = false)
            @test solve_status(native).publishable
            for r in results
                d=solve_diagnostics(r)
                @test d.objective≈native.objective atol=1e-6 rtol=1e-7
                @test d.dof==native.observability.redundancy==3n
                @test norm(r.evaluation.constraints)<1e-7
                for (node, k) in s.free_state_map
                    v=complex(
                        native.bus[node[1]][node[2]]["vr"],
                        native.bus[node[1]][node[2]]["vi"],
                    )
                    @test r.evaluation.voltage[s.node_index[node]]≈v atol=1e-5
                end
            end
            @test maximum(abs, results[1].state-results[2].state)<1e-6
            ipopt_diagnostics=solve_diagnostics(results[3])
            @test ipopt_diagnostics.optimizer_jacobian_entries <
                  length(initial)*length(results[3].constraint_rows)
        end
    end
end

@testset "Review: independent acceptance of an iteration-limited primal" begin
    case=review_feeder(1)
    s=compile_state_estimator(case.net, case.meters)
    p=SEParameters(s, case.meters)
    initial=initial_state_estimator(s, p)
    optimum=solve_sparse_state_estimator(s, p, initial)
    options=(max_iter = 0, tol = 1e-16, dual_inf_tol = 1e-16, acceptable_tol = 1e-16)
    accepted=solve_ipopt_state_estimator(s, p, optimum.state; solver_options = options)
    @test accepted.termination_status == "ITERATION_LIMIT"
    @test solve_status(accepted).publishable
    @test accepted.state ≈ optimum.state atol=1e-10
    rejected=solve_ipopt_state_estimator(s, p, initial; solver_options = options)
    @test rejected.termination_status == "ITERATION_LIMIT"
    @test !solve_status(rejected).publishable
end

@testset "Review: overlapping paths and a known device share one neutral" begin
    case=review_feeder(1)
    meters=Measurement[]
    for (terminal, reference, P, Q) in
        (("a", "b", -1000.0, -100.0), ("b", "n", -1000.0, -100.0), ("c", "n", -1000.0, -100.0))
        for (kind, value) in ((:pinj, P), (:qinj, Q))
            push!(meters, Measurement(; kind, bus = "b1", terminal, reference, value, sigma = 10.0))
        end
    end
    push!(
        meters,
        Measurement(
            kind = :vmag,
            bus = "b1",
            terminal = "a",
            reference = "n",
            value = 229.0,
            sigma = 1.0,
        ),
    )
    s=compile_state_estimator(case.net, meters)
    p=SEParameters(s, meters)
    @test length(injection_state_indices(s))==3
    native=solve_state_estimation(case.net, meters; per_unit = false)
    for f in REVIEW_SE_SOLVERS
        r=f(s, p, initial_state_estimator(s, p))
        @test solve_status(r).publishable
        @test solve_diagnostics(r).objective≈native.objective atol=1e-6
        @test solve_diagnostics(r).dof==solve_diagnostics(native).dof==1
    end
    # Keep the source/current sign distinction on the SAME bus: known positive
    # consumption and unknown negative customer power both draw from the source.
    device=ExactDeviceEquation(
        ConstantCurrentDevice([TerminalConnection(("b1", "a"), ("b1", "n"))], ComplexF64[2]),
    )
    others=filter(m->m.terminal!="a", meters)
    s=compile_state_estimator(case.net, others; exact_devices = [device])
    p=SEParameters(s, others; exact_devices = [device])
    for f in REVIEW_SE_SOLVERS
        r=f(s, p, initial_state_estimator(s, p))
        @test solve_status(r).publishable
        @test norm(r.evaluation.constraints)<1e-7
        @test norm(r.evaluation.residual)<1e-5
        @test real(r.evaluation.current[s.node_index[("b1", "a")]])≈-2 atol=1e-7
    end
end

@testset "Review: merged-current formulation is defined at zero voltage" begin
    case=review_feeder(1)
    net=case.net
    net["bus"]["alias"]=deepcopy(net["bus"]["b1"])
    net["switch"]=Dict(
        "bond"=>Dict(
            "bus_from"=>"b1",
            "bus_to"=>"alias",
            "terminal_map_from"=>["a", "b", "c", "n"],
            "terminal_map_to"=>["a", "b", "c", "n"],
            "status"=>"closed",
        ),
    )
    meters=[
        Measurement(; kind, bus, terminal = "a", reference = "n", value, sigma = 10.0) for
        (bus, power) in (("b1", 1000+200im), ("alias", -1000-200im)) for
        (kind, value) in ((:pinj, real(power)), (:qinj, imag(power)))
    ]
    s=compile_state_estimator(net, meters)
    p=SEParameters(s, meters)
    @test length(injection_state_indices(s))==2
    x=zeros(length(s.nominal_state))
    e=evaluate_state_estimator(s, p, x)
    @test all(isfinite, e.residual)
    @test all(isfinite, residual_jacobian(s, p, x))
    @test sum(abs2, e.residual)>1
end

@testset "Review: known and unknown customer on the same phase and neutral" begin
    net=Dict{String,Any}(
        "bus"=>Dict(
            "s"=>Dict("terminal_names"=>["p", "n"], "perfectly_grounded_terminals"=>["n"]),
            "b"=>Dict("terminal_names"=>["p", "n"]),
        ),
        "voltage_source"=>Dict(
            "s"=>Dict("bus"=>"s", "terminal_map"=>["p"], "v_magnitude"=>[230.0], "v_angle"=>[0.0]),
        ),
        "line"=>Dict(
            "l"=>Dict(
                "bus_from"=>"s",
                "bus_to"=>"b",
                "terminal_map_from"=>["p", "n"],
                "terminal_map_to"=>["p", "n"],
                "R_series_1_1"=>0.5,
                "R_series_2_2"=>0.5,
            ),
        ),
    )
    connection=TerminalConnection(("b", "p"), ("b", "n"))
    known=ExactDeviceEquation(ConstantCurrentDevice([connection], ComplexF64[2]))
    meters=[
        Measurement(; kind, bus = "b", terminal = "p", reference = "n", value, sigma = 10.0) for
        (kind, value) in ((:pinj, -1760.0), (:qinj, 0.0))
    ]
    # An explicit path declares an ADDITIONAL unknown customer. Without it,
    # matching telemetry would meter the known device itself.
    s=compile_state_estimator(
        net,
        meters;
        exact_devices = [known],
        injection_connections = [connection],
    )
    p=SEParameters(s, meters; exact_devices = [known])
    for f in REVIEW_SE_SOLVERS
        r=f(s, p, initial_state_estimator(s, p))
        @test solve_status(r).publishable
        @test r.evaluation.voltage[s.node_index[("b", "p")]]≈225 atol=1e-7
        @test r.evaluation.voltage[s.node_index[("b", "n")]]≈5 atol=1e-7
        @test r.evaluation.device_current[s.node_index[("b", "p")]]≈2 atol=1e-10
        @test r.state[injection_state_indices(s)[("b", "p")].real]≈-8 atol=1e-7
    end
    @test_throws SEUnsupportedNetwork compile_state_estimator(
        net,
        meters;
        exact_devices = [known],
        injection_connections = [connection],
        formulation = :voltage,
    )
end

@testset "Review: missing conductor declaration and backend options" begin
    case=review_feeder(1)
    zero=[("b1", t) for t in ("a", "b", "c", "n")]
    voltage=[
        Measurement(
            kind = :vmag,
            bus = "b1",
            terminal = "a",
            reference = "n",
            value = 230.0,
            sigma = 1.0,
        ),
    ]
    native=@test_logs min_level=Logging.Warn solve_state_estimation(
        case.net,
        voltage;
        zero_injection = zero,
    )
    @test solve_status(native).publishable
    earth=[
        Measurement(
            kind = :pinj,
            bus = "b1",
            terminal = "a",
            reference = nothing,
            value = -1000.0,
            sigma = 10.0,
        ),
    ]
    @test any(f->f.code==:undeclared_injection, state_estimator_preflight(case.net, earth).findings)
    s=compile_state_estimator(case.net, case.meters)
    p=SEParameters(s, case.meters)
    series=solve_time_series_state_estimator(
        s,
        [p, deepcopy(p)],
        initial_state_estimator(s, p);
        solver = :ipopt,
        solver_kwargs = (dense = (initial_radius = 0.5,), ipopt = (max_iterations = 30,)),
    )
    @test series.completed_snapshots==2
    @test_throws ArgumentError solve_time_series_state_estimator(
        s,
        [p],
        initial_state_estimator(s, p);
        solver = :ipopt,
        initial_radius = 0.5,
    )
    @test solve_state_estimation(s, p, initial_state_estimator(s, p)) isa
          SparseConstrainedStateEstimationResult
end

@testset "Review: fixed phasor versus voltage-following constant current" begin
    case=review_feeder(1)
    connection=TerminalConnection(("b1", "a"), ("b1", "n"))
    for reference in (:earth, :voltage)
        device=ExactDeviceEquation(
            ConstantCurrentDevice([connection], ComplexF64[3-im]; angle_reference = reference),
        )
        s=compile_state_estimator(case.net; exact_devices = [device])
        p=SEParameters(s; exact_devices = [device])
        x=initial_state_estimator(s, p)
        nf=length(s.free_state_map)
        k=s.free_state_map[("b1", "a")]
        x[k]=200cos(0.4)
        x[nf+k]=200sin(0.4)
        e=evaluate_state_estimator(s, p, x)
        i=e.device_current[s.node_index[("b1", "a")]]
        @test i≈(reference===:earth ? 3-im : (3-im)*cis(0.4)) atol=1e-10
    end
    coefficient=3-im
    drop=(0.02+0.005im)*coefficient
    magnitude=sqrt(230^2-imag(drop)^2)-real(drop)
    theta=-angle(magnitude+drop)
    current=coefficient*cis(theta)
    device=ExactDeviceEquation(
        ConstantCurrentDevice([connection], [coefficient]; angle_reference = :voltage),
    )
    s=compile_state_estimator(
        case.net;
        exact_devices = [device],
        zero_injection = [("b1", "b"), ("b1", "c")],
    )
    p=SEParameters(s; exact_devices = [device])
    for f in REVIEW_SE_SOLVERS
        r=f(s, p, initial_state_estimator(s, p))
        @test solve_status(r).publishable
        a=r.evaluation.voltage[s.node_index[("b1", "a")]]
        n=r.evaluation.voltage[s.node_index[("b1", "n")]]
        @test a-n≈magnitude*cis(theta) atol=1e-7
        @test r.evaluation.device_current[s.node_index[("b1", "a")]]≈current atol=1e-8
    end
    @test_throws ArgumentError ConstantCurrentDevice(
        [connection],
        [coefficient];
        angle_reference = :invalid,
    )
end

@testset "Review: nonlinear merged-meter statistical calibration (500 draws)" begin
    # Two customers on the same conductor, unequal channel variances, a repeat
    # reading, and voltage magnitude make this nonlinear and overdetermined.
    net=Dict{String,Any}(
        "bus"=>Dict(
            "s"=>Dict("terminal_names"=>["p"]),
            "b"=>Dict("terminal_names"=>["p"]),
            "z"=>Dict("terminal_names"=>["p"]),
        ),
        "voltage_source"=>Dict(
            "s"=>Dict("bus"=>"s", "terminal_map"=>["p"], "v_magnitude"=>[230.0], "v_angle"=>[0.0]),
        ),
        "line"=>Dict(
            "l"=>Dict(
                "bus_from"=>"s",
                "bus_to"=>"b",
                "terminal_map_from"=>["p"],
                "terminal_map_to"=>["p"],
                "R_series_1_1"=>0.5,
            ),
        ),
        "switch"=>Dict(
            "alias"=>Dict(
                "bus_from"=>"b",
                "bus_to"=>"z",
                "terminal_map_from"=>["p"],
                "terminal_map_to"=>["p"],
                "status"=>"closed",
            ),
        ),
    )
    voltage=225+0.5im
    total=voltage*conj(2(voltage-230))
    first=-1000-100im
    meters=Measurement[]
    for (bus, power, sigma) in (("b", first, 10.0), ("z", total-first, 15.0))
        for (kind, value) in ((:pinj, real(power)), (:qinj, imag(power)))
            push!(
                meters,
                Measurement(; kind, bus, terminal = "p", reference = nothing, value, sigma),
            )
        end
    end
    push!(
        meters,
        Measurement(
            kind = :pinj,
            bus = "b",
            terminal = "p",
            reference = nothing,
            value = real(first),
            sigma = 9.0,
        ),
        Measurement(
            kind = :vmag,
            bus = "b",
            terminal = "p",
            reference = nothing,
            value = abs(voltage),
            sigma = 0.1,
        ),
        Measurement(
            kind = :vi,
            bus = "b",
            terminal = "p",
            reference = nothing,
            value = imag(voltage),
            sigma = 0.05,
        ),
    )
    # Profiled and explicit-current states represent the same likelihood.
    s=compile_state_estimator(net, meters; formulation = :voltage)
    p=SEParameters(s, meters)
    augmented=compile_state_estimator(net, meters)
    pa=SEParameters(augmented, meters)
    truth=[real(voltage), imag(voltage)]
    expected=selected_state_covariance(s, p, truth, [1, 2])
    means=copy(p.measurement_values)
    rng=MersenneTwister(8941)
    estimates=zeros(2, 500)
    statistics=zeros(500)
    for draw = 1:500
        p.measurement_values.=means .+ p.covariance_values .* randn(rng, length(means))
        r=solve_compiled_state_estimator(s, p, truth; optimality_tolerance = 1e-6)
        @test solve_status(r).publishable
        d=solve_diagnostics(r)
        @test d.dof==3
        @test d.nuisance_dimension==2
        estimates[:, draw]=r.state
        statistics[draw]=d.objective
        if draw<=10
            pa.measurement_values.=p.measurement_values
            for f in REVIEW_SE_SOLVERS
                explicit=f(
                    augmented,
                    pa,
                    initial_state_estimator(augmented, pa);
                    optimality_tolerance = 1e-6,
                )
                @test solve_status(explicit).publishable
                @test explicit.state[1:2]≈r.state atol=1e-5
                @test solve_diagnostics(explicit).objective≈d.objective atol=1e-6
                @test solve_diagnostics(explicit).dof==d.dof
            end
        end
    end
    center=sum(estimates; dims = 2)/500
    empirical=(estimates .- center)*(estimates .- center)'/499
    @info "Nonlinear merged-meter calibration" draws=500 dof=3 mean_J=sum(statistics)/500 relative_covariance_error=norm(
        empirical-expected,
    )/norm(expected)
    @test abs(sum(statistics)/500-3)<0.5 # >4 standard errors of chi-square mean
    @test norm(empirical-expected)/norm(expected)<0.2
    @test all(abs.(vec(center)-truth) .< 4sqrt.(diag(expected)/500))
end
