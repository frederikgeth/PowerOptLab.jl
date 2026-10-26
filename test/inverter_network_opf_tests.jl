include(joinpath(@__DIR__,"..","scripts","formulations","inverter_network_opf.jl"))

@testset "Public solve status for composed JuMP models" begin
    unstarted = solve_status(Model())
    @test !unstarted.has_primal
    @test !unstarted.publishable
    @test unstarted.termination_status == "OPTIMIZE_NOT_CALLED"
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model,x)
    @constraint(model,x >= 1)
    @objective(model,Min,x)
    optimize!(model)
    @test solve_status(model).publishable
    @test solve_status(model).feasible
    @test value(x) ≈ 1. atol=1e-6
    @constraint(model,x <= 0)
    optimize!(model)
    @test !solve_status(model).publishable
end

@testset "Full inverter OPF: matched controls and independent network replay" begin
    r = inverter_full_opf_validation()
    @test length(r.case.buses) == 7
    @test length(r.case.plants) == 3
    @test isempty(get(r.case.net,"ibr",Dict()))
    @test r.matched.variables >= 200
    @test r.agreement.p_error_W < .2
    @test r.agreement.q_error_var < .1
    @test r.agreement.voltage_error_V < .001
    @test r.agreement.converter_current_error_A < .001
    @test r.agreement.branch_current_error_A < .002
    for name in (:forward,:matched,:unmanaged,:dispatch,:benchmark)
        solution,audit = getproperty(r,name),getproperty(r.audits,name)
        @test solution.status.publishable
        @test audit.replay.converged
        @test audit.replay.iteration <= 30
        @test audit.voltage_error_V < 1e-5
        @test audit.current_error_A < 1e-5
        @test audit.kcl_error_A < 1e-5
        @test audit.kvl_error_V < 1e-6
        @test audit.power_balance_error_VA < .01
        @test real(audit.line_power_VA) > 0
        @test audit.neutral_max_V > 1.
        @test audit.phase_neutral_min_V > 220.
        for (i,device) in enumerate(solution.devices)
            plant = r.case.plants[i]
            current = device.converter_terminal.phase_current
            @test abs(sum(current)) < 1e-6
            @test device.plant.i_negative < 1e-6
            @test maximum(abs,current) <= plant.i_max+1e-6
            @test hypot(device.plant.p_conv,device.plant.q_conv) <= plant.s_max+.01
            @test device.plant.p_dc <= solution.requests[i].p_available+.001
            @test device.plant.p_dc >= -.001
            @test device.plant.p_filter_loss ≈ plant.r_filter*sum(abs2,current) atol=.001
            @test device.plant.p_loss ≈ plant.p_loss_fixed+plant.c_loss*sum(abs2,current) atol=.001
            @test device.plant.p_dc ≈ device.plant.p_conv+device.plant.p_loss+device.plant.p_cap_loss atol=.001
            terminal_power = sum(device.converter_terminal.phase_voltage .* conj.(current))
            @test terminal_power ≈ complex(device.plant.p_conv,device.plant.q_conv) atol=.001
        end
    end
    for (forward,matched) in zip(r.forward.devices,r.matched.devices)
        @test forward.control.current_scale > .9998
        @test forward.control.power_scale > .9998
        @test matched.operation.converter_apparent_power_margin_VA > 1000.
        @test abs(matched.operation.dc_availability_margin_W) < .001
        @test abs(matched.operation.q_constraint_residual_var) < .001
        @test abs(matched.plant.q_poc) > 100.
    end

    # Same high availability as dispatch mode, but autonomous local controls
    # do not impose the new feeder-wide phase-to-neutral voltage requirement.
    @test [x.p_available for x in r.unmanaged.requests] == [x.p_available for x in r.dispatch.requests]
    @test r.audits.unmanaged.phase_neutral_max_V > r.dispatch.voltage_max+3.
    @test r.audits.dispatch.phase_neutral_max_V ≈ r.dispatch.voltage_max atol=1e-5
    @test r.audits.dispatch.feeder_max_A ≈ r.dispatch.feeder_current_max atol=1e-6
    @test r.audits.benchmark.phase_neutral_max_V < r.dispatch.voltage_max-.1
    @test r.audits.benchmark.feeder_max_A < r.dispatch.feeder_current_max-1.
    @test r.dispatch.objective > r.benchmark.objective+.1
    @test r.dispatch.objective ≈ sum(r.dispatch.weights[i]*r.dispatch.devices[i].plant.p_poc/10000-
        .05*(r.dispatch.devices[i].plant.q_poc/r.case.plants[i].s_max)^2 for i in 1:3) atol=1e-9
    @test minimum(d.operation.converter_apparent_power_margin_VA for d in r.dispatch.devices) ≈ 0. atol=.01
    @test maximum(abs(abs(v[k])-abs(v[k]-v[4])) for v in r.dispatch.voltages[2:end],k in 1:3) > .5

    # Change model coordinates, retaining the same physical network and objective.
    for mode in (:matched,:dispatch)
        other = solve_inverter_network(r.case;mode,s_base=1e4)
        reference = getproperty(r,mode)
        @test other.status.publishable
        @test other.objective ≈ reference.objective atol=1e-7
        for i in 1:3
            @test other.devices[i].plant.p_poc ≈ reference.devices[i].plant.p_poc atol=.01
            @test other.devices[i].plant.q_poc ≈ reference.devices[i].plant.q_poc atol=.01
        end
        @test maximum(maximum(abs.(other.voltages[b]-reference.voltages[b])) for b in 1:7) < 1e-5
        @test audit_inverter_network(r.case,other).voltage_error_V < 1e-5
    end

    # Independently written phase-neutral inequalities agree with the native
    # vpn_min/vpn_max implementation, rather than confusing it with v_min/v_max.
    hook = solve_inverter_network(r.case;mode=:dispatch,voltage_bounds=:hook)
    @test hook.status.publishable
    @test hook.objective ≈ r.dispatch.objective atol=1e-7
    for i in 1:3
        @test hook.devices[i].plant.p_poc ≈ r.dispatch.devices[i].plant.p_poc atol=.01
        @test hook.devices[i].plant.q_poc ≈ r.dispatch.devices[i].plant.q_poc atol=.01
    end
    @test maximum(maximum(abs.(hook.voltages[b]-r.dispatch.voltages[b])) for b in 1:7) < 1e-5
    @test audit_inverter_network(r.case,hook).phase_neutral_max_V ≈ 244. atol=1e-5

    failed = solve_inverter_network(r.case;mode=:dispatch,max_iter=0)
    @test !failed.status.publishable
    @test isnan(failed.objective)
    @test all(v -> all(isnan,real.(v)),failed.voltages)
    @test all(v -> all(isnan,real.(v)),failed.branch_current)
    @test all(d -> isnan(d.plant.p_poc),failed.devices)
    @test audit_inverter_network(r.case,failed) === nothing
    @test !replay_inverter_network(r.case,r.audits.dispatch.dispatch;max_iter=1).converged
    @test_throws ArgumentError replay_inverter_network(r.case,ComplexF64[])
    @test_throws ArgumentError replay_inverter_network(r.case,fill(NaN+0im,3))
    @test_throws ArgumentError solve_inverter_network(r.case;mode=:bogus)
    @test_throws ArgumentError solve_inverter_network(r.case;s_base=0.)
    @test_throws ArgumentError solve_inverter_network(r.case;feeder_current_max=0.)
    @test_throws ArgumentError solve_inverter_network(r.case;voltage_bounds=:wrong)
    @test_throws ArgumentError solve_inverter_network(r.case;fixed_dispatch=fill(0im,3))
    @test_throws ArgumentError solve_inverter_network(r.case;mode=:dispatch,fixed_dispatch=[0im])
end
