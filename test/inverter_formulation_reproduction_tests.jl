include(joinpath(@__DIR__,"..","scripts","formulations","ibr_reproductions.jl"))

@testset "IBR forward/optimisation reproductions" begin
    rows = ibr_formulation_reproductions()
    @test length(rows.feeder) == 3
    for row in rows.feeder
        @test row.forward.solve.publishable && row.constrained.solve.publishable
        @test row.forward.control.current_scale > .9998
        @test row.forward.control.power_scale > .9998
        @test row.p_error_W < .2
        @test row.q_error_var < .2
        @test row.current_error_A < .001
        @test row.voltage_error_V < .001
        @test row.constrained.operation.dc_lower_margin_W >= -.03
        @test abs(row.constrained.operation.q_constraint_residual_var) < .03
        @test minimum(row.constrained.operation.converter_current_margin_A) > 5.
    end
    for row in rows.proportional
        @test row.forward.solve.publishable && row.optimisation.status.publishable
        @test row.optimisation.scale ≈ row.analytic_scale atol=1e-7
        @test row.p_error_W < .4
        @test row.q_error_var < .2
        p,q = row.optimisation.result.plant.p_poc,row.optimisation.result.plant.q_poc
        @test p ≈ 6000row.analytic_scale atol=.03
        @test q ≈ 2640row.analytic_scale atol=.03
        @test q/p ≈ .44 atol=1e-8
        @test minimum(row.optimisation.result.operation.converter_current_margin_A) >= -1e-6
    end
    for row in rows.sequence_priority
        @test row.optimisation.status.publishable
        @test row.optimisation.scale ≈ row.analytic_scale atol=1e-7
        @test 0 < row.optimisation.scale < 1
        @test row.current_error_A < 1e-6
        current = row.optimisation.result.converter_terminal.current_sequence
        reserved_index = row.priority === :power ? 2 : 3
        reserved_command = row.priority === :power ? row.exact.i1_request : row.exact.i2_request
        @test current[reserved_index] ≈ reserved_command atol=1e-6
        @test maximum(row.optimisation.result.plant.i_mag) ≈ 6. atol=1e-6
        @test row.optimisation.result.plant.i_zero < 1e-7
    end
end

@testset "IBR P-Q tutorial uses independent capability oracles" begin
    rows = ibr_pq_tutorial()
    @test length(rows) == 8
    for row in rows
        @test row.result.solve.publishable
        @test row.result.plant.p_poc ≈ row.expected_p_W atol=.03
        @test row.result.plant.q_poc ≈ row.expected_q_var atol=.1
        @test hypot(row.result.plant.p_poc,row.result.plant.q_poc) ≈ row.capacity_VA atol=.03
        @test minimum(row.result.operation.converter_current_margin_A) >= -1e-6
        @test row.result.operation.dc_lower_margin_W >= -.03
    end
end

@testset "IBR unbalanced network reproductions and intentional differences" begin
    rows = ibr_network_reproductions()
    a = cis(2pi/3)
    seq(i) = ((i[1]+i[2]+i[3])/3,(i[1]+a*i[2]+a^2*i[3])/3,
        (i[1]+a^2*i[2]+a*i[3])/3)
    @test length(rows.pairs) == 4
    for row in rows.pairs
        @test row.forward.solve.publishable && row.constrained.solve.publishable
        if row.expected === :agreement
            @test row.p_error_W < .2
            @test row.q_error_var < .2
            @test row.current_error_A < .001
            @test row.voltage_error_V < .001
            @test row.forward.control.current_scale > .9998
            @test row.forward.control.power_scale > .9998
        else
            @test row.expected === :different_q_location
            @test row.q_error_var > 1000.
            @test row.constrained.plant.q_poc ≈ 600. atol=.03
        end
        # Independent complex KVL, capacitor admittance and power-location audits.
        for result in (row.forward,row.constrained)
            u = hasproperty(result,:control) ? result.control.phase_voltage : result.operation.phase_voltage
            ic,ig = result.converter_terminal.phase_current,result.grid_phase_current
            @test collect(u) ≈ collect(_repro_phases(230.0+0im,8cis(.2)))+.05collect(ig) atol=1e-5
            @test sum(u .* conj.(ig)) ≈ complex(result.plant.p_poc,result.plant.q_poc) atol=.03
            @test sum(result.converter_terminal.phase_voltage .* conj.(ic)) ≈
                complex(result.plant.p_conv,result.plant.q_conv) atol=.03
            @test result.plant.p_dc ≈ result.plant.p_conv+result.plant.p_loss+result.plant.p_cap_loss atol=.03
            target = row.target === :converter ? ic : ig
            @test abs(seq(target)[1]) < 1e-6
            @test abs(seq(target)[3]) < 1e-6
            if row.id !== :lcl_poc_q_mismatch
                @test imag(sum(u .* conj.(target))) ≈ 600. atol=.1
            end
            if row.plant.c_filter_mid > 0
                ym = inv(complex(row.plant.r_filter_damping,-1/(2pi*row.plant.f*row.plant.c_filter_mid)))
                vmid = u .+ complex(row.plant.r_filter_grid,row.plant.x_filter_grid).*ig
                @test collect(ic .- ig) ≈ collect(ym .* vmid) atol=1e-6
                @test abs(seq(row.target === :converter ? ig : ic)[3]) > .01
            end
        end
        if row.id === :lcl_poc_q_mismatch
            u = row.forward.control.phase_voltage
            @test imag(sum(u .* conj.(row.forward.converter_terminal.phase_current))) ≈ 600. atol=.1
        end
    end
    saturated = rows.saturated
    @test saturated.forward.solve.publishable && saturated.constrained.solve.publishable
    @test saturated.constrained.plant.q_poc ≈ 2640. atol=.03
    @test saturated.forward.plant.q_poc < saturated.constrained.plant.q_poc-15.
    @test saturated.constrained.plant.p_poc < saturated.forward.plant.p_poc-5.
    @test saturated.forward.control.current_scale < .999
    @test maximum(saturated.constrained.plant.i_mag) ≈ saturated.plant.i_max atol=1e-6
    u1 = seq(saturated.constrained.operation.phase_voltage)[2]
    @test hypot(saturated.constrained.plant.p_poc,saturated.constrained.plant.q_poc) ≈
        3abs(u1)*saturated.plant.i_max atol=.03
    @test saturated.constrained.operation.converter_apparent_power_margin_VA >= -.03
    @test minimum(saturated.constrained.operation.converter_current_margin_A) >= -1e-6
    @test maximum(saturated.forward.plant.i_mag) <= saturated.plant.i_max+1e-6

    service = rows.impossible_service
    @test service.minimum_required_S_VA > service.plant.s_max
    @test service.forward.solve.publishable
    @test hypot(service.forward.plant.p_conv,service.forward.plant.q_conv) <= service.plant.s_max+.03
    @test !service.constrained.solve.publishable
    @test isnan(service.constrained.plant.q_poc)
    @test isnan(service.constrained.operation.q_target_var)

    zero = rows.zero_sequence
    @test abs(zero.required_shunt_i0_A) > .02
    @test !zero.balanced.solve.publishable
    @test isnan(zero.balanced.plant.p_poc)
    @test zero.released.solve.publishable
    @test abs(seq(zero.released.converter_terminal.phase_current)[1]) < 1e-6
    ym = inv(complex(zero.plant.r_filter_damping,-1/(2pi*zero.plant.f*zero.plant.c_filter_mid)))
    zg = complex(zero.plant.r_filter_grid,zero.plant.x_filter_grid)
    @test seq(zero.released.grid_phase_current)[1] ≈ -ym*zero.u0/(1+ym*zg) atol=1e-6
end
