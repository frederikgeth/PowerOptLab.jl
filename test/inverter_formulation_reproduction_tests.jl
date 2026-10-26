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
