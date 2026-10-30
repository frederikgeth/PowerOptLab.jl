@testset "PWM closure: incomplete reserve cannot publish a fundamental optimum" begin
    inv=AdvancedInverter(id="closure",bus="poc",phase_terminals=["a","b","c"],neutral="n",
        topology=:THREE_LEG,s_max=20000.,i_max=50.,v_dc=750.,c_dc=.001,
        r_filter=.05,x_filter=.15,c_loss=.01,i_cap_max=.5,
        pwm_strategy=:CENTERED,f_sw=10000.)
    options=("tol"=>1e-8,"bound_relax_factor"=>0.)
    incomplete=solve_advanced_inverter(inv_grid3_bal(),inv;per_unit=true,s_base=1e5,
        pwm_max_iterations=1,solver_options=options)
    @test incomplete.termination_status in ("OPTIMAL","LOCALLY_SOLVED")
    @test incomplete.solve.feasible && incomplete.solve.optimal
    @test !incomplete.solve.publishable
    @test incomplete.pwm_iterations==1
    @test incomplete.i_cap_switching_reserved==0.
    @test incomplete.i_cap_switching>10inv.i_cap_max
    failed=Set(c.name for c in physical_feasibility(incomplete).checks if !c.passed)
    @test "pwm_reserve_closure" in failed
    @test "pwm_actual_capacitor_thermal_current" in failed

    # Identical hardware and source, with enough iterations for the declared
    # model to close. The carrier budget curtails dispatch instead of silently
    # accepting an undersized allowance.
    complete=solve_advanced_inverter(inv_grid3_bal(),inv;per_unit=true,s_base=1e5,
        solver_options=options)
    @test complete.solve.publishable
    @test physical_feasibility(complete).passed
    @test 1<complete.pwm_iterations<=16
    @test complete.p_poc < incomplete.p_poc/5
    actual=only(c for c in physical_feasibility(complete).checks if c.name=="pwm_actual_capacitor_thermal_current")
    @test actual.value <= inv.i_cap_max+actual.tolerance
    @test all(c.passed for c in physical_feasibility(complete).checks)
end

@testset "PWM closure: a nonfinite harmonic network cannot publish" begin
    fs=10000.;cap=.001
    # A lossless source branch exactly antiresonates with the DC link at the
    # first retained carrier harmonic. The independent audit cannot certify it.
    inv=AdvancedInverter(id="resonance",bus="poc",phase_terminals=["a","b","c"],neutral="n",
        topology=:THREE_LEG,s_max=6000.,i_max=15.,v_dc=750.,c_dc=cap,
        r_filter=.05,x_filter=.15,c_loss=.01,i_cap_max=10.,
        pwm_strategy=:CENTERED,f_sw=fs,pwm_dc_source_r=0.,
        pwm_dc_source_l=1/((2pi*fs)^2*cap))
    r=solve_advanced_inverter(inv_grid3_bal(),inv;per_unit=true,s_base=1e5)
    @test r.termination_status in ("OPTIMAL","LOCALLY_SOLVED")
    @test !r.solve.publishable
    @test !physical_feasibility(r).passed
    @test !isfinite(r.i_cap_switching)
    @test "pwm_capacitor_prediction" in Set(c.name for c in physical_feasibility(r).checks if !c.passed)
end

@testset "PWM closure tolerance cannot relax a small thermal rating" begin
    inv=AdvancedInverter(id="small_thermal",bus="poc",phase_terminals=["a","b","c"],neutral="n",
        topology=:THREE_LEG,s_max=20000.,i_max=50.,v_dc=750.,c_dc=.001,
        r_filter=.05,x_filter=.15,c_loss=.01,i_cap_max=.001,
        pwm_strategy=:CENTERED,f_sw=10000.)
    r=solve_advanced_inverter(inv_grid3_bal(),inv;per_unit=true,s_base=1e5,
        solver_options=("tol"=>1e-8,"bound_relax_factor"=>0.))
    @test r.solve.publishable
    check=only(c for c in physical_feasibility(r).checks if c.name=="pwm_actual_capacitor_thermal_current")
    @test check.value <= inv.i_cap_max+check.tolerance
    @test physical_feasibility(r).passed
end
