_pwm_alloc_bal(x)=x.*ComplexF64[1,cis(-2pi/3),cis(2pi/3)]
_pwm_alloc_inv(;kwargs...)=AdvancedInverter(;id="pv",bus="poc",topology=:THREE_LEG,
    phase_terminals=["a","b","c"],neutral=nothing,s_max=6000.,i_max=10.,
    v_dc=750.,c_dc=300e-6,i_cap_max=1.,dv2_max=7.5,power_limit_location=:poc,
    r_filter=.13225,x_filter=1.3225,p_loss_fixed=15.,a_loss=2.,c_loss=.3,
    pwm_strategy=:CENTERED,f_sw=10000.,kwargs...)

@testset "PWM allocation: analytic bound encloses independent carrier dwell integral" begin
    for four in (false,true),angle in (0.,.37,1.1),ratio in (0.,.01,.03)
        inv=_pwm_alloc_inv(topology=four ? :FOUR_LEG : :THREE_LEG,
            neutral=four ? "n" : nothing,In_max=four ? 10. : nothing)
        u=_pwm_alloc_bal(230cis(angle)).+ratio*230cis(.4).*[1,cis(2pi/3),cis(-2pi/3)]
        i=_pwm_alloc_bal(2-.7im).+.3cis(.2).*[1,cis(2pi/3),cis(-2pi/3)]
        four && (i.+=.1+.2im)
        s2=sum(u.*i);d=im*s2/(4pi*50*300e-6*750)
        bound=PowerOptLab._allocation_pwm_variance(inv,real.(u),imag.(u),real.(i),imag.(i),
            real(d),imag(d),real(sum(u.*conj.(i))),real(s2),imag(s2);headroom=0.)
        exact=PowerOptLab._allocation_pwm_carrier(inv,u,i,d;angles=8192)
        @test isfinite(bound) && bound>=0
        @test exact.rms^2<=bound+1e-5
        @test exact.modulation_margin>0
    end
    inv=_pwm_alloc_inv()
    @test PowerOptLab._allocation_pwm_variance(inv,zeros(3),zeros(3),zeros(3),zeros(3),0.,0.,0.,0.,0.)==0.
end

@testset "PWM allocation: fixed sizing, staged reactive priority and coordinate invariance" begin
    inv=_pwm_alloc_inv();v=_pwm_alloc_bal(230.)
    req=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    reference=nothing
    for base in (1e4,1e5,1e6)
        r=allocate_inverter_services(inv,req;grid_voltage=v,reactive_request=0.,s_base=base)
        @test solve_status(r).publishable
        @test 500<r.allocation.plant.p_poc<1500
        @test r.allocation.plant.i_cap_thermal<=1.00001
        @test all(solve_status(s).publishable for (_,s) in r.stages)
        checks=physical_feasibility(r).checks
        @test all(c.passed for c in checks)
        @test any(c.name=="allocation_pwm_actual_capacitor_current" for c in checks)
        if reference!==nothing
            @test r.allocation.plant.p_poc ≈ reference atol=.01
        end
        reference=r.allocation.plant.p_poc
    end
    reactive=allocate_inverter_services(inv,req;grid_voltage=v,reactive_request=-1800.)
    @test solve_status(reactive).publishable
    @test -1800<reactive.allocation.plant.q_poc < -1000
    @test abs(reactive.allocation.plant.q_poc-first(reactive.stages).second.plant.q_poc)<.00012
    four=_pwm_alloc_inv(topology=:FOUR_LEG,neutral="n",In_max=10.)
    phase=allocate_phase_inverter_services(four,req;grid_voltage=v,reactive_request=fill(-600.,3))
    @test solve_status(phase).publishable
    @test 0<phase.reactive_scale<1
    @test phase.allocation.plant.q_poc ≈ reactive.allocation.plant.q_poc atol=.05
    @test all(solve_status(s).publishable for (_,s) in phase.stages)
    @test abs(phase.allocation.plant.q_poc-first(phase.stages).second.plant.q_poc)<.00012
    @test inv.i_cap_max==four.i_cap_max==1.
    @test_throws ArgumentError ConstrainedInverter(_pwm_alloc_inv(pwm_strategy=:SPWM),InverterOperatingConstraints())
    @test_throws ArgumentError ConstrainedInverter(_pwm_alloc_inv(pwm_dc_source_r=.1),InverterOperatingConstraints())
    @test_throws ArgumentError ConstrainedInverter(inv,InverterOperatingConstraints();pwm_thermal_headroom=-.1)
end
