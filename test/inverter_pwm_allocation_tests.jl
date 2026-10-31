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
    @test PowerOptLab._allocation_pwm_variance(inv,zeros(3),zeros(3),zeros(3),zeros(3),0.,0.,0.,0.,0.)>0.
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
    # This legacy fixture has zero neutral inductance. Its all-duty coupled
    # thermal reserve alone exceeds 1 A; a failed bound-constrained solve must
    # remain unresolved rather than being reported as a physical impossibility.
    @test PowerOptLab._carrier_thermal_bound(four).variance_per_V2*four.v_dc^2>four.i_cap_max^2
    @test !solve_status(phase).publishable
    @test phase.current_command===nothing
    @test phase.plant_feasibility===:unresolved
    @test inv.i_cap_max==four.i_cap_max==1.
    @test_throws ArgumentError ConstrainedInverter(_pwm_alloc_inv(pwm_strategy=:SPWM),InverterOperatingConstraints())
    @test_throws ArgumentError ConstrainedInverter(_pwm_alloc_inv(pwm_dc_source_r=.1),InverterOperatingConstraints())
    @test_throws ArgumentError ConstrainedInverter(inv,InverterOperatingConstraints();pwm_thermal_headroom=-.1)
end

@testset "Coupled carrier thermal bound: independent gates, damping and fixed sizing" begin
    for four in (false,true),resistance in (.13225,132.25)
        inv=_pwm_alloc_inv(topology=four ? :FOUR_LEG : :THREE_LEG,
            neutral=four ? "n" : nothing,r_filter=resistance,
            r_filter_neutral=four ? resistance : 0.,x_filter_neutral=four ? 1.3225 : 0.)
        geometry=PowerOptLab._switching_envelope_geometry(inv)
        K=four ? [geometry.gate_map;geometry.neutral_gate_map'] : geometry.gate_map
        bound=PowerOptLab._carrier_thermal_bound(inv)
        for angle in (0.,.4,1.2),spread in (.1,.6,.95)
            refs=cos.(angle .+ [0.,-2pi/3,2pi/3]);four && (refs=[refs;0.])
            d=.5 .+ spread.*(refs.-(maximum(refs)+minimum(refs))/2)./(maximum(refs)-minimum(refs))
            exact=inv.v_dc^2/inv.f_sw^2*PowerOptLab._carrier_ac_bridge_variance(d,K)
            @test exact<=inv.v_dc^2*bound.variance_per_V2
            # Independent periodic passive R/L orbit, evaluated against a
            # triangular comparator. Integrate g'i, not frozen gate currents.
            rate=geometry.damping_rate;T=1/inv.f_sw
            edges=sort!(unique!([0.;d./2;(1 .- d./2);1.]))
            state=zeros(length(d));segments=[]
            for j in 1:length(edges)-1
                a,z=edges[j:j+1];dt=(z-a)*T;t=(a+z)/2
                gate=Float64.(2d.-1 .>= (t<.5 ? -1+4t : 3-4t))
                force=inv.v_dc.*(K*(gate.-d));push!(segments,(dt,gate,force))
                state=state.*exp(-rate*dt).+force.*(-expm1(-rate*dt)/rate)
            end
            state./=-expm1(-rate*T);initial=copy(state);first=0.;second=0.
            for (dt,gate,force) in segments
                for j in 1:256
                    t=(j-.5)*dt/256
                    ripple=state.*exp(-rate*t).+force.*(-expm1(-rate*t)/rate)
                    bridge=sum(gate.*ripple)
                    first+=bridge*dt/256;second+=bridge^2*dt/256
                end
                state=state.*exp(-rate*dt).+force.*(-expm1(-rate*dt)/rate)
            end
            rms=sqrt(max(0.,second/T-(first/T)^2))
            @test state ≈ initial atol=1e-10
            @test rms<=sqrt(exact)+inv.v_dc*bound.damping_rms_per_V+1e-7
        end
    end
    for four in (false,true)
        inv=_pwm_alloc_inv(topology=four ? :FOUR_LEG : :THREE_LEG,
            neutral=four ? "n" : nothing,r_filter_neutral=four ? .13225 : 0.,
            x_filter_neutral=four ? 1.3225 : 0.)
        idle=PowerOptLab._allocation_pwm_carrier(inv,_pwm_alloc_bal(230.),zeros(ComplexF64,3),0im)
        @test idle.frozen_rms==0.
        @test .2<idle.ac_lossless_rms<.5
        request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
        r=four ? allocate_phase_inverter_services(inv,request;
                grid_voltage=_pwm_alloc_bal(230.),reactive_request=fill(-600.,3)) :
            allocate_inverter_services(inv,request;
                grid_voltage=_pwm_alloc_bal(230.),reactive_request=-1800.)
        @test solve_status(r).publishable
        @test r.allocation.plant.i_cap_thermal<=1.00001
        @test inv.i_cap_max==1. && inv.v_dc==750. && inv.c_dc==300e-6
    end
    @test_throws r"LCL" ConstrainedInverter(_pwm_alloc_inv(c_filter_mid=10e-6,x_filter_grid=.1),InverterOperatingConstraints())
end
