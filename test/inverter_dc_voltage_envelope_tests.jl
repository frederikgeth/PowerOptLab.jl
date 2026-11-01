_dc_guard_bal(x)=x .* ComplexF64[1,cis(-2pi/3),cis(2pi/3)]
_dc_guard_inv(;kwargs...)=_envelope_inv(;i_max=10.,i_cap_max=1.,
    dc_energy_model=:nonlinear,v_dc_min=735.,v_dc_max=765.,
    pwm_current_envelope=true,pwm_dc_voltage_envelope=true,
    p_loss_fixed=15.,a_loss=2.,c_loss=.3,kwargs...)

@testset "DC carrier voltage: independent coupled carrier charge integration" begin
    for four in (false,true),angle in (0.,.7,1.3),modulation in (.4,.9),amplitude in (0.,3.)
        inv=_dc_guard_inv(topology=four ? :FOUR_LEG : :THREE_LEG,neutral=four ? "n" : nothing,
            In_max=four ? 10. : nothing,r_filter_neutral=four ? .13225 : 0.,
            x_filter_neutral=four ? 1.3225 : 0.)
        u=_dc_guard_bal(230modulation) .+ 6cis(.4) .* [1,cis(2pi/3),cis(-2pi/3)]
        i=_dc_guard_bal(amplitude*cis(.8)) .+ .1amplitude*cis(.2) .* [1,cis(2pi/3),cis(-2pi/3)]
        four && (i.+=.05amplitude*cis(.3))
        voltage=sqrt(2) .* real.(u .* cis(angle));fund=sqrt(2) .* real.(i .* cis(angle))
        reference=four ? [voltage;0.] : voltage
        duties=.5 .+ (reference .- (maximum(reference)+minimum(reference))/2) ./ 750.
        n=8192;h=1/(inv.f_sw*n);current=zeros(3,n);gates=zeros(length(duties),n)
        state=zeros(3);L=inv.x_filter/(2pi*inv.f);Ln=inv.x_filter_neutral/(2pi*inv.f)
        for k in 1:n
            tau=(k-.5)/n;carrier=tau<.5 ? -1+4tau : 3-4tau
            gate=Float64.(2duties .- 1 .>=carrier);gates[:,k]=gate
            error=750 .* (gate .- duties)
            force=four ? error[1:3] .- error[4] : error .- sum(error)/3
            slope=four ? (force .- Ln/(L+3Ln)*sum(force)) ./ L : force ./ L
            current[:,k]=state .+ h/2 .* slope;state.+=h .* slope
        end
        current.-=sum(current;dims=2) ./ n
        dc=zeros(n)
        for k in 1:n
            phase=fund .+ current[:,k]
            legs=four ? [phase;-sum(phase)] : phase
            dc[k]=sum(gates[:,k] .* legs)
        end
        dc.-=sum(dc)/n
        charge=cumsum(-h .* dc) ./ inv.c_dc;charge.-=sum(charge)/n
        parts=PowerOptLab._dc_switching_voltage_bound(inv,real.(i),imag.(i))
        @test maximum(abs,charge)<=parts.total+1e-5
        @test parts.ac_variation>0
        @test parts.total==parts.frozen+parts.ac_variation
        if amplitude==0.
            @test parts.frozen==0.
            @test maximum(abs,charge)>1e-4 # AC ripple still charges the link.
        end
    end
end

@testset "DC carrier voltage: full-cycle rail enclosure including reserve" begin
    for four in (false,true),angle in (0.,.4,1.3),factor in (1.,1.25,1.35)
        inv=_dc_guard_inv(topology=four ? :FOUR_LEG : :THREE_LEG,neutral=four ? "n" : nothing,
            In_max=four ? 10. : nothing,r_filter_neutral=four ? .13225 : 0.,
            x_filter_neutral=four ? 1.3225 : 0.)
        u=_dc_guard_bal(230factor*cis(angle))
        i=_dc_guard_bal(3-.5im) .+ .4cis(.7) .* [1,cis(2pi/3),cis(-2pi/3)]
        four && (i.+=.1cis(.2))
        e=PowerOptLab._dc_energy_state(inv,sum(u .* i))
        reserve=PowerOptLab._dc_switching_voltage_bound(inv,real.(i),imag.(i)).total
        r=continuous_rail_margin(inv,u,e.voltage_2omega;energy=e,switching_voltage=reserve)
        minimum_gap=Inf
        for theta in range(0,2pi;length=16385)
            bus=sqrt(e.center_V2+2inv.v_dc*real(e.energy_coordinate*cis(2theta)))-reserve
            v=sqrt(2) .* real.(u .* cis(theta))
            needed=maximum(abs(v[a]-v[b]) for (a,b) in ((1,2),(2,3),(3,1)))
            four && (needed=max(needed,maximum(abs,v)))
            minimum_gap=min(minimum_gap,inv.m_max*bus-needed)
        end
        @test r.converged
        @test r.lower_V<=minimum_gap+1e-8
        @test r.upper_V>=minimum_gap-5e-5
        @test 0<=r.upper_V-r.lower_V<=1.0001e-6
        @test (r.lower_V>0)==(factor<1.35)
        @test_throws ArgumentError continuous_rail_margin(inv,u,e.voltage_2omega;energy=e)
    end
end

@testset "DC carrier voltage: allocation preserves service within complete voltage bounds" begin
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    v=_dc_guard_bal(230.)
    reference=nothing
    for base in (1e4,1e5,1e6)
        inv=_dc_guard_inv(v_dc_min=747.,v_dc_max=753.)
        r=allocate_inverter_services(inv,request;grid_voltage=v,reactive_request=-300.,
            active_fraction=.1,negative_sequence_request=1.0+0im,s_base=base)
        @test solve_status(r).publishable
        @test all(solve_status(s).publishable for (_,s) in r.stages)
        @test r.allocation.plant.dc_minimum_bound_V>=747 .- 1e-5
        @test r.allocation.plant.dc_maximum_bound_V<=753 .+ 1e-5
        @test r.allocation.plant.dc_switching_voltage_bound_V>1.
        d=dc_link_energy_diagnostics(inv;converter_voltage=collect(r.allocation.converter_terminal.phase_voltage),
            converter_current=collect(r.allocation.converter_terminal.phase_current))
        @test d.combined_minimum_V ≈ r.allocation.plant.dc_minimum_bound_V atol=1e-9
        @test d.combined_maximum_V ≈ r.allocation.plant.dc_maximum_bound_V atol=1e-9
        @test min(d.combined_minimum_V-747,753-d.combined_maximum_V)<1e-4
        reference===nothing || @test r.allocation.plant.p_poc ≈ reference atol=.001
        reference=r.allocation.plant.p_poc
    end
    old=_dc_guard_inv(v_dc_min=747.,v_dc_max=753.,pwm_dc_voltage_envelope=false)
    legacy=allocate_inverter_services(old,request;grid_voltage=v,reactive_request=-300.,
        active_fraction=.1,negative_sequence_request=1.0+0im)
    @test solve_status(legacy).publishable
    e=dc_link_energy_diagnostics(old;converter_voltage=collect(legacy.allocation.converter_terminal.phase_voltage),
        converter_current=collect(legacy.allocation.converter_terminal.phase_current))
    carrier=PowerOptLab._dc_switching_voltage_bound(old,real.(legacy.allocation.converter_terminal.phase_current),
        imag.(legacy.allocation.converter_terminal.phase_current)).total
    @test e.minimum_V>=747 .- 1e-5
    @test e.minimum_V-carrier<746.5 # Retained LF-only engineering counterexample.
    four=_dc_guard_inv(topology=:FOUR_LEG,neutral="n",In_max=10.,r_filter_neutral=.13225,
        x_filter_neutral=1.3225,v_dc_min=747.,v_dc_max=753.)
    r=allocate_phase_inverter_services(four,request;grid_voltage=v,reactive_request=[-600.,0.,0.],active_fraction=fill(.1,3))
    @test solve_status(r).publishable
    @test 0<r.reactive_scale<1
    @test r.allocation.plant.dc_minimum_bound_V>=747 .- 1e-5
    @test r.allocation.plant.dc_maximum_bound_V<=753 .+ 1e-5
    @test all(c.passed for c in physical_feasibility(r).checks)
end

@testset "DC carrier voltage: reject missing peak/source evidence" begin
    @test_throws r"lower and upper DC" validate_device(_dc_guard_inv(v_dc_min=nothing))
    @test_throws r"AC current envelope" validate_device(_dc_guard_inv(pwm_current_envelope=false))
    @test_throws r"RMS-only" validate_device(_dc_guard_inv(i_sw=.1))
    @test_throws r"RMS-only" validate_device(_dc_guard_inv(pwm_ac_converter_reserve=(.1,.1,.1)))
    @test_throws r"open carrier-frequency source" validate_device(_dc_guard_inv(pwm_dc_source_r=.1))
end
