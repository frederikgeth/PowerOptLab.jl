_current_limit_bal(x)=x.*ComplexF64[1,cis(-2pi/3),cis(2pi/3)]

@testset "Declared phase/grid peaks: independent sinusoidal current oracle" begin
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    for grid in (false,true),base in (1e4,1e6)
        inv=_envelope_inv(pwm_strategy=:NONE,f_sw=nothing,i_max=10.,
            i_peak_max=grid ? nothing : 2.,i_grid_peak_max=grid ? 2. : nothing,
            p_loss_fixed=15.,a_loss=2.,c_loss=.3)
        r=allocate_inverter_services(inv,request;grid_voltage=_current_limit_bal(230.),
            reactive_request=0.,s_base=base)
        @test solve_status(r).publishable
        @test maximum(r.allocation.plant.i_mag) ≈ sqrt(2) atol=1e-5
        @test maximum(r.allocation.plant.i_total_peak_bound) ≈ 2. atol=1e-5
        currents=collect(r.allocation.converter_terminal.phase_current)
        peak=maximum(abs(sqrt(2)*real(z*cis(t))) for z in currents for t in range(0,2pi;length=8193))
        @test peak<=2.00001
        checks=physical_feasibility(r).checks
        name=grid ? "grid_peak_current_1" : "converter_peak_current_1"
        @test only(c for c in checks if c.name==name).passed
    end
end

@testset "Declared neutral peak limits per-phase reactive service" begin
    inv=_envelope_inv(pwm_strategy=:NONE,f_sw=nothing,topology=:FOUR_LEG,neutral="n",
        i_max=10.,In_max=10.,In_peak_max=1.,r_filter_neutral=.13225,x_filter_neutral=1.3225,
        p_loss_fixed=15.,a_loss=2.,c_loss=.3)
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    r=allocate_phase_inverter_services(inv,request;grid_voltage=_current_limit_bal(230.),
        reactive_request=[-600.,0.,0.])
    @test solve_status(r).publishable
    @test 0<r.reactive_scale<1
    @test sqrt(2)*r.allocation.plant.i_neutral ≈ 1. atol=1e-5
    @test only(c for c in physical_feasibility(r).checks if c.name=="neutral_peak_current").passed
    @test all(solve_status(s).publishable for (_,s) in r.stages)
end

@testset "Carrier RMS and peak envelopes constrain the shared plant" begin
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    for four in (false,true),peak in (false,true)
        inv=_envelope_inv(topology=four ? :FOUR_LEG : :THREE_LEG,neutral=four ? "n" : nothing,
            r_filter_neutral=four ? .13225 : 0.,x_filter_neutral=four ? 1.3225 : 0.,
            i_max=peak ? 10. : 3.7,In_max=four ? 10. : nothing,i_peak_max=peak ? 8. : nothing,
            pwm_current_envelope=true,dc_energy_model=:nonlinear,v_dc_min=700.,v_dc_max=800.,
            i_cap_max=20.,
            p_loss_fixed=15.,a_loss=2.,c_loss=.3)
        r=four ? allocate_phase_inverter_services(inv,request;grid_voltage=_current_limit_bal(230.),reactive_request=zeros(3)) :
            allocate_inverter_services(inv,request;grid_voltage=_current_limit_bal(230.),reactive_request=0.)
        @test solve_status(r).publishable
        @test r.allocation.plant.pwm_current_envelope
        @test minimum(r.allocation.plant.i_converter_ripple_bound)>3.
        @test maximum(r.allocation.plant.i_total_rms_bound)<=inv.i_max+1e-5
        if peak
            @test maximum(r.allocation.plant.i_total_peak_bound) ≈ 8. atol=1e-5
        else
            @test maximum(r.allocation.plant.i_total_rms_bound) ≈ 3.7 atol=1e-5
        end
        @test all(c.passed for c in physical_feasibility(r).checks)
    end
end

@testset "Peak limits cannot silently omit switching or unknown waveforms" begin
    @test_throws r"PWM peak-current limits" validate_device(_envelope_inv(i_peak_max=10.,i_cap_max=20.))
    @test_throws r"current envelope requires a declared v_dc_max" validate_device(_envelope_inv(pwm_current_envelope=true,i_cap_max=20.))
    @test_throws ArgumentError validate_device(_envelope_inv(pwm_strategy=:NONE,f_sw=nothing,i_peak_max=-1.))
    @test_throws ArgumentError validate_device(_envelope_inv(pwm_strategy=:NONE,f_sw=nothing,
        i_max=10.,i_peak_max=15.,pwm_ac_converter_reserve=(1.,1.,1.)))
end

@testset "Unmodeled AC switching is not reported as a certified bound" begin
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    inv=_envelope_inv(i_max=10.,i_cap_max=20.)
    r=allocate_inverter_services(inv,request;grid_voltage=_current_limit_bal(230.),reactive_request=0.)
    @test solve_status(r).publishable
    @test r.allocation.plant.i_total_rms_bound===nothing
    @test r.allocation.plant.i_total_peak_bound===nothing
end
