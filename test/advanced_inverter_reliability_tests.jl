@testset "Advanced inverter reliability: full PWM publication contract" begin
    common = (id="review", bus="poc", phase_terminals=["a","b","c"],
        topology=:FOUR_LEG, s_max=20e3, i_max=40.0, In_max=40.0,
        v_dc=700.0, c_dc=1.1e-3, i_cap_max=10.0,
        r_filter=0.05, x_filter=0.15, pwm_strategy=:CENTERED,
        f_sw=10e3, pwm_fundamental_samples=24, pwm_carrier_samples=64)
    kw = (per_unit=true, s_base=100e3)
    incomplete = solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; common...); kw..., pwm_max_iterations=1)
    @test incomplete.inner_solve.publishable
    @test !solve_status(incomplete).publishable
    @test incomplete.termination_status == "PWM_ITERATION_LIMIT"
    @test incomplete.pwm_status == :ITERATION_LIMIT
    @test incomplete.i_cap_switching > 10
    @test incomplete.i_cap_switching_reserved == 0
    @test incomplete.pwm_reserve_margin < -10
    @test isnan(incomplete.p_poc)
    @test all(isnan, incomplete.i_mag)
    @test solve_diagnostics(incomplete).inner_solve.publishable

    invalid = solve_advanced_inverter(inv_grid3_bal(), AdvancedInverter(;
        merge(common, (v_dc=600.0, pwm_strategy=:SPWM))...); kw...)
    @test invalid.inner_solve.publishable
    @test !solve_status(invalid).publishable
    @test invalid.pwm_status == :AUDIT_FAILED
    @test invalid.pwm_modulation_margin < -1
    @test isnan(invalid.p_poc)

    singular = solve_advanced_inverter(inv_grid3_bal(), AdvancedInverter(;
        merge(common, (pwm_dc_source_r=0.0,
            pwm_dc_source_l=1/((2pi*3common.f_sw)^2*common.c_dc),
            pwm_dc_harmonics=1))...); kw...)
    @test singular.inner_solve.publishable
    @test !solve_status(singular).publishable
    @test singular.pwm_status == :AUDIT_FAILED
    @test !isfinite(singular.i_cap_switching)

    # A rotated balanced waveform fits the coarse hull and carrier grids while
    # exceeding the rail between them. Only the independent dense audit sees it.
    rotated = inv_grid3_src(mags=[230.,230.,230.],
        angs=[pi/12, pi/12-2pi/3, pi/12+2pi/3])
    missed_peak = solve_advanced_inverter(rotated, AdvancedInverter(;
        merge(common, (v_dc=560.0, r_filter=0.0, x_filter=0.0,
            i_cap_max=1000.0, n_samples=4, pwm_fundamental_samples=12))...); kw...)
    @test missed_peak.inner_solve.publishable
    @test isfinite(missed_peak.i_cap_switching)
    @test missed_peak.pwm_modulation_margin >= 0
    @test missed_peak.switching_margin < -1
    @test missed_peak.pwm_status == :AUDIT_FAILED
    @test !solve_status(missed_peak).publishable

    complete = solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; common...); kw...)
    @test solve_status(complete).publishable
    @test complete.pwm_status == :CONVERGED
    @test complete.pwm_reserve_margin >= -0.02
    @test isfinite(complete.p_poc)
    @test complete.switching_margin >= -7e-4
    # Audit checks cannot be bypassed merely because no corresponding AC rating
    # was declared. Corrupt individual diagnostics on an otherwise valid point.
    replace_result(r; kwargs...) = InverterResult((get(kwargs, n, getfield(r, n))
        for n in fieldnames(InverterResult))...)
    inv = AdvancedInverter(; common...)
    @test !PowerOptLab._pwm_audit_valid(
        replace_result(complete; switching_margin=-1.0), inv)
    @test !PowerOptLab._pwm_audit_valid(
        replace_result(complete; pwm_dc_network_margin=0.0), inv)
    acinv = AdvancedInverter(; merge(common, (pwm_ac_ripple=true,
        pwm_ac_harmonics=32))...)
    @test !PowerOptLab._pwm_audit_valid(
        replace_result(complete; i_grid_switching_rms=[NaN,0.,0.]), acinv)

    inner_failed = solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; common...); kw..., solver_options=("max_iter"=>0,))
    @test !solve_status(inner_failed).publishable
    @test !inner_failed.inner_solve.publishable
    @test inner_failed.pwm_status == :INNER_SOLVE_FAILED
end

@testset "Advanced inverter reliability: normalized physical ratings" begin
    for sb in (1e5, 1e6, 1e8)
        r = solve_advanced_inverter(inv_grid(),
            AdvancedInverter(id="rating", bus="poc", s_max=5000.0);
            per_unit=true, s_base=sb)
        @test solve_status(r).publishable
        @test hypot(r.p_conv,r.q_conv) <= 5000.1
        @test r.p_poc > 4990
        @test r.pwm_status == :NOT_REQUESTED
        limited = solve_advanced_inverter(inv_grid(),
            AdvancedInverter(id="current", bus="poc", s_max=5000.0, i_max=10.0);
            per_unit=true, s_base=sb)
        @test solve_status(limited).publishable
        @test maximum(limited.i_mag) <= 10.001
        @test maximum(limited.i_mag) > 9.99
    end
end

@testset "Advanced inverter reliability: bounded omitted DC spectrum" begin
    a = cis(2pi/3)
    U = 0.3/sqrt(2)*ComplexF64[1,a^-1,a]
    currents = 1/sqrt(2)*ComplexF64[1,a^-1,a]
    common = (id="spectral", bus="poc", phase_terminals=["a","b","c"],
        topology=:FOUR_LEG, s_max=10.0, In_max=10.0, v_dc=1.0, c_dc=1.0,
        i_cap_max=10.0, f=0.01, f_sw=1.0, pwm_strategy=:CENTERED,
        pwm_fundamental_samples=24, pwm_carrier_samples=128,
        pwm_dc_source_r=1e-5, pwm_dc_source_l=1/(6pi)^2)
    audit(nh; overrides...) = PowerOptLab._pwm_ripple_audit(
        AdvancedInverter(; merge(common, (;pwm_dc_harmonics=nh, overrides...))...),
        real.(U), imag.(U), real.(currents), imag.(currents))
    coarse, fine = audit(1), audit(64)
    @test fine.i_rms > 100
    @test coarse.i_rms >= fine.i_rms
    @test coarse.dc_network_margin < 1e-3
    @test audit(2).i_rms >= fine.i_rms
    @test !isfinite(audit(1; pwm_dc_source_r=0.0).i_rms)
    # Independently evaluate the rational gain over many omitted harmonics;
    # include a resonance between integer harmonics and finite source damping.
    for r in (0.0, 1e-5, 0.1), l in (0.001, 0.023), nh in (1, 8)
        bound = PowerOptLab._dc_tail_bound(1.0, 1.0, r, l, nh)
        actual = maximum(begin
            w = 2pi*h
            (w^2*r^2+w^4*l^2)/((1-w^2*l)^2+w^2*r^2)
        end for h in nh+1:1000)
        @test bound.gain_sq + 1e-10*max(1,actual) >= actual
    end
    @test PowerOptLab._dc_tail_bound(1.,1.,0.01,0.,1).gain_sq == 1
end

@testset "Advanced inverter reliability: unequal-bank 2ω rail charge" begin
    inv = AdvancedInverter(id="split", bus="poc", phase_terminals=["a","b","c"],
        topology=:SPLIT_DC, s_max=20e3, In_max=40.0, v_dc=800.0,
        c_dc=1e-3, c_dc_upper=1e-3, c_dc_lower=3e-3)
    # At θ=0: equal series charge makes Vu=615 V and Vl=205 V.
    # An instantaneous -208 V command exceeds the physical lower rail by 3 V;
    # the previous equal split of d(t) incorrectly gave +2 V headroom.
    ure = [-208/sqrt(2),0.,0.]; uim=zeros(3)
    margin = PowerOptLab._switching_margin(inv,ure,uim,20.,0.,0.,0.,-200.)
    @test margin ≈ -3.0 atol=1e-8
    cu,cl = inv.c_dc_upper, inv.c_dc_lower
    for th in range(0,2pi;length=121)
        d=20cos(2th)
        n=-200+PowerOptLab._split_midpoint_bus_factor(inv)*d
        @test (800+d)/2-n ≈ cl/(cu+cl)*(800+d)
        @test (800+d)/2+n ≈ cu/(cu+cl)*(800+d)
    end
    pwm = AdvancedInverter(; merge(NamedTuple{fieldnames(AdvancedInverter)}(
        Tuple(getfield(inv,n) for n in fieldnames(AdvancedInverter))),
        (pwm_strategy=:SPWM, f_sw=10e3, i_cap_max=20., pwm_ac_ripple=true,
         x_filter=0.15))...)
    dc = PowerOptLab._pwm_ripple_audit(pwm,ure,uim,ones(3),zeros(3),20.,0.,0.,0.,-200.)
    ac = PowerOptLab._pwm_ac_ripple_audit(pwm,ure,uim,20.,0.,0.,0.,-200.)
    @test dc.modulation_margin ≈ -3 atol=1e-8
    @test !isfinite(dc.i_rms)
    @test all(isnan,ac.converter_rms)
end

@testset "Advanced inverter reliability: three-leg conductor dimensions" begin
    common=(id="three",bus="poc",phase_terminals=["a","b","c"],
        topology=:THREE_LEG,s_max=20e3,i_max=40.,v_dc=700.,c_dc=1e-3,
        i_cap_max=30.,pwm_strategy=:CENTERED,f_sw=10e3,pwm_ac_ripple=true,
        pwm_fundamental_samples=24,pwm_carrier_samples=64,pwm_ac_harmonics=32)
    a=cis(2pi/3); U=150*ComplexF64[1,a^-1,a]
    audit(inv)=PowerOptLab._pwm_ac_ripple_audit(inv,real.(U),imag.(U))
    grounded=AdvancedInverter(;common...,x_filter=0.15)
    floating=AdvancedInverter(;common...,neutral=nothing,x_filter=0.15)
    @test validate_device(grounded) === nothing
    @test validate_device(floating) === nothing
    rg,rf=audit(grounded),audit(floating)
    @test rg.converter_rms ≈ rf.converter_rms rtol=1e-12
    @test rg.neutral_rms < 1e-10
    X=[0.15 0.01 0.0 0.02;0.01 0.20 0.02 0.01;
       0.0 0.02 0.25 0.03;0.02 0.01 0.03 0.10]
    coupled=AdvancedInverter(;common...,x_filter_matrix=X,
        x_filter_grid=0.1,c_filter_mid=20e-6,r_filter_damping=1.0)
    @test validate_device(coupled) === nothing
    rc=audit(coupled)
    @test all(isfinite,rc.grid_rms)
    @test maximum(rc.converter_rms)>0
    @test maximum(rc.shunt_rms)>0
    @test rc.neutral_rms<1e-9
end
