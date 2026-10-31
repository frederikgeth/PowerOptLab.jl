_envelope_inv(;kwargs...)=AdvancedInverter(;id="envelope",bus="poc",s_max=6000.,
    phase_terminals=["a","b","c"],neutral=nothing,topology=:THREE_LEG,
    v_dc=750.,c_dc=300e-6,pwm_strategy=:CENTERED,f_sw=10000.,
    r_filter=.13225,x_filter=1.3225,kwargs...)

@testset "Carrier events: analytic square-wave orbit and dimensional scaling" begin
    inv=_envelope_inv();g=PowerOptLab._switching_envelope_geometry(inv)
    d=[.5,.5,0.];i=[2.,-1.,-1.]
    r=PowerOptLab._carrier_event_ripple(d,750.,i,g,1e-4,300e-6)
    slope=750g.gate_map*[.5,.5,0.]
    # Each phase is a centered triangular wave with amplitude |slope|*T/4.
    peak=abs.(slope) .* 1e-4/4
    @test r.phase_peak ≈ peak atol=1e-12
    @test r.phase_rms ≈ peak ./ sqrt(3) atol=1e-12
    @test r.phase_pp ≈ 2peak atol=1e-12
    @test r.neutral_peak==r.neutral_rms==0.
    @test maximum(abs,r.closure)<1e-11
    @test r.dc_max ≈ .5*1e-4/(4*300e-6) atol=1e-12
    @test r.dc_min ≈ -r.dc_max atol=1e-12
    scaled=PowerOptLab._carrier_event_ripple(d,1500.,2i,g,2e-4,600e-6)
    @test scaled.phase_rms ≈ 4r.phase_rms
    @test scaled.dc_max ≈ 2r.dc_max
    @test_throws ArgumentError PowerOptLab._carrier_event_ripple([1.1,.2,.3],750.,i,g,1e-4,300e-6)
    @test_throws ArgumentError PowerOptLab._carrier_event_ripple(d,750.,ones(3),g,1e-4,300e-6)
end

@testset "Carrier events: independent sampled integration and all-duty envelopes" begin
    for four in (false,true),angle in (0.,.41,1.19),modulation in (.15,.55,.9),unbalance in (0.,.03),power_angle in (-pi/2,0.,pi/2)
        inv=_envelope_inv(topology=four ? :FOUR_LEG : :THREE_LEG,neutral=four ? "n" : nothing,
            r_filter_neutral=four ? .13225 : 0.,x_filter_neutral=four ? 1.3225 : 0.)
        geometry=PowerOptLab._switching_envelope_geometry(inv)
        bounds=PowerOptLab._switching_current_bound_coefficients(inv)
        theta=angle .+ [0.,-2pi/3,2pi/3]
        u=modulation .* cos.(theta) .+ unbalance .* cos.(-theta .+ .3)
        i=3 .* cos.(theta .- power_angle) .+ .4 .* cos.(-theta .+ 1.)
        four && (i.+=.3cos(angle+.2))
        u=four ? [u;0.] : u
        currents=four ? [i;-sum(i)] : i
        duties=.5 .+ .5 .* (u .- (maximum(u)+minimum(u))/2)
        exact=PowerOptLab._carrier_event_ripple(duties,750.,currents,geometry,1e-4,300e-6)
        @test all(exact.phase_rms .<=750bounds.phase_rms .+ 1e-10)
        @test all(exact.phase_peak .<=750bounds.phase_peak .+ 1e-10)
        @test exact.neutral_rms<=750bounds.neutral_rms+1e-10
        @test exact.neutral_peak<=750bounds.neutral_peak+1e-10
        @test max(abs(exact.dc_min),abs(exact.dc_max))<=sum(abs,currents)/(8inv.f_sw*inv.c_dc)+1e-10
        @test maximum(abs,exact.closure)<1e-10
        # Independent midpoint gate sampling and trapezoidal integration.
        n=16384;h=1/(inv.f_sw*n);wave=zeros(5,n);state=zeros(5)
        rows=[geometry.gate_map;geometry.neutral_gate_map']
        for k in 1:n
            tau=(k-.5)/n
            carrier=tau<.5 ? -1+4tau : 3-4tau
            gates=Float64.(2duties .- 1 .>=carrier)
            e=gates .- duties
            slope=[750 .* (rows*e);-sum(e .* currents)/inv.c_dc]
            wave[:,k]=state .+ h/2 .* slope
            state.+=h .* slope
        end
        centered=wave .- sum(wave;dims=2) ./ n
        @test sqrt.(vec(sum(abs2,centered;dims=2)) ./ n)[1:3] ≈ exact.phase_rms atol=.002
        @test maximum(centered[5,:]) ≈ exact.dc_max atol=5e-5
    end
end

@testset "All-duty bounds: independent periodic damped RL orbits" begin
    for four in (false,true),resistance in (.13225,132.25)
        inv=_envelope_inv(topology=four ? :FOUR_LEG : :THREE_LEG,neutral=four ? "n" : nothing,
            r_filter=resistance,r_filter_neutral=four ? resistance : 0.,
            x_filter_neutral=four ? 1.3225 : 0.)
        g=PowerOptLab._switching_envelope_geometry(inv)
        b=PowerOptLab._switching_current_bound_coefficients(inv)
        duties=four ? [.08,.91,.39,.62] : [.08,.91,.39]
        # This oracle's carrier starts at -1 (the event implementation starts
        # at +1), so its equivalent gate pulses straddle the period boundary.
        edges=sort!(unique!([0.;duties ./ 2;1 .- duties ./ 2;1.]))
        rows=[g.gate_map;g.neutral_gate_map'];rate=g.damping_rate;period=1e-4
        orbit=zeros(4);segments=[]
        for k in 1:length(edges)-1
            a,z=edges[k:k+1];h=(z-a)*period
            # Deliberately evaluate the triangular carrier rather than events.
            tau=(a+z)/2;carrier=tau<.5 ? -1+4tau : 3-4tau
            gate=Float64.(2duties .- 1 .>=carrier)
            force=750rows*(gate .- duties)
            push!(segments,(h,force))
            orbit=orbit .* exp(-rate*h) .+ force .* (-expm1(-rate*h)/rate)
        end
        # Periodic boundary condition for x' + rate*x = force(t).
        initial=orbit ./ (-expm1(-rate*period));state=copy(initial)
        squares=zeros(4);peak=zeros(4);integral=zeros(4)
        for (h,force) in segments
            for j in 1:128
                t=(j-.5)*h/128
                value=state .* exp(-rate*t) .+ force .* (-expm1(-rate*t)/rate)
                squares.+=value.^2 .* (h/128);integral.+=value .* (h/128)
                peak=max.(peak,abs.(value))
            end
            state=state .* exp(-rate*h) .+ force .* (-expm1(-rate*h)/rate)
        end
        @test state ≈ initial atol=1e-10
        @test maximum(abs,integral ./ period)<2e-5
        @test all(sqrt.(squares ./ period) .<=[750b.phase_rms;750b.neutral_rms] .+ 1e-9)
        @test all(peak .<=[750b.phase_peak;750b.neutral_peak] .+ 1e-9)
    end
end

@testset "Envelope geometry rejects unsupported dynamics" begin
    @test_throws ArgumentError PowerOptLab._switching_envelope_geometry(_envelope_inv(x_filter=0.))
    @test_throws ArgumentError PowerOptLab._switching_envelope_geometry(_envelope_inv(pwm_strategy=:SPWM))
    @test_throws ArgumentError PowerOptLab._switching_envelope_geometry(_envelope_inv(c_filter_mid=10e-6,x_filter_grid=.1))
    bad=_envelope_inv(topology=:FOUR_LEG,neutral="n",r_filter_neutral=.1,x_filter_neutral=.2)
    @test_throws ArgumentError PowerOptLab._switching_envelope_geometry(bad)
end
