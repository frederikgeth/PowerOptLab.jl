# Carrier-period engineering envelopes for a monolithic, inductive bridge.
# These are separate from the calibrated DC thermal allocation allowance.
const _PWM_DC_VOLTAGE_NORM_EPS = 1e-3 # A; one-sided norm bias, not a voltage tolerance

function _switching_envelope_geometry(inv)
    inv.topology in (:THREE_LEG,:FOUR_LEG) || throw(ArgumentError(
        "switching engineering envelope requires a monolithic three/four-leg bridge"))
    inv.pwm_strategy===:CENTERED || throw(ArgumentError(
        "switching engineering envelope requires centered shared-carrier PWM"))
    _lcl_active(inv) && throw(ArgumentError(
        "inductive switching envelope does not cover an LCL resonance"))
    inv.f_sw!==nothing && isfinite(inv.f_sw) && inv.f_sw>0 ||
        throw(ArgumentError("positive finite switching frequency required"))
    R,X=_series_filter_matrices(inv,3,1.;side=:converter)
    if inv.topology===:THREE_LEG
        Q=[1/sqrt(2) 1/sqrt(6);-1/sqrt(2) 1/sqrt(6);0. -2/sqrt(6)]
        # A declared terminal neutral is absent from the three-wire current
        # path, even when the primitive includes that conductor.
        Rr=Q'*R[1:3,1:3]*Q;Lr=Q'*X[1:3,1:3]*Q/(2pi*inv.f)
    else
        size(R)==(4,4) || throw(ArgumentError("four-leg envelope requires a physical neutral conductor"))
        B=[Matrix{Float64}(I,3,3);-ones(1,3)]
        Rr=B'*R*B;Lr=B'*X*B/(2pi*inv.f)
    end
    isapprox(Lr,Lr';atol=1e-12,rtol=1e-12) &&
        isapprox(Rr,Rr';atol=1e-12,rtol=1e-12) || throw(ArgumentError("symmetric filter matrices required"))
    Lr=(Lr+Lr')/2;Rr=(Rr+Rr')/2
    minimum(eigvals(Symmetric(Lr)))>1e-12 || throw(ArgumentError(
        "switching envelope requires positive inductance in every active mode"))
    # A scalar passive R/L makes every physical output satisfy the same
    # first-order equation. General damping matrices require a separate norm
    # enclosure; silently treating them as scalar would invalidate peak bounds.
    rate=sum(Rr.*Lr)/sum(abs2,Lr)
    rate>=-1e-10 && norm(Rr-rate*Lr)<=1e-10*max(1.,norm(Rr)) ||
        throw(ArgumentError("switching envelope requires proportional passive R and L matrices"))
    gate_map=inv.topology===:THREE_LEG ? Q*(Lr\Q') :
        Lr\[Matrix{Float64}(I,3,3) -ones(3)]
    neutral=inv.topology===:THREE_LEG ? zeros(size(gate_map,2)) : -vec(sum(gate_map;dims=1))
    (gate_map=gate_map,neutral_gate_map=neutral,damping_rate=max(0.,rate))
end

"""
    _switching_current_bound_coefficients(inv)

Carrier-period phase/neutral ripple bounds per volt of maximum DC voltage.
For a periodic zero-mean forcing of range W, its primitive has range <= W*T/4.
Popoviciu's variance bound then gives RMS <= W*T/8. Proportional passive R/L
cannot increase that RMS; its peak is bounded by the primitive's whole range.
The bounds cover all gate states and duties, with no Fourier truncation or
empirical thermal factor. The carrier-period model freezes slow bus/reference
quantities; it is not a bound on transients or on errors in that separation.
"""
function _switching_current_bound_coefficients(inv)
    geometry=_switching_envelope_geometry(inv)
    widths=vec(sum(abs,geometry.gate_map;dims=2))
    neutral=sum(abs,geometry.neutral_gate_map)
    (phase_rms=widths/(8inv.f_sw),phase_peak=widths/(4inv.f_sw),
        neutral_rms=neutral/(8inv.f_sw),neutral_peak=neutral/(4inv.f_sw))
end

function _inverter_current_envelope(inv)
    phase=zeros(3);neutral=0.;peak=zeros(3);npeak=0.
    if inv.pwm_current_envelope
        bounds=_switching_current_bound_coefficients(inv)
        phase=bounds.phase_rms.*inv.v_dc_max;neutral=bounds.neutral_rms*inv.v_dc_max
        peak=bounds.phase_peak.*inv.v_dc_max;npeak=bounds.neutral_peak*inv.v_dc_max
    end
    # Expert reserves may occupy the same frequency band: sum their RMS
    # bounds by Minkowski, rather than assuming quadrature independence.
    (converter_rms=collect(inv.pwm_ac_converter_reserve).+phase,
        grid_rms=collect(inv.pwm_ac_grid_reserve).+phase,
        neutral_rms=inv.pwm_ac_neutral_reserve+neutral,
        phase_peak=peak,neutral_peak=npeak,
        rms_resolved=inv.pwm_strategy===:NONE || inv.pwm_current_envelope,
        peak_resolved=(inv.pwm_strategy===:NONE || inv.pwm_current_envelope) &&
            all(iszero,[inv.pwm_ac_converter_reserve...;inv.pwm_ac_grid_reserve...;inv.pwm_ac_neutral_reserve]))
end

function _validate_switching_current_limits(inv)
    peaks=(inv.i_peak_max,inv.i_grid_peak_max,inv.In_peak_max)
    for limit in peaks
        limit===nothing || isfinite(limit) && limit>0 ||
            throw(ArgumentError("declared peak-current limits must be positive and finite"))
    end
    inv.topology===:SINGLE_PHASE && inv.In_peak_max!==nothing &&
        throw(ArgumentError("neutral peak-current limit requires a three-phase topology"))
    if any(!isnothing,peaks)
        if inv.pwm_strategy!==:NONE && !inv.pwm_current_envelope
            throw(ArgumentError("PWM peak-current limits require pwm_current_envelope=true"))
        end
        any(!iszero,[inv.pwm_ac_converter_reserve...;inv.pwm_ac_grid_reserve...;inv.pwm_ac_neutral_reserve]) &&
            throw(ArgumentError("manual RMS reserves do not specify a peak waveform; cannot certify a peak limit"))
    end
    if inv.pwm_current_envelope
        inv.v_dc_max!==nothing || throw(ArgumentError("current envelope requires a declared v_dc_max"))
        !inv.pwm_ac_ripple || throw(ArgumentError("choose current envelope or harmonic AC closure, not both"))
        _switching_envelope_geometry(inv)
    end
    reserve=_inverter_current_envelope(inv)
    for (limit,allowance) in ((inv.i_peak_max,maximum(reserve.phase_peak)),
            (inv.i_grid_peak_max,maximum(reserve.phase_peak)),(inv.In_peak_max,reserve.neutral_peak))
        limit===nothing || limit>allowance || throw(ArgumentError("carrier peak bound exhausts the declared peak-current limit"))
    end
    if inv.pwm_dc_voltage_envelope
        inv.pwm_current_envelope && inv.dc_energy_model===:nonlinear ||
            throw(ArgumentError("DC voltage envelope requires nonlinear energy and the AC current envelope"))
        inv.v_dc_min!==nothing && inv.v_dc_max!==nothing ||
            throw(ArgumentError("DC voltage envelope requires declared lower and upper DC voltage limits"))
        inv.pwm_dc_source_r===nothing && inv.pwm_dc_source_l==0 ||
            throw(ArgumentError("DC voltage envelope requires an open carrier-frequency source"))
        reserve.peak_resolved && iszero(inv.i_sw) && iszero(inv.pwm_current_factor) ||
            throw(ArgumentError("RMS-only AC/DC allowances cannot establish the carrier voltage peak"))
    end
    nothing
end

"""
    _dc_switching_voltage_bound(inv, ir, ii; norm_epsilon=0)

Uniform carrier-voltage reserve, independent of the thermal calibration factor.
All leg currents sum to zero, so the gate-current range for frozen fundamental
currents is half their absolute sum. The primitive range bound W*T/4 gives
`sqrt(2)*sum(abs(I_leg))/(8*f_sw*C)` over the line cycle. Bounded AC ripple adds
a gate-current variation range no greater than the sum of leg peak bounds,
giving a separate `sum(peak_leg)/(4*f_sw*C)` term. The latter also covers the
mean-subtracted part of the gate/ripple-current product without assuming that
product is frozen or orthogonal to the fundamental carrier term.

This is a periodic carrier approximation, not an all-time transient enclosure.
At constant nonnegative input power, the DC source's incremental voltage term
is dissipative. The waveform mean and slow energy balance remain those of the
declared mean-loss model. NLP norms use a positive, one-sided smoothing bias;
the independent SI reconstruction uses exact magnitudes.
"""
function _dc_switching_voltage_bound(inv,ir,ii;norm_epsilon=0.)
    magnitudes=[sqrt(ir[k]^2+ii[k]^2+norm_epsilon^2) for k in 1:3]
    inv.topology===:FOUR_LEG && push!(magnitudes,
        sqrt(sum(ir)^2+sum(ii)^2+norm_epsilon^2))
    bounds=_inverter_current_envelope(inv)
    frozen=sqrt(2)*sum(magnitudes)/(8inv.f_sw*inv.c_dc)
    ac_variation=(sum(bounds.phase_peak)+bounds.neutral_peak)/(4inv.f_sw*inv.c_dc)
    (frozen=frozen,ac_variation=ac_variation,total=frozen+ac_variation)
end

# Integrate exact centered-PWM events, not sampled gate states or harmonics.
# AC is the lossless inductive carrier response. DC uses frozen leg currents.
# Neither result is a complete coupled switched-circuit solution.
function _carrier_event_ripple(duties,bus,leg_current,geometry,period,capacitance)
    length(duties)==length(leg_current)==size(geometry.gate_map,2) ||
        throw(DimensionMismatch("carrier leg dimensions differ"))
    all(isfinite,duties) && all(d->0<=d<=1,duties) ||
        throw(ArgumentError("carrier duties must lie in [0,1]; no clipping"))
    all(isfinite,leg_current) && isfinite(bus) && bus>0 &&
        isfinite(period) && period>0 && isfinite(capacitance) && capacitance>0 ||
        throw(ArgumentError("finite carrier data and positive bus, period and capacitance required"))
    abs(sum(leg_current))<=1e-10*max(1.,sum(abs,leg_current)) ||
        throw(ArgumentError("all bridge-leg currents must sum to zero"))
    rows=[geometry.gate_map;geometry.neutral_gate_map']
    x=zeros(5);first=zeros(5);second=zeros(5);lo=zeros(5);hi=zeros(5)
    events=sort!(unique!([0.;(1 .-duties)./2;(1 .+duties)./2;1.]))
    for j in 1:length(events)-1
        left,right=events[j:j+1];h=(right-left)*period
        middle=(right+left)/2
        gate=Float64.((middle .>= (1 .-duties)./2).&(middle .< (1 .+duties)./2))
        error=gate.-duties
        slope=[bus.*(rows*error);-sum(error.*leg_current)/capacitance]
        first .+= x .* h .+ slope .* (h^2/2)
        second .+= x.^2 .* h .+ x .* slope .* h^2 .+ slope.^2 .* (h^3/3)
        x .+= slope.*h
        lo=min.(lo,x);hi=max.(hi,x)
    end
    mean=first./period
    rms=sqrt.(max.(0.,second./period.-mean.^2))
    (phase_rms=rms[1:3],neutral_rms=rms[4],dc_rms=rms[5],
        phase_peak=max.(abs.(lo[1:3].-mean[1:3]),abs.(hi[1:3].-mean[1:3])),
        neutral_peak=max(abs(lo[4]-mean[4]),abs(hi[4]-mean[4])),
        dc_min=lo[5]-mean[5],dc_max=hi[5]-mean[5],
        phase_pp=hi[1:3].-lo[1:3],closure=x)
end
