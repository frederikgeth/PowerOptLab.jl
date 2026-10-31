# Shared-carrier DC-current variance, embedded in local allocation objectives.
# No point-calibrated current factor is transferred between dispatches.
const _PWM_ALLOCATION_HEADROOM = 0.10
const _PWM_ALLOCATION_VOLTAGE_EPS = 1e-3 # V; explicit analytic error enclosure
const _PWM_ALLOCATION_BUS_EPS = 1e-2 # V; conservative ripple enclosure at balance
const _PWM_SERVICE_PRESERVATION_VA = 1e-4 # SI Q/P budget, not a hardware tolerance

function _lock_pwm_service!(m,coordinate,value,base)
    # An exact tangent lock at the thermal maximum can violate constraint
    # qualification. Preserve service within a declared 100 µVA budget.
    JuMP.set_lower_bound(coordinate,max(0.,value-_PWM_SERVICE_PRESERVATION_VA/base))
    JuMP.set_upper_bound(coordinate,value)
end

function _validate_allocation_pwm(inv,headroom)
    isfinite(headroom) && headroom>=0 || throw(ArgumentError("PWM thermal headroom must be finite and nonnegative"))
    inv.pwm_strategy===:NONE && return
    inv.pwm_strategy===:CENTERED || throw(ArgumentError("embedded PWM allocation currently requires :CENTERED"))
    inv.topology in (:THREE_LEG,:FOUR_LEG) || throw(ArgumentError("embedded PWM allocation requires a monolithic three/four-leg link"))
    inv.pwm_dc_source_r===nothing && iszero(inv.pwm_dc_source_l) ||
        throw(ArgumentError("embedded PWM allocation requires an open carrier-frequency DC source"))
    !inv.pwm_ac_ripple || throw(ArgumentError("automatic AC-ripple reserve closure is separate from embedded DC thermal allocation"))
    _switching_envelope_geometry(inv)
end

# For RMS phasors z=Uj-Uk, let a_j(t)=sqrt(2)Re(Ij exp(jwt)). Then
# <a_j a_k |sqrt(2)Re(z exp(jwt))|> = 4sqrt(2)/(3pi) *
# (|z| Re(Ij conj(Ik)) + Re(Ij conj(z))*Re(Ik conj(z))/|z|).
# Regularizing |z| to R=sqrt(|z|²+eps²) changes this by at most
# K*eps*(|Ij|²+|Ik|²). Adding that bound after negation encloses M2.
function _allocation_pwm_variance(inv,ur,ui,ir,ii,dr,di,pconv,s2r,s2i;headroom=_PWM_ALLOCATION_HEADROOM,voltage_bounds=nothing)
    er=collect(ur);ei=collect(ui);ar=collect(ir);ai=collect(ii)
    if inv.topology===:FOUR_LEG
        push!(er,0.);push!(ei,0.);push!(ar,-sum(ir));push!(ai,-sum(ii))
    end
    epsilon=_PWM_ALLOCATION_VOLTAGE_EPS;K=4sqrt(2)/(3pi)
    second=0.
    for j in 1:length(ar)-1,k in j+1:length(ar)
        zr=er[j]-er[k];zi=ei[j]-ei[k]
        radius=sqrt(zr^2+zi^2+epsilon^2)
        dot=ar[j]*ar[k]+ai[j]*ai[k]
        pj=ar[j]*zr+ai[j]*zi;pk=ar[k]*zr+ai[k]*zi
        second+=K*(-radius*dot-pj*pk/radius+
            epsilon*(ar[j]^2+ai[j]^2+ar[k]^2+ai[k]^2))
    end
    ripple=sqrt(dr^2+di^2+_PWM_ALLOCATION_BUS_EPS^2)
    # M2 is nonnegative. Bound its denominator downwards and the subtracted
    # mean-square denominator upwards over the complete line cycle.
    lo,hi=voltage_bounds===nothing ? (inv.v_dc-ripple,inv.v_dc+ripple) : voltage_bounds
    upper=second/lo-(pconv^2+(s2r^2+s2i^2)/2)/hi^2
    _allocation_coupled_pwm_variance(inv,(1+headroom)^2*upper,hi)
end

# Independent numerical integration of exact carrier dwell probabilities.
# Common-mode injection cancels because all bridge-leg currents sum to zero.
function _allocation_pwm_carrier(inv,u,i,d;angles=2048,energy=nothing)
    total=0.;ac_total=0.;bus_max=0.;margin=Inf
    geometry=_switching_envelope_geometry(inv)
    K=inv.topology===:THREE_LEG ? geometry.gate_map :
        [geometry.gate_map;geometry.neutral_gate_map']
    for k in 0:angles-1
        th=2pi*(k+.5)/angles;rot=cis(th)
        v=energy===nothing ? inv.v_dc+real(d*rot^2) :
            sqrt(max(0.,energy.center_V2+2inv.v_dc*real(energy.energy_coordinate*rot^2)))
        v>0 || return (rms=NaN,modulation_margin=-Inf)
        voltage=sqrt(2).*real.(u.*rot);current=sqrt(2).*real.(i.*rot)
        if inv.topology===:FOUR_LEG
            voltage=[voltage;0.];current=[current;-sum(current)]
        end
        common=(maximum(voltage)+minimum(voltage))/2
        duties=.5 .+(voltage.-common)./v
        margin=min(margin,inv.m_max*v-(maximum(voltage)-minimum(voltage)))
        all(x->0<=x<=1,duties) || return (rms=NaN,modulation_margin=margin)
        average=sum(duties.*current)
        moment=sum(min(duties[j],duties[l])*current[j]*current[l]
            for j in eachindex(current),l in eachindex(current))
        total+=max(0.,moment-average^2)
        ac_total+=v^2/inv.f_sw^2*_carrier_ac_bridge_variance(duties,K)
        bus_max=max(bus_max,v)
    end
    bound=_carrier_thermal_bound(inv)
    (rms=sqrt((total+ac_total)/angles)+bound.damping_rms_per_V*bus_max,
        frozen_rms=sqrt(total/angles),ac_lossless_rms=sqrt(ac_total/angles),
        damping_reserve=bound.damping_rms_per_V*bus_max,modulation_margin=margin)
end

function _allocation_pwm_checks(inv,u,i,d,i2sq,base_switching;headroom=_PWM_ALLOCATION_HEADROOM,energy=nothing)
    carrier=_allocation_pwm_carrier(inv,u,i,d;energy)
    coarse=_allocation_pwm_carrier(inv,u,i,d;angles=1024,energy)
    wn,w2,ws=inv.cap_thermal_weights
    actual=sqrt(max(0.,w2*i2sq+ws*(base_switching+carrier.rms^2)))
    checks=InverterLimitCheck[
        _inverter_limit_check("allocation_pwm_carrier_prediction",:A,carrier.rms,0.,Inf),
        _inverter_limit_check("allocation_pwm_frozen_current",:A,carrier.frozen_rms,0.,Inf),
        _inverter_limit_check("allocation_pwm_ac_induced_current",:A,carrier.ac_lossless_rms,0.,Inf),
        _inverter_limit_check("allocation_pwm_damping_reserve",:A,carrier.damping_reserve,0.,Inf),
        _inverter_limit_check("allocation_pwm_modulation",:V,carrier.modulation_margin,0.,Inf),
        _inverter_limit_check("allocation_pwm_quadrature",:A,abs(carrier.rms-coarse.rms),0.,1e-4)]
    inv.i_cap_max===nothing || push!(checks,
        _inverter_limit_check("allocation_pwm_actual_capacitor_current",:A,actual,0.,inv.i_cap_max))
    checks
end
