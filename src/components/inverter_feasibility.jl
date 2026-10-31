# Independent SI checks of the declared fundamental-frequency inverter model.
# This is a physical-model certificate, not a device/grid-code certification.

"""Physical acceptance tolerances in SI, independent of the solver's NLP tolerance."""
Base.@kwdef struct InverterFeasibilityTolerances
    current_A::Float64 = 1e-5
    voltage_V::Float64 = 1e-5
    power_VA::Float64 = 0.01
    charge_C::Float64 = 1e-9
    relative::Float64 = 1e-8
    rail_enclosure_V::Float64 = 1e-6
end

"""One independently recomputed bound, equality or closure check; positive margin is feasible."""
struct InverterLimitCheck
    name::String
    unit::Symbol
    value::Float64
    lower::Float64
    upper::Float64
    tolerance::Float64
    margin::Float64
    passed::Bool
end

"""Independent SI-model audit. An unchecked report never passes. See `checks` for failures."""
struct InverterFeasibilityReport
    checked::Bool
    passed::Bool
    checks::Vector{InverterLimitCheck}
    rail::NamedTuple
end

"""
    physical_feasibility(result)

Return the independent inverter SI-model report (a dictionary for a fleet).
Solver success and physical acceptance are separate evidence. A rejected
candidate can retain diagnostic quantities; only `solve_status(result).publishable`
permits using it as an accepted result. The certificate covers the declared
fundamental model and calibrated reserves. Standalone PWM results additionally
check finite carrier predictions, reserve closure and actual thermal/total RMS
ratings. Network feasibility, unmodeled harmonics, regional compliance and
dynamic converter stability require separate evidence.
"""
function physical_feasibility end

_unchecked_inverter_feasibility() = InverterFeasibilityReport(false,false,
    InverterLimitCheck[], (lower_V=NaN,upper_V=NaN,angle_rad=NaN,
        branch="unchecked",converged=false,evaluations=0))

function _validate_feasibility_tolerances(t)
    all(x -> isfinite(x) && x >= 0, (t.current_A,t.voltage_V,t.power_VA,
        t.charge_C,t.relative)) && isfinite(t.rail_enclosure_V) && t.rail_enclosure_V > 0 ||
        throw(ArgumentError("physical tolerances must be finite and nonnegative; rail enclosure must be positive"))
end

function _inverter_limit_check(name,unit,value,lower,upper,t=InverterFeasibilityTolerances())
    base=unit===:A ? t.current_A : unit===:V ? t.voltage_V : unit===:C ? t.charge_C : t.power_VA
    scale=maximum(abs(x) for x in (lower,upper) if isfinite(x);init=0.)
    tolerance=base+t.relative*scale
    margin=min(value-lower,upper-value)
    InverterLimitCheck(name,unit,value,lower,upper,tolerance,margin,
        isfinite(value) && margin>=-tolerance)
end

function _extend_inverter_feasibility(report,checks)
    report.checked || return report
    InverterFeasibilityReport(true,report.passed&&all(c.passed for c in checks),
        [report.checks;checks],report.rail)
end

# Enclose min f(theta), theta in [0,2pi], for a degree-two trigonometric
# polynomial. Taylor's theorem gives a global interval lower bound from
# |f''| <= hypot(c1,s1)+4hypot(c2,s2). Refinement targets only intervals that
# can improve the best sampled upper bound. No stationary-root or sample-grid
# completeness assumption is used. Float64 roundoff is padded explicitly;
# this is not an interval-arithmetic proof of the hardware waveform.
function _trig_minimum(c0,c1,s1,c2,s2; atol=1e-6,max_intervals=100000)
    all(isfinite,(c0,c1,s1,c2,s2)) || return (
        lower=-Inf,upper=Inf,angle=NaN,converged=false,evaluations=0)
    isfinite(atol) && atol>0 || throw(ArgumentError("atol must be positive and finite"))
    max_intervals >= 16 || throw(ArgumentError("max_intervals must be at least 16"))
    curvature=hypot(c1,s1)+4hypot(c2,s2)
    pad=64eps(Float64)*max(1.,abs(c0)+abs(c1)+abs(s1)+abs(c2)+abs(s2))
    f(t)=c0+c1*cos(t)+s1*sin(t)+c2*cos(2t)+s2*sin(2t)
    if !isfinite(curvature) || !isfinite(pad) || 2pad>atol
        return (lower=c0-hypot(c1,s1)-hypot(c2,s2)-pad,
            upper=f(0.)+pad,angle=0.,converged=false,evaluations=1)
    end
    df(t)=-c1*sin(t)+s1*cos(t)-2c2*sin(2t)+2s2*cos(2t)
    curvature==0 && return (lower=c0-pad,upper=c0+pad,angle=0.,
        converged=2pad<=atol,evaluations=1)
    intervals=Tuple{Float64,Float64,Float64}[]
    best,angle=Inf,0.
    evaluations=0
    function interval(lo,hi)
        mid=(lo+hi)/2; half=(hi-lo)/2
        value=f(mid); evaluations+=1
        if value < best
            best,angle=value,mid
        end
        lower=value-abs(df(mid))*half-curvature*half^2/2-pad
        (lo,hi,lower)
    end
    for k in 0:15
        push!(intervals,interval(2pi*k/16,2pi*(k+1)/16))
    end
    while true
        lower,index=findmin(x[3] for x in intervals)
        best+pad-lower <= atol && return (
            lower=lower,upper=best+pad,angle=angle,converged=true,evaluations=evaluations)
        (length(intervals)>=max_intervals || evaluations>=2max_intervals) && return (
            lower=lower,upper=best+pad,angle=angle,converged=false,evaluations=evaluations)
        lo,hi,_=intervals[index]
        mid=(lo+hi)/2
        left,right=interval(lo,mid),interval(mid,hi)
        intervals[index]=left
        push!(intervals,right)
        # Discard intervals already proven unable to contain a better minimum.
        # Retaining bounds within atol also keeps the final enclosure explicit.
        filter!(x->x[3]<=best+pad,intervals)
    end
end

"""
    continuous_rail_margin(inv, voltage, ripple=0im;
        midpoint=0im, midpoint_mean=0, bus_fundamental=0im, atol=1e-6)

Enclose the minimum instantaneous rail headroom over the full electrical cycle.
Voltage/midpoint/bus-fundamental phasors are RMS volts; `ripple` is the peak
2ω DC-bus phasor. The result gives conservative `lower_V`, sampled `upper_V`,
a worst-branch angle and `converged`. It covers the declared degree-two averaged
rail waveform, including unequal split-link banks. It does not cover unmodeled
harmonics or switching transients. SINGLE_PHASE has no three-phase rail model.
"""
function continuous_rail_margin(inv,voltage,ripple=0im;
        midpoint=0im,midpoint_mean=0.,bus_fundamental=0im,atol=1e-6,energy=nothing)
    if inv.dc_energy_model===:nonlinear
        energy===nothing && throw(ArgumentError("nonlinear rail checks require capacitor-energy evidence; use inverter_physical_feasibility or supply energy"))
        return _dc_energy_rail_margin(inv,voltage,energy)
    end
    inv.topology===:SINGLE_PHASE && return (lower_V=Inf,upper_V=Inf,angle_rad=NaN,
        branch="not_applicable",converged=true,evaluations=0)
    length(voltage)==3 || throw(DimensionMismatch("rail check requires three phasors"))
    rows=Tuple{String,NTuple{5,Float64}}[]
    m=inv.m_max; dr,di=reim(ripple)
    if inv.topology in (:THREE_LEG,:FOUR_LEG)
        for (a,b) in ((1,2),(2,3),(3,1)), sg in (-1.,1.)
            v=voltage[a]-voltage[b]
            push!(rows,("line_$(a)_$(b)_$(sg)",
                (m*inv.v_dc,-sg*sqrt(2)*real(v),sg*sqrt(2)*imag(v),m*dr,-m*di)))
        end
        if inv.topology===:FOUR_LEG
            for k in 1:3,sg in (-1.,1.)
                v=voltage[k]
                push!(rows,("phase_$(k)_$(sg)",
                    (m*inv.v_dc,-sg*sqrt(2)*real(v),sg*sqrt(2)*imag(v),m*dr,-m*di)))
            end
        end
    elseif inv.topology===:SPLIT_DC
        ratio=_split_midpoint_2w_ratio(inv)
        for k in 1:3,sg in (-1.,1.)
            v=voltage[k]+midpoint
            push!(rows,("half_$(k)_$(sg)",
                (m*inv.v_dc/2-sg*midpoint_mean,
                 m/sqrt(2)*real(bus_fundamental)-sg*sqrt(2)*real(v),
                 -m/sqrt(2)*imag(bus_fundamental)+sg*sqrt(2)*imag(v),
                 (m/2-sg*ratio)*dr,(-m/2+sg*ratio)*di)))
        end
    else
        throw(ArgumentError("unsupported rail topology $(inv.topology)"))
    end
    bounds=[_trig_minimum(c...;atol=atol) for (_,c) in rows]
    index=argmin([b.upper for b in bounds])
    return (lower_V=minimum(b.lower for b in bounds),upper_V=bounds[index].upper,
        angle_rad=bounds[index].angle,branch=rows[index][1],
        converged=all(b.converged for b in bounds),evaluations=sum(b.evaluations for b in bounds))
end

"""
    inverter_physical_feasibility(inv; converter_voltage, converter_current,
        grid_voltage, grid_current=converter_current, filter_voltage=grid_voltage,
        midpoint_charge=0, dc_min=nothing, dc_max=nothing,
        tolerances=InverterFeasibilityTolerances())

Recompute engineering bounds and filter/topology equalities from SI RMS phasors,
without using the NLP's squared-current, ripple or power auxiliaries. Includes
continuous rail feasibility, source-dependent ripple, simultaneous capacitor
thermal budgets and the exact fitted converter-loss law. `midpoint_charge` is
signed coulombs. Optional `dc_min`/`dc_max` are watts; wrappers supply request/source bounds.
The optional `reported` values audit consistency of the NLP's SI quantities.
"""
function inverter_physical_feasibility(inv; converter_voltage,converter_current,
        grid_voltage,grid_current=converter_current,filter_voltage=grid_voltage,
        midpoint_charge=0.,reported=(;),dc_min=nothing,dc_max=nothing,
        tolerances=InverterFeasibilityTolerances(),pwm_thermal_headroom=nothing)
    t=tolerances; _validate_feasibility_tolerances(t)
    validate_device(inv)
    n=length(inv.phase_terminals)
    all(length(x)==n for x in (converter_voltage,converter_current,grid_voltage,grid_current,filter_voltage)) ||
        throw(DimensionMismatch("audit phasors must match inverter phase count"))
    u,ic,v,ig,vm=map(x->ComplexF64.(x),
        (converter_voltage,converter_current,grid_voltage,grid_current,filter_voltage))
    checks=InverterLimitCheck[]
    function check(name,unit,value,lower,upper)
        push!(checks,_inverter_limit_check(name,unit,value,lower,upper,t))
    end
    upper(name,unit,value,rating)=rating===nothing ? nothing : check(name,unit,value,-Inf,rating)
    zero(name,unit,value)=check(name,unit,abs(value),0.,0.)
    current_envelope=_inverter_current_envelope(inv)
    for k in 1:n
        upper("converter_current_$k",:A,hypot(abs(ic[k]),current_envelope.converter_rms[k]),inv.i_max)
        upper("grid_current_$k",:A,hypot(abs(ig[k]),current_envelope.grid_rms[k]),inv.i_grid_max)
        upper("converter_peak_current_$k",:A,sqrt(2)*abs(ic[k])+current_envelope.phase_peak[k],inv.i_peak_max)
        upper("grid_peak_current_$k",:A,sqrt(2)*abs(ig[k])+current_envelope.phase_peak[k],inv.i_grid_peak_max)
    end
    ineutral=-sum(ic)
    if inv.topology!==:SINGLE_PHASE
        upper("neutral_current",:A,hypot(abs(ineutral),current_envelope.neutral_rms),inv.In_max)
        upper("neutral_peak_current",:A,sqrt(2)*abs(ineutral)+current_envelope.neutral_peak,inv.In_peak_max)
        a=cis(2pi/3)
        i0=sum(ic)/3; i1=(ic[1]+a*ic[2]+a^2*ic[3])/3; i2=(ic[1]+a^2*ic[2]+a*ic[3])/3
        for (name,value,rating) in (("zero",i0,inv.i_zero_max),("positive",i1,inv.i_positive_max),("negative",i2,inv.i_negative_max))
            upper("$(name)_sequence_current",:A,abs(value),rating)
        end
        inv.topology===:THREE_LEG && zero("three_leg_zero_current",:A,ineutral)
    end
    # SI conductor-domain KVL, including mutual and neutral-conductor drops.
    for (side,current,drop) in ((:converter,ic,u-vm),(:grid,ig,vm-v))
        R,X=_series_filter_matrices(inv,n,1.;side=side)
        conductor=inv.neutral===nothing ? current : [current;-sum(current)]
        primitive=complex.(R,X)*conductor
        expected=inv.neutral===nothing ? primitive : primitive[1:n].-primitive[end]
        for k in 1:n
            zero("$(side)_filter_kvl_$k",:V,drop[k]-expected[k])
        end
    end
    gm,bm=_mid_filter_admittance(inv)
    for k in 1:n
        zero("filter_kcl_$k",:A,ic[k]-ig[k]-complex(gm,bm)*vm[k])
    end
    if inv.grid_forming && n==3
        a=cis(2pi/3)
        u1=(u[1]+a*u[2]+a^2*u[3])/3; u2=(u[1]+a^2*u[2]+a*u[3])/3
        zero("grid_forming_negative_voltage",:V,u2)
        inv.topology!==:THREE_LEG && zero("grid_forming_zero_voltage",:V,sum(u)/3)
        inv.v_int_min===nothing || check("internal_positive_voltage_min",:V,abs(u1),inv.v_int_min,Inf)
        upper("internal_positive_voltage_max",:V,abs(u1),inv.v_int_max)
    else
        for k in 1:n
            inv.v_int_min===nothing || check("internal_voltage_min_$k",:V,abs(u[k]),inv.v_int_min,Inf)
            upper("internal_voltage_max_$k",:V,abs(u[k]),inv.v_int_max)
        end
    end
    if inv.topology===:SINGLE_PHASE && inv.modulation_max!==nothing && inv.v_dc!==nothing
        for k in 1:n
            upper("single_phase_modulation_$k",:V,abs(u[k]),inv.modulation_max*inv.v_dc/sqrt(3))
        end
    end
    sconv=sum(u.*conj.(ic))
    spoc=sum(v.*conj.(ig))+im*inv.b_filter_shunt*sum(abs2,v)
    upper("$(inv.power_limit_location)_apparent_power",:VA,
        abs(inv.power_limit_location===:converter ? sconv : spoc),inv.s_max)
    ripple=sum(u.*ic)
    upper("single_phase_power_ripple",:VA,abs(ripple),inv.p_ripple_max)
    loss_currents=inv.topology===:FOUR_LEG ? [ic;ineutral] : ic
    loss=inv.p_loss_fixed+inv.a_loss*sum(abs,loss_currents)+inv.c_loss*sum(abs2,loss_currents)
    caploss=0.; d=0im; midpoint=0im; bus=0im; mean=0.;energy=nothing
    if inv.topology!==:SINGLE_PHASE
        d=-ripple*_dc_ripple_factor(inv)/inv.v_dc
        if inv.dc_energy_model===:nonlinear
            energy=_dc_energy_state(inv,ripple)
            d=energy.voltage_2omega
            check("dc_energy_quadrature",:V,energy.quadrature_error_V,0.,1e-6)
            check("dc_energy_positive_voltage",:V,energy.minimum_V,0.,Inf)
            check("dc_energy_maximum_voltage",:V,energy.maximum_V,0.,Inf)
            inv.v_dc_min===nothing || check("dc_voltage_min",:V,energy.minimum_V,inv.v_dc_min,Inf)
            upper("dc_voltage_max",:V,energy.maximum_V,inv.v_dc_max)
        end
        upper("dc_voltage_ripple",:V,abs(d),inv.dv2_max)
        # I_C,2w = j*2w*C_eq*D: phasor D is peak, current rating is RMS.
        i2sq=abs2(2*(2pi*inv.f)*_dc_link_capacitance(inv)*d)/2
        energy===nothing || (i2sq=energy.capacitor_lf_sq)
        if energy!==nothing || pwm_thermal_headroom!==nothing
            check("capacitor_low_frequency_current",:A,sqrt(i2sq),0.,Inf)
        end
        wn,w2,ws=inv.cap_thermal_weights
        switching=inv.i_sw^2+inv.pwm_current_factor^2*sum(abs2,loss_currents)
        if pwm_thermal_headroom!==nothing
            _validate_allocation_pwm(inv,pwm_thermal_headroom)
            append!(checks,_allocation_pwm_checks(inv,u,ic,d,i2sq,switching;headroom=pwm_thermal_headroom,energy))
            padded=energy===nothing ? nothing : begin
                delta=2inv.v_dc*sqrt(abs2(energy.energy_coordinate)+_PWM_ALLOCATION_BUS_EPS^2)
                (sqrt(max(0.,energy.center_V2-delta)),sqrt(energy.center_V2+delta))
            end
            switching+=max(0.,_allocation_pwm_variance(inv,real.(u),imag.(u),real.(ic),imag.(ic),
                real(d),imag(d),real(sconv),real(ripple),imag(ripple);headroom=pwm_thermal_headroom,voltage_bounds=padded))
        end
        if inv.topology===:SPLIT_DC
            share=_split_fundamental_coefficients(inv)
            midpoint=share.midpoint*(-ineutral); bus=share.bus*(-ineutral)
            cu,cl=_split_capacitances(inv)
            mean=inv.v_dc*(cu-cl)/(2*(cu+cl))+2midpoint_charge/(cu+cl)
            upper("midpoint_fundamental_ripple",:V,abs(midpoint),inv.dv_mid_max)
            upper("midpoint_mean",:V,abs(mean),inv.v_mid_mean_max)
            upper("midpoint_balance_charge",:C,abs(midpoint_charge),inv.q_mid_balance_max)
            inv.q_mid_balance_max===nothing && zero("no_midpoint_charge_actuator",:C,midpoint_charge)
            thermal=[w2*i2sq+wn*weight*abs2(ineutral)+ws*switching for weight in (share.upper_weight,share.lower_weight)]
            for (j,label,rating) in ((1,"upper",inv.i_cap_upper_max),(2,"lower",inv.i_cap_lower_max))
                upper("$(label)_capacitor_current",:A,sqrt(thermal[j]),inv.i_cap_max)
                upper("$(label)_capacitor_own_rating",:A,sqrt(thermal[j]),rating)
            end
            caploss=inv.esr_dc_upper*thermal[1]+inv.esr_dc_lower*thermal[2]
        else
            thermal=w2*i2sq+ws*switching
            upper("capacitor_current",:A,sqrt(thermal),inv.i_cap_max)
            caploss=inv.esr_dc*thermal
        end
    end
    derived=(p_conv=real(sconv),q_conv=imag(sconv),p_poc=real(spoc),q_poc=imag(spoc),
        p_dc=real(sconv)+loss+caploss,dc_ripple=d,midpoint_ripple=midpoint,midpoint_mean=mean)
    upper("dc_source_power_max",:VA,derived.p_dc,dc_max)
    dc_min===nothing || check("dc_source_power_min",:VA,derived.p_dc,dc_min,Inf)
    for name in keys(reported)
        unit=name in (:dc_ripple,:midpoint_ripple,:midpoint_mean) ? :V : :VA
        zero("reported_$(name)_residual",unit,getproperty(reported,name)-getproperty(derived,name))
    end
    rail=continuous_rail_margin(inv,u,d;midpoint=midpoint,midpoint_mean=mean,
        bus_fundamental=bus,atol=t.rail_enclosure_V,energy)
    if inv.topology!==:SINGLE_PHASE
        check("continuous_rail_headroom",:V,rail.lower_V,0.,Inf)
    end
    InverterFeasibilityReport(true,all(c.passed for c in checks)&&rail.converged,checks,rail)
end

# The solver outcome remains visible; physical rejection changes publication,
# never the recorded MOI termination/primal statuses or optimality claim.
function _physical_solve_status(status::SolveStatus,report::InverterFeasibilityReport)
    accepted=report.checked && report.passed
    SolveStatus(status.termination_status,status.primal_status,status.has_primal,
        status.feasible,status.optimal,status.publishable && accepted)
end
