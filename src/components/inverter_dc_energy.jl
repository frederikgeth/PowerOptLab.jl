# Constant input power, fixed mean converter loss, monolithic capacitor energy.
function _validate_dc_energy_model(inv)
    inv.dc_energy_model in (:linearized,:nonlinear) || throw(ArgumentError("unknown DC energy model"))
    if inv.dc_energy_model===:linearized
        inv.v_dc_min===nothing && inv.v_dc_max===nothing || throw(ArgumentError("explicit DC extrema bounds require dc_energy_model=:nonlinear"))
        return
    end
    inv.topology in (:THREE_LEG,:FOUR_LEG) || throw(ArgumentError("nonlinear DC energy currently requires a monolithic three/four-leg link"))
    inv.dc_source_r===nothing && inv.dc_source_l==0 ||
        throw(ArgumentError("nonlinear DC energy uses constant input power; a finite small-signal source branch is not supported"))
    inv.v_dc!==nothing && inv.c_dc!==nothing || throw(ArgumentError("nonlinear energy requires v_dc and c_dc"))
    for (label,value) in ((:v_dc_min,inv.v_dc_min),(:v_dc_max,inv.v_dc_max))
        value===nothing || (isfinite(value) && value>0) || throw(ArgumentError("$label must be positive and finite"))
    end
    something(inv.v_dc_min,0.)<=inv.v_dc<=something(inv.v_dc_max,Inf) || throw(ArgumentError("mean DC voltage must lie within its extrema bounds"))
end

function _dc_energy_state(inv,s2)
    v=inv.v_dc;c=inv.c_dc;omega=2pi*inv.f
    d=im*s2/(2omega*c*v);amplitude=2v*abs(d)
    lower=max(v^2,amplitude*(1+1e-12));upper=v^2+amplitude
    mean_at(center,n)=sum(sqrt(max(0.,center+amplitude*cos(2pi*k/n))) for k in 0:n-1)/n
    if mean_at(lower,256)>v+1e-8
        return (resolved=false,center_V2=NaN,energy_coordinate=d,voltage_2omega=ComplexF64(NaN,NaN),
            minimum_V=NaN,maximum_V=NaN,capacitor_lf_sq=NaN,quadrature_error_V=Inf)
    end
    for _ in 1:52
        mid=(lower+upper)/2
        mean_at(mid,256)<v ? (lower=mid) : (upper=mid)
    end
    center=(lower+upper)/2
    actual=amplitude==0 ? 0im : (2sum(sqrt(center+amplitude*cos(2pi*k/1024))*cos(2pi*k/1024) for k in 0:1023)/1024)*(d/abs(d))
    rad=sqrt(max(0.,center^2-amplitude^2))
    (resolved=true,center_V2=center,energy_coordinate=d,voltage_2omega=actual,
        minimum_V=sqrt(max(0.,center-amplitude)),maximum_V=sqrt(center+amplitude),
        capacitor_lf_sq=abs2(s2)/(center+rad),quadrature_error_V=abs(mean_at(center,1024)-v))
end

"""
    dc_link_energy_diagnostics(inv; converter_voltage, converter_current)

Nonlinear monolithic DC-capacitor diagnostics at constant input power and fixed
mean converter loss. The energy offset is solved so the cycle mean voltage is
exactly `v_dc`, rather than setting mean energy to `C*v_dc^2/2`. Returns actual
voltage extrema, its peak 2ω coefficient and total low-frequency capacitor RMS.
Not a finite-source, split-bank, switching, protection or stability certificate.
"""
function dc_link_energy_diagnostics(inv;converter_voltage,converter_current)
    validate_device(inv)
    inv.topology in (:THREE_LEG,:FOUR_LEG) || throw(ArgumentError("monolithic three/four-leg link required"))
    inv.dc_source_r===nothing && iszero(inv.dc_source_l) || throw(ArgumentError("constant-power DC source required"))
    length(converter_voltage)==length(converter_current)==3 || throw(DimensionMismatch("three phasors required"))
    all(isfinite,[converter_voltage;converter_current]) || throw(ArgumentError("finite converter phasors required"))
    s2=sum(converter_voltage.*converter_current);state=_dc_energy_state(inv,s2)
    carrier=inv.pwm_dc_voltage_envelope ? _dc_switching_voltage_bound(inv,
        real.(converter_current),imag.(converter_current)).total : nothing
    merge(state,(mean_voltage_V=inv.v_dc,mean_energy_J=inv.c_dc*state.center_V2/2,
        switching_voltage_bound_V=carrier,
        combined_minimum_V=carrier===nothing ? nothing : state.minimum_V-carrier,
        combined_maximum_V=carrier===nothing ? nothing : state.maximum_V+carrier,
        energy_2omega_peak_J=inv.c_dc*inv.v_dc*state.energy_coordinate,
        capacitor_lf_rms_A=sqrt(state.capacitor_lf_sq)))
end

function _stamp_dc_energy!(m,inv,er,ei;switching_reserve=0.)
    v=inv.v_dc;n=64
    center=JuMP.@variable(m,lower_bound=.5,upper_bound=3.,start=1.,base_name="dc_energy_center_$(inv.id)")
    radius=JuMP.@variable(m,lower_bound=.01,start=1.,base_name="dc_energy_radius_$(inv.id)")
    JuMP.@constraint(m,radius^2==center^2-4*(er^2+ei^2)/v^2)
    samples=JuMP.@variable(m,[1:n],lower_bound=.01,start=1.,base_name="dc_voltage_$(inv.id)")
    for k in 1:n
        th=2pi*(k-1)/n
        JuMP.@constraint(m,samples[k]^2==center+2/v*(er*cos(2th)-ei*sin(2th)))
    end
    JuMP.@constraint(m,sum(samples)/n==1.)
    real2=JuMP.@expression(m,2v/n*sum(samples[k]*cos(4pi*(k-1)/n) for k in 1:n))
    imag2=JuMP.@expression(m,-2v/n*sum(samples[k]*sin(4pi*(k-1)/n) for k in 1:n))
    if inv.v_dc_min!==nothing && !inv.pwm_dc_voltage_envelope
        bound=((inv.v_dc_min+switching_reserve)/v)^2
        JuMP.@constraint(m,center>=bound)
        JuMP.@constraint(m,4*(er^2+ei^2)/v^2<=(center-bound)^2)
    end
    if inv.v_dc_max!==nothing && !inv.pwm_dc_voltage_envelope
        bound=((inv.v_dc_max-switching_reserve)/v)^2
        JuMP.@constraint(m,center<=bound)
        JuMP.@constraint(m,4*(er^2+ei^2)/v^2<=(bound-center)^2)
    end
    low=high=nothing
    if inv.pwm_strategy!==:NONE
        low=JuMP.@variable(m,lower_bound=.01,start=1.,base_name="dc_pwm_lower_$(inv.id)")
        high=JuMP.@variable(m,lower_bound=.01,start=1.,base_name="dc_pwm_upper_$(inv.id)")
        padding=sqrt(er^2+ei^2+_PWM_ALLOCATION_BUS_EPS^2)
        JuMP.@constraint(m,low^2==center-2padding/v)
        JuMP.@constraint(m,high^2==center+2padding/v)
        if inv.pwm_dc_voltage_envelope
            JuMP.@constraint(m,low-switching_reserve/v>=.01)
            low=low-switching_reserve/v;high=high+switching_reserve/v
            # Use the positive voltage lifts directly. Squared small energy
            # gaps become nearly dependent at symmetric, tight voltage bands.
            # The existing 0.01 V energy-norm padding is a one-sided enclosure.
            lower_scale=max((v-inv.v_dc_min)/v,1e-3)
            upper_scale=max((inv.v_dc_max-v)/v,1e-3)
            JuMP.@constraint(m,(low-inv.v_dc_min/v)/lower_scale>=0.)
            JuMP.@constraint(m,(inv.v_dc_max/v-high)/upper_scale>=0.)
        end
    end
    (center=center,radius=radius,er=er,ei=ei,real2=real2,imag2=imag2,pwm_min=low,pwm_max=high,
        switching_reserve=switching_reserve)
end

function _stamp_dc_energy_rail!(m,inv,ur,ui,energy,theta)
    c,s=cos(theta),sin(theta)
    if inv.pwm_dc_voltage_envelope
        # A positive lift avoids square-root domain failures at infeasible NLP
        # trial points. These rows also apply to adaptive separating cuts.
        bus=JuMP.@variable(m,lower_bound=.01,start=1.,base_name="dc_carrier_rail_$(inv.id)")
        JuMP.@constraint(m,bus^2==energy.center+
            2/inv.v_dc*(energy.er*cos(2theta)-energy.ei*sin(2theta)))
        rail=inv.m_max*(bus-energy.switching_reserve/inv.v_dc)
        for (a,b) in ((1,2),(2,3),(3,1)),sg in (-1.,1.)
            value=sg*sqrt(2)*((ur[a]-ur[b])*c-(ui[a]-ui[b])*s)/inv.v_dc
            JuMP.@constraint(m,value<=rail)
        end
        if inv.topology===:FOUR_LEG
            for k in 1:3,sg in (-1.,1.)
                JuMP.@constraint(m,sg*sqrt(2)*(ur[k]*c-ui[k]*s)/inv.v_dc<=rail)
            end
        end
        return
    end
    squared=inv.m_max^2*(inv.v_dc^2*energy.center+
        2inv.v_dc*(energy.er*cos(2theta)-energy.ei*sin(2theta)))
    for (a,b) in ((1,2),(2,3),(3,1))
        value=sqrt(2)*((ur[a]-ur[b])*c-(ui[a]-ui[b])*s)
        JuMP.@constraint(m,value^2/inv.v_dc^2<=squared/inv.v_dc^2)
    end
    if inv.topology===:FOUR_LEG
        for k in 1:3
            value=sqrt(2)*(ur[k]*c-ui[k]*s)
            JuMP.@constraint(m,value^2/inv.v_dc^2<=squared/inv.v_dc^2)
        end
    end
end

function _dc_energy_rail_margin(inv,u,state;switching_voltage=0.,atol=1e-6)
    state.resolved || return (lower_V=-Inf,upper_V=Inf,angle_rad=NaN,branch="dc_energy_unresolved",converged=false,evaluations=0)
    if inv.pwm_dc_voltage_envelope
        return _dc_carrier_rail_margin(inv,u,state,switching_voltage;atol)
    end
    choices=[("line_$(a)_$(b)",u[a]-u[b]) for (a,b) in ((1,2),(2,3),(3,1))]
    inv.topology===:FOUR_LEG && append!(choices,[("phase_$k",u[k]) for k in 1:3])
    values=NamedTuple[];d=state.energy_coordinate;m=inv.m_max
    for (name,z) in choices
        constant=m^2*state.center_V2-abs2(z)
        cosine=m^2*2inv.v_dc*real(d)-(real(z)^2-imag(z)^2)
        sine=-m^2*2inv.v_dc*imag(d)+2real(z)*imag(z)
        norm=hypot(cosine,sine);pad=128eps(Float64)*max(1.,abs(constant)+norm)
        squared_lower=constant-norm-pad
        theta=(atan(sine,cosine)+pi)/2
        bus=sqrt(state.center_V2+2inv.v_dc*real(d*cis(2theta)))
        upper=m*bus-abs(sqrt(2)*real(z*cis(theta)))
        denom=squared_lower>=0 ? m*state.maximum_V+sqrt(2)*abs(z) : m*state.minimum_V
        push!(values,(lower=squared_lower/denom,upper=upper+pad/max(1.,denom),theta=theta,name=name))
    end
    index=argmin([r.upper for r in values])
    (lower_V=minimum(r.lower for r in values),upper_V=values[index].upper,
        angle_rad=values[index].theta,branch=values[index].name,converged=true,evaluations=length(values))
end

# A Taylor enclosure for a smooth periodic function with a global second-
# derivative bound. Unlike a dense sample, this returns both sides of min(f).
function _smooth_periodic_minimum(f,df,curvature,pad;atol=1e-6,max_intervals=100000)
    isfinite(atol) && atol>0 || throw(ArgumentError("positive finite rail enclosure tolerance required"))
    all(isfinite,(curvature,pad)) && curvature>=0 && pad>=0 && 2pad<atol ||
        return (lower=-Inf,upper=Inf,angle=NaN,converged=false,evaluations=0)
    intervals=Tuple{Float64,Float64,Float64}[];best=Inf;angle=0.;evaluations=0
    function interval(lo,hi)
        mid=(lo+hi)/2;half=(hi-lo)/2;value=f(mid);evaluations+=1
        if value<best
            best,angle=value,mid
        end
        (lo,hi,value-abs(df(mid))*half-curvature*half^2/2-pad)
    end
    for k in 0:15
        push!(intervals,interval(2pi*k/16,2pi*(k+1)/16))
    end
    while true
        lower,index=findmin(x[3] for x in intervals)
        best+pad-lower<=atol && return (
            lower=lower,upper=best+pad,angle=angle,converged=true,evaluations=evaluations)
        evaluations>=2max_intervals && return (
            lower=lower,upper=best+pad,angle=angle,converged=false,evaluations=evaluations)
        lo,hi,_=intervals[index];mid=(lo+hi)/2
        left,right=interval(lo,mid),interval(mid,hi)
        intervals[index]=left;push!(intervals,right)
        filter!(x->x[3]<=best+pad,intervals)
    end
end

function _dc_carrier_rail_margin(inv,u,state,reserve;atol=1e-6)
    isfinite(reserve) && reserve>=0 && state.minimum_V>0 ||
        return (lower_V=-Inf,upper_V=Inf,angle_rad=NaN,branch="dc_carrier_unresolved",converged=false,evaluations=0)
    choices=[("line_$(a)_$(b)",u[a]-u[b]) for (a,b) in ((1,2),(2,3),(3,1))]
    inv.topology===:FOUR_LEG && append!(choices,[("phase_$k",u[k]) for k in 1:3])
    dr,di=reim(state.energy_coordinate);B=2inv.v_dc*dr;C=-2inv.v_dc*di;H=hypot(B,C)
    results=[]
    for (name,z) in choices
        energy(t)=state.center_V2+B*cos(2t)+C*sin(2t)
        value(t)=inv.m_max*(sqrt(energy(t))-reserve)-sqrt(2)*real(z*cis(t))
        gradient(t)=inv.m_max*(-B*sin(2t)+C*cos(2t))/sqrt(energy(t))+
            sqrt(2)*(real(z)*sin(t)+imag(z)*cos(t))
        curvature=inv.m_max*(2H/state.minimum_V+H^2/state.minimum_V^3)+sqrt(2)*abs(z)
        pad=128eps(Float64)*max(1.,inv.m_max*(state.maximum_V+reserve)+sqrt(2)*abs(z))
        # E(t+pi)=E(t) and voltage(t+pi)=-voltage(t): one signed branch over
        # a complete cycle covers both absolute-voltage signs exactly.
        bound=_smooth_periodic_minimum(value,gradient,curvature,pad;atol)
        push!(results,(name=name,bound=bound))
    end
    index=argmin([r.bound.upper for r in results])
    (lower_V=minimum(r.bound.lower for r in results),upper_V=results[index].bound.upper,
        angle_rad=results[index].bound.angle,branch=results[index].name,
        converged=all(r.bound.converged for r in results),evaluations=sum(r.bound.evaluations for r in results))
end
