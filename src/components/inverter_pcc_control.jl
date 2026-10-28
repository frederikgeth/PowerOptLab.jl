# Idealised terminal-power controller. The policy is defined in physical PCC
# coordinates first; the JuMP graph below approximates that same policy.
# All complex algebra uses rectangular pairs so the geometry is shared by
# numerical replay and the symbolic stamp, without complex JuMP expressions.
_pcadd(a, b) = (a[1]+b[1], a[2]+b[2])
_pcsub(a, b) = (a[1]-b[1], a[2]-b[2])
_pcscale(t, a) = (t*a[1], t*a[2])
_pcabs2(a) = a[1]^2+a[2]^2
_pcdot(a, b) = a[1]*b[1]+a[2]*b[2]
_pcpair(z::Complex) = (real(z), imag(z))
_pcsum(a) = reduce(_pcadd, a; init=(0.0, 0.0))
_pcpower(v, i) = _pcsum(_cmul(v[k][1], v[k][2], i[k][1], -i[k][2]) for k in 1:3)
_pcphadd(a, b) = ntuple(k -> _pcadd(a[k], b[k]), 3)
_pcphscale(t, a) = ntuple(k -> _pcscale(t, a[k]), 3)
_pcsequences(a) = _sequence_pair([x[1] for x in a], [x[2] for x in a])

# Positive-sequence grid current delivering complex power s. The terminal
# controller's domain excludes zero positive-sequence voltage; unlike the
# standalone service law this inverse must not attenuate shunt compensation.
function _pcc_i1(s, u1)
    all(_is_identically_zero, s) && return (0., 0.)
    _pcscale(1/(3*_pcabs2(u1)), _cmul(s[1], -s[2], u1[1], u1[2]))
end

function _pcc_geometry(voltage, shunt, b, p_import, p_service, q, i2)
    _, u1, _ = _pcsequences(voltage)
    sh0, _, _ = _pcsequences(shunt)
    # THREE_LEG requires I0_converter=0. Its grid arm may carry the opposite
    # of the LCL capacitor zero sequence; this is not a zero-grid-I0 policy.
    grid0 = ntuple(_ -> _pcscale(-1, sh0), 3)
    shunt_power = (0.0, b*sum(_pcabs2, voltage))
    compensation = _pcsub((p_import, 0.0),
        _pcadd(_pcpower(voltage, grid0), shunt_power))
    base = _pcphadd(grid0, _phase_pairs(_pcc_i1(compensation, u1), (0., 0.)))
    pdir = _phase_pairs(_pcc_i1((p_service, 0.), u1), (0., 0.))
    qdir = _phase_pairs(_pcc_i1((0., q), u1), (0., 0.))
    negative = _phase_pairs((0., 0.), i2)
    # Negative-sequence support changes phase currents, not requested net P/Q.
    ncomp = _pcc_i1(_pcscale(-1, _pcpower(voltage, negative)), u1)
    ndir = _pcphadd(negative, _phase_pairs(ncomp, (0., 0.)))
    return base, pdir, qdir, ndir
end

# Largest t in [0,1] satisfying |a+t*b| <= limit, from a feasible anchor a.
# This keeps the complex cross term. The triangle inequality is not used.
function _pcc_segment_exact(a, b, limit, epsilon)
    if iszero(limit)
        iszero(_pcabs2(a)) || throw(ArgumentError(
            "PCC compensation exceeds zero fundamental-current headroom"))
        return iszero(_pcabs2(b)) ? 1.0 : 0.0
    end
    h = 1-_pcabs2(a)/limit^2
    h >= -1e-8 || throw(ArgumentError(
        "PCC controller compensation/priority anchor exceeds a declared limit"))
    h = max(h, 0.0)
    d = _pcabs2(b)/limit^2
    iszero(d) && return 1.0
    c = _pcdot(a, b)/limit^2
    root = sqrt(max(c^2+d*h, 0.0))
    t = c >= 0 ? (iszero(c+root) ? 0.0 : h/(c+root)) : (root-c)/d
    return clamp(t, 0.0, 1.0)
end

function _pcc_segment_smooth!(m, a, b, limit, epsilon)
    if iszero(limit)
        for value in a
            _is_identically_zero(value) || JuMP.@constraint(m, value == 0)
        end
        return all(_is_identically_zero, b) ? 1.0 : 0.0
    end
    all(_is_identically_zero, b) && return 1.0
    h = _defined_control_expression!(m, 1-_pcabs2(a)/limit^2)
    # An inverter unable to sustain its compensated zero-service anchor is
    # outside this controller's operating domain, not silently re-dispatched.
    JuMP.@constraint(m, h >= 0)
    d = _defined_control_expression!(m, _pcabs2(b)/limit^2; start=.25)
    c = _defined_control_expression!(m, _pcdot(a, b)/limit^2; start=0.)
    # Use a slightly stricter quadratic, then smoothly clip its positive
    # root to one. Positive coefficients remain separated from zero even
    # when a requested direction vanishes at the solved point. The root's
    # nonnegative-numerator constraints select the outer boundary and retain
    # inward motion from a saturated anchor.
    delta = epsilon^2
    headroom = h
    coefficient = d+delta
    projection = c+delta
    root = JuMP.@variable(m, lower_bound=0., start=.5)
    # A plain squared-root equality loses physical accuracy at a PQ-axis
    # endpoint: an O(tol) squared residual permits O(sqrt(tol)) output.
    # This positive denominator changes only the residual scaling, keeping
    # both the root and the requested smoothing law unchanged.
    JuMP.@constraint(m,
        (root^2-projection^2-coefficient*headroom)/(root+epsilon) == 0)
    JuMP.@constraint(m, root >= projection)
    JuMP.@constraint(m, root >= -projection)
    numerator = _defined_control_expression!(m, root-projection; start=.5)
    denominator = _smooth_max_expression!(m, coefficient, numerator, delta)
    return _defined_control_expression!(m, numerator/denominator; start=.9)
end

function _pcc_allocate(limiter, base, pdir, qdir, ndir, vconv, shunt,
                       power_request, inverter, ib, sb, segment, combine, peps, ieps)
    p_import, p_service, q = power_request
    function scale_at(anchor, direction, power_anchor, power_direction)
        t = all(_is_identically_zero, power_direction) ? 1.0 :
            segment(power_anchor, power_direction, inverter.s_max/sb, peps)
        for (rating, reserves, offset) in (
                (inverter.i_max, inverter.pwm_ac_converter_reserve, shunt),
                (inverter.i_grid_max, inverter.pwm_ac_grid_reserve, nothing))
            rating === nothing && continue
            for k in 1:3
                a = offset === nothing ? anchor[k] : _pcadd(anchor[k], offset[k])
                limit = sqrt(max(rating^2-reserves[k]^2, 0.0))/ib
                t = combine(t, segment(a, direction[k], limit, ieps))
            end
        end
        if inverter.dv2_max !== nothing && !iszero(_dc_ripple_factor(inverter))
            conv = _pcphadd(anchor, shunt)
            a = _pcsum(_cmul(vconv[k]..., conv[k]...) for k in 1:3)
            b = _pcsum(_cmul(vconv[k]..., direction[k]...) for k in 1:3)
            limit = inverter.dv2_max*inverter.v_dc/abs(_dc_ripple_factor(inverter))/sb
            t = combine(t, segment(a, b, limit, peps))
        end
        return t
    end
    # Allocate services in the declared order, retaining already allocated
    # services unchanged. Proportional common-scale mode scales all together.
    priority = limiter.pq_priority
    if limiter isa CommonScaleLimiter && priority === :proportional
        direction = _pcphadd(_pcphadd(pdir, qdir), ndir)
        t = scale_at(base, direction, (p_import,0.), (p_service,q))
        return _pcphadd(base, _pcphscale(t, direction)), t, t, t
    end
    anchor = base
    power_anchor = (p_import,0.)
    tp = tq = tn = 1.0
    balance_first = limiter isa SequencePriorityLimiter && limiter.priority === :balance
    if balance_first
        tn = scale_at(anchor, ndir, power_anchor, (0.,0.))
        anchor = _pcphadd(anchor, _pcphscale(tn, ndir))
    end
    if priority === :proportional
        direction = _pcphadd(pdir, qdir)
        tp = tq = scale_at(anchor, direction, power_anchor, (p_service,q))
        anchor = _pcphadd(anchor, _pcphscale(tp, direction))
        power_anchor = _pcadd(power_anchor, _pcscale(tp,(p_service,q)))
    else
        first, second = priority === :watt ? (pdir, qdir) : (qdir, pdir)
        first_power, second_power = priority === :watt ?
            ((p_service,0.), (0.,q)) : ((0.,q), (p_service,0.))
        t1 = scale_at(anchor, first, power_anchor, first_power)
        anchor = _pcphadd(anchor, _pcphscale(t1, first))
        power_anchor = _pcadd(power_anchor, _pcscale(t1,first_power))
        t2 = scale_at(anchor, second, power_anchor, second_power)
        power_anchor = _pcadd(power_anchor, _pcscale(t2,second_power))
        anchor = _pcphadd(anchor, _pcphscale(t2, second))
        tp, tq = priority === :watt ? (t1, t2) : (t2, t1)
    end
    if !balance_first
        tn = scale_at(anchor, ndir, power_anchor, (0.,0.))
        anchor = _pcphadd(anchor, _pcphscale(tn, ndir))
    end
    return anchor, tp, tq, tn
end

function _evaluate_pcc_control_exact(controller, request, inverter,
        sensed, voltage, vconv, shunt, loss)
    u0, u1, u2 = _phasor_sequences(sensed)
    abs(u1) >= controller.power_voltage_floor || throw(ArgumentError(
        "PCC terminal controller requires |U1| >= power_voltage_floor"))
    fraction, q, vmin, vmax = _positive_exact(controller.positive, u1, abs.(sensed), request)
    p = _active_power_exact(controller.positive, request, fraction, loss)
    p_import = min(request.p_available-loss, 0.)
    abs(p_import) <= inverter.s_max || throw(ArgumentError(
        "PCC loss import exceeds the terminal apparent-power rating"))
    p_service = p-p_import
    i1request = conj(complex(p_service, q))/(3conj(u1))
    i2v, i2r, i2request, eta = _unbalance_exact(controller.unbalance, u1, u2, i1request)
    v, vc, ish = _pcpair.(voltage), _pcpair.(vconv), _pcpair.(shunt)
    base, pd, qd, nd = _pcc_geometry(v, ish, inverter.b_filter_shunt,
        p_import, p_service, q, _pcpair(i2request))
    grid, tp, tq, tn = _pcc_allocate(controller.limiter, base, pd, qd, nd,
        vc, ish, (p_import,p_service,q), inverter, 1., 1., _pcc_segment_exact, min, 0., 0.)
    target = controller.current_target isa ConverterCurrentTarget ? _pcphadd(grid, ish) : grid
    phase = ntuple(k -> ComplexF64(target[k]...), 3)
    _, i1, i2 = _phasor_sequences(phase)
    raw_norm = hypot(p_service, q)
    service_scale = iszero(raw_norm) ? 1. : hypot(tp*p_service, tq*q)/raw_norm
    return InverterControlResult(sensed, (u0,u1,u2), vmin,vmax,eta,
        p, q, service_scale, i1request,i2v,i2r,i2request, service_scale,
        i1,i2,phase,3u1*conj(i1),sum(sensed .* conj.(phase)),3(u1*i2+u2*i1))
end

function _stamp_pcc_control!(ctx, controller, request, inverter, handles)
    m = _opf_model(ctx)
    sb, vb, ib = handles.sb, handles.vb, handles.ib
    ratings = InverterControlRatings(inverter, controller.current_target)
    sensed = ntuple(k -> _dv(ctx, inverter.bus, inverter.phase_terminals[k],
        sensing_terminal(controller, inverter)), 3)
    voltage = ntuple(k -> (handles.vrpoc[k], handles.vipoc[k]), 3)
    vconv = ntuple(k -> (handles.vrint[k], handles.viint[k]), 3)
    shunt = ntuple(k -> (handles.cri[k]-handles.gri[k], handles.cii[k]-handles.gii[k]), 3)
    u0,u1,u2 = _pcsequences(sensed)
    JuMP.@constraint(m, _pcabs2(u1) >= (controller.power_voltage_floor/vb)^2)
    vstart = _control_voltage_start(controller.positive, vb)
    magnitudes = [_implicit_sqrt!(m, _pcabs2(v); start=vstart, scale=vstart) for v in sensed]
    fraction,q,vmin,vmax = _positive_smooth(ctx, controller.positive, magnitudes,u1,request,sb,vb)
    smax = inverter.s_max/sb
    seps = _power_epsilon(controller.limiter, ratings)/sb
    loss = _defined_control_expression!(m, handles.p_dc-handles.p_poc; scale=smax, start=0.)
    available = request.p_available/sb-loss
    p = _volt_watt_basis(controller.positive) === :available ?
        available-(1-fraction)*_smooth_positive_expression!(m,available,seps) :
        _smooth_min_expression!(m,request.p_rated/sb*fraction,available,seps)
    p_import = _defined_control_expression!(m,
        _smooth_negative_expression!(m,available,_IMPORT_WIDTH_FRACTION*smax); scale=smax,start=0.)
    p_service = _defined_control_expression!(m,p-p_import; scale=smax)
    i1request = _pcc_i1((p_service,q),u1)
    i2v,i2r,i2request,eta = _unbalance_smooth(ctx,controller.unbalance,u1,u2,i1request,vb,ib,vstart)
    base,pd,qd,nd = _pcc_geometry(voltage,shunt,inverter.b_filter_shunt*vb^2/sb,
        p_import,p_service,q,i2request)
    peps = seps/smax
    ieps = _current_epsilon(controller.limiter,ratings)/ratings.i_scale
    grid,tp,tq,tn = _pcc_allocate(controller.limiter,base,pd,qd,nd,vconv,shunt,(p_import,p_service,q),
        inverter,ib,sb,(a,b,l,e)->_pcc_segment_smooth!(m,a,b,l,e),
        (a,b)->_combine_safe_scale_implicit!(m,a,b,min(peps,ieps)),peps,ieps)
    target = controller.current_target isa ConverterCurrentTarget ? _pcphadd(grid,shunt) : grid
    _,i1,i2 = _pcsequences(target)
    # Stamp the same physical equalities for both actuator representations.
    # LCL KCL makes grid=command equivalent to converter=command+shunt.
    _,command1,command2 = _pcsequences(grid)
    _,actual1,actual2 = _sequence_pair(collect(handles.gri),collect(handles.gii))
    # Converter I0=0 and LCL KCL already supply the zero-sequence equation.
    command_scale = ratings.i_scale/ib
    for (actual,command) in ((actual1,command1),(actual2,command2)), j in 1:2
        JuMP.@constraint(m,(actual[j]-command[j])/command_scale == 0)
    end
    service_scale = _defined_control_expression!(m,
        sqrt((tp*p_service)^2+(tq*q)^2+seps^2)/sqrt(p_service^2+q^2+seps^2))
    sp = _pcscale(3,_cmul(u1[1],u1[2],i1[1],-i1[2]))
    ripple = _pcscale(3,_pcadd(_cmul(u1...,i2...),_cmul(u2...,i1...)))
    return _ControlHandles(sensed,(u0,u1,u2),vmin,vmax,eta,p,q,service_scale,
        i1request,i2v,i2r,i2request,service_scale,i1,i2,target,sp,
        _pcpower(sensed,target),ripple,sb,vb,ib,controller,request,ratings)
end
