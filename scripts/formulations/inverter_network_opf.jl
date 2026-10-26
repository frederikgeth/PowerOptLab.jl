# Full network composition plus an independent four-wire radial replay.
# Run: julia --project=. scripts/formulations/inverter_network_opf.jl
using PowerOptLab, JuMP, Ipopt, LinearAlgebra
import BMOPFTools, JSON3

const _NOPF_PHASES = ("a","b","c","n")
const _NOPF_A = cis(2pi/3)

function inverter_network_case()
    buses = ["source"; ["b$i" for i in 1:6]]
    parents = [1,2,2,3,4,4] # each child 2:7, in topological order
    length_factor = [.7,.9,.8,1.1,.9,1.2]
    r = fill(.018,4,4); x = fill(.035,4,4)
    for k in 1:4
        r[k,k] = k == 4 ? .32 : .23
        x[k,k] = k == 4 ? .18 : .15
    end
    impedances = [f*complex.(r,x) for f in length_factor]
    loads = ComplexF64[
        0 0 0;
        400+100im 300+60im 500+120im;
        700+160im 400+80im 300+75im;
        300+60im 850+180im 400+80im;
        500+100im 250+60im 750+170im;
        900+200im 450+100im 250+50im;
        200+40im 700+160im 450+100im]
    source = ComplexF64[241,242cis(-2pi/3),240cis(2pi/3),0]
    data = Dict{String,Any}("bus"=>Dict{String,Any}(),"linecode"=>Dict{String,Any}(),
        "line"=>Dict{String,Any}(),"load"=>Dict{String,Any}(),
        "voltage_source"=>Dict("grid"=>Dict("bus"=>"source","terminal_map"=>["a","b","c"],
            "v_magnitude"=>abs.(source[1:3]),"v_angle"=>angle.(source[1:3]))))
    for (i,bus) in enumerate(buses)
        # Native v_min/v_max are terminal-to-ground; the example adds explicit
        # phase-to-neutral service limits in its hook below.
        data["bus"][bus] = Dict("terminal_names"=>collect(_NOPF_PHASES),
            "v_min"=>fill(200.,3),"v_max"=>fill(280.,3))
        i == 1 && (data["bus"][bus]["perfectly_grounded_terminals"]=["n"])
        for k in 1:3
            iszero(loads[i,k]) && continue
            data["load"]["load_$(i)_$(k)"] = Dict("bus"=>bus,
                "terminal_map"=>[_NOPF_PHASES[k],"n"],"configuration"=>"SINGLE_PHASE",
                "p_nom"=>[real(loads[i,k])],"q_nom"=>[imag(loads[i,k])])
        end
    end
    for e in 1:6
        code = Dict{String,Any}()
        for i in 1:4,j in 1:4
            code["R_series_$(i)_$(j)"]=real(impedances[e][i,j])
            code["X_series_$(i)_$(j)"]=imag(impedances[e][i,j])
        end
        data["linecode"]["z$e"] = code
        data["line"]["line$e"] = Dict("bus_from"=>buses[parents[e]],"bus_to"=>buses[e+1],
            "terminal_map_from"=>collect(_NOPF_PHASES),"terminal_map_to"=>collect(_NOPF_PHASES),
            "linecode"=>"z$e","length"=>1.)
    end
    net = BMOPFTools.parse_bmopf(JSON3.write(data);from_string=true)
    plants = [AdvancedInverter(id="pv$i",bus=bus,topology=:THREE_LEG,
        phase_terminals=["a","b","c"],neutral="n",s_max=rating,i_max=rating/500.,
        v_dc=800.,c_dc=.002,r_filter=.05,x_filter=.15,p_loss_fixed=10.0+5i,c_loss=.03)
        for (i,bus,rating) in ((1,"b2",6000.),(2,"b4",6000.),(3,"b6",5000.))]
    vv = PiecewiseLinearLaw([230.,250.],[.3,-.3];formulation=LocalC2Formulation(.01))
    policies = [AverageVoltageVoltVarWatt(volt_var=vv),
        PositiveSequenceVoltVarWatt(volt_var=vv),
        WorstPhaseVoltVarWatt(volt_var=vv,conflict_policy=:net)]
    availability = [3500.,4500.,3000.]
    return (;net,buses,parents,impedances,loads,source,plants,policies,availability,
        inverter_bus=[findfirst(==(p.bus),buses) for p in plants])
end

_nopf_v(ctx,bus,terminal,component) = BMOPFTools.opf_object(ctx,
    BMOPFTools.opf_bus_voltage_key(bus,terminal;component))

"""Compose all inverters in one network model. `:forward` prescribes controls;
`:matched` maximises export under matching voltage laws; `:dispatch` frees Q
and applies coordinated weights/network limits. No native IBR ports are present.
"""
function solve_inverter_network(case; mode=:matched,s_base=1e6,voltage_max=244.,
        feeder_current_max=12.,max_iter=1000,fixed_dispatch=nothing,voltage_bounds=:native)
    mode in (:forward,:matched,:dispatch) || throw(ArgumentError("unknown network mode"))
    voltage_bounds in (:native,:hook) || throw(ArgumentError("unknown voltage-bound formulation"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError("s_base must be positive"))
    isfinite(voltage_max) && voltage_max > 0 || throw(ArgumentError("voltage_max must be positive"))
    isfinite(feeder_current_max) && feeder_current_max > 0 || throw(ArgumentError("current rating must be positive"))
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    if fixed_dispatch !== nothing
        mode === :dispatch || throw(ArgumentError("fixed dispatch requires dispatch mode"))
        length(fixed_dispatch)==length(case.plants) && all(isfinite,fixed_dispatch) ||
            throw(ArgumentError("fixed dispatch must contain one finite complex P/Q value per inverter"))
    end
    requests = [InverterControlRequest(p_available=mode === :dispatch ? inv.s_max : case.availability[i],
        p_rated=inv.s_max,q_scale=inv.s_max) for (i,inv) in enumerate(case.plants)]
    devices = [mode === :forward ? ControlledDevice(inv,SequenceController(case.policies[i];
            limiter=CommonScaleLimiter(pq_priority=:var),power_voltage_floor=1e-3)) :
        ConstrainedInverter(inv,mode === :matched ? VoltVarWattConstraints(case.policies[i]) :
            InverterOperatingConstraints(p_min=0.,q_min=-.4inv.s_max,q_max=.4inv.s_max);
            dc_source=DCSourceBounds(p_min=0.)) for (i,inv) in enumerate(case.plants)]
    handles = Any[]
    weights = [1.,1.25,.9]
    model_net = deepcopy(case.net)
    if mode === :dispatch && voltage_bounds === :native
        for bus in case.buses[2:end]
            model_net["bus"][bus]["vpn_min"] = fill(220.,3)
            model_net["bus"][bus]["vpn_max"] = fill(voltage_max,3)
        end
    end
    hook! = ctx -> begin
        m = BMOPFTools.opf_model(ctx)
        for i in eachindex(devices)
            push!(handles,stamp_device!(ctx,devices[i];request=requests[i]))
        end
        if mode === :dispatch
            for bus in (voltage_bounds === :hook ? case.buses[2:end] : String[])
                scale = BMOPFTools.opf_coordinate_bases(ctx,bus).voltage
                nr,ni = _nopf_v(ctx,bus,"n",:real),_nopf_v(ctx,bus,"n",:imag)
                for phase in _NOPF_PHASES[1:3]
                    vr,vi = _nopf_v(ctx,bus,phase,:real)-nr,_nopf_v(ctx,bus,phase,:imag)-ni
                    @constraint(m,(vr^2+vi^2)/(voltage_max/scale)^2 <= 1)
                    @constraint(m,(vr^2+vi^2)/(220.0/scale)^2 >= 1)
                end
            end
            base = BMOPFTools.opf_coordinate_bases(ctx,"source").current
            for k in 1:4
                ir = BMOPFTools.opf_object(ctx,BMOPFTools.opf_line_current_key("line1",k))
                ii = BMOPFTools.opf_object(ctx,BMOPFTools.opf_line_current_key("line1",k;component=:imag))
                @constraint(m,(ir^2+ii^2)/(feeder_current_max/base)^2 <= 1)
            end
        end
        h = [inverter_handles(devices[i],handles[i]) for i in eachindex(devices)]
        if fixed_dispatch !== nothing
            for i in eachindex(h)
                rating = case.plants[i].s_max
                @constraint(m,h[i].p_poc/(rating/h[i].sb) == real(fixed_dispatch[i])/rating)
                @constraint(m,h[i].q_poc/(rating/h[i].sb) == imag(fixed_dispatch[i])/rating)
            end
        end
        if mode === :forward
            @objective(m,Min,sum(x.p_loss+x.p_cap_loss for x in h))
        elseif mode === :matched
            @objective(m,Max,sum(x.p_poc*x.sb/10000 for x in h))
        else
            @objective(m,Max,sum(weights[i]*h[i].p_poc*h[i].sb/10000-
                .05*(h[i].q_poc*h[i].sb/case.plants[i].s_max)^2 for i in eachindex(h)))
        end
    end
    ctx = BMOPFTools.build_opf_model(model_net;per_unit=true,s_base=Float64(s_base),
        add_objective=false,model_hook! = hook!,optimizer=Ipopt.Optimizer)
    model = BMOPFTools.opf_model(ctx)
    set_optimizer_attribute(model,"tol",1e-9)
    set_optimizer_attribute(model,"max_iter",max_iter)
    BMOPFTools.optimize_opf!(ctx)
    status = solve_status(model)
    physical(x,base) = status.publishable ? value(x)*base : NaN
    voltages = [ComplexF64[complex(physical(_nopf_v(ctx,b,t,:real),
        BMOPFTools.opf_coordinate_bases(ctx,b).voltage),physical(_nopf_v(ctx,b,t,:imag),
        BMOPFTools.opf_coordinate_bases(ctx,b).voltage)) for t in _NOPF_PHASES] for b in case.buses]
    branch_current = [ComplexF64[complex(
        physical(BMOPFTools.opf_object(ctx,BMOPFTools.opf_line_current_key("line$e",k)),
            BMOPFTools.opf_coordinate_bases(ctx,case.buses[case.parents[e]]).current),
        physical(BMOPFTools.opf_object(ctx,BMOPFTools.opf_line_current_key("line$e",k;component=:imag)),
            BMOPFTools.opf_coordinate_bases(ctx,case.buses[case.parents[e]]).current)) for k in 1:4] for e in 1:6]
    results = [extract_device(devices[i],handles[i],status) for i in eachindex(devices)]
    return (;mode,status,devices=results,voltages,branch_current,requests,weights,
        voltage_max,feeder_current_max,s_base,fixed_dispatch,voltage_bounds,variables=num_variables(model),
        objective=status.publishable ? objective_value(model) : NaN)
end

"""Independent SI demand currents from constant-PQ loads and balanced positive-
sequence inverter POC dispatch. This oracle has no JuMP or controller evaluator.
"""
function _nopf_demand(case,voltages,dispatch)
    demand = [zeros(ComplexF64,4) for _ in case.buses]
    for b in eachindex(case.buses)
        u = voltages[b][1:3] .- voltages[b][4]
        for k in 1:3
            current = conj(case.loads[b,k]/u[k])
            demand[b][k] += current
            demand[b][4] -= current
        end
    end
    for (i,b) in enumerate(case.inverter_bus)
        u = voltages[b][1:3] .- voltages[b][4]
        u1 = (u[1]+_NOPF_A*u[2]+_NOPF_A^2*u[3])/3
        i1 = conj(dispatch[i]/(3u1))
        demand[b][1:3] .-= [i1,_NOPF_A^2*i1,_NOPF_A*i1]
    end
    return demand
end

function _nopf_backward(case,demand)
    accumulated = deepcopy(demand)
    currents = [zeros(ComplexF64,4) for _ in 1:6]
    for e in 6:-1:1
        currents[e] = copy(accumulated[e+1])
        accumulated[case.parents[e]] .+= currents[e]
    end
    return currents
end

function replay_inverter_network(case,dispatch;tol=1e-10,max_iter=200)
    length(dispatch)==length(case.plants) || throw(ArgumentError("one dispatch per inverter required"))
    all(isfinite,dispatch) || throw(ArgumentError("dispatch must be finite"))
    isfinite(tol) && tol > 0 && max_iter > 0 || throw(ArgumentError("invalid replay tolerance/iterations"))
    voltages = [copy(case.source) for _ in case.buses]
    error = Inf
    for iteration in 1:max_iter
        currents = _nopf_backward(case,_nopf_demand(case,voltages,dispatch))
        next = [copy(case.source) for _ in case.buses]
        for e in 1:6
            next[e+1] = next[case.parents[e]]-case.impedances[e]*currents[e]
        end
        error = maximum(maximum(abs.(next[b]-voltages[b])) for b in eachindex(voltages))
        voltages = next
        if error < tol
            currents = _nopf_backward(case,_nopf_demand(case,voltages,dispatch))
            return (;converged=true,iteration,error,voltages,branch_current=currents)
        end
    end
    return (;converged=false,iteration=max_iter,error,voltages,
        branch_current=_nopf_backward(case,_nopf_demand(case,voltages,dispatch)))
end

function audit_inverter_network(case,solution)
    solution.status.publishable || return nothing
    dispatch = [complex(r.plant.p_poc,r.plant.q_poc) for r in solution.devices]
    replay = replay_inverter_network(case,dispatch)
    demand = _nopf_demand(case,solution.voltages,dispatch)
    kcl = [copy(x) for x in demand]
    kvl = Float64[]
    for e in 1:6
        kcl[e+1] .-= solution.branch_current[e]
        kcl[case.parents[e]] .+= solution.branch_current[e]
        push!(kvl,maximum(abs.(solution.voltages[e+1]-solution.voltages[case.parents[e]]+
            case.impedances[e]*solution.branch_current[e])))
    end
    source_power = sum(case.source .* conj.(solution.branch_current[1]))
    line_power = sum(dot(solution.branch_current[e],case.impedances[e]*solution.branch_current[e]) for e in 1:6)
    phase_v = [abs.(v[1:3] .- v[4]) for v in solution.voltages[2:end]]
    return (;replay,dispatch,
        voltage_error_V=maximum(maximum(abs.(solution.voltages[b]-replay.voltages[b])) for b in eachindex(case.buses)),
        current_error_A=maximum(maximum(abs.(solution.branch_current[e]-replay.branch_current[e])) for e in 1:6),
        kcl_error_A=maximum(maximum(abs.(x)) for x in kcl[2:end]),kvl_error_V=maximum(kvl),
        source_power_VA=source_power,line_power_VA=line_power,
        power_balance_error_VA=abs(source_power+sum(dispatch)-sum(case.loads)-line_power),
        phase_neutral_max_V=maximum(maximum,phase_v),phase_neutral_min_V=minimum(minimum,phase_v),
        neutral_max_V=maximum(abs(v[4]) for v in solution.voltages),
        feeder_max_A=maximum(abs,solution.branch_current[1]))
end

function inverter_full_opf_validation(;s_base=1e6)
    case = inverter_network_case()
    forward = solve_inverter_network(case;mode=:forward,s_base)
    matched = solve_inverter_network(case;mode=:matched,s_base)
    high_case = merge(case,(availability=[p.s_max for p in case.plants],))
    unmanaged = solve_inverter_network(high_case;mode=:forward,s_base)
    dispatch = solve_inverter_network(case;mode=:dispatch,s_base)
    benchmark = solve_inverter_network(case;mode=:dispatch,s_base,
        fixed_dispatch=fill(2000.0+0im,3))
    agreement = (
        p_error_W=maximum(abs(forward.devices[i].plant.p_poc-matched.devices[i].plant.p_poc) for i in 1:3),
        q_error_var=maximum(abs(forward.devices[i].plant.q_poc-matched.devices[i].plant.q_poc) for i in 1:3),
        voltage_error_V=maximum(maximum(abs.(forward.voltages[b]-matched.voltages[b])) for b in 1:7),
        converter_current_error_A=maximum(maximum(abs.(forward.devices[i].converter_terminal.phase_current .-
            matched.devices[i].converter_terminal.phase_current)) for i in 1:3),
        branch_current_error_A=maximum(maximum(abs.(forward.branch_current[e]-matched.branch_current[e])) for e in 1:6))
    return (;case,forward,matched,unmanaged,dispatch,benchmark,agreement,
        audits=(forward=audit_inverter_network(case,forward),matched=audit_inverter_network(case,matched),
            unmanaged=audit_inverter_network(case,unmanaged),
            dispatch=audit_inverter_network(case,dispatch),benchmark=audit_inverter_network(case,benchmark)))
end

function print_inverter_network_validation(r)
    println("Matched-control agreement: ",r.agreement)
    for name in (:forward,:matched,:unmanaged,:dispatch,:benchmark)
        solution,audit = getproperty(r,name),getproperty(r.audits,name)
        println(name," status: ",solution.status)
        println(name," P/Q: ",[(x.plant.p_poc,x.plant.q_poc) for x in solution.devices])
        println(name," replay: ",audit === nothing ? nothing :
            (;voltage_error_V=audit.voltage_error_V,current_error_A=audit.current_error_A,
              kcl_error_A=audit.kcl_error_A,kvl_error_V=audit.kvl_error_V,
              power_balance_error_VA=audit.power_balance_error_VA,
              phase_neutral_max_V=audit.phase_neutral_max_V,feeder_max_A=audit.feeder_max_A,
              neutral_max_V=audit.neutral_max_V))
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    print_inverter_network_validation(inverter_full_opf_validation())
end
