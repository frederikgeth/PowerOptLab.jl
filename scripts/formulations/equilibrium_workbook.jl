# Teaching models, not package APIs or calibrated feeder/battery models.
# Run: julia --project=. scripts/formulations/equilibrium_workbook.jl
using JuMP, Ipopt

function wb_parameters(; available=1., resistance=.1, gain=12., source_v=1., reference_v=1.)
    p = (;available=Float64(available),resistance=Float64(resistance),gain=Float64(gain),
        source_v=Float64(source_v),reference_v=Float64(reference_v))
    all(isfinite,values(p)) || throw(ArgumentError("workbook parameters must be finite"))
    p.available >= 0 && p.resistance >= 0 && p.gain >= 0 ||
        throw(ArgumentError("available power, resistance and gain must be nonnegative"))
    p.source_v > 0 && p.reference_v > 0 || throw(ArgumentError("voltages must be positive"))
    return p
end

wb_voltage(p, power) = p.source_v+p.resistance*power
wb_raw(p, voltage) = p.available-p.gain*(voltage-p.reference_v)
wb_command(p, voltage) = clamp(wb_raw(p,voltage),0.,p.available)
wb_equilibrium_oracle(p) = clamp(
    (p.available-p.gain*(p.source_v-p.reference_v))/(1+p.gain*p.resistance),0.,p.available)

function wb_iteration(p; relaxation=1., steps::Integer=20, initial=0.)
    isfinite(relaxation) && 0 < relaxation <= 1 || throw(ArgumentError("relaxation must be in (0,1]"))
    steps >= 0 || throw(ArgumentError("steps must be nonnegative"))
    isfinite(initial) && 0 <= initial <= p.available || throw(ArgumentError("initial power is outside source bounds"))
    power = [Float64(initial)]
    for _ in 1:steps
        command = wb_command(p,wb_voltage(p,last(power)))
        push!(power,(1-relaxation)*last(power)+relaxation*command)
    end
    return (;power,voltage=wb_voltage.(Ref(p),power),
        residual=[x-wb_command(p,wb_voltage(p,x)) for x in power])
end

function _wb_model()
    m = Model(Ipopt.Optimizer)
    set_silent(m)
    set_optimizer_attribute(m,"tol",1e-10)
    set_optimizer_attribute(m,"bound_relax_factor",0.)
    return m
end

function _wb_solved(m)
    termination_status(m) in (JuMP.MOI.OPTIMAL,JuMP.MOI.LOCALLY_SOLVED) && primal_status(m) == JuMP.MOI.FEASIBLE_POINT
end

"""Enumerate the three GUARDED affine clip regions. This transparent teaching
oracle is practical for this tiny example; it is not a scalable OPF algorithm.
`relationship=:capability` replaces p=command with p<=command, retaining guards.
"""
function wb_solve(p; relationship=:controller,objective=:feasibility,target=.2,voltage_max=nothing)
    relationship in (:controller,:capability) || throw(ArgumentError("unknown relationship"))
    objective in (:feasibility,:max_p,:min_p,:tracking) || throw(ArgumentError("unknown objective"))
    isfinite(target) || throw(ArgumentError("tracking target must be finite"))
    voltage_max === nothing || (isfinite(voltage_max) && voltage_max > 0) ||
        throw(ArgumentError("voltage_max must be finite and positive"))
    candidates = NamedTuple[]
    for branch in (:full,:slope,:off)
        m = _wb_model()
        @variable(m,0 <= power <= p.available)
        @variable(m,voltage >= 0,start=p.source_v)
        @constraint(m,voltage == p.source_v+p.resistance*power)
        raw = @expression(m,p.available-p.gain*(voltage-p.reference_v))
        cap = if branch === :full
            @constraint(m,raw >= p.available)
            p.available
        elseif branch === :slope
            @constraint(m,raw >= 0)
            @constraint(m,raw <= p.available)
            raw
        else
            @constraint(m,raw <= 0)
            0.
        end
        if relationship === :controller
            @constraint(m,power == cap)
        else
            @constraint(m,power <= cap)
        end
        voltage_max === nothing || @constraint(m,voltage <= voltage_max)
        if objective === :max_p
            @objective(m,Max,power)
        elseif objective === :min_p
            @objective(m,Min,power)
        elseif objective === :tracking
            @objective(m,Min,(power-target)^2)
        else
            @objective(m,Min,0.)
        end
        optimize!(m)
        solved = _wb_solved(m)
        push!(candidates,(;branch,solved,status=string(termination_status(m)),
            power=solved ? value(power) : NaN,voltage=solved ? value(voltage) : NaN))
    end
    feasible = filter(x -> x.solved,candidates)
    # Closed guards overlap at a clip knot. Deduplicate equal operating points,
    # rather than counting the same boundary equilibrium as multiple solutions.
    solutions = NamedTuple[]
    for c in feasible
        any(s -> abs(s.power-c.power)<1e-7,solutions) || push!(solutions,c)
    end
    score(c) = objective === :max_p ? -c.power : objective === :min_p ? c.power :
        objective === :tracking ? (c.power-target)^2 : 0.
    selected = isempty(solutions) ? nothing : solutions[argmin(score.(solutions))]
    return (;relationship,objective,target,voltage_max,candidates,solutions,selected)
end

"""Explain declared equalities and signed inequality margins in the toy model.
`controller_gap` is a comparison, not a constraint residual in capability mode.
"""
function wb_report(p, solved; tolerance=1e-7)
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError("tolerance must be positive"))
    s = solved.selected
    s === nothing && return nothing
    command = wb_command(p,s.voltage)
    equalities = [(name=:network,residual=s.voltage-wb_voltage(p,s.power),unit=:voltage_pu)]
    if solved.relationship === :controller
        push!(equalities,(name=:controller,residual=s.power-command,unit=:power_pu))
    end
    margins = [(name=:source_lower,margin=s.power,unit=:power_pu),
        (name=:source_upper,margin=p.available-s.power,unit=:power_pu),
        (name=:voltage_lower,margin=s.voltage,unit=:voltage_pu)]
    raw = wb_raw(p,s.voltage)
    if s.branch === :full
        push!(margins,(name=:full_guard,margin=raw-p.available,unit=:power_pu))
    elseif s.branch === :slope
        push!(margins,(name=:slope_lower_guard,margin=raw,unit=:power_pu))
        push!(margins,(name=:slope_upper_guard,margin=p.available-raw,unit=:power_pu))
    else
        push!(margins,(name=:off_guard,margin=-raw,unit=:power_pu))
    end
    if solved.relationship === :capability
        push!(margins,(name=:command_ceiling,margin=command-s.power,unit=:power_pu))
    end
    solved.voltage_max === nothing || push!(margins,
        (name=:voltage_upper,margin=solved.voltage_max-s.voltage,unit=:voltage_pu))
    return (;branch=s.branch,available_p=p.available,raw_p=raw,command_p=command,
        delivered_p=s.power,controller_gap=s.power-command,equalities,margins,
        binding=[x.name for x in margins if abs(x.margin)<=tolerance])
end

"""Physical-time teaching model: ideal signed battery power (kW, positive
discharge), energy (kWh), duration (h), and a toy algebraic voltage each period.
Known prices and a fixed terminal energy are explicit planning assumptions.
"""
function wb_battery_schedule(; voltage_max=1.03)
    isfinite(voltage_max) && voltage_max > 0 || throw(ArgumentError("voltage_max must be positive"))
    prices, duration = [1.,5.,3.],ones(3)
    m = _wb_model()
    @variable(m,-1 <= power[1:3] <= 1)
    @variable(m,0 <= energy[1:4] <= 2)
    @variable(m,voltage[1:3])
    @constraint(m,energy[1] == 1)
    @constraint(m,energy[4] == 1)
    for t in 1:3
        @constraint(m,energy[t+1] == energy[t]-duration[t]*power[t])
        @constraint(m,voltage[t] == 1+.03power[t])
        @constraint(m,voltage[t] <= voltage_max)
    end
    @objective(m,Max,sum(prices[t]*duration[t]*power[t] for t in 1:3))
    optimize!(m)
    solved = _wb_solved(m)
    pp = solved ? value.(power) : fill(NaN,3)
    ee = solved ? value.(energy) : fill(NaN,4)
    vv = solved ? value.(voltage) : fill(NaN,3)
    replay = accumulate(-,duration .* pp;init=1.)
    return (;solved,status=string(termination_status(m)),power=pp,energy=ee,voltage=vv,
        replay_energy=vcat(1.,replay),prices,duration,voltage_max,
        revenue=sum(prices .* duration .* pp))
end

"""A real memory update: upper threshold sets curtailment, lower resets it,
and the open band retains the preceding mode. Equality at either threshold
triggers its corresponding transition. No optimisation encoding is implied.
"""
function wb_hysteresis(voltage; initial=false,lower=.98,upper=1.02)
    isfinite(lower) && isfinite(upper) && 0 < lower < upper || throw(ArgumentError("invalid hysteresis band"))
    all(isfinite,voltage) || throw(ArgumentError("voltage measurements must be finite"))
    mode = Bool(initial)
    modes = Bool[]
    for v in voltage
        if v >= upper
            mode = true
        elseif v <= lower
            mode = false
        end
        push!(modes,mode)
    end
    return modes
end

function equilibrium_workbook()
    p = wb_parameters()
    controller = wb_solve(p)
    cap_max = wb_solve(p;relationship=:capability,objective=:max_p)
    cap_min = wb_solve(p;relationship=:capability,objective=:min_p)
    cap_track = wb_solve(p;relationship=:capability,objective=:tracking,target=.2)
    voltage_controller = wb_solve(p;voltage_max=1.03)
    voltage_dispatch = wb_solve(p;relationship=:capability,objective=:max_p,voltage_max=1.03)
    return (;parameters=p,controller,cap_max,cap_min,cap_track,voltage_controller,voltage_dispatch,
        direct=wb_iteration(p),relaxed=wb_iteration(p;relaxation=.5),
        battery=wb_battery_schedule(),restricted_battery=wb_battery_schedule(voltage_max=1.015))
end

if abspath(PROGRAM_FILE) == @__FILE__
    w = equilibrium_workbook()
    println("Controller equilibrium: ",w.controller.selected)
    println("Direct iteration's last two powers: ",w.direct.power[end-1:end])
    println("Relaxed iteration residual: ",last(w.relaxed.residual))
    for name in (:cap_max,:cap_min,:cap_track,:voltage_dispatch)
        println(name,": ",wb_report(w.parameters,getproperty(w,name)))
    end
    println("Strict controller with voltage ceiling: ",w.voltage_controller.selected)
    println("Battery power: ",w.battery.power,"; replay energy: ",w.battery.replay_energy)
    println("Restricted battery power: ",w.restricted_battery.power)
end
