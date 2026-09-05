"""
    VoltVarWattIntent(; volt_watt=nothing, volt_var=nothing,
        sensing=:worst_phase, volt_watt_basis=:available,
        conflict_policy=:dominant, conflict_epsilon=0.01,
        worst_phase_watt_guard=true)

Encoding-independent positive-sequence control intent. Curves are immutable
`PWLFunction`s: inputs are physical volts, volt-watt outputs are fractions in
`[0,1]`, and volt-var outputs multiply the request's `q_scale` (positive injection).
Both are non-increasing and extend flat beyond their endpoint knots. No smoothing
family, solver or model-coordinate scale is stored here.

Sensing is `:worst_phase`, `:average_voltage`, or `:positive_sequence`. Conflict
settings apply to worst-phase sensing; the optional watt guard applies to
positive-sequence sensing. `conflict_epsilon` is part of the existing firmware
blend semantics, including exact replay; it is not a solver approximation knob.
"""
struct VoltVarWattIntent
    volt_watt::Union{Nothing,PWLFunction}
    volt_var::Union{Nothing,PWLFunction}
    sensing::Symbol
    volt_watt_basis::Symbol
    conflict_policy::Symbol
    conflict_epsilon::Float64
    worst_phase_watt_guard::Bool
    function VoltVarWattIntent(;volt_watt=nothing,volt_var=nothing,
            sensing::Symbol=:worst_phase,volt_watt_basis::Symbol=:available,
            conflict_policy::Symbol=:dominant,conflict_epsilon::Real=.01,
            worst_phase_watt_guard::Bool=true)
        sensing in (:worst_phase,:average_voltage,:positive_sequence) ||
            throw(ArgumentError("Unknown voltage sensing policy"))
        _validate_volt_watt_basis(volt_watt_basis)
        conflict_policy in (:dominant,:net,:low_voltage,:high_voltage) ||
            throw(ArgumentError("Unknown voltage-conflict policy"))
        curves = map((volt_watt,volt_var)) do curve
            curve === nothing && return nothing
            curve isa PWLFunction || throw(ArgumentError("Intent curves must be PWLFunction objects"))
            curve.input_unit in (:unitless,:V) && curve.output_unit in (:unitless,:pu) ||
                throw(ArgumentError("Intent curves use volts and dimensionless fractions; convert units explicitly"))
            all(<=(0),curve.slopes) || throw(ArgumentError("Control curves must be non-increasing"))
            PWLFunction(curve.breakpoints,curve.values;input_unit=:V,output_unit=:pu)
        end
        curves[1] === nothing || all(y -> 0<=y<=1,curves[1].values) ||
            throw(ArgumentError("Volt-watt fractions must lie in [0,1]"))
        new(curves...,sensing,volt_watt_basis,conflict_policy,
            _positive_width(conflict_epsilon),worst_phase_watt_guard)
    end
end

"""
    VoltVarWattEncoding(; volt_watt=nothing, volt_var=nothing, extrema_epsilon=0.05)

Numerical choices kept separate from `VoltVarWattIntent`. Assign an explicit
formulation to every present curve. `extrema_epsilon` is the physical voltage
width of the established algebraic phase extrema/guard in the full smooth policy.
Graph formulations can be used for individual bounded curve ports; complete IBR
policy lowering currently accepts smooth formulations only.
"""
struct VoltVarWattEncoding
    volt_watt::Union{Nothing,AbstractPWLFormulation}
    volt_var::Union{Nothing,AbstractPWLFormulation}
    extrema_epsilon::Float64
    function VoltVarWattEncoding(;volt_watt=nothing,volt_var=nothing,extrema_epsilon::Real=.05)
        new(volt_watt,volt_var,_positive_width(extrema_epsilon))
    end
end

"""
    lower_positive_policy(intent, encoding)

Lower semantic volt-var/watt intent to the existing smooth IBR policy without
building a JuMP model. Reuse the result when building many snapshots. The result
fits `SequenceController` and its exact/smooth evaluators and staged builder.
Non-smooth graph encodings throw `UnsupportedFormulation` here: support for an
individual PWL graph does not imply support for the whole AC controller.
"""
function lower_positive_policy(intent::VoltVarWattIntent,encoding::VoltVarWattEncoding)
    for role in (:volt_watt,:volt_var)
        curve,rep = getproperty(intent,role),getproperty(encoding,role)
        curve === nothing && continue
        rep isa AbstractPWLSmoothing || throw(UnsupportedFormulation(
            "Full IBR policy lowering requires a smooth $role encoding; graph lowering is available for individual bounded curve ports"))
    end
    laws = map((:volt_watt,:volt_var)) do role
        f = getproperty(intent,role)
        f === nothing ? nothing : PiecewiseLinearLaw(collect(f.breakpoints),collect(f.values);
            formulation=getproperty(encoding,role))
    end
    shared = (volt_watt=laws[1],volt_var=laws[2],volt_watt_basis=intent.volt_watt_basis)
    intent.sensing == :worst_phase && return WorstPhaseVoltVarWatt(;shared...,
        conflict_policy=intent.conflict_policy,conflict_epsilon=intent.conflict_epsilon,
        extrema_epsilon=encoding.extrema_epsilon)
    intent.sensing == :average_voltage && return AverageVoltageVoltVarWatt(;shared...,
        extrema_epsilon=encoding.extrema_epsilon)
    return PositiveSequenceVoltVarWatt(;shared...,
        worst_phase_watt_guard=intent.worst_phase_watt_guard,guard_epsilon=encoding.extrema_epsilon)
end

function _control_port_curve(intent,role,domain)
    role in (:volt_watt,:volt_var) || throw(ArgumentError("Choose :volt_watt or :volt_var"))
    f = getproperty(intent,role)
    f === nothing && throw(ArgumentError("The intent has no $role curve"))
    length(domain)==2 || throw(ArgumentError("Domain needs two voltage endpoints"))
    lo,hi = Float64.(domain)
    all(isfinite,(lo,hi)) && lo<hi || throw(ArgumentError("Need finite increasing voltage endpoints"))
    (lo,hi)==(first(f.breakpoints),last(f.breakpoints)) && return f
    knots = [lo; [k for k in f.breakpoints if lo<k<hi]; hi]
    values = [primitive_value(f,k;domain_policy=:flat_extension) for k in knots]
    return PWLFunction(knots,values;input_unit=:V,output_unit=:pu)
end

"""
    formulate_control_curve!(model_or_context, intent, role, sensed_voltage,
        encoding; domain, input_scale=1, output_scale=1, reuse=true)

Lower one `:volt_watt` or `:volt_var` curve port at an already selected voltage.
The explicit finite `domain=(lower_V,upper_V)` restricts the port's physical
voltage and may extend into the intent's flat tails. `encoding` is a
`VoltVarWattEncoding`; its selected role may be smooth, MPCC, exact MIP or hull.
Return a `PWLFormulationHandle` in fraction output units. Sensing, P/Q conversion,
conflict resolution and plant coupling are not added by this curve-port operation.

With `reuse=true`, identical calls at the same `VariableRef` return the same
handle and rows. Caches belong to one model; they are never shared across models.
Mutable/general input expressions are not cached. `reuse=false` forces a fresh
output/constraint block. Domain, role, curve data, encoding and scales key reuse.
"""
function formulate_control_curve!(target,intent::VoltVarWattIntent,role::Symbol,x,
        encoding::VoltVarWattEncoding;domain,input_scale::Real=1,output_scale::Real=1,
        reuse::Bool=true)
    role in (:volt_watt,:volt_var) || throw(ArgumentError("Choose :volt_watt or :volt_var"))
    f = getproperty(intent,role)
    f === nothing && throw(ArgumentError("The intent has no $role curve"))
    length(domain)==2 || throw(ArgumentError("Domain needs two voltage endpoints"))
    bounds = Tuple(Float64.(domain))
    all(isfinite,bounds) && bounds[1]<bounds[2] ||
        throw(ArgumentError("Need finite increasing voltage endpoints"))
    rep = getproperty(encoding,role)
    rep === nothing && throw(ArgumentError("Supply a $role encoding"))
    # Preserve the input recorded for auditing if a caller subsequently edits
    # an affine/quadratic expression to build a different port.
    x = x isa Union{AffExpr,QuadExpr} ? copy(x) : x
    si,so = _pwl_scale(input_scale),_pwl_scale(output_scale)
    model = target isa JuMP.Model ? target : BMOPFTools.opf_model(target)
    build = function()
        if rep isa AbstractPWLSmoothing
            # Restrict the voltage, not the curve before smoothing: cutting a
            # sloped segment would introduce artificial endpoint hinges.
            @constraint(model,bounds[1]/si <= x <= bounds[2]/si)
            expression = smooth_pwl_expression(target,f,x,rep;input_scale=si,output_scale=so)
            if all(==(first(f.values)),f.values)
                y = first(f.values)/so
            else
                y = @variable(model,base_name="control_curve_output")
                @constraint(model,y == expression)
            end
            return PWLFormulationHandle(model,f,rep,x,y,si,so,si,
                Tuple{VariableRef,VariableRef}[],bounds)
        end
        return formulate_pwl!(target,_control_port_curve(intent,role,bounds),x,rep;
            input_scale=si,output_scale=so)
    end
    reuse && x isa VariableRef || return build()
    owner_model(x) === model || throw(ArgumentError("Voltage variable belongs to another model"))
    cache = get!(model.ext,:PowerOptLabControlCurveHandles) do
        Dict{Any,PWLFormulationHandle}()
    end
    return get!(build,cache,(target,role,f.breakpoints,f.values,bounds,x,rep,si,so))
end

"""
    lowering_statistics(model_or_context)

Inspect model structure and PowerOptLab cache sizes: variables, constraints,
PWL/hinge scalar operators and shared control-curve blocks. Counts describe
construction, not solve cost or a performance benchmark. Upstream registration
caches are not counted in the PowerOptLab operator fields.
"""
function lowering_statistics(target)
    m = target isa JuMP.Model ? target : BMOPFTools.opf_model(target)
    count_cache(key) = length(get(m.ext,key,()))
    return (variables=num_variables(m),
        constraints=num_constraints(m;count_variable_in_set_constraints=true),
        pwl_operators=count_cache(:PowerOptLabPWLOperators),
        hinge_operators=count_cache(:PowerOptLabHingeOperators),
        shared_curve_blocks=count_cache(:PowerOptLabControlCurveHandles))
end
