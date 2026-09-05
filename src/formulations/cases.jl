"""
    resistive_equilibria(curve, source_voltage, resistance; atol=1e-9)

Enumerate canonical equilibria of `V = source_voltage + resistance * curve(V)`
on the curve's bounded domain by solving each affine segment. Return isolated
`points` and degenerate `intervals`; uniqueness is not assumed. `atol` is an
absolute voltage tolerance for endpoint inclusion/deduplication. Continuum
detection uses exact zero coefficients in floating-point arithmetic; nearly
parallel segments remain isolated, potentially ill-conditioned intersections.
This is a scalar electrical reference, not a general AC power-flow solver.
"""
function resistive_equilibria(f::PWLFunction,source_voltage::Real,resistance::Real;
                              atol::Real=1e-9)
    isfinite(source_voltage) && isfinite(resistance) && resistance>=0 ||
        throw(ArgumentError("Need finite source voltage and nonnegative resistance"))
    isfinite(atol) && atol>=0 || throw(ArgumentError("atol must be finite and nonnegative"))
    points,intervals = NamedTuple[],NamedTuple[]
    for j in 1:length(f.breakpoints)-1
        a,b = f.breakpoints[j:j+1]
        slope = (f.values[j+1]-f.values[j])/(b-a)
        intercept = f.values[j]-slope*a
        denominator = 1-resistance*slope
        numerator = source_voltage+resistance*intercept
        if iszero(denominator)
            iszero(numerator) && push!(intervals,(voltage_lower=a,voltage_upper=b,
                current_lower=f.values[j],current_upper=f.values[j+1]))
        else
            v = numerator/denominator
            if a-atol <= v <= b+atol
                v = clamp(v,a,b)
                all(p -> abs(p.voltage-v)>atol,points) &&
                    push!(points,(voltage=v,current=primitive_value(f,v)))
            end
        end
    end
    filter!(p -> !any(i -> i.voltage_lower-atol <= p.voltage <= i.voltage_upper+atol,
                     intervals),points)
    return (points=points,intervals=intervals)
end

"""
    resistive_control_case(curve; source_voltage, resistance, id="resistive")

Create a configurable scalar electrical case with a bounded PWL current law.
Configuration keys override `source_voltage`, `resistance`, `input_scale`,
`output_scale`, `start_input`, and `objective` (`:max_current` or `:zero`). Values
are physical volts, amperes and ohms except the normalized model coordinates.
The method representation is either a formulation or a callable `(curve, config)`
returning one. Metrics include canonical equilibria and physical equation errors.
"""
function resistive_control_case(f::PWLFunction;source_voltage::Real,resistance::Real,
                               id="resistive")
    build = function(representation,config)
        vs,r = get(config,:source_voltage,source_voltage),get(config,:resistance,resistance)
        reference = resistive_equilibria(f,vs,r)
        si,so = _pwl_scale(get(config,:input_scale,1.)),_pwl_scale(get(config,:output_scale,1.))
        start = get(config,:start_input,(first(f.breakpoints)+last(f.breakpoints))/2)
        rep = representation isa AbstractPWLFormulation ? representation : representation(f,config)
        objective = get(config,:objective,:max_current)
        objective in (:max_current,:zero) || throw(ArgumentError("Unknown scalar objective"))
        model = Model()
        @variable(model,x,start=start/si)
        h = formulate_pwl!(model,f,x,rep;input_scale=si,output_scale=so)
        h.output isa VariableRef && set_start_value(h.output,
            primitive_value(f,start;domain_policy=:flat_extension)/so)
        @constraint(model,x == vs/si+r*so/si*h.output)
        @objective(model,Max,objective == :max_current ? h.output : 0.)
        metrics = () -> begin
            v,i = value(x)*si,value(h.output)*so
            (voltage_V=v,current_A=i,power_W=v*i,
             electrical_residual_V=v-vs-r*i,reference=reference)
        end
        return (model=model,observations=[h],metrics=metrics)
    end
    return FormulationCase(id,build;metadata=(model="scalar resistive export",
        source_voltage_V=source_voltage,resistance_ohm=resistance))
end

"""
    controlled_inverter_case(network, device_builder, request;
                             id="controlled_inverter", metrics=nothing)

Adapt the existing single-inverter staged model to formulation experiments.
`device_builder(representation, config)` returns a `ControlledDevice`; `network`
is a dictionary or `config -> dictionary`, and `request` is a request or callback.
Each network is copied. Configuration supports `s_base` and `selection_objective`
(`:loss` or `:zero`). Supply additional model changes through a custom case.

Candidate metrics contain physical POC voltage/current phasors and a same-voltage
exact/smooth controller target gap **before plant capability backoff**. They are
diagnostics even when the solver does not succeed; they do not certify equilibrium
or hardware feasibility. Optional `metrics(ctx, handles, device, request)` supplies
additional named/dictionary data under `custom`. Optimizers belong to the method.
"""
function controlled_inverter_case(network,device_builder,request;id="controlled_inverter",
                                 metrics=nothing)
    build = function(rep,config)
        net = deepcopy(network isa AbstractDict ? network : network(config))
        device = device_builder(rep,config)
        req = request isa InverterControlRequest ? request : request(config)
        validate_device(device,(net,);periods=1)
        objective = get(config,:selection_objective,:loss)
        objective in (:loss,:zero) || throw(ArgumentError("Unknown inverter selection objective"))
        handles = Ref{_ControlledHandles}()
        hook! = ctx -> begin
            handles[] = stamp_device!(ctx,device;request=req)
            @objective(_opf_model(ctx),Min,objective == :loss ?
                handles[].plant.p_loss+handles[].plant.p_cap_loss : 0.)
        end
        ctx = build_opf_model(net;per_unit=true,s_base=Float64(get(config,:s_base,1e6)),
            add_objective=false,model_hook! = hook!,optimizer=nothing,verbose=false)
        enforce_kcl!(ctx)
        candidate_metrics = () -> begin
            h = handles[].control
            voltage = [complex(value(p[1]),value(p[2]))*h.vb for p in h.phase_voltage]
            current = [complex(value(p[1]),value(p[2]))*h.ib for p in h.phase_current]
            measurement = InverterControlMeasurement(voltage)
            ratings = InverterControlRatings(device.device,device.controller.current_target)
            exact = evaluate_exact(device.controller,measurement,req,ratings)
            smooth = evaluate_smooth(device.controller,measurement,req,ratings)
            (phase_voltage_V=voltage,phase_current_A=current,
             requested_active_power_W=value(h.p_request)*h.sb,
             requested_reactive_power_var=value(h.q_request)*h.sb,
             pre_capability_target_gap_A=maximum(abs,exact.phase_current.-smooth.phase_current),
             custom=metrics === nothing ? NamedTuple() : metrics(ctx,handles[],device,req))
        end
        return (model=_opf_model(ctx),metrics=candidate_metrics)
    end
    return FormulationCase(id,build;metadata=(model="controlled inverter",per_unit=true))
end
