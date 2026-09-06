# A thin composition of the native staged engine and generator component stamps.

function _gg_replacements(net,devices,replacements)
    ids=[d.id for d in devices]
    allunique(ids) || throw(ArgumentError("generator ids must be unique"))
    all(k -> k in ids,keys(replacements)) || throw(ArgumentError("replacement has unknown device id"))
    allunique(collect(values(replacements))) || throw(ArgumentError("a native component can be replaced only once"))
    for (id,(family,nativeid)) in replacements
        family in (:generator,:voltage_source) || throw(ArgumentError("replacement family must be :generator or :voltage_source"))
        native=get(get(net,string(family),Dict()),nativeid,nothing)
        native===nothing && throw(ArgumentError("native replacement target not found"))
        d=only(filter(d -> d.id==id,devices))
        native["bus"]==d.bus || throw(ArgumentError("replacement bus must match native bus"))
    end
end

function _gg_retire_native!(ctx,family,id)
    # Select the public ledger's already-declared current keys rather than
    # reproducing the engine's terminal-arity heuristics.
    keymaker=family==:generator ? BMOPFTools.opf_generator_current_key : BMOPFTools.opf_voltage_source_current_key
    exemplar=keymaker(id,1)
    for key in BMOPFTools.opf_object_keys(ctx)
        key.kind==exemplar.kind && key.family in (exemplar.family,keymaker(id,1;component=:imag).family) || continue
        key.index isa Tuple && first(key.index)==id || continue
        v=BMOPFTools.opf_object(ctx,key)
        JuMP.fix(v,0.0;force=true); JuMP.set_start_value(v,0.0)
    end
end

"""
    build_generator_model(net, devices; replacements=Dict(), objective=:cost,
        generator_hook! = (ctx, handles)->nothing, kwargs...)

Compose `GeneralizedGenerator`/`SourceGenerator` stamps with BMOPFTools' staged
engine, returning `(context, handles, devices, replacements)`. Does not enforce
network KCL or solve. `replacements` maps custom id to `(family, native_id)`, with
family `:generator` or `:voltage_source`; native physics is replaced through
`OpfDeviceBuilder` and unused native currents are fixed to zero. New devices are
stamped through the model hook. `objective=:cost` sums native costs plus each
device's PCC cost, `:loss` minimizes custom series/ground losses, and
`:feasibility` sets zero. The optional hook can override the objective and add
problem-specific constraints/reference choices after all devices are available.
Other keywords are passed to `build_opf_model` (including scaling/optimizer).
"""
function build_generator_model(net,devices::AbstractVector;
        replacements=Dict{String,Tuple{Symbol,String}}(),objective::Symbol=:cost,
        generator_hook! = (ctx,handles)->nothing,kwargs...)
    all(d -> d isa _GeneralizedDevice,devices) || throw(ArgumentError("expected generator devices"))
    objective in (:cost,:loss,:feasibility) || throw(ArgumentError("unknown generator objective"))
    _gg_replacements(net,devices,replacements)
    foreach(d -> validate_device(d,(net,)),devices)
    handles=Dict{String,Any}()
    builders=Dict{Tuple{Symbol,String},BMOPFTools.OpfDeviceBuilder}()
    for d in devices
        haskey(replacements,d.id) || continue
        family,id=replacements[d.id]
        builders[(family,id)]=BMOPFTools.OpfDeviceBuilder(:PowerOptLab, (ctx,ids) -> begin
            _gg_retire_native!(ctx,family,id)
            handles[d.id]=stamp_device!(ctx,d)
            ctx
        end)
    end
    spec=BMOPFTools.OpfBuildSpec(component_builders=builders)
    hook! = ctx -> begin
        for d in devices
            haskey(handles,d.id) || (handles[d.id]=stamp_device!(ctx,d))
        end
        m=_opf_model(ctx)
        if objective==:cost
            JuMP.@objective(m,Min,BMOPFTools.generation_cost(ctx)+sum(h.cost for h in values(handles);init=0.0))
        elseif objective==:loss
            JuMP.@objective(m,Min,sum((h.series_loss+h.ground_loss)*(h.bases.power/1000) for h in values(handles);init=0.0))
        else
            JuMP.@objective(m,Min,0.0)
        end
        generator_hook!(ctx,handles)
    end
    ctx=BMOPFTools.initialize_opf_model(net;build_spec=spec,kwargs...)
    BMOPFTools.set_opf_start_values!(ctx)
    # The pinned engine skips bus limits at terminals listed in native source
    # metadata, even when their physics has a custom owner. During the public
    # limits stage, present only the remaining native sources. Restore the owned
    # working dictionary before the device stage (and on any failure).
    working=BMOPFTools.opf_network(ctx)
    sources=get(working,"voltage_source",Dict{String,Any}())
    replaced_sources=Set(id for (family,id) in values(replacements) if family==:voltage_source)
    if isempty(replaced_sources)
        BMOPFTools.add_opf_operational_limits!(ctx)
    else
        working["voltage_source"]=Dict(id=>s for (id,s) in sources if !(id in replaced_sources))
        try
            BMOPFTools.add_opf_operational_limits!(ctx)
        finally
            working["voltage_source"]=sources
        end
    end
    BMOPFTools.add_opf_device_constraints!(ctx)
    hook!(ctx)
    (context=ctx,handles=handles,devices=collect(devices),replacements=copy(replacements))
end

"""
    GeneratorOPFResult

Result of `solve_generator_opf`: authoritative custom `devices` in SI, native
`network` results with replaced placeholder records removed, `objective`,
`build_manifest`, and normalized `solve` status. The native `custom_injection`
power ledger includes custom PCC injections exactly once. It is not a dynamic
stability or global-optimality certificate.
"""
struct GeneratorOPFResult <: AbstractSolveResult
    devices::Dict{String,GeneratorResult}
    network::Dict{String,Any}
    objective::Float64
    build_manifest::BMOPFTools.OpfBuildManifest
    solve::SolveStatus
end
solve_status(r::GeneratorOPFResult)=r.solve
solve_diagnostics(r::GeneratorOPFResult)=(objective=r.objective,
    device_power_balance=Dict(id=>d.power_balance_error for (id,d) in r.devices))

"""
    solve_generator_opf(net, devices; solver_options=(), kwargs...)

Build with `build_generator_model`, enforce native network KCL, solve, and
extract generator/native results. All builder keywords are supported. Numerical
measurements are published only for fully feasible optimal/locally solved
outcomes under the package's standard status contract.
"""
function solve_generator_opf(net,devices::AbstractVector;solver_options=(),kwargs...)
    built=build_generator_model(net,devices;kwargs...)
    ctx=built.context; m=_opf_model(ctx)
    _set_solver_options!(m,solver_options)
    BMOPFTools.enforce_kcl!(ctx); JuMP.optimize!(m)
    outcome=_solve_outcome(m); status=SolveStatus(outcome)
    results=Dict(d.id=>extract_device(d,built.handles[d.id],status) for d in devices)
    network=_extract_result(ctx,outcome)
    for (family,id) in values(built.replacements)
        delete!(get(network,string(family),Dict()),id)
    end
    prior=get(network,"custom_injection",Dict("p"=>0.0,"q"=>0.0))
    network["custom_injection"]=Dict("p"=>get(prior,"p",0.0)+sum(d.p for d in values(results);init=0.0),
                                     "q"=>get(prior,"q",0.0)+sum(d.q for d in values(results);init=0.0))
    GeneratorOPFResult(results,network,status.publishable ? JuMP.objective_value(m) : NaN,
        BMOPFTools.opf_build_manifest(ctx),status)
end
