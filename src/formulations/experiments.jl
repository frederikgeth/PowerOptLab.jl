"""
    FormulationCase(id, build; metadata=NamedTuple())

A researcher-defined case. `build(representation, configuration)` returns a named
record with `model::JuMP.Model`, optional `observations` (graph or relation handles belonging to that model), and optional
`metrics` (a zero-argument callback returning a dictionary/named tuple). Each run
must build a fresh model. Cases may represent primitives, inverters, feeders or DOEs.
"""
struct FormulationCase{B,M}
    id::String
    build::B
    metadata::M
end
FormulationCase(id,build;metadata=NamedTuple()) = FormulationCase(string(id),build,metadata)

"""
    FormulationMethod(id, representation, optimizer;
                      options=[], configure! = identity, metadata=NamedTuple())

Representation and external solver configuration. The representation may be a
custom object or a named tuple selecting multiple primitives. `configure!(model)`
can install bridges before the optimizer is attached. Solver options and method
metadata are caller-owned; no continuation, retry or tolerance defaults are added.
"""
struct FormulationMethod{R,O,S,C,M}
    id::String
    representation::R
    optimizer::O
    options::S
    configure!::C
    metadata::M
end
FormulationMethod(id,r,optimizer;options=Pair{String,Any}[],configure! = identity,
                  metadata=NamedTuple()) =
    FormulationMethod(string(id),r,optimizer,options,configure!,metadata)

"""`UnsupportedFormulation(reason)` can be thrown by a case to record an unsupported combination without treating it as solver failure."""
struct UnsupportedFormulation <: Exception
    reason::String
end
Base.showerror(io::IO,e::UnsupportedFormulation) = print(io,e.reason)

"""
    run_formulation_experiment(cases, methods;
        configurations=[NamedTuple()], solve=true, assess=nothing,
        on_result=identity, on_error=:record)

Run the caller's case/method/configuration product with one fresh model per entry.
Raw solver outcomes, candidate audits, custom metrics and optional `assess(row)`
results remain distinct. Non-success candidates are available for diagnostics;
only strict MOI success sets `strict_solver_success`. No acceptance tolerance is
imposed. `on_error=:throw` stops on ordinary errors; interrupts always propagate.
`solve=false` records model structure without requiring a solver.

Callbacks are trusted user code. Timing includes first-use compilation and is not
by itself a controlled performance benchmark. Configurations should be immutable
or treated as read-only by callbacks.
"""
function run_formulation_experiment(cases,methods;configurations=[NamedTuple()],
        solve::Bool=true,assess=nothing,on_result=identity,on_error::Symbol=:record)
    on_error in (:record,:throw) || throw(ArgumentError("Choose :record or :throw"))
    rows = Dict{String,Any}[]
    # Materialize once so generator configurations can be reused across methods.
    configs = collect(configurations)
    method_list = collect(methods)
    for case in cases, method in method_list, (index,config) in enumerate(configs)
        row = Dict{String,Any}("case"=>case.id,"method"=>method.id,
            "configuration_index"=>index,"configuration"=>config,
            "case_metadata"=>case.metadata,"method_metadata"=>method.metadata,
            "solver_options"=>method.options,"run_status"=>"building")
        try
            started = time_ns()
            payload = case.build(method.representation,config)
            model = payload.model
            model isa JuMP.Model || throw(ArgumentError("Case must return a JuMP model"))
            method.configure!(model)
            method.optimizer === nothing || set_optimizer(model,method.optimizer)
            _set_solver_options!(model,method.options)
            row["build_seconds"] = (time_ns()-started)/1e9
            row["variables"] = num_variables(model)
            row["constraints"] = num_constraints(model;count_variable_in_set_constraints=true)
            row["run_status"] = "built"
            if solve
                row["run_status"] = "solving"
                row["solve_seconds"] = @elapsed optimize!(model)
                outcome = _solve_outcome(model)
                row["termination_status"] = string(outcome.termination_status)
                row["primal_status"] = string(outcome.primal_status)
                row["strict_solver_success"] = _publishable(outcome)
                certificate = outcome.primal_status in (MOI.INFEASIBILITY_CERTIFICATE,
                    MOI.NEARLY_INFEASIBILITY_CERTIFICATE,MOI.REDUCTION_CERTIFICATE,
                    MOI.NEARLY_REDUCTION_CERTIFICATE)
                row["primal_is_certificate"] = certificate
                row["candidate_available"] = outcome.has_primal && !certificate
                row["result_count"] = outcome.result_count
                # Optional backend diagnostics must not erase a returned primal
                # candidate if an external package's attribute getter fails.
                detail_errors = Dict{String,String}()
                for (key,getter) in (("raw_solver_status",raw_status),("dual_status",dual_status))
                    try
                        row[key] = string(getter(model))
                    catch error
                        error isa InterruptException && rethrow()
                        row[key] = nothing
                        detail_errors[key] = sprint(showerror,error)
                    end
                end
                row["solver_detail_errors"] = detail_errors
                row["run_status"] = "finished"
                if row["candidate_available"]
                    row["observations"] = [begin
                        h.model === model || throw(ArgumentError(
                            "Observation handle belongs to a different model"))
                        _formulation_observation(h)
                    end for h in get(payload,:observations,())]
                    row["metrics"] = haskey(payload,:metrics) ? payload.metrics() : NamedTuple()
                end
            end
            assess === nothing || (row["assessment"] = assess(row))
        catch error
            error isa InterruptException && rethrow()
            if error isa UnsupportedFormulation
                row["run_status"] = "unsupported"
                row["error"] = error.reason
            else
                on_error == :throw && rethrow()
                row["error_stage"] = row["run_status"]
                row["run_status"] = "error"
                row["error_type"] = string(typeof(error))
                row["error"] = sprint(showerror,error)
            end
        end
        push!(rows,row)
        on_result(row)
    end
    return rows
end

# Export only plain research data. Unsupported objects fail explicitly instead
# of embedding opaque solver/model representations that cannot be replayed.
_research_data(x::Union{AbstractString,Bool,Integer,AbstractFloat}) = x
_research_data(x::Symbol) = string(x)
_research_data(x::Enum) = string(x)
_research_data(::Nothing) = Dict("kind"=>"nothing")
_research_data(::Missing) = Dict("kind"=>"missing")
_research_data(x::Complex) = Dict("real"=>real(x),"imaginary"=>imag(x))
_research_data(x::Pair) = Dict("key"=>_research_data(first(x)),"value"=>_research_data(last(x)))
_research_data(x::Union{AbstractVector,Tuple}) = [_research_data(v) for v in x]
_research_data(x::Union{NamedTuple,AbstractDict}) = Dict(string(k)=>_research_data(v) for (k,v) in pairs(x))
_research_data(x) = throw(ArgumentError("Cannot export $(typeof(x)); convert custom metrics/configuration to plain data"))

"""
    write_formulation_results(path, rows; metadata=NamedTuple(), sources=[])

Write a versioned TOML result bundle with Julia/platform/core package versions and
SHA-256 hashes of caller-selected source files. The bundle does not automatically
capture a complete dependency lock or arbitrary callback code: supply a manifest,
configuration script and relevant source paths when reproducibility requires them.
Custom configuration/metrics must be plain data. Existing files are overwritten.
"""
function write_formulation_results(path::AbstractString,rows;metadata=NamedTuple(),sources=[])
    fingerprints = Dict(abspath(file)=>bytes2hex(sha256(read(file))) for file in sources)
    data = Dict("schema_version"=>1,"created_at_utc"=>string(now(UTC)),
        "runtime"=>Dict("julia"=>string(VERSION),"kernel"=>string(Sys.KERNEL),
            "architecture"=>string(Sys.ARCH),"PowerOptLab"=>string(Base.pkgversion(@__MODULE__)),
            "JuMP"=>string(Base.pkgversion(JuMP)),"BMOPFTools"=>string(Base.pkgversion(BMOPFTools))),
        "source_sha256"=>fingerprints,"metadata"=>_research_data(metadata),
        "runs"=>_research_data(rows))
    # Complete conversion before opening the file, so invalid metrics do not
    # truncate an existing result bundle.
    mkpath(dirname(abspath(path)))
    open(path,"w") do io
        TOML.print(io,data;sorted=true)
    end
    return abspath(path)
end

# Keep graph records compatible while allowing occurrence-specific relation audits.
function _formulation_observation(h::PWLFormulationHandle)
    merge(audit_pwl(h),
        (observation_kind=:graph,contract=formulation_contract(h.curve,h.formulation),
         formulation_type=string(typeof(h.formulation)),
         input_scale=h.input_scale,output_scale=h.output_scale,domain=h.domain,
         complementarity_scale=h.complementarity_scale,
         curve=(breakpoints=h.curve.breakpoints,values=h.curve.values,
                input_unit=h.curve.input_unit,output_unit=h.curve.output_unit)))
end
