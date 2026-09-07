# Convert vendored Australian networks and validate determined power flows.
using PowerOptLab, BMOPFTools, Ipopt, JSON3
using SHA
const ROOT = normpath(joinpath(@__DIR__, "..", "..", "studies", "australian_pv"))
const OPTIONS = ("tol"=>1e-8, "max_iter"=>3000, "bound_relax_factor"=>0.0)
writejson(path, value) = open(io -> JSON3.pretty(io, value), path, "w")
phasor(t) = complex(t["vr"], t["vi"])
function validate_pf(net)
    pu = solve_pf(net; per_unit=true, solver_options=OPTIONS)
    si = solve_pf(net; per_unit=false, solver_options=OPTIONS)
    for r in (pu, si)
        @assert r["termination_status"] in ("LOCALLY_SOLVED", "OPTIMAL")
        @assert r["feasible"] && r["is_power_flow"]
        @assert all(isfinite(phasor(t)) for b in values(r["bus"]) for t in values(b))
    end
    delta = maximum(abs(phasor(t)-phasor(si["bus"][b][p]))
        for (b,ts) in pu["bus"] for (p,t) in ts)
    @assert delta < 0.01 "SI/per-unit phasor disagreement exceeds 0.01 V"
    pu, Dict("per_unit_status"=>pu["termination_status"],
        "SI_status"=>si["termination_status"], "max_SI_per_unit_difference_V"=>delta,
        "bus_count"=>length(net["bus"]), "customer_count"=>length(get(net,"load",Dict())),
        "losses"=>pu["losses"])
end
function prepare()
    mkpath(joinpath(ROOT,"networks")); mkpath(joinpath(ROOT,"results"))
    report = Dict{String,Any}()
    for name in ("LV3_55bus", "MV21_328bus")
        net = from_dss(joinpath(ROOT,"source",name,"Master.dss"); name=name)
        if name == "LV3_55bus"
            # Source Master intends 1 kW/customer via BatchEdit, which this
            # importer does not apply. Set both P and Q, preserving default PF=.88.
            for load in values(net["load"])
                load["p_nom"] = [1000.0]
                load["q_nom"] = [1000.0*tan(acos(0.88))]
            end
        end
        path = joinpath(ROOT,"networks",name*".bmopf.json")
        write_bmopf(net,path)
        # Validate the serialized/reparsed deliverable, not just memory state.
        parsed = parse_bmopf(path)
        r,audit = validate_pf(parsed)
        audit["sha256"] = bytes2hex(sha256(read(path)))
        audit["conversion_warnings"] = get(get(net,"_meta",Dict()),"powerio_warnings",[])
        report[name] = audit
        writejson(joinpath(ROOT,"results",name*"_pf.json"),r)
        println(name," validation: ",audit["per_unit_status"],", SI/PU difference ",audit["max_SI_per_unit_difference_V"]," V")
    end
    writejson(joinpath(ROOT,"results","validation.json"), report)
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    prepare()
end
