# Rated-capability gate for the study fixture (roadmap R1/R2).
#
# A study inverter declared as 6 kVA must be able to deliver that nameplate at
# the AS/NZS 4777.2 Australia A volt-var extremes. Each corner fixes balanced POC
# voltages, requests the corner's reactive power and lets active power rise to
# the apparent-power or phase-current limit. A corner passes only if:
#   - the solve and physical report pass;
#   - the requested Q is delivered (within 1 var);
#   - the limiting quantity is POC apparent power or phase current, never an
#     internal budget (capacitor, peak, DC voltage, modulation).
# A fixture that fails this gate cannot represent certifiable hardware and must
# not be used for a feeder comparison.
include("engineering_records.jl")

# (POC phase-neutral voltage V, Q as a fraction of rated S). Australia A:
# -60 % at 258 V, +44 % at 207 V; ±60 % at 230 V covers the 0.8 PF range.
const RATED_CAPABILITY_CORNERS=((230.,0.),(258.,-.60),(207.,.44),(230.,-.60),(230.,.60))
const RATED_CAPABILITY_ALLOWED=("poc_apparent_power","converter_current_","grid_current_")

function rated_capability(inv;corners=RATED_CAPABILITY_CORNERS,solver_options=ENGINEERING_OPTIONS)
    a=cis(-2pi/3);rows=NamedTuple[]
    for (V,qf) in corners
        S=inv.s_max;q=qf*S
        request=InverterControlRequest(p_available=1.5S,p_rated=S,q_scale=S)
        voltage=(V+0im,V*a,V*conj(a))
        result=inv.topology===:THREE_LEG ?
            allocate_inverter_services(inv,request;grid_voltage=voltage,reactive_request=q,solver_options) :
            allocate_phase_inverter_services(inv,request;grid_voltage=voltage,reactive_request=(q/3,q/3,q/3),solver_options)
        status=solve_status(result);report=physical_feasibility(result)
        plant=result.allocation.plant
        limits=engineering_active_limits(report)
        internal=[l for l in limits if !any(startswith(l,p) for p in RATED_CAPABILITY_ALLOWED)]
        value(name)=(i=findfirst(c->c.name==name,report.checks);i===nothing ? NaN : report.checks[i].value)
        accepted=status.publishable && report.passed && abs(plant.q_poc-q)<=1. &&
            !isempty(limits) && isempty(internal)
        push!(rows,(topology=String(inv.topology),voltage_V=V,q_fraction=qf,q_request_var=q,
            p_poc_W=plant.p_poc,q_poc_var=plant.q_poc,s_poc_VA=hypot(plant.p_poc,plant.q_poc),
            capacitor_A=value("capacitor_current"),capacitor_limit_A=inv.i_cap_max,
            limiting=join(sort!(unique(limits)),"|"),internal_limits=join(sort!(unique(internal)),"|"),
            accepted=accepted))
    end
    rows
end

"""Run the gate on one plant of each topology in a benchmark fleet."""
function fleet_rated_capability(mixed)
    seen=Set{Symbol}();rows=NamedTuple[]
    for id in sort!(collect(keys(mixed.spec.devices)))
        inv=mixed.spec.devices[id].device
        inv.topology in seen && continue
        push!(seen,inv.topology);append!(rows,rated_capability(inv))
    end
    rows
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    limits=get(ARGS,1,"declared")=="legacy" ? ENGINEERING_LEGACY_LIMITS : ENGINEERING_OPERATING_LIMITS
    condition=first(engineering_conditions())
    net=voltage_aware_network(parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json")),condition)
    rows=fleet_rated_capability(engineering_fleet(net,condition;limits))
    for r in rows
        println(rpad(r.topology,10),rpad(r.voltage_V,7),rpad(r.q_fraction,6)," S=",round(r.s_poc_VA)," VA",
            " Icap=",round(r.capacitor_A,digits=2),"/",r.capacitor_limit_A," limiting=",r.limiting,
            r.accepted ? "  PASS" : "  FAIL")
    end
    all(r.accepted for r in rows) || exit(1)
end
