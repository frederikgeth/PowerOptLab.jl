module ControlLoweringExample
using PowerOptLab, JuMP, Ipopt

# Canonical physical intent, reusable across encodings and independent models.
const INTENT = VoltVarWattIntent(
    volt_watt=PWLFunction([220.,240.,250.,270.],[1.,1.,0.,0.]),
    volt_var=PWLFunction([210.,220.,240.,250.],[.5,0.,0.,-.5]),
    sensing=:average_voltage,volt_watt_basis=:available)

function curve_case(;intent=INTENT,role=:volt_watt)
    FormulationCase("$(role) at a sensed voltage",(rep,config) -> begin
        m = Model()
        si = get(config,:input_scale,230.)
        so = get(config,:output_scale,1.)
        voltage = get(config,:voltage_V,245.)
        @variable(m,v,start=get(config,:start_voltage_V,voltage)/si)
        @constraint(m,v == voltage/si)
        encoding = role == :volt_watt ? VoltVarWattEncoding(volt_watt=rep) :
            VoltVarWattEncoding(volt_var=rep)
        h = formulate_control_curve!(m,intent,role,v,encoding;
            domain=get(config,:domain_V,(210.,280.)),input_scale=si,output_scale=so)
        # Maximization exposes the difference between the graph and its hull.
        @objective(m,Max,h.output)
        (model=m,observations=[h],metrics=() -> (
            fraction=value(h.output)*so,
            exact_fraction=primitive_value(getproperty(intent,role),voltage;
                domain_policy=:flat_extension),
            lowering=lowering_statistics(m)))
    end;metadata=(sensing="voltage supplied to curve port",role=role,
        scope="curve only; no AC network, sensing or capability equations"))
end

function smooth_methods(optimizer=Ipopt.Optimizer;solver_name="Ipopt",fraction_error=1e-3)
    FormulationMethod[FormulationMethod("$family / $solver_name",
        smoothing_for_error(INTENT.volt_watt,family,fraction_error),optimizer;
        configure! = set_silent,options=(tol=1e-9,),
        metadata=(optimizer=solver_name,fraction_error=fraction_error))
        for family in (SoftplusFormulation,LocalC2Formulation,AlgebraicFormulation)]
end

function run_core()
    methods = smooth_methods()
    push!(methods,FormulationMethod("hull / Ipopt",PWLConvexHull(),Ipopt.Optimizer;
        configure! = set_silent,options=(tol=1e-9,)))
    run_formulation_experiment([curve_case()],methods;on_error=:throw)
end

# Structural comparison, with no fragile elapsed-time pass/fail threshold.
function reuse_demo(;ports=10,reuse=true)
    m = Model(); @variable(m,v)
    e = VoltVarWattEncoding(volt_watt=LocalC2Formulation(.1))
    for _ in 1:ports
        formulate_control_curve!(m,INTENT,:volt_watt,v,e;domain=(210.,280.),reuse)
    end
    lowering_statistics(m)
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    using .ControlLoweringExample, PowerOptLab
    rows = ControlLoweringExample.run_core()
    for r in rows
        println(r["method"],": ",r["termination_status"],"; fraction=",
            r["metrics"].fraction,"; exact graph error=",only(r["observations"]).exact_graph_error)
    end
    println("Repeated port, shared: ",ControlLoweringExample.reuse_demo())
    println("Repeated port, fresh:  ",ControlLoweringExample.reuse_demo(reuse=false))
    if haskey(ENV,"POL_FORMULATION_RESULTS")
        write_formulation_results(ENV["POL_FORMULATION_RESULTS"],rows;sources=[@__FILE__])
    end
end
