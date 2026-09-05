module BoundedRelationExample
using PowerOptLab,JuMP,Ipopt
const INTENT=VoltVarWattIntent(volt_watt=PWLFunction([220.,240.,250.,270.],[1.,1.,0.,0.]))
function case()
    FormulationCase("bounded volt-watt limit",(rep,c)->begin
        m=Model()
        @variable(m,c.lower_V/230<=v<=c.upper_V/230,start=c.voltage_V/230)
        @variable(m,p>=0,start=.5)
        @constraint(m,v==c.voltage_V/230)
        h=formulate_control_relation!(m,INTENT,:volt_watt,v,p;
            relation=c.relation,formulation=rep,input_scale=230.,specialize=c.specialize,
            domain=(c.lower_V,c.upper_V))
        @objective(m,Max,p)
        (model=m,observations=[h])
    end)
end
function run(optimizer=Ipopt.Optimizer;options=(tol=1e-9,),graphs=false)
    reps = graphs ? (:auto,PWLConvexHull(),ExactPWLGraph()) :
        (:auto,PWLConvexHull(),SoftplusFormulation(.1),LocalC2Formulation(.1))
    methods=FormulationMethod[FormulationMethod(string(r),r,optimizer;
        configure! = set_silent,options) for r in reps]
    configurations=[(lower_V=lo,upper_V=hi,voltage_V=245.,relation=relation,specialize=s)
        for (lo,hi) in ((220.,250.),(220.,270.),(242.,248.)),
            relation in (:upper,:equal),s in (true,false)]
    # :auto is an exact-specialization request, so only run its specialized rows.
    rows=Dict{String,Any}[]
    for method in methods
        configs=method.representation===:auto ? filter(c->c.specialize,vec(configurations)) : vec(configurations)
        append!(rows,run_formulation_experiment([case()],[method];configurations=configs,on_error=:throw))
    end
    rows
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    using .BoundedRelationExample,PowerOptLab
    rows=BoundedRelationExample.run()
    for r in rows
        println(r["method"]," ",r["configuration"],": ",r["run_status"]," ",get(r,"observations",nothing))
    end
    haskey(ENV,"POL_FORMULATION_RESULTS") && write_formulation_results(
        ENV["POL_FORMULATION_RESULTS"],rows;sources=[@__FILE__])
end
