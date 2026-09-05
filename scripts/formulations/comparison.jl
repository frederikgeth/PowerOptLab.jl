module PWLComparison
using PowerOptLab, JuMP, Ipopt

# Exact real-current, resistive feeder: V = Vs + R*I, P = V*I.
# This is a deliberately small electrical reference, not a three-phase AC model.
const CURVE = PWLFunction([220.,240.,250.,270.],[20.,20.,0.,0.];
    input_unit=:V,output_unit=:A)
const ERROR_BUDGET_A = .01

function representations()
    # Signed hinge masses are +2 and -2 A/V for this particular curve.
    return (softplus=SoftplusFormulation(ERROR_BUDGET_A/(2log(2.))),
        local_c2=LocalC2Formulation(ERROR_BUDGET_A/(2*3/16)),
        complementarity=ComplementarityGraph(), exact_graph=ExactPWLGraph(),
        hull=PWLConvexHull())
end

"""Independent segment-by-segment solution of V=Vs+R*f(V), in physical units."""
function reference_equilibrium(; source_voltage=230., resistance=1.)
    xs,ys = CURVE.breakpoints,CURVE.values
    candidates = Float64[]
    for j in 1:length(xs)-1
        slope = (ys[j+1]-ys[j])/(xs[j+1]-xs[j])
        intercept = ys[j]-slope*xs[j]
        voltage = (source_voltage+resistance*intercept)/(1-resistance*slope)
        if xs[j]-1e-10 <= voltage <= xs[j+1]+1e-10
            all(v -> abs(v-voltage)>1e-9,candidates) && push!(candidates,voltage)
        end
    end
    length(candidates) == 1 || error("Reference needs exactly one equilibrium in the domain")
    voltage = only(candidates)
    return (voltage=voltage,current=(voltage-source_voltage)/resistance)
end

function run_case(representation; optimizer=Ipopt.Optimizer,
                  prepare! = identity, solver_options=["tol"=>1e-9],
                  source_voltage=230., resistance=1.,
                  voltage_scale=230.,current_scale=20.,start_voltage=245.)
    isfinite(resistance) && resistance > 0 || error("Resistance must be finite and positive")
    isfinite(source_voltage) || error("Source voltage must be finite")
    reference = reference_equilibrium(;source_voltage,resistance)
    m = Model()
    @variable(m,v,start=start_voltage/voltage_scale)
    h = formulate_pwl!(m,CURVE,v,representation;
        input_scale=voltage_scale,output_scale=current_scale)
    if h.output isa VariableRef
        set_start_value(h.output,primitive_value(CURVE,start_voltage)/current_scale)
    end
    for ((positive,negative),(_,k)) in zip(h.pairs,
            filter(ck -> first(CURVE.breakpoints)<ck[2]<last(CURVE.breakpoints),
                   PowerOptLab._pwl_hinges(CURVE)))
        # Independent physical starts, with both complementarity sides interior.
        set_start_value(positive,max(start_voltage-k,0.)/voltage_scale + .01)
        set_start_value(negative,max(k-start_voltage,0.)/voltage_scale + .01)
    end
    @constraint(m,v == source_voltage/voltage_scale +
        resistance*current_scale/voltage_scale*h.output)
    # Linear objective is shared by NLP, MPCC, MIP and LP. The canonical
    # equilibrium is unique; the hull can optimize over additional fake states.
    @objective(m,Max,h.output)
    prepare!(m)
    set_optimizer(m,optimizer)
    set_silent(m)
    for (key,value) in solver_options
        set_optimizer_attribute(m,key,value)
    end
    elapsed = @elapsed optimize!(m)
    contract = formulation_contract(CURVE,representation)
    row = Dict{String,Any}(
        "representation"=>string(nameof(typeof(representation))),
        "semantics"=>string(contract.semantics),
        "source_voltage_V"=>source_voltage,"resistance_ohm"=>resistance,
        "solver_options"=>Dict(string(k)=>v for (k,v) in solver_options),
        "voltage_scale"=>voltage_scale,
        "current_scale"=>current_scale,"start_voltage_V"=>start_voltage,
        "termination_status"=>string(termination_status(m)),
        "primal_status"=>string(primal_status(m)),
        "strict_solver_success"=>termination_status(m) in (MOI.OPTIMAL,MOI.LOCALLY_SOLVED) &&
            primal_status(m)==MOI.FEASIBLE_POINT,
        "elapsed_seconds"=>elapsed,"candidate_available"=>has_values(m),
        "reference_voltage_V"=>reference.voltage,"reference_current_A"=>reference.current)
    if has_values(m)
        a = audit_pwl(h)
        row["voltage_V"] = a.input
        row["current_A"] = a.output
        row["power_W"] = a.input*a.output
        row["electrical_residual_V"] = abs(a.input-source_voltage-resistance*a.output)
        row["domain_violation_V"] = a.domain_violation
        row["current_limit_violation_A"] = max(-a.output,a.output-20.,0.)
        row["exact_graph_error_A"] = a.exact_graph_error
        row["reference_current_error_A"] = a.output-reference.current
        row["complementarity_minimum_V"] = a.complementarity_minimum
        row["complementarity_product_V2"] = a.complementarity_product
        a.surrogate_equation_error === nothing ||
            (row["surrogate_equation_error_A"] = a.surrogate_equation_error)
        row["canonical_equations_satisfied"] = abs(a.exact_graph_error)<=1e-6 &&
            row["electrical_residual_V"]<=1e-6 && a.domain_violation<=1e-6 &&
            row["current_limit_violation_A"]<=1e-6
    end
    contract.width === nothing || (row["smoothing_width_V"] = contract.width)
    if contract.error_lower !== nothing
        row["error_lower_A"] = contract.error_lower
        row["error_upper_A"] = contract.error_upper
    end
    # No MPCC stationarity certificate is inferred from an MOI termination flag.
    row["stationarity_certificate"] = "not_independently_assessed"
    return row
end

function basic_comparison()
    r = representations()
    return [run_case(rep) for rep in (r.softplus,r.local_c2,r.hull)]
end
end # module

if abspath(PROGRAM_FILE) == @__FILE__
    foreach(println,PWLComparison.basic_comparison())
end
