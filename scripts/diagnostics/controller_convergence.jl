# Run in an optional diagnostic environment; see controller_convergence.md.
using PowerOptLab, BMOPFTools, JuMP, Ipopt, MadNLP, LinearAlgebra, TOML, Pkg
const POL = PowerOptLab
const MOI = JuMP.MOI
include(joinpath(pkgdir(PowerOptLab), "test", "fixtures.jl"))
LinearAlgebra.BLAS.set_num_threads(1)

function fixture(name)
    if name == "lcl_grid"
        inverter = AdvancedInverter(id="target", bus="poc", phase_terminals=["a","b","c"],
            neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40., i_grid_max=35.,
            v_dc=750., c_dc=1.1e-3, r_filter=.02, x_filter=.06,
            r_filter_grid=.03, x_filter_grid=.09, c_filter_mid=30e-6, r_filter_damping=.5)
        controller = SequenceController(PositiveSequenceVoltVarWatt(); current_target=:grid)
        return inv_grid3_bal(), ControlledDevice(inverter, controller), InverterControlRequest(p_available=9e3, q_scale=0.)
    end
    limit = name == "ripple_0.7" ? .7 : name == "ripple_0.1" ? .1 : error("Unknown fixture $name")
    inverter = AdvancedInverter(id="ripple_limit", bus="poc", phase_terminals=["a","b","c"],
        neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40., v_dc=700.,
        c_dc=1.1e-3, dv2_max=limit, r_filter=.05, x_filter=.15)
    return inv_grid3_unbal(), ControlledDevice(inverter, SequenceController(AverageVoltageVoltVarWatt())),
        InverterControlRequest(p_available=12e3, q_scale=8e3)
end

# Recognize precisely a*x^2 == 0. No approximate coefficient threshold is used.
function zero_squares(model)
    found = []
    for c in all_constraints(model; include_variable_in_set_constraints=false)
        obj = constraint_object(c)
        f = obj.func
        obj.set isa MOI.EqualTo || continue
        f isa JuMP.QuadExpr || continue
        constant(f) == obj.set.value || continue
        all(iszero(a) for (a, _) in linear_terms(f)) || continue
        terms = [(a,x,y) for (a,x,y) in quad_terms(f) if !iszero(a)]
        length(terms) == 1 || continue
        a,x,y = only(terms)
        x == y || continue
        push!(found, (c,x,a))
    end
    return found
end

function build_case(name, optimizer; remove_zero_squares=false)
    net, device, request = fixture(name)
    h = Ref{Any}()
    hook! = ctx -> begin
        h[] = stamp_device!(ctx, device; request=request)
        @objective(POL._opf_model(ctx), Min, h[].plant.p_loss + h[].plant.p_cap_loss)
    end
    ctx = build_opf_model(net; per_unit=true, s_base=1e6, add_objective=false,
        model_hook! = hook!, optimizer=optimizer, verbose=true)
    enforce_kcl!(ctx)
    model = POL._opf_model(ctx)
    roots = zero_squares(model)
    if remove_zero_squares
        for (c,x,_) in roots
            delete(model, c)
            fix(x, 0.; force=true)
            set_start_value(x, 0.)
        end
    end
    return model, h[], device, request, roots
end

# These are explicitly raw candidate diagnostics, even for a failed solve.
# Never fabricate a publishable SolveStatus to extract a failed candidate.
function candidate_metrics(h, device, request)
    p, c = h.plant, h.control
    pair(z, scale) = complex(value(z[1]),value(z[2])) * scale
    voltage = ntuple(k -> pair(c.phase_voltage[k], c.vb), 3)
    command = ntuple(k -> pair(c.phase_current[k], c.ib), 3)
    converter_voltage = ntuple(k -> complex(value(p.vrint[k]),value(p.viint[k])) * p.vb, 3)
    converter_current = ntuple(k -> complex(value(p.cri[k]),value(p.cii[k])) * p.ib, 3)
    grid_current = ntuple(k -> complex(value(p.gri[k]),value(p.gii[k])) * p.ib, 3)
    exact = evaluate_exact(device.controller, InverterControlMeasurement(collect(voltage)), request,
        InverterControlRatings(device.device, device.controller.current_target))
    exact = POL._apply_plant_capability_exact(exact, converter_voltage, converter_current,
        grid_current, device.device, device.controller.current_target)
    target = device.controller.current_target isa POL.GridCurrentTarget ? grid_current : converter_current
    inv = device.device
    ripple = abs(sum(converter_voltage .* converter_current)) / (2 * 2pi * inv.f * inv.c_dc * inv.v_dc)
    return Dict("candidate_dv2_V" => ripple,
        "candidate_current_scale" => value(c.current_scale),
        "candidate_exact_smooth_A" => maximum(abs, exact.phase_current .- command),
        "candidate_target_command_A" => maximum(abs, target .- command),
        "candidate_converter_current_A" => maximum(abs, converter_current))
end

function run_case(name, solver, variant, ablation)
    solver in ("ipopt", "madnlp") || error("Unknown solver $solver")
    variant in ("default", "unscaled", "unscaled_adaptive", "unscaled_unrelaxed") || error("Unknown variant $variant")
    optimizer = solver == "ipopt" ? Ipopt.Optimizer : MadNLP.Optimizer
    model, h, device, request, roots = build_case(name, optimizer; remove_zero_squares=ablation)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "tol", 1e-8)
    variant == "unscaled_unrelaxed" && set_optimizer_attribute(model, "bound_relax_factor", 0.)
    if solver == "ipopt"
        variant != "default" && set_optimizer_attribute(model, "nlp_scaling_method", "none")
        variant == "unscaled_adaptive" && set_optimizer_attribute(model, "mu_strategy", "adaptive")
    else
        variant != "default" && set_optimizer_attribute(model, "nlp_scaling", false)
        variant == "unscaled_adaptive" && set_optimizer_attribute(model, "barrier", MadNLP.QualityFunctionUpdate(1e-8, 10.))
    end
    label = "$name/$solver/$variant/zero_squares_removed=$ablation"
    println("CASE ",label); flush(stdout)
    optimize!(model)
    row = Dict{String,Any}("fixture"=>name, "solver"=>solver, "variant"=>variant,
        "zero_squares_removed"=>ablation, "zero_square_count"=>length(roots),
        "termination"=>string(termination_status(model)), "primal_status"=>string(primal_status(model)),
        "publishable"=>POL.SolveStatus(POL._solve_outcome(model)).publishable,
        "iterations"=>MOI.get(model, MOI.BarrierIterations()),
        "variables"=>num_variables(model))
    row["zero_square_coefficients"] = [a for (_,_,a) in roots]
    if solver == "madnlp"
        opt = unsafe_backend(model)
        row["solver_primal_infeasibility"] = opt.result.primal_feas
        row["solver_dual_infeasibility"] = opt.result.dual_feas
        row["solver_complementarity"] = opt.solver.inf_compl
    end
    if has_values(model)
        row["objective"] = objective_value(model)
        report = primal_feasibility_report(model; atol=0.)
        row["max_original_row_violation"] = maximum(values(report); init=0.)
        row["max_zero_root_value"] = maximum((abs(value(x)) for (_,x,_) in roots); init=0.)
        merge!(row, candidate_metrics(h,device,request))
    end
    println("RESULT ",row); flush(stdout)
    return row
end

function main()
    names = split(get(ENV,"POL_DIAG_CASES","ripple_0.7,ripple_0.1,lcl_grid"), ',')
    solvers = split(get(ENV,"POL_DIAG_SOLVERS","ipopt,madnlp"), ',')
    variants = split(get(ENV,"POL_DIAG_VARIANTS","default,unscaled,unscaled_adaptive,unscaled_unrelaxed"), ',')
    ablations = get(ENV,"POL_DIAG_ABLATIONS","both") == "both" ? (false,true) : (false,)
    packages = Dict(info.name => Dict("version"=>string(info.version), "revision"=>something(info.git_revision,""))
        for info in values(Pkg.dependencies()) if info.name in ("BMOPFTools","JuMP","Ipopt","MadNLP","Ipopt_jll","MUMPS_seq_jll"))
    output = Dict("julia"=>string(VERSION), "machine"=>Sys.MACHINE, "packages"=>packages,
        "poweroptlab_revision"=>readchomp(`git -C $(pkgdir(PowerOptLab)) rev-parse HEAD`),
        "blas"=>sprint(show, BLAS.get_config()), "runs"=>Any[])
    path = get(ENV,"POL_DIAG_OUTPUT","controller_convergence_results.toml")
    for name in names, solver in solvers, variant in variants, ablation in ablations
        push!(output["runs"], run_case(name,solver,variant,ablation))
        open(path,"w") do io
            TOML.print(io, output; sorted=true)
        end
    end
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
