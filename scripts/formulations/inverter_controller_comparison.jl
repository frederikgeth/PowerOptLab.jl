# Run with the environment created by scripts/formulations/setup.jl.
# One fresh model per method/case. Failures are recorded, never promoted.
using Test, PowerOptLab, JuMP, Ipopt, MadNLP, TOML
using MathOptComplements, NLPModelsJuMP, CCOpt
import Pkg
include(joinpath(@__DIR__,"..","..","test","fixtures.jl"))

function comparison_device(encoding, config)
    watt = PiecewiseLinearLaw([230.,240.,250.],[1.,1.,.2];smoothing_epsilon=.05)
    var = PiecewiseLinearLaw([210.,220.,240.,250.],[.3,0.,0.,-.3];smoothing_epsilon=.05)
    gain = PiecewiseLinearLaw([0.,.01,.10],[0.,0.,.08];smoothing_epsilon=1e-4)
    elaborate = config.case == "unbalanced_droop"
    positive = elaborate ? WorstPhaseVoltVarWatt(volt_watt=watt,volt_var=var,
        conflict_policy=:net,extrema_epsilon=1e-3,conflict_epsilon=1e-6) :
        config.case == "balanced_curve" ? AverageVoltageVoltVarWatt(volt_watt=watt) :
        AverageVoltageVoltVarWatt()
    unbalance = elaborate ? NegativeSequenceAdmittanceDroop(gain;
        impedance_angle=0.,ripple_blend=.5,voltage_floor=1.) : NoUnbalanceControl()
    lcl = config.case == "lcl_grid_target"
    controller = SequenceController(positive,unbalance,CommonScaleLimiter();
        encoding, current_target=lcl ? :grid : :converter)
    plant = AdvancedInverter(id="comparison",bus="poc",phase_terminals=["a","b","c"],
        neutral="n",topology=:THREE_LEG,
        s_max=config.case == "apparent_power_saturation" ? 8e3 : 20e3,
        i_max=40.,i_grid_max=lcl ? 35. : nothing,
        v_dc=lcl ? 750. : 700.,c_dc=1.1e-3,
        r_filter=lcl ? .02 : .05,x_filter=lcl ? .06 : .15,
        r_filter_grid=lcl ? .03 : 0.,x_filter_grid=lcl ? .09 : 0.,
        c_filter_mid=lcl ? 30e-6 : 0.,r_filter_damping=lcl ? .5 : 0.,m_max=.96)
    return ControlledDevice(plant,controller)
end

function comparison_metrics(ctx,handles,device,request)
    plant,control = handles.plant,handles.control
    phasors(re,im,scale) = [complex(value(r),value(i))*scale for (r,i) in zip(re,im)]
    voltage = [complex(value(p[1]),value(p[2]))*control.vb for p in control.phase_voltage]
    command = [complex(value(p[1]),value(p[2]))*control.ib for p in control.phase_current]
    vc = phasors(plant.vrint,plant.viint,plant.vb)
    ic = phasors(plant.cri,plant.cii,plant.ib)
    ig = phasors(plant.gri,plant.gii,plant.ib)
    exact = evaluate_exact(device.controller,InverterControlMeasurement(voltage),request,
        InverterControlRatings(device.device,device.controller.current_target))
    exact = PowerOptLab._apply_plant_capability_exact(exact,Tuple(vc),Tuple(ic),Tuple(ig),
        device.device,device.controller.current_target)
    max_equality = 0.; max_inequality = 0.
    worst_equality = ""
    model = PowerOptLab._opf_model(ctx)
    for (F,S) in list_of_constraint_types(model), ref in all_constraints(model,F,S)
        c = constraint_object(ref)
        if c.set isa MOI.EqualTo
            residual = abs(value(c.func)-c.set.value)
            if residual > max_equality
                max_equality = residual; worst_equality = string(ref)
            end
        elseif c.set isa MOI.LessThan
            max_inequality = max(max_inequality,value(c.func)-c.set.upper)
        elseif c.set isa MOI.GreaterThan
            max_inequality = max(max_inequality,c.set.lower-value(c.func))
        end
    end
    return (exact_current_replay_error_A=maximum(abs,command.-exact.phase_current),
        converter_current_margin_A=device.device.i_max-maximum(abs,ic),
        converter_apparent_power_margin_VA=device.device.s_max-abs(sum(vc.*conj.(ic))),
        max_equality_residual_raw_units=max_equality,
        max_inequality_violation_raw_units=max_inequality,worst_equality=worst_equality)
end

configs = [(case=name,s_base=1e6,selection_objective=:loss) for name in (
    "balanced_constant","balanced_curve","unbalanced_droop",
    "apparent_power_saturation","lcl_grid_target")]
case = controlled_inverter_case(
    c -> c.case == "unbalanced_droop" ? inv_grid3_unbal() : inv_grid3_bal(),
    comparison_device,
    c -> InverterControlRequest(p_available=c.case == "lcl_grid_target" ? 9e3 : 12e3,q_scale=8e3);
    id="controller_encoding_comparison",metrics=comparison_metrics)
methods = [
    FormulationMethod("smooth / Ipopt",:smooth,Ipopt.Optimizer;
        options=(tol=1e-8,max_iter=500),configure! = set_silent),
    FormulationMethod("smooth / MadNLP",:smooth,MadNLP.Optimizer;
        options=(tol=1e-8,max_iter=500),configure! = set_silent),
    FormulationMethod("MPCC / CCOpt 0.1.0",ComplementarityGraph(),CCOpt.Optimizer;
        options=(tol=1e-8,max_iter=500),
        configure! = m -> (MathOptComplements.Bridges.add_all_bridges(m);set_silent(m))),
    FormulationMethod("MPCC / CCOpt tighter relaxation",ComplementarityGraph(),CCOpt.Optimizer;
        options=(tol=1e-9,max_iter=500,bound_relax_factor=1e-12,
            relaxation_update=CCOpt.ProportionalRelaxationUpdate(sigma_mu_ratio=1e-4,sigma_min=1e-14)),
        configure! = m -> (MathOptComplements.Bridges.add_all_bridges(m);set_silent(m))),
    # Paired tests: only bound_relax_factor changes relative to each method above.
    FormulationMethod("smooth / MadNLP zero bound relaxation",:smooth,MadNLP.Optimizer;
        options=(tol=1e-8,max_iter=500,bound_relax_factor=0.0),configure! = set_silent),
    FormulationMethod("MPCC / CCOpt zero bound relaxation",ComplementarityGraph(),CCOpt.Optimizer;
        options=(tol=1e-8,max_iter=500,bound_relax_factor=0.0),
        configure! = m -> (MathOptComplements.Bridges.add_all_bridges(m);set_silent(m))),
    FormulationMethod("MPCC / CCOpt tighter relaxation with zero bounds",ComplementarityGraph(),CCOpt.Optimizer;
        options=(tol=1e-9,max_iter=500,bound_relax_factor=0.0,
            relaxation_update=CCOpt.ProportionalRelaxationUpdate(sigma_mu_ratio=1e-4,sigma_min=1e-14)),
        configure! = m -> (MathOptComplements.Bridges.add_all_bridges(m);set_silent(m))),
]
rows = run_formulation_experiment([case],methods;configurations=configs,on_error=:throw,
    on_result = r -> begin
        if hasproperty(r["solver_options"],:relaxation_update)
            update = r["solver_options"].relaxation_update
            description = (type="CCOpt.ProportionalRelaxationUpdate",
                sigma_mu_ratio=update.sigma_mu_ratio,sigma_mu_exp=update.sigma_mu_exp,
                sigma_min=update.sigma_min,monotone=update.monotone)
            r["solver_options"] = merge(r["solver_options"],(relaxation_update=description,))
        end
        metrics = get(r,"metrics",(;))
        gap = hasproperty(metrics,:custom) ? metrics.custom.exact_current_replay_error_A : NaN
        println(r["configuration"].case," | ",r["method"]," | ",r["termination_status"]," | replay A = ",gap)
        flush(stdout)
    end)
output = get(ENV,"POL_FORMULATION_RESULTS",joinpath(tempdir(),"controller-encoding-comparison.toml"))
versions = Dict(info.name=>string(info.version) for info in values(Pkg.dependencies()) if info.version !== nothing)
write_formulation_results(output,rows;
    sources=[@__FILE__,joinpath(@__DIR__,"..","..","src","components","inverter_controls.jl"),
        joinpath(@__DIR__,"..","..","src","formulations","selectors.jl"),
        joinpath(dirname(Base.active_project()),"Manifest.toml")],
    metadata=(purpose="Matched controller encoding failure report",stationarity="unassessed",
        package_versions=versions,acceptance="Raw solver outcomes and candidate audits are separate"))
@testset "Matched controller comparison records" begin
    @test length(rows)==35
    @test all(r -> r["run_status"]=="finished",rows)
    @test all(r -> r["strict_solver_success"],rows[1:5])
    for r in rows
        r["candidate_available"] || continue
        @test isfinite(r["metrics"].custom.exact_current_replay_error_A)
    end
end
println("Results: ",output)
