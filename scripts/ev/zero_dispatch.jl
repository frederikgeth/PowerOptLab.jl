# Bounded sensitivity study of the analytically known pc=pd=0 MPCC solution.
# Run in the same optional environment as ccopt.jl.
using PowerOptLab, JuMP, Ipopt, CCOpt, MathOptComplements, NLPModelsJuMP, TOML, SHA, Pkg, Test
include(joinpath(@__DIR__, "..", "..", "test", "fixtures.jl"))

# Solver-status publication and this experiment's physical acceptance are
# deliberately independent. A successful solver can miss a tighter SI budget.
function assess_zero_dispatch(publishable, pc, pd)
    within_budget = all(isfinite, (pc,pd)) &&
        max(abs(pc),abs(pd)) < 0.01 && abs(0.9pc-pd/0.9) < 0.01
    accepted = publishable && within_budget
    (; within_budget, accepted, charge_w=accepted ? pc : NaN,
       discharge_w=accepted ? pd : NaN)
end

function zero_dispatch_case(; per_unit=true, s_base=1e6, initialization=:default,
        strategy=:proportional, tol=1e-8, floor=1e-12, respect_bounds=false,
        bound_relaxation=1e-8)
    reduced = strategy == :branch_reference
    ev = EVDevice(id="full", bus="bus1", p_charge_max=1e4, p_discharge_max=reduced ? 0.0 : 1e4,
        energy_max=4e4, energy_init=4e4, energy_final=4e4,
        eff_charge=0.9, eff_discharge=0.9, available=[true],
        departure_energy=4e4, operation=:complementarity)
    model_ref = Ref{JuMP.Model}()
    function configure(m)
        model_ref[] = m
        MathOptComplements.Bridges.add_all_bridges(m)
        if initialization == :zero
            # The engine chains this hook after its KCL guard. Only device
            # starts change; upstream voltage/current initialization is kept.
            JuMP.set_optimize_hook(m, function (model; kwargs...)
                for name in ("pc_full", "pd_full", "mode_charge_full", "mode_discharge_full",
                             "cr_full_1", "ci_full_1")
                    set_start_value(variable_by_name(model, name), 0.0)
                end
                optimize!(model; ignore_optimize_hook=true, kwargs...)
            end)
        end
    end
    update = strategy == :rolloff ? CCOpt.RolloffRelaxationUpdate(sigma_min=floor) :
        CCOpt.ProportionalRelaxationUpdate(sigma_min=floor)
    options = reduced ? (; tol, bound_relax_factor=bound_relaxation,max_iter=1000) :
        (; tol, relaxation_update=update, respect_comp_bounds=respect_bounds,
            bound_relax_factor=bound_relaxation,max_iter=1000)
    result = solve_multiperiod_opf([single_bus_net(src_cost=-0.1, pload=1e4)], [ev];
        per_unit, s_base, optimizer=reduced ? Ipopt.Optimizer : CCOpt.Optimizer,
        configure! = configure, solver_options=options)
    m = model_ref[]
    sb = per_unit ? s_base : 1.0
    pc, pd = [has_values(m) ? value(variable_by_name(m, n))*sb : NaN for n in ("pc_full", "pd_full")]
    energy_error = 0.9pc-pd/0.9
    assessment = assess_zero_dispatch(solve_status(result).publishable,pc,pd)
    Dict("per_unit"=>per_unit, "s_base"=>s_base, "initialization"=>string(initialization),
        "strategy"=>string(strategy), "tol"=>tol, "sigma_floor"=>floor,
        "optimizer"=>reduced ? "Ipopt" : "CCOpt", "max_iter"=>1000,
        "respect_comp_bounds"=>respect_bounds,"bound_relax_factor"=>bound_relaxation,
        "termination_status"=>result.termination_status,
        "publishable"=>solve_status(result).publishable,
        "raw_charge_w"=>pc, "raw_discharge_w"=>pd,
        "raw_energy_balance_wh"=>energy_error,
        "raw_mode_minimum"=>min(pc,pd)/1e4,
        "raw_zero_solution_error_w"=>max(abs(pc),abs(pd)),
        "within_physical_budget"=>assessment.within_budget,
        "study_accepted"=>assessment.accepted,
        "study_charge_w"=>assessment.charge_w,"study_discharge_w"=>assessment.discharge_w,
        "published_charge_w"=>only(result.dispatch["full"].p_charge),
        "published_discharge_w"=>only(result.dispatch["full"].p_discharge),
        "published_dispatch_is_nan"=>all(isnan,result.dispatch["full"].p_charge))
end

function zero_dispatch_study()
    rows = Dict{String,Any}[]
    cases = ((label="baseline",), (label="zero_start", initialization=:zero),
        (label="equipment_base", s_base=1e4),
        (label="tight_proportional", tol=1e-10, floor=1e-16),
        (label="rolloff", strategy=:rolloff, tol=1e-10, floor=1e-16),
        (label="respect_bounds", respect_bounds=true, tol=1e-10, floor=1e-16),
        (label="no_bound_relaxation", bound_relaxation=0.0, tol=1e-10, floor=1e-16),
        (label="analytic_branch_reference", strategy=:branch_reference, bound_relaxation=0.0, tol=1e-10))
    for pu in (true,false), case in cases
        args = Base.structdiff(case, (; label=case.label))
        row = zero_dispatch_case(; per_unit=pu, args...)
        row["case"] = case.label
        push!(rows,row)
        println(row)
        flush(stdout)
    end
    repo = normpath(joinpath(@__DIR__,"..",".."))
    sources = ("scripts/ev/zero_dispatch.jl", "src/components/devices.jl",
        "src/problems/multiperiod.jl", "src/interfaces.jl", "test/fixtures.jl", "Project.toml")
    output = get(ENV,"POL_EV_ZERO_RESULTS",joinpath(tempdir(),"pol-ev-zero-results.toml"))
    open(output,"w") do io
        TOML.print(io,Dict("julia_version"=>string(VERSION), "runs"=>rows,
            "source_sha256"=>Dict(p=>bytes2hex(sha256(read(joinpath(repo,p)))) for p in sources),
            "packages"=>Dict(x.name=>string(x.version) for x in values(Pkg.dependencies())
                if x.name in ("CCOpt","MadNLP","Ipopt","JuMP","MathOptComplements","NLPModelsJuMP","MathOptInterface")));sorted=true)
    end
    @testset "Zero-dispatch reference and publication contract" begin
        # Linux CI returned LOCALLY_SOLVED at these powers. Regardless of
        # platform/termination variability, it must fail our unchanged 0.01 W budget.
        linux_candidate = assess_zero_dispatch(true,0.03350078299065128,0.027135634223650056)
        @test !linux_candidate.within_budget
        @test !linux_candidate.accepted
        @test isnan(linux_candidate.charge_w) && isnan(linux_candidate.discharge_w)
        @test assess_zero_dispatch(true,0.0,0.0).accepted
        @test !assess_zero_dispatch(false,0.0,0.0).accepted
        @test !assess_zero_dispatch(true,NaN,0.0).accepted
        for row in rows
            if row["case"] == "analytic_branch_reference"
                @test row["publishable"]
                @test row["within_physical_budget"]
            end
            if row["publishable"]
                # Publication follows solver status, not a study-specific budget.
                @test row["published_charge_w"] == row["raw_charge_w"]
                @test row["published_discharge_w"] == row["raw_discharge_w"]
            else
                @test row["published_dispatch_is_nan"]
            end
            if !row["publishable"] || !row["within_physical_budget"]
                @test !row["study_accepted"]
                @test isnan(row["study_charge_w"]) && isnan(row["study_discharge_w"])
            else
                @test row["study_accepted"]
                @test row["study_charge_w"] == row["raw_charge_w"]
            end
        end
    end
    rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    zero_dispatch_study()
end
