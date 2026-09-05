using Test, PowerOptLab, JuMP, Ipopt, PiecewiseLinearOpt, HiGHS
using MathOptComplements, NLPModelsJuMP, CCOpt
using Pkg, TOML, SHA
include("comparison.jl")

rows = Dict{String,Any}[]
# First collect the entire experiment, then persist evidence before assertions.
# CCOpt is a characterization backend: success and exact-law agreement are not
# assumed. Required exact-graph accuracy is tested against the independent MIP.
r = PWLComparison.representations()
for vs in (210.,220.,230.,250.,260.), (si,so) in ((1.,1.),(230.,20.)), start in (230.,260.)
    for (rep,optimizer,prepare,options) in (
        (r.exact_graph,HiGHS.Optimizer,identity,Pair{String,Any}[]),
        (r.complementarity,CCOpt.Optimizer,MathOptComplements.Bridges.add_all_bridges,
         ["tol"=>1e-9]))
        push!(rows,PWLComparison.run_case(rep; optimizer,prepare! = prepare,
            solver_options=options,source_voltage=vs,
            voltage_scale=si,current_scale=so,start_voltage=start))
    end
end
push!(rows,PWLComparison.run_case(r.hull;optimizer=HiGHS.Optimizer,solver_options=[]))

if haskey(ENV,"POL_FORMULATION_RESULTS")
    packages = Dict(info.name=>string(info.version) for info in values(Pkg.dependencies())
        if info.name in ("CCOpt","MadNLP","MathOptComplements","NLPModelsJuMP",
                         "HiGHS","PiecewiseLinearOpt","JuMP","MathOptInterface","DiffOpt"))
    bmopf = only(i for i in values(Pkg.dependencies()) if i.name == "BMOPFTools")
    repo = normpath(joinpath(@__DIR__,"..",".."))
    files = ("src/formulations/primitives.jl","src/formulations/jump.jl",
        "ext/PowerOptLabPiecewiseLinearOptExt.jl","scripts/formulations/comparison.jl",
        "scripts/formulations/optional_tests.jl","scripts/formulations/setup.jl")
    fingerprints = Dict(file=>bytes2hex(sha256(read(joinpath(repo,file)))) for file in files)
    open(ENV["POL_FORMULATION_RESULTS"],"w") do io
        TOML.print(io,Dict("julia_version"=>string(VERSION),"kernel"=>string(Sys.KERNEL),
            "architecture"=>string(Sys.ARCH),"packages"=>packages,"source_sha256"=>fingerprints,
            "bmopftools_revision"=>bmopf.git_revision,"runs"=>rows);sorted=true)
    end
end

@testset "External exact graph and CCOpt integration" begin
    for row in rows
        @testset "$(row["representation"]) Vs=$(row["source_voltage_V"]) scale=$(row["voltage_scale"]) start=$(row["start_voltage_V"])" begin
            @test row["candidate_available"]
            if row["candidate_available"]
                @test row["electrical_residual_V"] <= 1e-6
                @test row["domain_violation_V"] <= 1e-6
                if row["representation"] == "ExactPWLGraph"
                    @test row["strict_solver_success"]
                    @test row["canonical_equations_satisfied"]
                    @test abs(row["reference_current_error_A"]) <= 1e-6
                elseif row["representation"] == "ComplementarityGraph"
                    # Verify graph encoding against the independent reference
                    # where this external solver meets microunit accuracy. Kinks
                    # and normalized-coordinate runs remain reported experiments.
                    if row["voltage_scale"] == 1. && row["source_voltage_V"] != 250.
                        @test row["canonical_equations_satisfied"]
                        @test abs(row["reference_current_error_A"]) <= 1e-6
                    end
                    @test row["stationarity_certificate"] == "not_independently_assessed"
                else
                    @test row["strict_solver_success"]
                    @test row["semantics"] == "outer_relaxation"
                    @test !row["canonical_equations_satisfied"]
                    @test row["current_A"] ≈ 16. atol=1e-7
                    @test row["exact_graph_error_A"] ≈ 8. atol=1e-7
                end
            end
        end
    end
end
for kind in ("ExactPWLGraph","ComplementarityGraph","PWLConvexHull")
    subset = filter(r -> r["representation"]==kind,rows)
    println(kind,": ",count(r -> r["strict_solver_success"],subset),"/",length(subset),
        " strict solver successes; ",count(r -> get(r,"canonical_equations_satisfied",false),subset),
        "/",length(subset)," canonical-equation checks within 1e-6 physical units")
end

# Exercise the configurable runner with external graph and MPCC backends. These
# are integration checks, not a new convergence/normalization campaign.
configurable_case = resistive_control_case(PWLComparison.CURVE;
    source_voltage=230.,resistance=1.)
configurable_methods = [
    FormulationMethod("exact / HiGHS",ExactPWLGraph(),HiGHS.Optimizer;configure! = set_silent),
    FormulationMethod("hull / HiGHS",PWLConvexHull(),HiGHS.Optimizer;configure! = set_silent),
    FormulationMethod("MPCC / CCOpt",ComplementarityGraph(scale=1.),CCOpt.Optimizer;
        configure! = m -> (MathOptComplements.Bridges.add_all_bridges(m);set_silent(m)),
        options=(tol=1e-9,))]
configurable_rows = run_formulation_experiment([configurable_case],configurable_methods;
    configurations=[(input_scale=230.,output_scale=20.,start_input=245.)])
if haskey(ENV,"POL_FORMULATION_RESULTS")
    write_formulation_results(splitext(ENV["POL_FORMULATION_RESULTS"])[1]*"-configurable.toml",
        configurable_rows;sources=[@__FILE__],
        metadata=(purpose="optional backend API integration",stationarity="unassessed"))
end
@testset "Configurable external backend integration" begin
    @test all(r -> r["run_status"]=="finished",configurable_rows)
    @test configurable_rows[1]["strict_solver_success"]
    @test abs(only(configurable_rows[1]["observations"]).exact_graph_error)<1e-6
    @test configurable_rows[2]["strict_solver_success"]
    @test only(configurable_rows[2]["observations"]).exact_graph_error ≈ 8. atol=1e-6
    @test configurable_rows[3]["candidate_available"]
    @test only(configurable_rows[3]["observations"]).complementarity_scale==1.
end
