# Run in the optional environment created by scripts/formulations/setup.jl.
# No integer solver is involved. The complete AC OPF is passed to CCOpt.
using PowerOptLab, JuMP, Test, TOML, Pkg, SHA
using CCOpt, MathOptComplements, NLPModelsJuMP
include(joinpath(@__DIR__,"..","..","test","fixtures.jl"))
include(joinpath(@__DIR__,"..","..","examples","ev_charging.jl"))

rows = Dict{String,Any}[]
@testset "EV / CCOpt exact charge-discharge modes" begin
    for pu in (true,false)
        full=EVDevice(id="full",bus="bus1",p_charge_max=10e3,p_discharge_max=10e3,
            energy_max=40e3,energy_init=40e3,eff_charge=0.9,eff_discharge=0.9,
            available=[true],departure_energy=40e3,energy_final=40e3,
            operation=:complementarity)
        model_ref=Ref{Any}()
        r=solve_multiperiod_opf([single_bus_net(src_cost=-0.1,pload=10e3)],[full];
            per_unit=pu,optimizer=CCOpt.Optimizer,
            configure! = m -> (model_ref[]=m; MathOptComplements.Bridges.add_all_bridges(m)),
            solver_options=(tol=1e-8,relaxation_update=CCOpt.ProportionalRelaxationUpdate(sigma_min=1e-12),max_iter=1000,))
        d=r.dispatch["full"]
        m=model_ref[]
        sb=pu ? 1e6 : 1.0
        rawpc=has_values(m) ? value(variable_by_name(m,"pc_full"))*sb : NaN
        rawpd=has_values(m) ? value(variable_by_name(m,"pd_full"))*sb : NaN
        push!(rows,Dict("case"=>"full_battery", "per_unit"=>pu,
            "termination_status"=>r.termination_status,
            "strict_solver_success"=>solve_status(r).publishable,
            "raw_candidate_charge_w"=>rawpc,"raw_candidate_discharge_w"=>rawpd,
            "raw_candidate_energy_balance_wh"=>0.9rawpc-rawpd/0.9,
            "raw_candidate_mode_minimum"=>min(rawpc,rawpd)/1e4,
            "stationarity_certificate"=>"not_independently_assessed"))
        # This biactive, zero-dispatch MPCC is a solver characterization case.
        # A local failure does not disprove feasibility: pc=pd=0 is an analytic
        # solution. Never turn a rejected raw candidate into published dispatch.
        if solve_status(r).publishable
            @test maximum(abs.(d.diagnostics.energy_balance_wh)) < 0.1
            @test maximum(d.diagnostics.complementarity_minimum) < 1e-4
            @test maximum(d.p_charge) < 2.0
            @test maximum(d.p_discharge) < 2.0
        else
            @test all(isnan,d.p_charge)
            @test all(isnan,d.energy_wh)
        end
        println(rows[end])
        flush(stdout)
    end
    # The converter must change mode across intervals and replenish initial
    # energy: initial-energy depletion cannot explain the expensive discharge.
    ev=EVDevice(id="cycle",bus="bus1",p_charge_max=10e3,p_discharge_max=10e3,
        energy_max=40e3,energy_init=20e3,eff_charge=0.9,eff_discharge=0.9,
        available=[true,true],departure_energy=20e3,energy_final=20e3,
        operation=:complementarity)
    r=solve_multiperiod_opf([single_bus_net(src_cost=p,pload=10e3) for p in (0.05,0.25)],
        [ev];optimizer=CCOpt.Optimizer,
        configure! = MathOptComplements.Bridges.add_all_bridges,solver_options=(tol=1e-8,relaxation_update=CCOpt.ProportionalRelaxationUpdate(sigma_min=1e-12),max_iter=1000,))
    @test solve_status(r).publishable
    d=r.dispatch["cycle"]
    @test d.p_charge[1] ≈ 10e3 atol=2.0
    @test d.p_discharge[2] ≈ 8.1e3 atol=2.0
    @test maximum(d.diagnostics.complementarity_minimum) < 1e-4
    @test maximum(abs.(d.diagnostics.energy_balance_wh)) < 0.1
    push!(rows,Dict("case"=>"replenished_cycle","per_unit"=>true,
        "termination_status"=>r.termination_status,
        "strict_solver_success"=>solve_status(r).publishable,
        "charge_w"=>d.p_charge,"discharge_w"=>d.p_discharge,
        "energy_wh"=>d.energy_wh,
        "mode_minimum"=>d.diagnostics.complementarity_minimum,
        "stationarity_certificate"=>"not_independently_assessed"))
    curve=PWLFunction([0.0,0.5,1.0],[7000.0,2000.0,0.0];input_unit=:unitless,output_unit=:W)
    car=EV(id="car",energy_max=40e3,onboard_charge_max=7000.0,charge_acceptance=curve)
    outlet=EVSE(id="outlet",bus="poc",s_max=7000.0,i_max=32.0)
    visit=ChargingSession(id="acceptance",ev=car,evse=outlet,available=[true],
        energy_init=10e3,departure_energy=10e3,acceptance_formulation=ComplementarityGraph())
    r=solve_multiperiod_opf([EVChargingExamples.grid(-0.1)],[visit];optimizer=CCOpt.Optimizer,
        configure! = MathOptComplements.Bridges.add_all_bridges,
        solver_options=(tol=1e-7,relaxation_update=CCOpt.ProportionalRelaxationUpdate(sigma_min=1e-12),max_iter=1000,))
    @test solve_status(r).publishable
    d=r.dispatch["acceptance"]
    @test d.p_charge[1] ≈ 3600.0 atol=0.1
    @test maximum(d.diagnostics.acceptance_violation_w) < 0.1
    push!(rows,Dict("case"=>"nonconcave_acceptance","per_unit"=>true,
        "nlp_tolerance"=>1e-7,
        "termination_status"=>r.termination_status,
        "strict_solver_success"=>solve_status(r).publishable,
        "charge_w"=>d.p_charge,"energy_wh"=>d.energy_wh,
        "acceptance_violation_w"=>d.diagnostics.acceptance_violation_w,
        "stationarity_certificate"=>"not_independently_assessed"))
end
if haskey(ENV,"POL_EV_RESULTS")
    packages=Dict(x.name=>string(x.version) for x in values(Pkg.dependencies())
        if x.name in ("CCOpt","MadNLP","MathOptComplements","NLPModelsJuMP","JuMP","MathOptInterface"))
    repo=normpath(joinpath(@__DIR__,"..",".."))
    sources=("src/components/devices.jl","src/components/evse.jl",
        "src/problems/multiperiod.jl","src/interfaces.jl","scripts/ev/ccopt.jl",
        "test/fixtures.jl","examples/ev_charging.jl","Project.toml")
    hashes=Dict(p=>bytes2hex(sha256(read(joinpath(repo,p)))) for p in sources)
    pin=TOML.parsefile(joinpath(repo,"Project.toml"))["sources"]["BMOPFTools"]["rev"]
    open(ENV["POL_EV_RESULTS"],"w") do io
        TOML.print(io,Dict("julia_version"=>string(VERSION),"packages"=>packages,
            "source_sha256"=>hashes,"bmopftools_revision"=>pin,
            "solver_options"=>Dict("tol"=>1e-8,"relaxation_update.sigma_min"=>1e-12,"max_iter"=>1000),
            "runs"=>rows);sorted=true)
    end
end
