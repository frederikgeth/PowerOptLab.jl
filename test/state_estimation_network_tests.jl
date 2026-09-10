using Test, PowerOptLab, LinearAlgebra, SparseArrays, JuMP, Ipopt, Random
include("state_estimation_feeder_fixtures.jl")

@testset "SE network support contract" begin
    net = _se_net_autotransformer()
    report = state_estimator_preflight(net)
    @test report.supported
    @test only(filter(e->e.element=="transformer/single_phase_autotransformer/reg",report.elements)).tap == 1.05
    @test any(o->o.element=="load/ld",report.omitted)
    @test length(report.sources) == 1
    s = compile_state_estimator(net)
    @test s.preflight.supported
    p = SEParameters(s)
    x = initial_state_estimator(s,p)
    k = s.free_state_map[("reg","1")]
    @test x[k] ≈ 2520. atol=.1
    p.fixed_voltages[s.node_index[("src","1")]] = 2300. + 0im
    @test initial_state_estimator(s,p)[k] ≈ 2415. atol=.1
    for mutate in (
        n -> empty!(n["transformer"]["single_phase_autotransformer"]["reg"]),
        n -> (n["transformer"]["single_phase_autotransformer"]["reg"]["tap_ratio"] = 0.),
        n -> (n["transformer"]["single_phase_autotransformer"]["reg"]["terminal_map_to"] = ["missing"]),
        n -> (n["transformer"]["single_phase_autotransformer"]["reg"]["regulator_type"] = "C"),
        n -> (n["mystery_device"] = Dict("x"=>Dict("bus"=>"reg"))))
        bad=deepcopy(net); mutate(bad)
        @test !state_estimator_preflight(bad).supported
        @test_throws SEUnsupportedNetwork compile_state_estimator(bad)
    end
    @test !state_estimator_preflight(net; zero_injection=[("absent","1")]).supported
    @test !state_estimator_preflight(net; zero_injection=[("src","1")]).supported
    @test !state_estimator_preflight(net, Any[nothing]).supported
    limited=deepcopy(net)
    limited["transformer"]["single_phase_autotransformer"]["reg"]["tap_min"]=.9
    @test any(o->o.field=="tap_min",state_estimator_preflight(limited).omitted)
    ideal=deepcopy(net); t=ideal["transformer"]["single_phase_autotransformer"]["reg"]
    for key in ("r_series_from","r_series_to","x_series_from","x_series_to");t[key]=0.;end
    @test any(f->f.code==:ideal_transformer,state_estimator_preflight(ideal).findings)
    @test_throws SEUnsupportedNetwork compile_state_estimator(ideal)
    bad=_se_net_center_tap();bad["transformer"]["center_tap"]["ct"]["terminal_map_to"]=["1","n"]
    @test !state_estimator_preflight(bad).supported # exporter would otherwise skip it
    for build in (_se_net_1ph_xfmr, _se_net_yd_xfmr, _se_net_dy_xfmr, _se_net_center_tap)
        bad=build()
        for group in values(bad["transformer"]), transformer in values(group)
            for key in ("r_series_from","r_series_to","x_series_from","x_series_to"); transformer[key]=0.; end
        end
        @test !state_estimator_preflight(bad).supported
    end
    malformed=_se_net_autotransformer()
    malformed["shunt"]["empty"]=Dict{String,Any}("bus"=>"reg","terminal_map"=>["1"])
    @test !state_estimator_preflight(malformed).supported
    net2=_se_net_1ph_xfmr();ss=compile_state_estimator(net2);xx=initial_state_estimator(ss,SEParameters(ss))
    sc=PowerOptLab._se_voltage_scale(ss,SEParameters(ss),xx)
    @test all(isfinite,sc)
    @test maximum(sc)/minimum(sc) > 2 # separate voltage levels / neutral bus levels
end

# Independent instrument construction: terminal power comes from specified
# load currents and oracle voltages, not from the estimator's Ybus evaluator.
function _se_load_currents(net, voltage)
    currents=Dict{Tuple{String,String},ComplexF64}()
    for (_,load) in get(net,"load",Dict())
        bus=load["bus"]; ts=load["terminal_map"]; powers=complex.(load["p_nom"],load["q_nom"])
        config=get(load,"configuration","WYE")
        phases=filter(!=("n"),ts)
        pairs = config=="DELTA" ? [(phases[k],phases[mod1(k+1,length(phases))]) for k in eachindex(phases)] :
                config=="SINGLE_PHASE" ? [(ts[1],length(ts)>1 ? ts[2] : nothing)] :
                [(t,"n" in ts ? "n" : nothing) for t in phases]
        for ((a,b),S) in zip(pairs,powers)
            va=get(voltage,(bus,a),0im);vb=b===nothing ? 0im : get(voltage,(bus,b),0im)
            current=conj(S/(va-vb))
            currents[(bus,a)]=get(currents,(bus,a),0im)+current
            b===nothing || (currents[(bus,b)]=get(currents,(bus,b),0im)-current)
        end
    end
    currents
end

function _se_instrument_set(net, v)
    seed=compile_state_estimator(net); nf=length(seed.free_state_map)
    current=_se_load_currents(net,v)
    ms=Measurement[]; zi=Tuple{String,String}[]
    for (node,k) in sort!(collect(seed.free_state_map);by=last)
        bus,t=node; value=v[node]; level=max(abs(value),100.)
        push!(ms,Measurement(kind=:vr,bus=bus,terminal=t,reference=nothing,value=real(value),sigma=.002level),
                 Measurement(kind=:vi,bus=bus,terminal=t,reference=nothing,value=imag(value),sigma=.002level))
        if haskey(current,node) && abs(value)>1.
            S=-value*conj(current[node])
            push!(ms,Measurement(kind=:pinj,bus=bus,terminal=t,reference=nothing,value=real(S),sigma=max(.01abs(S),10.)),
                     Measurement(kind=:qinj,bus=bus,terminal=t,reference=nothing,value=imag(S),sigma=max(.01abs(S),10.)))
        elseif !haskey(current,node)
            push!(zi,node)
        end
    end
    ms,zi
end

# JuMP differentiates these scalar equations independently of residual_jacobian.
function _se_jump_reference(s,p,x0)
    model=Model(Ipopt.Optimizer);set_silent(model)
    set_optimizer_attribute(model,"tol",1e-7);set_optimizer_attribute(model,"max_iter",500)
    levels = 1 ./ PowerOptLab._se_voltage_scale(s,p,x0)
    @variable(model,y[j=1:length(x0)],start=0.)
    dx = levels .* y
    x = x0 + dx
    n=length(s.nodes); b=copy(s.fixed_voltage_state)
    for node in keys(s.reference_map)
        i=s.node_index[node];b[i]=real(p.fixed_voltages[i]);b[n+i]=imag(p.fixed_voltages[i])
    end
    u=(s.voltage_state_jacobian*x0+b)+s.voltage_state_jacobian*dx
    vf=complex.(b[1:n],b[n+1:end]);iff=s.passive_pattern*vf
    q=(s.current_state_jacobian*x0+vcat(real(iff),imag(iff)))+s.current_state_jacobian*dx
    for i in s.constraint_pattern
        current_level = max(norm(s.current_state_jacobian[i,:] .* levels),1.)
        @constraint(model,q[i]/current_level==0);@constraint(model,q[n+i]/current_level==0)
    end
    h=Any[]
    for m in s.measurement_pattern
        i,j=m.terminal,m.reference
        vr=u[i]-(j==0 ? 0 : u[j]);vi=u[n+i]-(j==0 ? 0 : u[n+j])
        push!(h,m.kind==:vr ? vr : m.kind==:vi ? vi : m.kind==:pinj ? vr*q[i]+vi*q[n+i] : vi*q[i]-vr*q[n+i])
    end
    @objective(model,Min,.5sum(((h[k]-p.measurement_values[k])/p.covariance_values[k])^2 for k in eachindex(h)))
    optimize!(model)
    is_solved_and_feasible(model) || @info "JuMP oracle status" termination_status(model) primal_status(model)
    @test is_solved_and_feasible(model)
    value.(x), objective_value(model)
end

if isdefined(Main,:OpenDSSDirect)
    function _se_dss_voltages(file; commands=String[])
        OpenDSSDirect.dss("Clear")
        OpenDSSDirect.dss("Redirect \"$(joinpath(@__DIR__,"data","state_estimation",file))\"")
        foreach(OpenDSSDirect.dss, commands)
        OpenDSSDirect.dss("Set MaxIterations=1000 algorithm=newton")
        OpenDSSDirect.dss("Solve")
        @test OpenDSSDirect.Solution.Converged()
        Dict(lowercase(name)=>ComplexF64(v) for (name,v) in zip(OpenDSSDirect.Circuit.AllNodeNames(),OpenDSSDirect.Circuit.AllBusVolts()))
    end
    @testset "SE finite transformers and fixed regulators: independent oracles" begin
        cases=[("single_phase",_se_net_1ph_xfmr,"pf_1ph_xfmr.dss"),
               ("single_phase_tap",_se_net_1ph_xfmr,"pf_1ph_xfmr.dss"),
               ("wye_delta",_se_net_yd_xfmr,"pf_yd_xfmr.dss"),
               ("delta_wye",_se_net_dy_xfmr,"pf_dy_xfmr.dss"),
               ("regulator",_se_net_autotransformer,"pf_autotransformer.dss"),
               ("open_delta",_se_net_open_delta_reg,"pf_open_delta_reg.dss"),
               ("center_tap",_se_net_center_tap,"center_tap.dss"),
               ("SE-IEEE13-finite",_se_ieee13_case,"")]
        for (name,build,file) in cases
            @testset "$name" begin
                net=build()
                commands=String[]
                if name == "single_phase_tap"
                    only(values(net["transformer"]["single_phase"]))["tap"] = 1.03
                    push!(commands,"Edit Transformer.t1 wdg=1 tap=1.03")
                end
                dss=isempty(file) ? _se_ieee13_opendss_voltages() : _se_dss_voltages(file;commands)
                mapping=Dict("a"=>"1","b"=>"2","c"=>"3","n"=>"4")
                v=Dict((String(bus),String(t))=>get(dss,lowercase("$bus.$(get(mapping,t,t))"),0im)
                       for (bus,b) in net["bus"] for t in b["terminal_names"])
                ms,zi=_se_instrument_set(net,v)
                s=compile_state_estimator(net,ms;zero_injection=zi);p=SEParameters(s,ms)
                for node in keys(s.reference_map);p.fixed_voltages[s.node_index[node]]=v[node];end
                x0=initial_state_estimator(s,p)
                nf=length(s.free_state_map);truth=zeros(2nf)
                for (node,k) in s.free_state_map;truth[k]=real(v[node]);truth[nf+k]=imag(v[node]);end
                voltage_levels = max.(hypot.(truth[1:nf],truth[nf+1:end]),100.)
                state_levels = vcat(voltage_levels,voltage_levels)
                rng=MersenneTwister(349)
                for noisy in (false,true)
                    p.measurement_values .= [m.value for m in ms]
                    noisy && (p.measurement_values .+= p.covariance_values .* randn(rng,length(ms)))
                    results=[solver(s,p,x0;max_iterations=250,initial_radius=.5,optimality_tolerance=1e-6) for solver in (solve_compiled_state_estimator,solve_sparse_state_estimator)]
                    for result in results
                        result.status == :converged_unique || @info "SE failed convergence" name noisy solver=typeof(result) status=result.status constraints=norm(result.evaluation.constraints) history=last(result.history)
                        @test result.status==:converged_unique
                        @test norm(result.evaluation.constraints) <= max(1e-7,1e-11PowerOptLab._se_current_scale(s,p))
                        @test maximum(abs.(result.state-truth)./state_levels) < (noisy ? .006 : .002)
                        @test observability_diagnostics(s,p,result.state).unobservable_dimension==0
                        covariance=selected_state_covariance(s,p,result.state,[1])
                        @test covariance[1,1]>=0 && isfinite(covariance[1,1])
                    end
                    @info "SE oracle metrics" name noisy voltage_error=maximum(abs.(results[1].state-truth)) relative_voltage_error=maximum(abs.(results[1].state-truth)./state_levels) objective=.5sum(abs2,results[1].evaluation.residual)
                    jx,jobj=_se_jump_reference(s,p,x0)
                    @test .5sum(abs2,results[1].evaluation.residual) ≈ jobj atol=1e-5 rtol=1e-5
                    @test maximum(abs.(jx-results[1].state)) < .01
                    @test results[1].state ≈ results[2].state atol=.01 rtol=1e-6
                end
            end
        end
    end
end
