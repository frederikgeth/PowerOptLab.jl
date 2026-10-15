using Test, PowerOptLab, LinearAlgebra, Random
using BMOPFTools: parse_bmopf

const CONSISTENCY_SE_SOLVERS = (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator)

# Independent two-conductor circuit. All reference predictions below use
# series Ohm's law and electrode KCL, never Ybus or a package evaluator.
function consistency_circuit(; voltage_gain=1., current_gain=1., angle=.37,
                               reversed=false, swapped=false, renamed=false)
    source=230cis(angle)*voltage_gain
    impedance=ComplexF64[.3+.12im .04+.015im; .04+.015im .6+.2im] .* (voltage_gain/current_gain)
    electrode=(.17-.07im)*(current_gain/voltage_gain)
    supply=renamed ? "z_supply" : "src"; load=renamed ? "a_customer" : "b"
    terms=swapped ? ["n","p"] : ["p","n"]
    perm=swapped ? [2,1] : [1,2]; z=impedance[perm,perm]
    code=Dict{String,Any}()
    for i in 1:2,j in 1:2
        code["R_series_$(i)_$(j)"]=real(z[i,j]);code["X_series_$(i)_$(j)"]=imag(z[i,j])
    end
    net=Dict{String,Any}(
        "bus"=>Dict(supply=>Dict("terminal_names"=>terms,"perfectly_grounded_terminals"=>["n"]),
                    load=>Dict("terminal_names"=>terms)),
        "voltage_source"=>Dict("source"=>Dict("bus"=>supply,"terminal_map"=>["p"],
            "v_magnitude"=>[abs(source)],"v_angle"=>[angle])),
        "linecode"=>Dict("code"=>code),
        "line"=>Dict("line"=>Dict("bus_from"=>reversed ? load : supply,"bus_to"=>reversed ? supply : load,
            "terminal_map_from"=>terms,"terminal_map_to"=>terms,"linecode"=>"code","length"=>1.)),
        "shunt"=>Dict("earth"=>Dict("bus"=>load,"terminal_map"=>["n"],
            "G_1_1"=>real(electrode),"B_1_1"=>imag(electrode))))
    injected=-(9-3im)*cis(angle)*current_gain
    admittance=inv(impedance)+Diagonal(ComplexF64[0,electrode])
    voltage=admittance \ (ComplexF64[injected,-injected]+impedance\ComplexF64[source,0])
    (;net,source,impedance,electrode,supply,load,voltage,voltage_gain,current_gain,reversed)
end

function consistency_affine_prediction(case,x)
    v=complex.(x[1:2],x[3:4])
    line=case.impedance \ (v-ComplexF64[case.source,0])
    shunt=case.electrode*v[2]
    net=line+ComplexF64[0,shunt]
    vcat(real.(v),imag.(v),real.(line),imag.(line),real.(net),imag.(net),real(shunt),imag(shunt))
end
function consistency_balance(case,x)
    v=complex.(x[1:2],x[3:4])
    current=sum(case.impedance \ (v-ComplexF64[case.source,0]))+case.electrode*v[2]
    [real(current),imag(current)]
end
function consistency_affine_operator(f,n;step=1.)
    b=f(zeros(n));A=hcat([(f(step.*Float64.(1:n .== k))-f(-step.*Float64.(1:n .== k)))./(2step) for k in 1:n]...)
    A,b
end
function consistency_meters(case; seed=823, noisy=true)
    truth=vcat(real.(case.voltage),imag.(case.voltage))
    sigmas=vcat([.2,.08,.13,.3].*case.voltage_gain,
                [.03,.07,.05,.04,.09,.02,.06,.08,.015,.025].*case.current_gain)
    values=consistency_affine_prediction(case,truth)
    noisy && (values+=sigmas.*randn(MersenneTwister(seed),length(values)))
    meters=Any[];row=1
    for kind in (:vr,:vi),terminal in ("p","n")
        push!(meters,Measurement(;kind,bus=case.load,terminal,reference=nothing,value=values[row],sigma=sigmas[row]));row+=1
    end
    for net in (false,true),kind in (:ire,:iim),terminal in ("p","n")
        meter=net ? NetInjectionMeasurement(;kind,bus=case.load,terminal,value=values[row],sigma=sigmas[row]) :
            BranchMeasurement(;kind,line="line",side=case.reversed ? :from : :to,terminal,value=values[row],sigma=sigmas[row])
        push!(meters,meter);row+=1
    end
    for kind in (:ire,:iim)
        push!(meters,ShuntMeasurement(;kind,shunt="earth",terminal="n",value=values[row],sigma=sigmas[row]));row+=1
    end
    meters,values,sigmas
end
function consistency_logical_order(s,load)
    nf=length(s.free_state_map);a=s.free_state_map[(load,"p")];b=s.free_state_map[(load,"n")]
    [a,b,nf+a,nf+b]
end

@testset "SE independent affine oracle: solution, covariance and representation" begin
    # Same statistical problem under different names, terminal order, branch
    # direction, reading order and SI voltage/current magnitudes.
    reference=nothing
    for config in ((;),(reversed=true,),(swapped=true,renamed=true,),
                   (voltage_gain=.001,current_gain=.01,),
                   (voltage_gain=1000.,current_gain=100.,))
        case=consistency_circuit(;config...)
        meters,values,sigmas=consistency_meters(case)
        H,b=consistency_affine_operator(x->consistency_affine_prediction(case,x),4;step=230case.voltage_gain)
        C,d=consistency_affine_operator(x->consistency_balance(case,x),4;step=230case.voltage_gain)
        # Analytic equality-constrained linear WLS via a particular solution
        # and null-space QR. This problem has a unique global voltage optimum.
        feasible=C\(-d);Z=nullspace(C)
        Hw=H./sigmas;target=(values-b)./sigmas
        expected=feasible+Z*((Hw*Z)\(target-Hw*feasible))
        covariance=Z*inv((Hw*Z)'*(Hw*Z))*Z'
        objective=.5sum(abs2,Hw*expected-target)
        perm=case.reversed || case.voltage_gain!=1 ? reverse(eachindex(meters)) : eachindex(meters)
        meters=meters[perm]
        s=compile_state_estimator(case.net,meters;injection_connections=[TerminalConnection((case.load,"p"),(case.load,"n"))])
        p=SEParameters(s,meters);order=consistency_logical_order(s,case.load)
        # Evaluate a point away from the answer, so a shared solver/model error
        # cannot pass just by returning the supplied initial state.
        trial=initial_state_estimator(s,p);trial[order]=vcat(real.(case.voltage),imag.(case.voltage)).+.03case.voltage_gain
        e=evaluate_state_estimator(s,p,trial)
        @test e.predicted ≈ (H*trial[order]+b)[perm] atol=1e-8max(case.voltage_gain,case.current_gain) rtol=1e-10
        @test Matrix(residual_jacobian(s,p,trial))[:,order] ≈ (H./sigmas)[perm,:] rtol=1e-9 atol=1e-9/case.voltage_gain
        for solver in CONSISTENCY_SE_SOLVERS
            r=solver(s,p,initial_state_estimator(s,p);optimality_tolerance=1e-7)
            @test r.status == :converged_unique
            @test r.state[order] ≈ expected atol=2e-6case.voltage_gain rtol=1e-8
            @test .5sum(abs2,r.evaluation.residual) ≈ objective atol=1e-7 rtol=1e-8
            @test norm(consistency_balance(case,r.state[order])) < 1e-7case.current_gain
            @test selected_state_covariance(s,p,r.state,order) ≈ covariance rtol=1e-8 atol=1e-10case.voltage_gain^2
        end
        normalized=(state=expected./case.voltage_gain,covariance=covariance./case.voltage_gain^2,objective=objective)
        if reference===nothing
            reference=normalized
        else
            @test normalized.state ≈ reference.state atol=1e-8
            @test normalized.covariance ≈ reference.covariance rtol=1e-9
            @test normalized.objective ≈ reference.objective atol=1e-8
        end
    end
end

@testset "SE nonlinear noisy fit: native SI/per-unit and all shared solvers" begin
    for seed in (27,912),angle in (-.41,.63)
        case=consistency_circuit(;angle)
        v=case.voltage;injected=(case.impedance\(v-ComplexF64[case.source,0]))[1]
        power=(v[1]-v[2])*conj(injected)
        raw=[real(v[1]),imag(v[1]),real(v[2]),imag(v[2]),real(power),imag(power)]
        sigma=[.2,.3,.08,.12,20.,15.]
        raw+=sigma.*randn(MersenneTwister(seed),6)
        meters=Measurement[]
        for (k,(kind,terminal,reference)) in enumerate(((:vr,"p",nothing),(:vi,"p",nothing),
                (:vr,"n",nothing),(:vi,"n",nothing),(:pinj,"p","n"),(:qinj,"p","n")))
            push!(meters,Measurement(;kind,bus=case.load,terminal,reference,value=raw[k],sigma=sigma[k]))
        end
        native=[solve_state_estimation(case.net,meters;per_unit,solver_options=("tol"=>1e-10,)) for per_unit in (false,true)]
        @test all(r->solve_status(r).publishable,native)
        voltages=[ComplexF64[complex(r.bus[case.load][t]["vr"],r.bus[case.load][t]["vi"]) for t in ("p","n")] for r in native]
        @test voltages[1] ≈ voltages[2] atol=1e-6
        @test native[1].objective ≈ native[2].objective atol=1e-7
        s=compile_state_estimator(case.net,meters);p=SEParameters(s,meters);order=consistency_logical_order(s,case.load)
        expected=vcat(real.(voltages[1]),imag.(voltages[1]))
        for solver in CONSISTENCY_SE_SOLVERS
            r=solver(s,p,initial_state_estimator(s,p);optimality_tolerance=1e-6)
            solve_status(r).publishable || @info "Noisy backend inconsistency" seed angle solver status=r.status diagnostics=solve_diagnostics(r)
            @test solve_status(r).publishable
            @test r.state[order] ≈ expected atol=1e-5
            @test sum(abs2,r.evaluation.residual) ≈ native[1].objective atol=1e-7 rtol=1e-7
            @test norm(consistency_balance(case,r.state[order]))<1e-7
        end
    end
end

function consistency_floating_case()
    net=parse_bmopf("""
    {"bus":{"s":{"terminal_names":["p","n"]},"b":{"terminal_names":["p","n"]}},
     "linecode":{"z":{"R_series_1_1":0.5,"R_series_2_2":0.5}},
     "line":{"l":{"bus_from":"s","bus_to":"b","terminal_map_from":["p","n"],
                   "terminal_map_to":["p","n"],"linecode":"z","length":1.}}}
    """;from_string=true)
    voltage=Dict(("s","p")=>233+2im,("s","n")=>3+2im,("b","p")=>228+3im,("b","n")=>8+im)
    meters=[Measurement(;kind,bus=node[1],terminal=node[2],reference=nothing,value,sigma=.2)
        for (node,v) in sort!(collect(voltage);by=first)
        for (kind,value) in ((:vr,real(v)),(:vi,imag(v)))]
    net,voltage,meters
end
function consistency_pack(s,voltage)
    nf=length(s.free_state_map);x=zeros(2nf)
    for (node,k) in s.free_state_map
        x[k]=real(voltage[node]);x[nf+k]=imag(voltage[node])
    end
    x
end

@testset "SE floating exact laws: redundant KCL retains the physical optimum" begin
    net,voltage,meters=consistency_floating_case()
    c=[TerminalConnection(("b","p"),("b","n"))]
    current=10-2im;power=(220+2im)*conj(current)
    for model in (ConstantCurrentDevice(c,[current]),ConstantPowerDevice(c,[power]),
                  ZIPDevice(c,[.5power],[.25current],[.25current/(220+2im)]))
        exact=[ExactDeviceEquation(model)]
        s=compile_state_estimator(net,meters;exact_devices=exact,
            injection_connections=[TerminalConnection(("s","p"),("s","n"))])
        p=SEParameters(s,meters;exact_devices=exact);truth=consistency_pack(s,voltage)
        initial=truth+.1cos.(eachindex(truth))
        @test size(constraint_jacobian(s,p,truth))==(6,8)
        @test rank(Matrix(constraint_jacobian(s,p,truth)))==4
        snapshot=deepcopy(p)
        for solver in CONSISTENCY_SE_SOLVERS
            r=solver(s,p,initial)
            @test r.status == :converged_unique
            @test r.state ≈ truth atol=1e-7
            @test norm(r.evaluation.constraints)<1e-8
            @test norm(r.evaluation.residual)<1e-6
            if solver===solve_ipopt_state_estimator
                diagnostics=solve_diagnostics(r)
                @test diagnostics.physical_constraint_count==6
                @test diagnostics.optimizer_constraint_count==4
                # Identities hold far from the initial Jacobian, including
                # changed nonlinear currents, not just at the fitted point.
                signature=PowerOptLab._se_ipopt_constraint_signature(s,p)
                rows=diagnostics.optimizer_constraint_rows
                for x in (truth,truth+.7sin.(eachindex(truth)))
                    value=evaluate_state_estimator(s,p,x).constraints
                    map=signature*pinv(signature[rows,:])
                    @test value ≈ map*value[rows] atol=1e-10
                end
            end
        end
        @test p.device_powers==snapshot.device_powers
        @test p.continuation_alpha==snapshot.continuation_alpha
        @test p.regularization_voltage==snapshot.regularization_voltage
    end
end

@testset "SE equation reduction: constants and singular initial Jacobians" begin
    net,voltage,meters=consistency_floating_case()
    connections=[TerminalConnection((bus,"p"),(bus,"n")) for bus in ("s","b")]
    for balanced in (true,false)
        exact=[ExactDeviceEquation(ConstantCurrentDevice(connections,[-10+2im+(balanced ? 0 : .1),10-2im]))]
        s=compile_state_estimator(net,meters;exact_devices=exact);p=SEParameters(s,meters;exact_devices=exact)
        initial=consistency_pack(s,voltage)
        for solver in CONSISTENCY_SE_SOLVERS
            r=solver(s,p,initial)
            @test solve_status(r).publishable==balanced
            if balanced
                @test norm(r.evaluation.constraints)<1e-8
                @test r.state ≈ initial atol=1e-7
            else
                @test norm(r.evaluation.constraints)>1e-3
            end
        end
    end
    # At the two-bus maximum-power point an exact equation has zero first
    # derivative, but is not an identity. Preserve both real/imag equations.
    grounded=parse_bmopf("""
    {"bus":{"s":{"terminal_names":["p"]},"b":{"terminal_names":["p"]}},
     "voltage_source":{"v":{"bus":"s","terminal_map":["p"],"v_magnitude":[230.],"v_angle":[0.]}},
     "linecode":{"z":{"R_series_1_1":1.}},
     "line":{"l":{"bus_from":"s","bus_to":"b","terminal_map_from":["p"],"terminal_map_to":["p"],"linecode":"z","length":1.}}}
    """;from_string=true)
    exact=[ExactDeviceEquation(ConstantPowerDevice([TerminalConnection(("b","p"),nothing)],ComplexF64[13225]))]
    s=compile_state_estimator(grounded;exact_devices=exact);p=SEParameters(s;exact_devices=exact)
    @test rank(Matrix(constraint_jacobian(s,p,[115.,0.])))==1
    @test length(PowerOptLab._se_ipopt_constraint_rows(s,p))==2
    p.regularization_voltage=1e-300
    @test length(PowerOptLab._se_ipopt_constraint_rows(s,p))==2
    # Weak grounding is physical information. The structural reduction must
    # not use the looser first-order observability threshold to discard it.
    weak,_,meters=consistency_floating_case()
    weak["shunt"]=Dict("earth"=>Dict("bus"=>"b","terminal_map"=>["n"],"G_1_1"=>1e-10))
    exact=[ExactDeviceEquation(ConstantCurrentDevice(connections,ComplexF64[0,0]))]
    sw=compile_state_estimator(weak,meters;exact_devices=exact);pw=SEParameters(sw,meters;exact_devices=exact)
    @test length(PowerOptLab._se_ipopt_constraint_rows(sw,pw))>=6 # full-system fallback is also safe
end

@testset "SE gauge freedom: compare observable quantities and uncertainty" begin
    net,voltage,_=consistency_floating_case()
    meters=Measurement[];means=Float64[];variances=Float64[];objective=0.
    for bus in ("s","b"),kind in (:vr,:vi)
        drop=voltage[(bus,"p")]-voltage[(bus,"n")]
        base=kind===:vr ? real(drop) : imag(drop)
        values=base.+[.2,-.35];sigma=[.3,.7];weights=1 ./ sigma.^2
        mean=sum(weights.*values)/sum(weights)
        push!(means,mean);push!(variances,1/sum(weights))
        objective+=.5sum(abs2,(values.-mean)./sigma)
        for (value,sd) in zip(values,sigma)
            push!(meters,Measurement(;kind,bus,terminal="p",reference="n",value,sigma=sd))
        end
    end
    cs=[TerminalConnection((bus,"p"),(bus,"n")) for bus in ("s","b")]
    s=compile_state_estimator(net,meters;injection_connections=cs)
    nf=length(s.free_state_map);J=zeros(4,2nf)
    for (row,(bus,kind)) in enumerate((("s",:vr),("s",:vi),("b",:vr),("b",:vi)))
        offset=kind===:vi ? nf : 0
        J[row,offset+s.free_state_map[(bus,"p")]]=1.
        J[row,offset+s.free_state_map[(bus,"n")]]=-1.
    end
    initial=consistency_pack(s,voltage)
    src_n=s.free_state_map[("s","n")]
    for anchored in (false,true)
        prior=anchored ? StatePrior([src_n,nf+src_n],[3.,2.],[.6,.8]) : nothing
        p=SEParameters(s,meters;prior)
        for solver in CONSISTENCY_SE_SOLVERS
            r=solver(s,p,initial)
            @test r.status == (anchored ? :converged_unique : :converged_underobserved)
            @test J*r.state ≈ means atol=1e-7
            @test .5sum(abs2,r.evaluation.residual) ≈ objective atol=1e-9
            @test derived_covariance(s,p,r.state,J) ≈ Diagonal(variances) atol=1e-10
            if anchored
                us=complex(means[1],means[2]);ub=complex(means[3],means[4]);sn=3+2im
                expected=Dict(("s","p")=>us+sn,("s","n")=>sn,
                    ("b","p")=>(us+2sn+ub)/2,("b","n")=>(us+2sn-ub)/2)
                @test r.state ≈ consistency_pack(s,expected) atol=1e-7
            else
                @test observability_diagnostics(s,p,r.state).unobservable_dimension==2
                @test_throws ArgumentError selected_state_covariance(s,p,r.state,[src_n])
                shifted=r.state+vcat(fill(17.,nf),fill(-9.,nf))
                @test evaluate_state_estimator(s,p,shifted).residual ≈ r.evaluation.residual atol=1e-10
                @test norm(evaluate_state_estimator(s,p,shifted).constraints)<1e-10
            end
        end
    end
end

if isdefined(Main,:OpenDSSDirect)
    @testset "SE three-winding transformer: independent OpenDSS port telemetry" begin
        OpenDSSDirect.dss("Clear")
        deck=joinpath(@__DIR__,"data","state_estimation","three_winding.dss")
        OpenDSSDirect.dss("Redirect \"$deck\"")
        @test OpenDSSDirect.Solution.Converged()
        volts=Dict(lowercase(node)=>ComplexF64(v) for (node,v) in
            zip(OpenDSSDirect.Circuit.AllNodeNames(),OpenDSSDirect.Circuit.AllBusVolts()))
        OpenDSSDirect.Circuit.SetActiveElement("Transformer.t")
        currents=ComplexF64.(OpenDSSDirect.CktElement.Currents())
        powers=1000 .* ComplexF64.(OpenDSSDirect.CktElement.Powers()) # kW/kvar -> W/var
        nc=OpenDSSDirect.CktElement.NumConductors()
        @test OpenDSSDirect.CktElement.NumTerminals()==3
        net=Dict{String,Any}(
            "bus"=>Dict(bus=>Dict("terminal_names"=>["a","n"],"perfectly_grounded_terminals"=>["n"]) for bus in ("h","l","t")),
            "voltage_source"=>Dict("s"=>Dict("bus"=>"h","terminal_map"=>["a"],"v_magnitude"=>[abs(volts["h.1"])],"v_angle"=>[angle(volts["h.1"])])),
            "transformer"=>Dict("n_winding"=>Dict("tx"=>Dict{String,Any}(
                "windings"=>[Dict{String,Any}("bus"=>bus,"terminal_map"=>["a","n"],"v_nom"=>v,"r_winding"=>r)
                    for (bus,v,r) in (("h",2400.,1.),("l",1200.,.25),("t",600.,.0625))],
                "x_sc"=>Dict("1_2"=>2.,"1_3"=>3.,"2_3"=>2.)))))
        meters=Any[]
        for (bus,S) in (("l",-10000-2000im),("t",-5000-1000im))
            for (kind,value,sigma) in ((:vr,real(volts["$bus.1"]),.1),(:vi,imag(volts["$bus.1"]),.1),(:pinj,real(S),20.),(:qinj,imag(S),20.))
                push!(meters,Measurement(;kind,bus,terminal="a",value,sigma))
            end
        end
        for port in 1:3
            row=(port-1)*nc+1
            for (kind,value,sigma) in ((:ire,real(currents[row]),.02),(:iim,imag(currents[row]),.02),
                                      (:pflow,real(powers[row]),20.),(:qflow,imag(powers[row]),20.))
                push!(meters,TransformerMeasurement(;kind,transformer="tx",variant="n_winding",port,terminal="a",value,sigma))
            end
        end
        s=compile_state_estimator(net,meters);p=SEParameters(s,meters)
        voltage=Dict((bus,"a")=>volts["$bus.1"] for bus in ("l","t"))
        truth=consistency_pack(s,voltage)
        e=evaluate_state_estimator(s,p,truth)
        @test norm(e.residual,Inf)<1e-4
        for solver in CONSISTENCY_SE_SOLVERS
            r=solver(s,p,initial_state_estimator(s,p);optimality_tolerance=1e-6)
            @test solve_status(r).publishable
            @test r.state ≈ truth atol=1e-4
            @test norm(r.evaluation.residual,Inf)<1e-4
        end
        for per_unit in (false,true)
            native=solve_state_estimation(net,[m for m in meters if m isa Measurement];per_unit)
            @test solve_status(native).publishable
            @test native.objective<1e-6
            for bus in ("l","t")
                @test complex(native.bus[bus]["a"]["vr"],native.bus[bus]["a"]["vi"]) ≈ volts["$bus.1"] atol=1e-4
            end
        end
    end
end

@testset "SE floating transformer: retain the turns-ratio current identity" begin
    net=Dict{String,Any}(
        "bus"=>Dict(bus=>Dict("terminal_names"=>["1","n"],"perfectly_grounded_terminals"=>["n"]) for bus in ("high","low")),
        "transformer"=>Dict("single_phase"=>Dict("tx"=>Dict{String,Any}(
            "bus_from"=>"high","bus_to"=>"low","terminal_map_from"=>["1","n"],"terminal_map_to"=>["1","n"],
            "v_nom_from"=>1000.,"v_nom_to"=>10.,"r_series_from"=>1.,"x_series_from"=>1.,"r_series_to"=>.0001))))
    currents=[ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection((bus,"1"),(bus,"n"))],[ComplexF64(current)])) for (bus,current) in (("high",-.01),("low",1.))]
    voltage=Dict(("high","1")=>1000+0im,("low","1")=>(1000-(2+im)*.01)/100)
    meters=[Measurement(;kind,bus=node[1],terminal=node[2],value,sigma=node[1]=="high" ? 1. : .01)
        for (node,v) in sort!(collect(voltage);by=first) for (kind,value) in ((:vr,real(v)),(:vi,imag(v)))]
    s=compile_state_estimator(net,meters;exact_devices=currents);p=SEParameters(s,meters;exact_devices=currents)
    truth=consistency_pack(s,voltage)
    rows,expansion=PowerOptLab._se_ipopt_constraint_basis(s,p)
    signature=PowerOptLab._se_ipopt_constraint_signature(s,p)
    @test length(rows)==2
    @test expansion*signature[rows,:] ≈ signature atol=1e-10
    hi=findfirst(==(s.node_index[("high","1")]),s.constraint_pattern)
    lo=findfirst(==(s.node_index[("low","1")]),s.constraint_pattern)
    @test expansion[2lo-1,:] ≈ -100expansion[2hi-1,:] atol=1e-10
    initial=truth+.001cos.(eachindex(truth))
    for solver in CONSISTENCY_SE_SOLVERS
        r=solver(s,p,initial)
        @test solve_status(r).publishable
        @test r.state ≈ truth atol=1e-7
        @test norm(r.evaluation.constraints)<1e-8
    end
end
