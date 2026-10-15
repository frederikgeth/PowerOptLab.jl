using Test, PowerOptLab, LinearAlgebra
import BMOPFTools

# Independent three-winding star circuit, not the package's ZB/Yprim stamp.
# Each core has a hidden common emf E; sum((u[k]/N[k]-E)/z[k]) == 0.
# All impedances below are referred to the high winding's coil base.
function tg_case(configurations; grounding=:floating, roll=-1, order=[1,2,3])
    buses=["h","l","t"]; nominal=[2400.,1200.,600.]
    z=ComplexF64[1+.8im,1+1.2im,1+2.2im]; turns=nominal./nominal[1]
    phases=["a","b","c"]
    terms=[c=="WYE" ? [phases;"n"] : copy(phases) for c in configurations]
    nodes=[(b,t) for (b,ts) in zip(buses,terms) for t in ts]
    other(k)=mod1(k+roll,3)
    function port_currents(v)
        current=Dict(node=>0.0im for node in nodes)
        for k in 1:3
            refs=[c=="WYE" ? "n" : phases[other(k)] for c in configurations]
            u=[v[(b,phases[k])]-v[(b,ref)] for (b,ref) in zip(buses,refs)]
            emf=sum(u./turns./z)/sum(1 ./z)
            coil=(u./turns.-emf)./z./turns
            for (b,ref,i) in zip(buses,refs,coil)
                current[(b,phases[k])]+=i;current[(b,ref)]-=i
            end
        end
        current
    end
    # Unbalanced current loads use real terminal-to-terminal return paths.
    load=Dict(node=>0.0im for node in nodes);connections=TerminalConnection[];amps=ComplexF64[]
    for w in 2:3, k in 1:3
        ref=configurations[w]=="WYE" ? "n" : phases[other(k)]
        i=([7-2im,4-1im,5-3im][k])*(w==2 ? 1. : .7)
        a=(buses[w],phases[k]);b=(buses[w],ref)
        push!(connections,TerminalConnection(a,b));push!(amps,i)
        load[a]+=i;load[b]-=i
    end
    electrodes=Dict{Tuple{String,String},ComplexF64}()
    if grounding===:impedance
        # A phase-earth path and a neutral-earth electrode permit a nonzero
        # neutral displacement and electrode current; neither is a fake anchor.
        for w in 2:3
            electrodes[(buses[w],"a")]=.003-.001im
            configurations[w]=="WYE" && (electrodes[(buses[w],"n")]=.12-.04im)
        end
    end
    fixed=Dict(("h",t)=>2400cis(-2pi*(k-1)/3) for (k,t) in enumerate(phases))
    configurations[1]=="WYE" && (fixed[("h","n")]=0im)
    if grounding===:solid
        for w in 2:3
            configurations[w]=="WYE" && (fixed[(buses[w],"n")]=0im)
        end
    end
    free=[node for node in nodes if !haskey(fixed,node)]
    function balance(x)
        v=merge(Dict(node=>ComplexF64(x[k]) for (k,node) in enumerate(free)),fixed)
        i=port_currents(v)
        [i[node]+load[node]+get(electrodes,node,0im)*v[node] for node in free]
    end
    b=balance(zeros(ComplexF64,length(free)))
    A=hcat([begin
        v=Dict(node=>ComplexF64(node==column) for node in nodes)
        i=port_currents(v)
        [i[node]+get(electrodes,node,0im)*v[node] for node in free]
    end for column in free]...)
    # Minimum norm chooses one gauge for the oracle only. Tests observe voltage
    # differences when a winding has no earth path, never that arbitrary gauge.
    x=-pinv(A)*b;voltage=merge(Dict(node=>x[k] for (k,node) in enumerate(free)),fixed)
    net=Dict{String,Any}("bus"=>Dict(b=>Dict{String,Any}("terminal_names"=>ts,
        "perfectly_grounded_terminals"=>[t for t in ts if haskey(fixed,(b,t)) && t=="n"])
        for (b,ts) in zip(buses,terms)),
        "voltage_source"=>Dict("s"=>Dict("bus"=>"h","terminal_map"=>phases,
            "v_magnitude"=>fill(2400.,3),"v_angle"=>[0.,-2pi/3,-4pi/3])),
        "shunt"=>Dict("e$k"=>Dict("bus"=>node[1],"terminal_map"=>[node[2]],
            "G_1_1"=>real(y),"B_1_1"=>imag(y)) for (k,(node,y)) in enumerate(sort!(collect(electrodes);by=first))))
    windings=[Dict{String,Any}("bus"=>buses[k],"terminal_map"=>terms[k],"configuration"=>configurations[k],
        "delta_roll"=>roll,"v_nom"=>nominal[k],"r_winding"=>real(z[k])*turns[k]^2) for k in order]
    # Re-refer physical pair reactances when another winding is first.
    xsc=Dict("$(i)_$(j)"=>(imag(z[order[i]])+imag(z[order[j]]))*turns[order[1]]^2 for i in 1:3 for j in i+1:3)
    net["transformer"]=Dict("n_winding"=>Dict("tx"=>Dict{String,Any}("windings"=>windings,"x_sc"=>xsc)))
    exact=[ExactDeviceEquation(ConstantCurrentDevice(connections,amps))]
    (;net,voltage,port_currents,electrodes,exact,balance,x,free,terms,configurations,order,turns,z)
end

function tg_meters(case)
    meters=Any[];phases=["a","b","c"];buses=["h","l","t"]
    current=case.port_currents(case.voltage)
    for w in 2:3
        bus=buses[w]
        # Delta line-line readings span its two differential complex states.
        for (terminal,reference) in (case.configurations[w]=="WYE" ? [(t,"n") for t in phases] : [("a","b"),("b","c")])
            v=case.voltage[(bus,terminal)]-case.voltage[(bus,reference)]
            for (kind,value) in ((:vr,real(v)),(:vi,imag(v)))
                push!(meters,Measurement(;kind,bus,terminal,reference,value,sigma=.1))
            end
        end
    end
    for (port,w) in enumerate(case.order),terminal in case.terms[w]
        bus=buses[w];i=current[(bus,terminal)]
        reference=case.configurations[w]=="WYE" ? "n" : phases[mod1(findfirst(==(terminal),phases)+1,3)]
        # Neutral power is metered only with a finite earth path; phase power
        # uses a local voltage reference that is invariant to a common shift.
        v=case.voltage[(bus,terminal)]-case.voltage[(bus,reference)]
        kinds=[(:ire,real(i),.02),(:iim,imag(i),.02)]
        if terminal!="n"
            append!(kinds,[(:pflow,real(v*conj(i)),20.),(:qflow,imag(v*conj(i)),20.)])
        elseif haskey(case.electrodes,(bus,"n"))
            reference=nothing;v=case.voltage[(bus,"n")]
            append!(kinds,[(:pflow,real(v*conj(i)),20.),(:qflow,imag(v*conj(i)),20.)])
        end
        for (kind,value,sigma) in kinds
            push!(meters,TransformerMeasurement(;kind,transformer="tx",variant="n_winding",port,terminal,reference,value,sigma))
        end
    end
    for (id,e) in case.net["shunt"]
        bus=e["bus"];terminal=only(e["terminal_map"]);v=case.voltage[(bus,terminal)]
        i=case.electrodes[(bus,terminal)]*v
        for (kind,value) in ((:vr,real(v)),(:vi,imag(v)),(:vmag,abs(v)))
            push!(meters,Measurement(;kind,bus,terminal,reference=nothing,value,sigma=.1))
        end
        for (kind,value,sigma) in ((:ire,real(i),.02),(:iim,imag(i),.02),(:pflow,real(v*conj(i)),20.),(:qflow,imag(v*conj(i)),20.))
            push!(meters,ShuntMeasurement(;kind,shunt=id,terminal,value,sigma))
        end
    end
    meters
end

@testset "SE independent mixed-winding and grounding circuit" begin
    for configurations in (["WYE","DELTA","WYE"],["DELTA","WYE","DELTA"],["WYE","WYE","DELTA"]),
        grounding in (:floating,:solid,:impedance), roll in (-1,1)
        @testset "$(join(configurations,'/')) $grounding roll=$roll" begin
            for order in ([1,2,3],[3,1,2])
                case=tg_case(configurations;grounding,roll,order)
                @test norm(case.balance(case.x),Inf)<1e-8
                # The independently specified resistive leakage is passive.
                i=case.port_currents(case.voltage)
                @test real(sum(case.voltage[n]*conj(i[n]) for n in keys(i)))>=0
                meters=tg_meters(case)
                s=compile_state_estimator(case.net,meters;exact_devices=case.exact)
                p=SEParameters(s,meters;exact_devices=case.exact)
                nf=length(s.free_state_map);truth=zeros(2nf)
                for (node,k) in s.free_state_map
                    truth[k]=real(case.voltage[node]);truth[nf+k]=imag(case.voltage[node])
                end
                e=evaluate_state_estimator(s,p,truth)
                @test norm(e.constraints,Inf)<1e-8
                @test norm(e.residual,Inf)<1e-7
                gauges=2count(w->configurations[w]=="DELTA" || grounding===:floating,2:3)
                grounding===:impedance && (gauges=0)
                @test observability_diagnostics(s,p,truth).unobservable_dimension==gauges
                @test size(unobservable_directions(s,p,truth),2)==gauges
                H=residual_jacobian(s,p,truth);C=constraint_jacobian(s,p,truth)
                legacy=PowerOptLab._se_observability_metrics(H,C)
                @test legacy.observable==(gauges==0)
                # Exact affine device laws already determine every differential
                # voltage. Its covariance is zero, even with an unresolved gauge.
                bus="l";ref=configurations[2]=="WYE" ? "n" : "b"
                J=zeros(1,2nf);J[1,s.free_state_map[(bus,"a")]]=1.
                haskey(s.free_state_map,(bus,ref)) && (J[1,s.free_state_map[(bus,ref)]]-=1.)
                @test norm(derived_covariance(s,p,truth,J))<1e-20
                if gauges>0
                    @test observability_diagnostics(s,p,truth).condition_number==Inf
                    w=findfirst(w->configurations[w]=="DELTA" || grounding===:floating,2:3)+1
                    @test_throws ArgumentError selected_state_covariance(s,p,truth,[s.free_state_map[(["h","l","t"][w],"a")]])
                end
                for solver in (solve_compiled_state_estimator,solve_sparse_state_estimator,solve_ipopt_state_estimator)
                    r=solver(s,p,truth+.01cos.(eachindex(truth));optimality_tolerance=1e-6)
                    @test solve_status(r).publishable
                    @test r.status==(gauges==0 ? :converged_unique : :converged_underobserved)
                    @test norm(r.evaluation.constraints,Inf)<1e-7
                    @test norm(r.evaluation.residual,Inf)<1e-5
                    gauges==0 && @test r.state≈truth atol=1e-5
                end
            end
        end
    end
end

@testset "SE floating winding: weak real reference information survives projection" begin
    case=tg_case(["WYE","DELTA","WYE"])
    meters=tg_meters(case);s=compile_state_estimator(case.net,meters;exact_devices=case.exact)
    nf=length(s.free_state_map);truth=zeros(2nf)
    for (node,k) in s.free_state_map
        truth[k]=real(case.voltage[node]);truth[nf+k]=imag(case.voltage[node])
    end
    indices=[s.free_state_map[(bus,"a")]+offset for bus in ("l","t") for offset in (0,nf)]
    prior=StatePrior(indices,truth[indices],fill(1e7,4))
    p=SEParameters(s,meters;exact_devices=case.exact,prior)
    @test observability_diagnostics(s,p,truth).unobservable_dimension==0
    @test observability_diagnostics(s,p,truth;include_prior=false).unobservable_dimension==4
    @test selected_state_covariance(s,p,truth,indices)≈1e14Matrix{Float64}(I,4,4) rtol=1e-5
end

@testset "SE mixed-winding preflight preserves the requested circuit" begin
    case=tg_case(["WYE","DELTA","WYE"])
    for mutate in (
        w->(w[2]["delta_roll"]=0),
        w->(w[2]["delta_roll"]=2),
        w->push!(w[2]["terminal_map"],"n"),
        w->(w[3]["r_winding"]=-.1),
        w->(w[3]["r_winding"]=Inf),
        w->(w[3]["connection"]="DELTA"),
        w->(w[3]["tap"]=1.03),
        w->(w[3]["tap_ratio"]=1.03),
        w->(w[3]["r_neutral"]=10.),
        w->(w[3]["x_neutral"]=2.))
        net=deepcopy(case.net);mutate(net["transformer"]["n_winding"]["tx"]["windings"])
        # A declared but unused bus neutral is valid; it is specifically the
        # DELTA winding map that must not pretend to connect that conductor.
        "n" in net["transformer"]["n_winding"]["tx"]["windings"][2]["terminal_map"] &&
            push!(net["bus"]["l"]["terminal_names"],"n")
        report=state_estimator_preflight(net)
        @test !report.supported
        @test any(f->f.code in (:invalid_connection,:invalid_impedance,:unsupported_winding_field),report.findings)
        @test_throws SEUnsupportedNetwork compile_state_estimator(net)
    end
    for field in ("r_neutral_from","x_neutral_from","r_neutral_to","x_neutral_to")
        net=deepcopy(case.net);net["transformer"]["n_winding"]["tx"][field]=10.
        @test any(f->f.code==:unsupported_winding_field,state_estimator_preflight(net).findings)
        @test_throws SEUnsupportedNetwork compile_state_estimator(net)
    end
end

@testset "SE winding return identities retain even weak physical electrodes" begin
    case=tg_case(["WYE","DELTA","WYE"])
    connections=[TerminalConnection((bus,terminal),reference) for (bus,ref) in (("l",nothing),("t","n"))
        for terminal in ("a","b","c") for reference in (ref===nothing ? nothing : (bus,ref),)]
    for g in (0.,1e-12,.1)
        net=deepcopy(case.net)
        g>0 && (net["shunt"]["weak"]=Dict("bus"=>"t","terminal_map"=>["n"],"G_1_1"=>g))
        s=compile_state_estimator(net;injection_connections=connections)
        @test size(s.return_current_pattern,1)==(g==0 ? 0 : 1)
        if g>0
            p=SEParameters(s)
            @test length(PowerOptLab._se_ipopt_constraint_rows(s,p))>=1
        end
    end
end

@testset "SE mixed-winding native and shared statistical models" begin
    for configurations in (["WYE","DELTA","WYE"],["DELTA","WYE","DELTA"],["WYE","WYE","DELTA"]), grounding in (:floating,:solid,:impedance)
        case=tg_case(configurations;grounding)
        meters=Measurement[];current=case.port_currents(case.voltage)
        for w in 2:3
            bus=["h","l","t"][w]
            for terminal in case.terms[w]
                v=case.voltage[(bus,terminal)]
                # Absolute voltage readings are explicit stochastic reference
                # information for this comparison; they are not perfect bonds.
                for (kind,value) in ((:vr,real(v)),(:vi,imag(v)))
                    push!(meters,Measurement(;kind,bus,terminal,reference=nothing,value,sigma=.1))
                end
                terminal=="n" && continue
                reference=configurations[w]=="WYE" ? "n" : nothing
                u=v-(reference===nothing ? 0im : case.voltage[(bus,reference)])
                i=current[(bus,terminal)]+get(case.electrodes,(bus,terminal),0im)*v
                power=u*conj(i)
                for (kind,value) in ((:pinj,real(power)),(:qinj,imag(power)))
                    push!(meters,Measurement(;kind,bus,terminal,reference,value,sigma=20.))
                end
            end
        end
        s=compile_state_estimator(case.net,meters);p=SEParameters(s,meters)
        nf=length(s.free_state_map);truth=zeros(2nf)
        for (node,k) in s.free_state_map
            truth[k]=real(case.voltage[node]);truth[nf+k]=imag(case.voltage[node])
        end
        @test norm(evaluate_state_estimator(s,p,truth).constraints,Inf)<1e-8
        for solver in (solve_compiled_state_estimator,solve_sparse_state_estimator,solve_ipopt_state_estimator)
            r=solver(s,p,initial_state_estimator(s,p);optimality_tolerance=1e-6)
            @test r.status==:converged_unique
            @test r.state≈truth atol=1e-4
        end
        for per_unit in (false,true)
            r=solve_state_estimation(case.net,meters;per_unit)
            @test solve_status(r).publishable
            @test r.objective<1e-6
            for (node,k) in s.free_state_map
                v=complex(r.bus[node[1]][node[2]]["vr"],r.bus[node[1]][node[2]]["vi"])
                @test v≈case.voltage[node] atol=1e-4
            end
        end
    end
end

if isdefined(Main,:OpenDSSDirect)
    @testset "SE mixed-winding external OpenDSS terminal equations" begin
        for cs in (["WYE","DELTA","WYE"],["DELTA","WYE","DELTA"],["WYE","WYE","DELTA"]), roll in (-1,1)
            case=tg_case(cs;grounding=:impedance,roll)
            OpenDSSDirect.dss("Clear")
            OpenDSSDirect.dss("Set DefaultBaseFreq=50")
            OpenDSSDirect.dss("New Circuit.mixed phases=3 bus1=h.1.2.3.0 basekv=4.156921938165306 model=ideal")
            OpenDSSDirect.dss("New Transformer.t phases=3 windings=3")
            # Reverse the conductor sequence for the opposite coil rotation.
            # This avoids depending on OpenDSS's leadlag handling of winding 3.
            permutation=roll==1 ? [1,2,3,4] : [1,3,2,4]
            for k in 1:3
                bus=["h","l","t"][k]
                kv=2400case.turns[k]*(cs[k]=="WYE" ? sqrt(3.) : 1.)/1000
                base=3*(2400case.turns[k])^2/100000 # per-coil ohms at 100 kVA
                percent_r=100real(case.z[k])*case.turns[k]^2/base
                OpenDSSDirect.dss("~ wdg=$k bus=$bus.$(join(permutation,'.')) conn=$(lowercase(cs[k])) kv=$kv kva=100 %r=$percent_r rneut=-1")
            end
            base=3*2400^2/100000
            leadlag=cs[1]=="DELTA" ? "lead" : "lag"
            OpenDSSDirect.dss("~ XHL=$(100*2/base) XHT=$(100*3/base) XLT=$(100*3.4/base) %noloadloss=0 %imag=0 ppm_antifloat=0 leadlag=$leadlag")
            OpenDSSDirect.dss("CalcVoltageBases")
            OpenDSSDirect.Circuit.SetActiveElement("Transformer.t")
            @test OpenDSSDirect.CktElement.NumTerminals()==3
            @test OpenDSSDirect.CktElement.NodeOrder()==repeat(permutation,3)
            # Evaluate external terminal equations at arbitrary operating
            # voltages, including displaced floating neutrals. No DSS gauge or
            # load-flow convergence is needed for this component comparison.
            nodes=[(["h","l","t"][w],["a","b","c","n"][j]) for w in 1:3 for j in permutation]
            v=ComplexF64[get(case.voltage,node,0im) for node in nodes]
            currents=OpenDSSDirect.CktElement.YPrim()*v
            expected=case.port_currents(case.voltage)
            @test currents≈[get(expected,node,0im) for node in nodes] atol=1e-8
            # Compare every represented port conductor with external telemetry.
            meters=Any[]
            for port in 1:3,terminal in case.terms[port]
                row=findfirst(==((["h","l","t"][port],terminal)),nodes)
                for (kind,value) in ((:ire,real(currents[row])),(:iim,imag(currents[row])))
                    push!(meters,TransformerMeasurement(;kind,variant="n_winding",transformer="tx",port,terminal,value,sigma=.02))
                end
            end
            s=compile_state_estimator(case.net,meters;exact_devices=case.exact)
            p=SEParameters(s,meters;exact_devices=case.exact);nf=length(s.free_state_map)
            truth=zeros(2nf)
            for (node,k) in s.free_state_map
                truth[k]=real(case.voltage[node]);truth[nf+k]=imag(case.voltage[node])
            end
            @test norm(evaluate_state_estimator(s,p,truth).residual,Inf)<1e-6
        end
    end
end
