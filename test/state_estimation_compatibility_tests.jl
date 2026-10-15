using Test, PowerOptLab, LinearAlgebra
using BMOPFTools: parse_bmopf
import BMOPFTools
isdefined(@__MODULE__, :_se_net_1ph_xfmr) || include("state_estimation_feeder_fixtures.jl")

const COMPAT_SE_SOLVERS = (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator)

function compat_ground_net(mode)
    grounded = mode === :solid ? ",\"perfectly_grounded_terminals\":[\"n\"]" : ""
    source_ground = mode === :multi_impedance ? "[]" : "[\"n\"]"
    net = parse_bmopf("""
    {"bus":{"src":{"terminal_names":["1","n"],"perfectly_grounded_terminals":$source_ground},
            "b":{"terminal_names":["1","n"]$grounded}},
     "voltage_source":{"s":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.],"v_angle":[0.]}},
     "linecode":{"lc":{"R_series_1_1":0.5,"R_series_2_2":0.5}},
     "line":{"l":{"bus_from":"src","bus_to":"b","terminal_map_from":["1","n"],
                    "terminal_map_to":["1","n"],"linecode":"lc","length":1.}}}
    """; from_string=true)
    if mode in (:impedance,:multi_impedance)
        net["shunt"] = Dict("earth" => Dict{String,Any}("bus"=>"b", "terminal_map"=>["n"], "G_1_1"=>.2,"B_1_1"=>-.05))
        if mode === :multi_impedance
            net["shunt"]["source_earth"] = Dict{String,Any}("bus"=>"src","terminal_map"=>["n"],"G_1_1"=>.1)
        end
    end
    net
end

function compat_state(s, voltage)
    x = zeros(2length(s.free_state_map)); nf=length(x)÷2
    for (node,k) in s.free_state_map
        x[k]=real(voltage[node]); x[nf+k]=imag(voltage[node])
    end
    x
end
function compat_fd(f,x)
    h=1e-5; eye=Matrix{Float64}(I,length(x),length(x))
    hcat([(f(x+h*eye[:,k])-f(x-h*eye[:,k]))/(2h) for k in eachindex(x)]...)
end

@testset "SE shared backends: ungrounded reference and gauge" begin
    net=compat_ground_net(:floating)
    empty!(net["voltage_source"])
    net["bus"]["src"]["perfectly_grounded_terminals"] = String[]
    v=Dict(("src","1")=>233+2im,("src","n")=>3+2im,
           ("b","1")=>228+3im,("b","n")=>8+im)
    cs=[TerminalConnection((bus,"1"),(bus,"n")) for bus in ("src","b")]
    absolute=[Measurement(kind=kind,bus=bus,terminal=t,reference=nothing,value=value,sigma=.2)
              for ((bus,t),voltage) in sort!(collect(v);by=first)
              for (kind,value) in ((:vr,real(voltage)),(:vi,imag(voltage)))]
    relative=[Measurement(kind=kind,bus=bus,reference="n",value=value,sigma=.2)
              for bus in ("src","b") for (kind,value) in
              ((:vr,real(v[(bus,"1")]-v[(bus,"n")])),(:vi,imag(v[(bus,"1")]-v[(bus,"n")])))]
    for (meters,expected) in ((absolute,:converged_unique),(relative,:converged_underobserved))
        s=compile_state_estimator(net,meters;injection_connections=cs);p=SEParameters(s,meters)
        truth=compat_state(s,v);x0=truth+.1cos.(eachindex(truth))
        @test norm(evaluate_state_estimator(s,p,truth).constraints)<1e-10
        @test any(f->f.code==:no_fixed_source,s.preflight.findings)
        for solver in COMPAT_SE_SOLVERS
            r=solver(s,p,x0)
            @test r.status == expected
            @test solve_status(r).publishable
            @test norm(r.evaluation.residual)<1e-5
            @test norm(r.evaluation.constraints)<1e-7
            expected === :converged_unique && @test r.state ≈ truth atol=1e-5
        end
    end
end

@testset "SE shared backends: grounding and neutral telemetry" begin
    for mode in (:floating, :impedance, :solid, :multi_impedance)
        @testset "$mode" begin
            net=compat_ground_net(mode)
            g=mode in (:impedance,:multi_impedance) ? .2-.05im : 0.
            ip=-10+2im; vp=230+.5ip
            vn=mode === :solid ? 0im : mode === :multi_impedance ? -ip/(g+2*.1/(2+.1)) : -ip/(2+g)
            source_neutral=mode === :multi_impedance ? 2vn/(2+.1) : 0im
            power=(vp-vn)*conj(ip)
            meters=Any[Measurement(kind=:pinj,bus="b",value=real(power),sigma=10.),
                       Measurement(kind=:qinj,bus="b",value=imag(power),sigma=10.),
                       Measurement(kind=:vr,bus="b",terminal="n",reference=nothing,value=real(vn),sigma=.2),
                       Measurement(kind=:vi,bus="b",terminal="n",reference=nothing,value=imag(vn),sigma=.2),
                       Measurement(kind=:vmag,bus="b",terminal="n",reference=nothing,value=abs(vn),sigma=.2)]
            # Neutral line current and finite-electrode current are separate
            # physical readings. An ideal-ground reaction has no Ybus row.
            line_neutral=2(vn-source_neutral)
            for (kind,value) in ((:ire,real(line_neutral)),(:iim,imag(line_neutral)),(:imag,abs(line_neutral)))
                push!(meters,BranchMeasurement(;kind,line="l",side=:to,terminal="n",value,sigma=.2))
            end
            # Phase-neutral branch power explicitly selects its voltage reference.
            push!(meters,BranchMeasurement(kind=:pflow,line="l",side=:to,terminal="1",reference="n",value=real(power),sigma=10.))
            if mode !== :solid
                net_current=-ip
                sn=vn*conj(net_current)
                for (kind,value) in ((:ire,real(net_current)),(:iim,imag(net_current)),(:imag,abs(net_current)),
                                     (:pinj,real(sn)),(:qinj,imag(sn)))
                    push!(meters,NetInjectionMeasurement(;kind,bus="b",terminal="n",value,sigma=.2))
                end
            else
                for kind in (:pinj,:qinj)
                    push!(meters,NetInjectionMeasurement(;kind,bus="b",terminal="n",value=0.,sigma=.2))
                end
                @test !state_estimator_preflight(net,[NetInjectionMeasurement(kind=:ire,bus="b",terminal="n",value=0.,sigma=1.)]).supported
            end
            if mode in (:impedance,:multi_impedance)
                ie=g*vn; se=vn*conj(ie)
                for (kind,value) in ((:ire,real(ie)),(:iim,imag(ie)),(:imag,abs(ie)),(:pflow,real(se)),(:qflow,imag(se)))
                    push!(meters,ShuntMeasurement(;kind,shunt="earth",terminal="n",value,sigma=.2))
                end
                @test line_neutral+ie ≈ -ip # independently derived current split
            end
            if mode === :multi_impedance
                for (kind,value) in ((:vr,real(source_neutral)),(:vi,imag(source_neutral)))
                    push!(meters,Measurement(;kind,bus="src",terminal="n",reference=nothing,value,sigma=.2))
                end
                source_electrode=.1source_neutral
                for (kind,value) in ((:ire,real(source_electrode)),(:iim,imag(source_electrode)))
                    push!(meters,ShuntMeasurement(;kind,shunt="source_earth",terminal="n",value,sigma=.2))
                end
                @test source_electrode+g*vn ≈ -ip
            end
            @test state_estimator_preflight(net,meters).supported
            zi=mode === :multi_impedance ? [("src","n")] : Tuple{String,String}[]
            s=compile_state_estimator(net,meters;zero_injection=zi)
            # Load-free starts have zero neutral V/I. Explicit smoothing below
            # instrument resolution makes nonconstant magnitude rows differentiable.
            smoothing=mode === :solid ? 0. : 1e-6
            p=SEParameters(s,meters;magnitude_epsilon=smoothing,current_epsilon=smoothing)
            truth=compat_state(s,Dict(("b","1")=>vp,("b","n")=>vn,("src","n")=>source_neutral))
            e=evaluate_state_estimator(s,p,truth)
            @test norm(e.residual,Inf)<1e-8
            @test norm(e.constraints,Inf)<1e-8
            @test Matrix(residual_jacobian(s,p,truth)) ≈ compat_fd(x->evaluate_state_estimator(s,p,x).residual,truth) atol=1e-6 rtol=1e-6
            for solver in COMPAT_SE_SOLVERS
                r=solver(s,p,initial_state_estimator(s,p);optimality_tolerance=1e-6)
                @test solve_status(r).publishable
                @test r.state ≈ truth atol=1e-4
                @test norm(r.evaluation.constraints,Inf)<1e-7
            end
            # The original native-engine WLS path must agree on these three
            # grounding regimes too, including constant zero neutral magnitude.
            customer=[m for m in meters if m isa Measurement]
            for per_unit in (false,true)
                if mode === :multi_impedance
                    @test_throws ArgumentError solve_state_estimation(net,customer;per_unit)
                    continue
                end
                legacy=solve_state_estimation(net,customer;per_unit)
                @test solve_status(legacy).publishable
                @test legacy.objective < 1e-6
                @test legacy.bus["b"]["n"]["vr"] ≈ real(vn) atol=1e-5
            end
        end
    end
end

@testset "SE functional measurements: derivatives, identity and stationarity" begin
    net=compat_ground_net(:floating)
    nodes=[("b","1"),("b","n"),("src","n")]
    f=z->((z[1]-z[2])^2+(z[4]-z[5])^2)/230
    custom=FunctionalMeasurement(name=:scaled_square_voltage,nodes=nodes,prediction=f,value=(220^2+2^2)/230,sigma=.2)
    meters=Any[Measurement(kind=:vr,bus="b",value=220.,sigma=.2),
               Measurement(kind=:vi,bus="b",value=2.,sigma=.2),custom]
    cs=[TerminalConnection(("b","1"),("b","n"))]
    s=compile_state_estimator(net,meters;injection_connections=cs);p=SEParameters(s,meters)
    x=compat_state(s,Dict(("b","1")=>225+im,("b","n")=>5-im))
    @test norm(evaluate_state_estimator(s,p,x).residual)<1e-8
    @test Matrix(residual_jacobian(s,p,x)) ≈ compat_fd(y->evaluate_state_estimator(s,p,y).residual,x) atol=1e-6
    for solver in COMPAT_SE_SOLVERS
        r=solver(s,p,initial_state_estimator(s,p))
        @test solve_status(r).publishable
        @test r.state ≈ x atol=1e-5
    end
    changed=copy(meters)
    changed[3]=FunctionalMeasurement(name=custom.name,nodes=nodes,prediction=f,value=custom.value+.1,sigma=.4)
    @test SEParameters(s,changed).measurement_values[end] == custom.value+.1
    changed[3]=FunctionalMeasurement(name=custom.name,nodes=reverse(nodes),prediction=f,value=custom.value,sigma=.2)
    @test_throws ArgumentError SEParameters(s,changed)
    bad=FunctionalMeasurement(name=:bad,nodes=nodes,prediction=z->NaN,value=0.,sigma=1.)
    sb=compile_state_estimator(net,[bad]);pb=SEParameters(sb,[bad])
    @test solve_ipopt_state_estimator(sb,pb,initial_state_estimator(sb,pb)).status == :invalid_initial_domain
    # A user-supplied gradient allows callbacks which accept Float64 only.
    manual(z::Vector{Float64})=z[1]^2/230
    gradient=z->[2z[1]/230,0.]
    mm=FunctionalMeasurement(name=:manual,nodes=[("b","1")],prediction=manual,gradient=gradient,value=225^2/230,sigma=1.)
    sm=compile_state_estimator(net,[mm]);pm=SEParameters(sm,[mm])
    @test Matrix(residual_jacobian(sm,pm,x)) ≈ compat_fd(y->evaluate_state_estimator(sm,pm,y).residual,x) atol=1e-6
    # Zero derivative of a nonlinear residual is not sufficient for publication.
    square=FunctionalMeasurement(name=:square,nodes=[("b","1")],prediction=z->sum(abs2,z),value=1.,sigma=1.)
    ss=compile_state_estimator(compat_ground_net(:solid),[square]);ps=SEParameters(ss,[square])
    for solver in COMPAT_SE_SOLVERS
        r=solver(ss,ps,[0.,0.])
        @test r.status == :stationary_not_minimum
        @test !solve_status(r).publishable
    end
    series=solve_time_series_state_estimator(s,[p,deepcopy(p)],initial_state_estimator(s,p);solver=:ipopt,previous_state_sigma=1.)
    @test solve_status(series).publishable
    @test series.completed_snapshots == 2
end

function compat_nwinding_net()
    Dict{String,Any}(
        "bus"=>Dict("h"=>Dict("terminal_names"=>["a","n"],"perfectly_grounded_terminals"=>["n"]),
                    "l"=>Dict("terminal_names"=>["a","n"]),"t"=>Dict("terminal_names"=>["a","n"])),
        "voltage_source"=>Dict("s"=>Dict("bus"=>"h","terminal_map"=>["a"],"v_magnitude"=>[2400.],"v_angle"=>[0.])),
        "transformer"=>Dict("n_winding"=>Dict("nw"=>Dict{String,Any}(
            "windings"=>[Dict{String,Any}("bus"=>"h","terminal_map"=>["a","n"],"v_nom"=>2400.,"r_winding"=>1.),
                         Dict{String,Any}("bus"=>"l","terminal_map"=>["a","n"],"v_nom"=>1200.,"r_winding"=>.25),
                         Dict{String,Any}("bus"=>"t","terminal_map"=>["a","n"],"v_nom"=>600.,"r_winding"=>.0625)],
            "x_sc"=>Dict("1_2"=>2.,"1_3"=>3.,"2_3"=>2.)))) )
end

function compat_native_case(net)
    net=deepcopy(net)
    for group in ("load","generator","ibr","dc_source")
        pop!(net,group,nothing)
    end
    limits=Set(vcat(collect(PowerOptLab._SE_BUS_LIMIT_FIELDS),collect(PowerOptLab._SE_BRANCH_LIMIT_FIELDS)))
    function strip!(item)
        if item isa AbstractDict
            for key in collect(keys(item))
                key in limits ? delete!(item,key) : strip!(item[key])
            end
        elseif item isa AbstractVector
            foreach(strip!,item)
        end
    end
    strip!(net)
    # Test the native transformer implementation under a common explicit
    # grounding regime. Other tests above cover floating/electrode neutrals.
    for bus in values(net["bus"])
        "n" in bus["terminal_names"] && (bus["perfectly_grounded_terminals"]=["n"])
    end
    s=compile_state_estimator(net);p=SEParameters(s)
    nominal=evaluate_state_estimator(s,p,initial_state_estimator(s,p)).voltage
    sources=Set(vs["bus"] for vs in values(net["voltage_source"]))
    currents=Dict{Tuple{String,String},ComplexF64}()
    total=zeros(ComplexF64,length(s.nodes))
    for (bus,b) in net["bus"], t in PowerOptLab._phase_terminals(b,"n")
        bus in sources && continue
        node=(bus,t);i=s.node_index[node]
        injection=-.1cis(angle(nominal[i]))
        currents[node]=injection;total[i]+=injection
    end
    free=sort!(collect(s.free_state_map);by=last)
    indices=[s.node_index[node] for (node,_) in free]
    fixed=copy(nominal);fixed[indices].=0
    voltage=copy(fixed)
    voltage[indices]=qr(s.passive_pattern[indices,indices]) \ (total[indices]-(s.passive_pattern*fixed)[indices])
    @test norm((s.passive_pattern*voltage-total)[indices],Inf)<1e-7
    meters=Measurement[]
    fixed_indices=Set(s.node_index[n] for n in keys(s.reference_map))
    for (node,current) in sort!(collect(currents);by=first)
        v=voltage[s.node_index[node]];power=v*conj(current)
        for (kind,value,sigma) in ((:vr,real(v),.1),(:vi,imag(v),.1),(:pinj,real(power),10.),(:qinj,imag(power),10.))
            kind in (:pinj,:qinj) && s.node_index[node] in fixed_indices && continue
            push!(meters,Measurement(;kind,bus=node[1],terminal=node[2],reference=nothing,value,sigma))
        end
    end
    net,meters,voltage,s
end

@testset "SE transformer telemetry: all finite exporter variants" begin
    builds=(_se_net_1ph_xfmr,_se_net_yd_xfmr,_se_net_dy_xfmr,_se_net_autotransformer,
            _se_net_open_delta_reg,_se_net_center_tap,compat_nwinding_net)
    for build in builds
        @testset "$(nameof(build))" begin
            net=build();seed=compile_state_estimator(net);p0=SEParameters(seed)
            x0=initial_state_estimator(seed,p0)
            truth=x0+[.02sin(k) for k in eachindex(x0)]
            v=evaluate_state_estimator(seed,p0,truth).voltage
            meters=Any[]
            for node in keys(seed.free_state_map), (kind,value) in ((:vr,real(v[seed.node_index[node]])),(:vi,imag(v[seed.node_index[node]])))
                push!(meters,Measurement(;kind,bus=node[1],terminal=node[2],reference=nothing,value,sigma=.2))
            end
            for (variant,group) in net["transformer"], (id,el) in group
                ports=variant=="n_winding" ? el["windings"] :
                    [Dict("bus"=>el["bus_from"],"terminal_map"=>el["terminal_map_from"]),
                     Dict("bus"=>el["bus_to"],"terminal_map"=>el["terminal_map_to"])]
                nodes,Y=variant=="n_winding" ? BMOPFTools.nwinding_yprim(el) : BMOPFTools.transformer_yprim(el,variant)
                volts=[i==0 ? 0im : v[i] for i in [seed.node_index[node] for node in nodes]]
                currents=Y*volts
                for (port,w) in enumerate(ports), terminal in w["terminal_map"]
                    if variant == "open_delta_regulator" && (w["bus"],terminal) in PowerOptLab._se_open_delta_bond(el)
                        rejected=TransformerMeasurement(;kind=:ire,transformer=id,variant,port,terminal,value=0.,sigma=1.)
                        @test !state_estimator_preflight(net,[rejected]).supported
                        continue
                    end
                    row=findfirst(==((w["bus"],terminal)),nodes)
                    row===nothing && continue
                    c=currents[row]; power=volts[row]*conj(c)
                    for (kind,value) in ((:ire,real(c)),(:iim,imag(c)),(:pflow,real(power)),(:qflow,imag(power)))
                        push!(meters,TransformerMeasurement(;kind,transformer=id,variant,port,terminal,value,sigma=10.))
                    end
                end
            end
            @test state_estimator_preflight(net,meters).supported
            s=compile_state_estimator(net,meters);p=SEParameters(s,meters)
            @test norm(evaluate_state_estimator(s,p,truth).residual,Inf)<1e-7
            @test Matrix(residual_jacobian(s,p,truth)) ≈ compat_fd(y->evaluate_state_estimator(s,p,y).residual,truth) atol=2e-5 rtol=2e-5
            for solver in COMPAT_SE_SOLVERS
                r=solver(s,p,x0;optimality_tolerance=1e-5,max_iterations=500)
                solve_status(r).publishable || @info "Transformer telemetry solve failed" fixture=nameof(build) solver status=r.status diagnostics=solve_diagnostics(r)
                @test solve_status(r).publishable
                @test norm(r.state-truth,Inf)<1e-4
                @test observability_diagnostics(s,p,r.state).unobservable_dimension==0
            end
            first_meter=first(m for m in meters if m isa TransformerMeasurement)
            wrong=copy(meters); k=findfirst(m->m isa TransformerMeasurement,wrong)
            wrong[k]=TransformerMeasurement(kind=first_meter.kind,transformer=first_meter.transformer,
                variant=first_meter.variant,port=99,terminal=first_meter.terminal,value=0.,sigma=1.)
            @test_throws ArgumentError SEParameters(s,wrong)
            @test !state_estimator_preflight(net,wrong).supported
            native,customer,expected,layout=compat_native_case(net)
            classic=solve_state_estimation(native,customer;per_unit=false,solver_options=("tol"=>1e-9,"max_iter"=>1000))
            @test solve_status(classic).publishable
            @test classic.objective < 1e-6
            for node in keys(layout.free_state_map)
                e=classic.bus[node[1]][node[2]]
                @test complex(e["vr"],e["vi"]) ≈ expected[layout.node_index[node]] atol=1e-5
            end
        end
    end
end

@testset "SE legacy compatibility: rectangular phasors and shared contract" begin
    net=compat_ground_net(:solid)
    ms=[Measurement(kind=:vr,bus="b",reference=nothing,value=225.,sigma=.1),
        Measurement(kind=:vi,bus="b",reference=nothing,value=1.,sigma=.1),
        Measurement(kind=:pinj,bus="b",reference=nothing,value=real((225+im)*conj(-10+2im)),sigma=10.),
        Measurement(kind=:qinj,bus="b",reference=nothing,value=imag((225+im)*conj(-10+2im)),sigma=10.)]
    for per_unit in (false,true)
        legacy=solve_state_estimation(net,ms;per_unit)
        @test solve_status(legacy).publishable
        @test legacy.objective < 1e-7
        @test legacy.observability.observable
    end
    s=compile_state_estimator(net,ms);p=SEParameters(s,ms)
    @test solve_status(solve_state_estimation(s,p,initial_state_estimator(s,p))).publishable
    @test_throws ArgumentError solve_ipopt_state_estimator(s,p,[230.,0.];solver_options=(hessian_approximation="exact",))
    @test_throws ArgumentError solve_ipopt_state_estimator(s,p,[230.,0.];constraint_tolerance=Inf)
end

@testset "SE open-delta bonds: connections and native consistency" begin
    for connection in ("ABBC","BCAC","CABA"), regulator_type in ("A","B")
        net=_se_net_open_delta_reg()
        el=only(values(net["transformer"]["open_delta_regulator"]))
        el["connection"]=connection;el["regulator_type"]=regulator_type
        native,meters,voltage,s=compat_native_case(net)
        a,b=PowerOptLab._se_open_delta_bond(el)
        @test s.node_index[a] == s.node_index[b]
        @test !haskey(s.free_state_map,b)
        r=solve_state_estimation(native,meters;per_unit=false)
        @test solve_status(r).publishable
        @test r.objective < 1e-6
        for node in keys(s.free_state_map)
            @test complex(r.bus[node[1]][node[2]]["vr"],r.bus[node[1]][node[2]]["vi"]) ≈ voltage[s.node_index[node]] atol=1e-5
        end
        @test !haskey(net,"switch") # compilation must not mutate the input
    end
end

@testset "SE net telemetry does not allocate merged customer power" begin
    net=compat_ground_net(:solid)
    net["bus"]["z"]=deepcopy(net["bus"]["b"])
    net["switch"]=Dict("closed"=>Dict("bus_from"=>"b","bus_to"=>"z",
        "terminal_map_from"=>["1","n"],"terminal_map_to"=>["1","n"],"status"=>"closed"))
    cs=[TerminalConnection((bus,"1"),(bus,"n")) for bus in ("b","z")]
    meters=Any[Measurement(kind=:pinj,bus="b",value=-1000.,sigma=10.),
        NetInjectionMeasurement(kind=:pinj,bus="z",value=-2250.,sigma=1.),
        NetInjectionMeasurement(kind=:qinj,bus="b",value=0.,sigma=1.)]
    s=compile_state_estimator(net,meters;injection_connections=cs);p=SEParameters(s,meters)
    # A missing customer's P leaves the individual P row uninformative about
    # voltage; the separate net meter still observes the total -2250 W.
    truth=[225.,0.]
    @test norm(evaluate_state_estimator(s,p,truth).residual)<1e-10
    @test norm(Matrix(residual_jacobian(s,p,truth))[1,:]) == 0.
    for solver in COMPAT_SE_SOLVERS
        r=solver(s,p,initial_state_estimator(s,p))
        @test solve_status(r).publishable
        @test r.state ≈ truth atol=1e-5
    end
end

@testset "SE shared Ipopt: fixed-state continuation remains internal" begin
    net=compat_ground_net(:solid)
    delete!(net,"line");delete!(net["bus"],"b")
    exact=ExactDeviceEquation(ConstantPowerDevice(
        [TerminalConnection(("src","1"),nothing)],ComplexF64[1000]))
    s=compile_state_estimator(net;exact_devices=[exact])
    p=SEParameters(s;exact_devices=[exact],continuation_alpha=0.)
    @test isempty(initial_state_estimator(s,p))
    r=solve_ipopt_state_estimator(s,p,Float64[])
    @test r.status == :converged_regularized
    @test !solve_status(r).publishable
end

@testset "SE shared backends: exact device laws and conflicting telemetry" begin
    net=compat_ground_net(:floating)
    c=[TerminalConnection(("b","1"),("b","n"))]
    current=10-2im;drop=220+2im;power=drop*conj(current)
    models=(ConstantCurrentDevice(c,[current]),ConstantPowerDevice(c,[power]),
        ZIPDevice(c,[.5power],[.25current],[.25current/drop]))
    # Deliberately inconsistent voltage data must not relax exact device KCL.
    meters=[Measurement(kind=:vr,bus="b",value=210.,sigma=.1)]
    for model in models
        exact=[ExactDeviceEquation(model)]
        s=compile_state_estimator(net,meters;exact_devices=exact)
        p=SEParameters(s,meters;exact_devices=exact)
        truth=compat_state(s,Dict(("b","1")=>225+im,("b","n")=>5-im))
        @test norm(evaluate_state_estimator(s,p,truth).constraints)<1e-9
        for solver in COMPAT_SE_SOLVERS
            r=solver(s,p,initial_state_estimator(s,p))
            solve_status(r).publishable || @info "Exact-device solve failed" model=typeof(model) solver status=r.status constraints=norm(r.evaluation.constraints)
            @test solve_status(r).publishable
            @test r.state ≈ truth atol=1e-5
            @test norm(r.evaluation.constraints)<1e-7
            @test norm(r.evaluation.residual)>90.
        end
        if model isa Union{ConstantPowerDevice,ZIPDevice}
            p.continuation_alpha=0.
            r=solve_ipopt_state_estimator(s,p,initial_state_estimator(s,p))
            @test r.status == :converged_regularized
            @test !solve_status(r).publishable
        end
    end
end

@testset "SE source-fixed customer identity and net telemetry" begin
    net=compat_ground_net(:solid)
    net["bus"]["alias"]=deepcopy(net["bus"]["src"])
    net["switch"]=Dict("closed"=>Dict("bus_from"=>"src","bus_to"=>"alias",
        "terminal_map_from"=>["1","n"],"terminal_map_to"=>["1","n"],"status"=>"closed"))
    load=[Measurement(kind=kind,bus="b",value=value,sigma=10.)
          for (kind,value) in ((:pinj,-2250.),(:qinj,-450.))]
    for bus in ("src","alias")
        bad=vcat(load,[Measurement(kind=kind,bus=bus,value=0.,sigma=10.) for kind in (:pinj,:qinj)])
        report=state_estimator_preflight(net,bad)
        @test !report.supported
        @test any(f->occursin("NetInjectionMeasurement",f.message),report.findings)
        @test_throws SEUnsupportedNetwork compile_state_estimator(net,bad)
        @test_throws SEUnsupportedNetwork solve_state_estimation(net,bad)
    end
    partial=deepcopy(net)
    push!(partial["bus"]["src"]["terminal_names"],"2")
    @test_throws ArgumentError solve_state_estimation(partial,load)
    # An unmetered source alias must not be mistaken for an uncovered customer.
    native=solve_state_estimation(net,load;per_unit=false)
    @test solve_status(native).publishable
    @test native.objective < 1e-6
    net_meters=[NetInjectionMeasurement(kind=kind,bus="alias",value=value,sigma=1.)
                for (kind,value) in ((:pinj,2300.),(:qinj,460.))]
    s=compile_state_estimator(net,net_meters;injection_connections=[TerminalConnection(("b","1"),("b","n"))])
    p=SEParameters(s,net_meters)
    for solver in COMPAT_SE_SOLVERS
        r=solver(s,p,initial_state_estimator(s,p))
        @test solve_status(r).publishable
        @test r.state ≈ [225.,1.] atol=1e-5
    end
end
