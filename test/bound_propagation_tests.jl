using Test
using JSON3
using JSONSchema
using BMOPFTools

_bq(f,id,q,index=();component=:magnitude)=QuantityKey(f,id,q,index;component)

function _bp_fixture(;load=3000.0)
    Dict{String,Any}(
        "terminal_conventions"=>Dict("phase"=>["a"],"neutral"=>["return"]),
        "bus"=>Dict(
            "src"=>Dict{String,Any}("terminal_names"=>["a","return"],"perfectly_grounded_terminals"=>["return"]),
            "dst"=>Dict{String,Any}("terminal_names"=>["a","return"],"perfectly_grounded_terminals"=>["return"],
                "vpn_min"=>[210.0],"vpn_max"=>[250.0])),
        "voltage_source"=>Dict("s"=>Dict{String,Any}("bus"=>"src","terminal_map"=>["a"],"v_magnitude"=>[240.0],"v_angle"=>[0.0])),
        "linecode"=>Dict("lc"=>Dict{String,Any}("R_series_1_1"=>0.1,"X_series_1_1"=>0.0)),
        "line"=>Dict("l"=>Dict{String,Any}("bus_from"=>"src","bus_to"=>"dst","terminal_map_from"=>["a"],
            "terminal_map_to"=>["a"],"linecode"=>"lc","length"=>1.0)),
        "load"=>Dict("d"=>Dict{String,Any}("bus"=>"dst","terminal_map"=>["a","return"],"configuration"=>"SINGLE_PHASE",
            "p_nom"=>[load],"q_nom"=>[4000.0])))
end

@testset "Independent feasible witnesses: coupled π line" begin
    # Binary-exact coefficients and phasors make this a rational oracle.
    # Construct a feasible state first, then check every relevant enclosure.
    vt=ComplexF64[10+2im,-1+1im]; series=ComplexF64[2-1im,-1+2im]
    Z=ComplexF64[1+1im 1im; 1im 2]
    Yf=ComplexF64[0.5 0.25im;0.25im 0.5]
    Yt=ComplexF64[0.25 0;0 0.25]
    vf=vt+Z*series; shf=Yf*vf; sht=Yt*vt
    iff=series+shf; it=-series+sht
    line=Dict{String,Any}("bus_from"=>"f","bus_to"=>"t","terminal_map_from"=>["a","b"],"terminal_map_to"=>["a","b"],
        "i_max"=>[ceil(max(abs(iff[k]),abs(it[k]))) for k in 1:2])
    for (matrix,rprefix,iprefix) in ((Z,"R_series_","X_series_"),(Yf,"G_from_","B_from_"),(Yt,"G_to_","B_to_"))
        for i in 1:2,j in 1:2
            line["$rprefix$(i)_$(j)"]=real(matrix[i,j]); line["$iprefix$(i)_$(j)"]=imag(matrix[i,j])
        end
    end
    buses=Dict{String,Any}(); gens=Dict{String,Any}()
    for (id,v,i) in (("f",vf,iff),("t",vt,it))
        buses[id]=Dict("terminal_names"=>["a","b"],"v_max"=>ceil.(abs.(v)))
        s=v.*conj.(i)
        gens[id]=Dict("bus"=>id,"terminal_map"=>["a","b"],"p_min"=>real.(s),"p_max"=>real.(s),"q_min"=>imag.(s),"q_max"=>imag.(s))
    end
    net=Dict{String,Any}("bus"=>buses,"line"=>Dict("l"=>line),"generator"=>gens)
    report=propagate_bounds(net)
    @test report.status!=:conflict
    function contains_complex(family,id,quantity,index,z)
        rc=quantity==:power ? :active : :real; ic=quantity==:power ? :reactive : :imag
        for (component,value) in ((rc,real(z)),(ic,imag(z)))
            d=bound(report,QuantityKey(family,id,quantity,index;component))
            @test d.lower<=value<=d.upper
        end
        d=bound(report,QuantityKey(family,id,quantity,index))
        @test Rational{BigInt}(d.lower)^2<=Rational{BigInt}(real(z))^2+Rational{BigInt}(imag(z))^2
        @test !isfinite(d.upper) || Rational{BigInt}(d.upper)^2>=Rational{BigInt}(real(z))^2+Rational{BigInt}(imag(z))^2
    end
    for (k,t) in enumerate(("a","b"))
        contains_complex(:bus,"f",:voltage,(t,),vf[k])
        contains_complex(:bus,"t",:voltage,(t,),vt[k])
        contains_complex(:line,"l",:current,("series",string(k)),series[k])
        for (side,v,i,sh) in (("from",vf,iff,shf),("to",vt,it,sht))
            contains_complex(:line,"l",:current,(side,t),i[k])
            contains_complex(:line,"l",:shunt_current,(side,t),sh[k])
            contains_complex(:line,"l",:power,(side,t),v[k]*conj(i[k]))
        end
    end
end

@testset "Both power-flow branches remain enclosed" begin
    # E=10 V, R=1 ohm, P=9 W: v*(10-v)=9 admits v=1 and v=9.
    net=_bp_fixture(load=9.0)
    net["voltage_source"]["s"]["v_magnitude"]=[10.0]
    net["linecode"]["lc"]["R_series_1_1"]=1.0
    net["load"]["d"]["q_nom"]=[0.0]
    net["bus"]["dst"]["vpn_min"]=[0.0]
    net["bus"]["dst"]["vpn_max"]=[10.0]
    r=propagate_bounds(net)
    @test r.status!=:conflict
    v=bound(r,_bq(:bus,"dst",:voltage,"a")); i=bound(r,_bq(:load,"d",:current,("coil","1")))
    @test v.lower<=1<=9<=v.upper
    @test i.lower<=1<=9<=i.upper
end

@testset "BMOPF schema and exchange round trip" begin
    net=_bp_fixture()
    net["meta"]=Dict("\$schema"=>"https://raw.githubusercontent.com/frederikgeth/bmopf-report/main/draft_schema_and_networks/draft_bmopf_schema.json")
    report=propagate_bounds(net)
    @test report.validation["jsonschema_valid"]
    filled,_=materialize_bounds(net,report)
    io=IOBuffer(); BMOPFTools.write_bmopf(filled,io)
    serialized=String(take!(io))
    raw=JSON3.read(serialized,Dict{String,Any})
    @test !haskey(raw,"_meta")
    @test !haskey(raw["bus"]["dst"],"neutral_terminal")
    findings=BMOPFTools.Finding[]
    @test BMOPFTools.schema_check(raw,findings)["jsonschema_valid"]
    @test propagate_bounds(raw).status!=:conflict
end

@testset "Outward interval arithmetic" begin
    P=PowerOptLab
    for (a,b) in ((0.1,0.2),(floatmin(Float64),0.1),(nextfloat(0.0),0.5),
                  (floatmax(Float64),2.0),(-floatmax(Float64),2.0),(-1.3,0.7))
        for (op,func) in ((+,P._badd),(*,P._bmul))
            d=func(P._bpoint(a),P._bpoint(b)); exact=op(Rational{BigInt}(a),Rational{BigInt}(b))
            @test !isfinite(d.lo) || Rational{BigInt}(d.lo)<=exact
            @test !isfinite(d.hi) || Rational{BigInt}(d.hi)>=exact
            @test !isnan(d.lo) && !isnan(d.hi)
        end
    end
    for x in (0.0,nextfloat(0.0),2.0,1e300)
        d=P._bsqrt(P._bpoint(x))
        @test Rational{BigInt}(d.lo)^2<=Rational{BigInt}(x)<=Rational{BigInt}(d.hi)^2
    end
    @test P._bmul(P._bpoint(0),P._BENTIRE)==P._bpoint(0)
    @test P._bdiv(P._bpoint(1),P._BInterval(-1,1))==P._BENTIRE
    # Enclosures must not depend on a caller's global BigFloat configuration.
    setprecision(BigFloat,2) do
        for q in (big(1)//big(10),-(big(1)//big(3)),big(1)//(big(1)<<1100),
                  (big(1)<<1100)//big(3),Rational{BigInt}(floatmax(Float64))+1)
            lo=P._bround(q,false); hi=P._bround(q,true)
            @test !isfinite(lo) || Rational{BigInt}(lo)<=q
            @test !isfinite(hi) || Rational{BigInt}(hi)>=q
        end
    end
end

@testset "Power inversion and evidence" begin
    P=PowerOptLab; e=P._BEngine(0.0,0.0)
    v=P._bcomplex!(e,:test,"t",:voltage,(),:V)
    i=P._bcomplex!(e,:test,"t",:current,(),:A)
    s=P._bcomplex!(e,:test,"t",:power,(),:VA)
    P._brelation!(e,:power,[v...,i...,s...],nothing,["/test"])
    for (id,val) in ((i[1],3),(i[2],4),(s[1],56),(s[2],-33))
        P._bdeclare!(e,id,val,val,"/test")
    end
    status,steps=P._brun!(e;max_steps=1000,time_limit=Inf)
    @test status==:fixed_point
    @test e.domains[v[1]].lo==e.domains[v[1]].hi==12
    @test e.domains[v[2]].lo==e.domains[v[2]].hi==5
    @test e.domains[v[3]].lo==e.domains[v[3]].hi==13
    @test all(all(p->p<n,proof.parents) for (n,proof) in enumerate(e.evidence))
end

@testset "BMOPF propagation and fill-in" begin
    net=_bp_fixture(); original=deepcopy(net)
    r=propagate_bounds(net)
    @test net==original
    @test r.status in (:fixed_point,:tolerance_stopped)
    i=bound(r,_bq(:load,"d",:current,("coil","1")))
    @test 20<=i.lower<=i.upper<=nextfloat(5000/210)
    line=bound(r,_bq(:line,"l",:current,("from","a")))
    @test line.upper<=nextfloat(5000/210)
    @test bound(r,_bq(:bus,"dst",:voltage,"return")).upper==0
    @test bound(r,_bq(:bus,"dst",:voltage,"a")).upper<=250
    @test !isempty(explain_bound(r,_bq(:load,"d",:current,("coil","1"))))
    @test all(all(p->p<n,proof.parents) for (n,proof) in enumerate(r.evidence))
    @test any(x->x.object=="/voltage_source/s" && x.status==:relaxed,coverage(r))
    @test any(x->x.object=="/linecode/lc" && x.status==:data,coverage(r))
    filled,patch=materialize_bounds(net,r)
    @test net==original
    @test !isempty(patch)
    @test filled["bus"]["dst"]["vpn_max"]==[250.0]
    @test filled["bus"]["dst"]["v_max"][1]<=250.0
    @test filled["line"]["l"]["i_max"][1]<=nextfloat(5000/210)
    reparsed=BMOPFTools.parse_bmopf(JSON3.write(filled);from_string=true)
    @test reparsed["line"]["l"]["i_max"]==filled["line"]["l"]["i_max"]
    @test propagate_bounds(filled).status!=:conflict
    changed=deepcopy(net); changed["load"]["d"]["p_nom"]=[3001.0]
    @test_throws ArgumentError materialize_bounds(changed,r)
    @test_throws ArgumentError materialize_bounds(net,r;mode=:replace)
    @test_throws ArgumentError explain_bound(r,first(keys(r.domains));endpoint=:wrong)
    @test_throws KeyError bound(r,_bq(:bus,"missing",:voltage,"a"))
    side=bounds_sidecar(r)
    @test JSON3.read(JSON3.write(side))["format"]=="PowerOptLab.bounds"
    schema=JSONSchema.Schema(JSON3.read(read(joinpath(@__DIR__,"..","data","schema","bounds.schema.json"),String)))
    @test JSONSchema.validate(schema,side)===nothing
    @test bounds_sidecar(propagate_bounds(JSON3.write(net)))==side
    reversed=Dict{String,Any}(k=>net[k] for k in reverse(sort(collect(keys(net)))))
    @test bounds_sidecar(propagate_bounds(reversed))==side
    limited=propagate_bounds(net;max_steps=0)
    @test limited.status==:budget_exhausted
    @test limited.steps==0
    @test propagate_bounds(net;time_limit=0).status==:budget_exhausted
    @test_throws ArgumentError propagate_bounds(net;rtol=-1)
end

@testset "Floating neutral and missing voltage lower bound" begin
    n=Dict{String,Any}("bus"=>Dict("b"=>Dict{String,Any}("terminal_names"=>["a","n"],"vpn_min"=>[210.0],"vpn_max"=>[250.0])))
    r=propagate_bounds(n)
    @test bound(r,_bq(:bus,"b",:connection_voltage,("a","n"))).upper==250
    @test isinf(bound(r,_bq(:bus,"b",:voltage,"a")).upper)
    @test isinf(bound(r,_bq(:bus,"b",:voltage,"n")).upper)
    @test isempty(explain_bound(r,_bq(:bus,"b",:voltage,"n")))
    out,_=materialize_bounds(n,r)
    @test !haskey(out["bus"]["b"],"v_max")
    @test occursin("null",JSON3.write(bounds_sidecar(r)))
    n["bus"]["b"]["vn_max"]=10.0
    r=propagate_bounds(n)
    d=bound(r,_bq(:bus,"b",:voltage,"a"))
    @test d.lower==200 && d.upper==260
    # No finite current bound may be inferred by dividing by nominal voltage.
    n["bus"]["b"]=Dict("terminal_names"=>["a","n"],"vpn_max"=>[250.0])
    n["load"]=Dict("d"=>Dict("bus"=>"b","terminal_map"=>["a","n"],"p_nom"=>[1.0],"q_nom"=>[0.0]))
    n["generator"]=Dict("g"=>Dict("bus"=>"b","terminal_map"=>["a","n"]))
    @test isinf(bound(propagate_bounds(n),_bq(:load,"d",:current,("coil","1"))).upper)
end

@testset "Unsupported devices preserve injections; conflicts are explicit" begin
    n=_bp_fixture(); n["load"]["d"]["model"]="constant_impedance"; n["load"]["d"]["v_nom"]=230.0
    r=propagate_bounds(n)
    @test any(x->x.object=="/load/d" && x.status==:relaxed,coverage(r))
    @test bound(r,_bq(:load,"d",:power,("coil","1");component=:active)).lower!=3000
    n=Dict{String,Any}("bus"=>Dict("b"=>Dict("terminal_names"=>["a"],"v_min"=>[1.0],"v_max"=>[2.0])),
        "load"=>Dict("d"=>Dict("bus"=>"b","terminal_map"=>["a"],"p_nom"=>[1.0],"q_nom"=>[0.0])))
    @test propagate_bounds(n).status==:conflict
    n["ibr"]=Dict("i"=>Dict("bus"=>"b","terminal_map"=>["a"]))
    r=propagate_bounds(n)
    @test r.status!=:conflict
    @test any(x->x.object=="/ibr/i" && x.status==:unsupported,coverage(r))
    n["bus"]["b"]["v_min"]=[3.0]
    conflict=propagate_bounds(n)
    @test conflict.status==:conflict
    @test !isempty(conflict.conflicts)
    @test_throws ArgumentError materialize_bounds(n,conflict)
    @test JSON3.read(JSON3.write(bounds_sidecar(conflict)))["status"]=="conflict"
end

@testset "Parallel sources, source exemptions and zero current" begin
    n=_bp_fixture(load=0.0); n["load"]["d"]["q_nom"]=[0.0]
    n["voltage_source"]["second"]=deepcopy(n["voltage_source"]["s"])
    r=propagate_bounds(n)
    @test isinf(bound(r,_bq(:voltage_source,"s",:current,("coil","1"))).upper)
    @test bound(r,_bq(:load,"d",:current,("coil","1"))).upper==0
    n["bus"]["src"]["v_max"]=[1.0] # Native builder skips source-fixed terminals.
    @test propagate_bounds(n).status!=:conflict
    @test bound(propagate_bounds(n),_bq(:bus,"src",:voltage,"a")).lower==240
    out,_=materialize_bounds(n,propagate_bounds(n);mode=:intersect)
    @test out["bus"]["src"]["v_max"]==[1.0]
    @test !haskey(out["bus"]["src"],"v_min")
    n=_bp_fixture(); n["line"]["l"]["i_max"]=[0.0]
    @test propagate_bounds(n).status==:conflict
end

@testset "Line shunts and rating precedence" begin
    n=Dict{String,Any}(
        "bus"=>Dict("f"=>Dict("terminal_names"=>["a"],"v_max"=>[100.0]),"t"=>Dict("terminal_names"=>["a"],"v_max"=>[100.0])),
        "linecode"=>Dict("lc"=>Dict("R_series_1_1"=>0.0,"G_from_1_1"=>0.05,"G_to_1_1"=>0.05,"i_max"=>[100.0])),
        "line"=>Dict("l"=>Dict("bus_from"=>"f","bus_to"=>"t","terminal_map_from"=>["a"],"terminal_map_to"=>["a"],"linecode"=>"lc","length"=>1.0)),
        "generator"=>Dict("g1"=>Dict("bus"=>"f","terminal_map"=>["a"]),"g2"=>Dict("bus"=>"t","terminal_map"=>["a"])))
    r=propagate_bounds(n)
    @test r.status!=:conflict
    @test 100<bound(r,_bq(:line,"l",:current,("series","1"))).upper<=nextfloat(105.0)
    @test bound(r,_bq(:line,"l",:current,("from","a"))).upper==100
    out,_=materialize_bounds(n,r)
    @test out["line"]["l"]["i_max"]==[100.0]
    n["line"]["l"]["i_max"]=[120.0] # Device field overrides linecode.
    @test bound(propagate_bounds(n),_bq(:line,"l",:current,("from","a"))).upper==120
    # Inline matrices are absolute, and do not scale with the line's length.
    n["line"]["l"]["R_series_1_1"]=0.0
    n["line"]["l"]["length"]=1000.0
    @test bound(propagate_bounds(n),_bq(:line,"l",:shunt_current,("from","a"))).upper==0
end

@testset "Time scope and terminal connections" begin
    n=_bp_fixture()
    n["time_series"]=Dict("p"=>Dict("values"=>[1.0,0.5]))
    n["load"]["d"]["time_series"]=Dict("p_nom"=>"p")
    r=propagate_bounds(n;scope=BoundScope(snapshot=2))
    @test bound(r,_bq(:load,"d",:power,("coil","1");component=:active)).lower==1500
    out,_=materialize_bounds(n,r)
    @test !haskey(out,"time_series")
    @test out["load"]["d"]["p_nom"]==[1500.0]
    @test_throws BoundsError propagate_bounds(n;scope=BoundScope(snapshot=3))
    @test_throws ArgumentError BoundScope(snapshot=0)
    n=Dict{String,Any}("bus"=>Dict("b"=>Dict("terminal_names"=>["a","b"],"vpp_min"=>[200.0],"vpp_max"=>[250.0])),
        "load"=>Dict("d"=>Dict("bus"=>"b","terminal_map"=>["a","b"],"configuration"=>"SINGLE_PHASE","p_nom"=>[2000.0],"q_nom"=>[0.0])),
        "generator"=>Dict("g"=>Dict("bus"=>"b","terminal_map"=>["a","b"],"configuration"=>"SINGLE_PHASE")))
    for cfg in ("SINGLE_PHASE","DELTA")
        n["load"]["d"]["configuration"]=cfg
        r=propagate_bounds(n)
        @test r.status!=:conflict
        d=bound(r,_bq(:load,"d",:current,("coil","1")))
        @test d.lower==8 && d.upper==10
    end
end
