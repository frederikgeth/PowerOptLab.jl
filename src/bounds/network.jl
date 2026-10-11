# Versioned, deliberately small BMOPF adapter. All unsupported devices retain
# free terminal injections; their absence must not create fictitious KCL bounds.
const _BOUND_ADAPTER_VERSION = "bmopf-72e6cec-static-bounds/v1"
_bp(s) = replace(replace(string(s),"~"=>"~0"),"/"=>"~1")
_bpath(parts...) = "/"*join(_bp.(parts),"/")
_bsorted(d) = sort!(collect(keys(d));by=string)

function _bcanonical(x)
    if x isa AbstractDict
        return "{"*join((JSON3.write(string(k))*":"*_bcanonical(x[k]) for k in _bsorted(x)),",")*"}"
    elseif x isa AbstractVector || x isa Tuple
        return "["*join(_bcanonical.(x),",")*"]"
    end
    String(JSON3.write(x))
end
function _bfingerprint(net,scope)
    clean=deepcopy(net)
    haskey(clean,"_meta") && delete!(clean["_meta"],"parsed_at")
    bytes2hex(sha256(_BOUND_ADAPTER_VERSION*"|"*string(scope.snapshot)*"|"*_bcanonical(clean)))
end

function _binput(input)
    if input isa AbstractDict
        # Use the same terminal/migration path for dictionaries and JSON.
        return BMOPFTools.parse_bmopf(String(JSON3.write(input));from_string=true)
    elseif input isa AbstractString
        # Raw JSON is recognized before isfile (long JSON can exceed path limits).
        startswith(lstrip(input),"{") && return BMOPFTools.parse_bmopf(String(input);from_string=true)
        return BMOPFTools.parse_bmopf(String(input))
    end
    throw(ArgumentError("Expected a BMOPF dictionary, JSON object string, or file path."))
end
function _bnumber(x,label;infinite=false)
    x isa Real && !(x isa Bool) || throw(ArgumentError("$label must be numeric"))
    f=Float64(x)
    (isnan(f) || (!infinite && !isfinite(f))) && throw(ArgumentError("$label must be finite"))
    f
end
function _bcollection(net,name)
    d=get(net,name,Dict{String,Any}())
    d isa AbstractDict || throw(ArgumentError("$name must be an object collection"))
    d
end
function _bcover!(e,object,status,reason)
    push!(e.objects,(object=object,status=status,reason=reason))
end
function _bdeclare!(e,i,lo,hi,path)
    a=_bnumber(lo,path;infinite=true); b=_bnumber(hi,path;infinite=true)
    _bupdate!(e,i,_BInterval(a,b),:declared,Int[],[path];force=true)
end
function _blimits!(e,ids,obj,field,path;lower=false,exact_length=false)
    v=get(obj,field,nothing); v===nothing && return
    v isa AbstractVector || throw(ArgumentError("$path/$field must be an array"))
    (length(v)<=length(ids) && (!exact_length || length(v)==length(ids))) ||
        throw(ArgumentError("$path/$field has the wrong number of entries"))
    for k in eachindex(v)
        val=_bnumber(v[k],"$path/$field/$k";infinite=true)
        e.keys[ids[k]].component==:magnitude && val<0 && throw(ArgumentError("Negative magnitude at $path/$field/$k"))
        _bdeclare!(e,ids[k],lower ? val : -Inf,lower ? Inf : val,"$path/$field/$(k-1)")
    end
end

mutable struct _BNetwork
    engine::_BEngine
    net::Dict{String,Any}
    voltages::Dict{Tuple{String,String},NTuple{3,Int}}
    kcl::Dict{Tuple{String,String},Vector{Tuple{NTuple{3,Int},Int}}}
    grounded::Set{Tuple{String,String}}
    fixed::Set{Tuple{String,String}}
end
function _bvoltage(c,bus,t)
    k=(string(bus),string(t))
    haskey(c.voltages,k) || throw(ArgumentError("Unknown bus terminal $k"))
    c.voltages[k]
end
function _bconnection!(c,family,id,index,bus,pos,neg,path)
    e=c.engine; vp=_bvoltage(c,bus,pos)
    v=_bcomplex!(e,family,id,:connection_voltage,index,:V)
    # Reuse a bus-level pair domain directly. Two separate equations involving
    # unbounded common-mode voltages cannot transfer the pair's magnitude bound
    # through scalar interval propagation alone.
    if family!=:bus && neg!==nothing
        for (a,b,sign) in ((pos,neg,1),(neg,pos,-1))
            key=QuantityKey(:bus,string(bus),:connection_voltage,(a,b))
            if haskey(e.lookup,key)
                z=(e.lookup[QuantityKey(:bus,string(bus),:connection_voltage,(a,b);component=:real)],
                   e.lookup[QuantityKey(:bus,string(bus),:connection_voltage,(a,b);component=:imag)],e.lookup[key])
                _brow!(e,[(v,1),(z,-sign)],[path])
                return v
            end
        end
    end
    terms=Any[(v,1),(vp,-1)]
    neg===nothing || push!(terms,(_bvoltage(c,bus,neg),1))
    _brow!(e,terms,[path]); v
end
function _binject!(c,bus,t,current,sign)
    t===nothing && return
    key=(string(bus),string(t))
    haskey(c.kcl,key) || throw(ArgumentError("Unknown injection terminal $key"))
    push!(c.kcl[key],(current,sign))
end
function _bpower!(c,family,id,index,v,i,path)
    s=_bcomplex!(c.engine,family,id,:power,index,:VA)
    _brelation!(c.engine,:power,[v...,i...,s...],nothing,[path]); s
end

function _bbuses!(c)
    e=c.engine; net=c.net
    for id in _bsorted(_bcollection(net,"bus"))
        b=net["bus"][id]; path=_bpath("bus",id)
        ts=string.(get(b,"terminal_names",String[]))
        !isempty(ts) && length(unique(ts))==length(ts) || throw(ArgumentError("$path needs distinct terminals"))
        for t in ts
            c.voltages[(id,t)]=_bcomplex!(e,:bus,id,:voltage,(t,),:V)
            c.kcl[(id,t)]=Tuple{NTuple{3,Int},Int}[]
        end
        for t in get(b,"perfectly_grounded_terminals",String[])
            _bvoltage(c,id,t); push!(c.grounded,(id,string(t)))
        end
    end
    # Matches BMOPFTools' native source-ground and bus-limit applicability.
    for sid in _bsorted(_bcollection(net,"voltage_source"))
        s=net["voltage_source"][sid]; b=string(s["bus"])
        haskey(net["bus"],b) || throw(ArgumentError("Unknown source bus $b"))
        get(s,"configuration","WYE") in ("WYE","SINGLE_PHASE") || continue
        for t in get(s,"terminal_map",String[])
            _bvoltage(c,b,t); push!(c.fixed,(b,string(t)))
        end
        nt=_neutral_terminal(net["bus"][b])
        nt===nothing || push!(c.grounded,(b,nt))
    end
    for (b,t) in sort!(collect(c.grounded))
        v=_bvoltage(c,b,t)
        _bdeclare!(e,v[3],0,0,_bpath("bus",b,"ground",t))
        i=_bcomplex!(e,:bus,b,:ground_current,(t,),:A)
        _binject!(c,b,t,i,1)
    end
    for id in _bsorted(_bcollection(net,"bus"))
        b=net["bus"][id]; path=_bpath("bus",id)
        ts=string.(b["terminal_names"]); nt=_neutral_terminal(b)
        phases=[t for t in ts if t!=nt]
        applicable(t)=!((id,t) in c.grounded || (id,t) in c.fixed)
        for field in ("v_min","v_max","vpn_min","vpn_max")
            val=get(b,field,nothing)
            val===nothing || (val isa AbstractVector && length(val)==length(phases)) ||
                throw(ArgumentError("$path/$field must have one entry per phase"))
        end
        for (k,t) in enumerate(phases)
            vm=_bvoltage(c,id,t)[3]
            if applicable(t)
                for (f,lower) in (("v_min",true),("v_max",false))
                    val=get(b,f,nothing); val===nothing && continue
                    v=_bnumber(val[k],"$path/$f";infinite=true)
                    v>=0 || throw(ArgumentError("Negative voltage limit"))
                    _bdeclare!(e,vm,lower ? v : 0,lower ? Inf : v,"$path/$f/$(k-1)")
                end
            end
            if nt!==nothing
                cv=_bconnection!(c,:bus,id,(t,nt),id,t,nt,path)
                if applicable(t)
                    for (f,lower) in (("vpn_min",true),("vpn_max",false))
                        val=get(b,f,nothing); val===nothing && continue
                        v=_bnumber(val[k],"$path/$f";infinite=true)
                        v>=0 || throw(ArgumentError("Negative connection voltage limit"))
                        _bdeclare!(e,cv[3],lower ? v : 0,lower ? Inf : v,"$path/$f/$(k-1)")
                    end
                end
            end
        end
        if nt!==nothing && haskey(b,"vn_max") && !((id,nt) in c.grounded)
            v=_bnumber(b["vn_max"],"$path/vn_max";infinite=true)
            v>=0 || throw(ArgumentError("Negative neutral voltage limit"))
            _bdeclare!(e,_bvoltage(c,id,nt)[3],0,v,"$path/vn_max")
        end
        pairs=[(phases[i],phases[j]) for i in eachindex(phases) for j in i+1:length(phases)]
        for f in ("vpp_min","vpp_max")
            val=get(b,f,nothing)
            val===nothing || (val isa AbstractVector && length(val)==length(pairs)) || throw(ArgumentError("$path/$f has wrong pair count"))
        end
        for (k,(t,u)) in enumerate(pairs)
            v=_bconnection!(c,:bus,id,(t,u),id,t,u,path)
            applicable(t) && applicable(u) || continue
            for (f,lower) in (("vpp_min",true),("vpp_max",false))
                val=get(b,f,nothing); val===nothing && continue
                x=_bnumber(val[k],"$path/$f";infinite=true)
                x>=0 || throw(ArgumentError("Negative pair voltage limit"))
                _bdeclare!(e,v[3],lower ? x : 0,lower ? Inf : x,"$path/$f/$(k-1)")
            end
        end
        omitted=[f for f in ("vpos_min","vpos_max","vneg_max","vzero_max","vuf_max","va_diff_min","va_diff_max") if haskey(b,f)]
        _bcover!(e,path,isempty(omitted) ? :supported : :relaxed,
            isempty(omitted) ? "Terminal, neutral and pair voltage bounds; native applicability." : "Omitted sequence/angle constraints: "*join(omitted,", "))
    end
end

# BMOPF matrix shorthand: direct entries take precedence; missing reciprocal
# entries are mirrored, other missing entries are zero. Keep coefficients exact.
function _bmatrix(d,prefix,n,path;scale=_brat(1))
    for key in keys(d)
        m=match(Regex("^"*prefix*raw"(\d+)_(\d+)$"),string(key))
        m===nothing && continue
        i,j=parse.(Int,m.captures)
        1<=i<=n && 1<=j<=n || throw(ArgumentError("$path/$key exceeds terminal-map dimension"))
    end
    [_brat(_bnumber(get(d,"$(prefix)$(i)_$(j)",get(d,"$(prefix)$(j)_$(i)",0.0)),path))*scale for i in 1:n,j in 1:n]
end
function _bmatrixdim(d,prefix)
    n=0
    for key in keys(d)
        m=match(Regex("^"*prefix*raw"(\d+)_(\d+)$"),string(key))
        m===nothing || (n=max(n,parse(Int,m[1]),parse(Int,m[2])))
    end
    n
end
function _bports(c,obj,path)
    f=string(obj["bus_from"]); t=string(obj["bus_to"])
    fm=string.(obj["terminal_map_from"]); tm=string.(obj["terminal_map_to"])
    !isempty(fm) && length(fm)==length(tm) || throw(ArgumentError("$path terminal maps must have equal nonzero lengths"))
    length(unique(fm))==length(fm) && length(unique(tm))==length(tm) || throw(ArgumentError("$path repeats a terminal"))
    for x in fm; _bvoltage(c,f,x); end
    for x in tm; _bvoltage(c,t,x); end
    (f,t,fm,tm)
end
function _blines!(c)
    e=c.engine; lcs=_bcollection(c.net,"linecode")
    for id in _bsorted(_bcollection(c.net,"line"))
        line=c.net["line"][id]; path=_bpath("line",id)
        f,t,fm,tm=_bports(c,line,path); n=length(fm)
        inline=any(k->occursin(r"^[RX]_series_\d+_\d+$",k),keys(line))
        lc=get(lcs,get(line,"linecode",""),nothing)
        data=inline ? line : lc
        if data===nothing || !any(k->occursin(r"^[RX]_series_\d+_\d+$",k),keys(data))
            _bunknown!(c,path,line,"No resolved series impedance; free terminal injections."); continue
        end
        len=inline ? 1.0 : _bnumber(get(line,"length",1.0),"$path/length")
        len>=0 || throw(ArgumentError("Negative line length"))
        src=inline ? path : _bpath("linecode",line["linecode"])
        sources=inline ? [path] : [path,src]
        dr=_bmatrixdim(data,"R_series_"); dx=_bmatrixdim(data,"X_series_")
        max(dr,dx)==n && dr in (0,n) && dx in (0,n) ||
            throw(ArgumentError("$path series matrix dimensions do not match its terminal maps"))
        R=_bmatrix(data,"R_series_",n,src;scale=_brat(len)); X=_bmatrix(data,"X_series_",n,src;scale=_brat(len))
        series=[_bcomplex!(e,:line,id,:current,("series",string(k)),:A) for k in 1:n]
        for k in 1:n
            terms=Any[(_bvoltage(c,f,fm[k]),1),(_bvoltage(c,t,tm[k]),-1)]
            append!(terms,[(series[j],-complex(R[k,j],X[k,j])) for j in 1:n])
            _brow!(e,terms,sources)
        end
        ratings=Dict{String,Any}()
        for fld in ("i_max","s_max")
            val=get(line,fld,lc===nothing ? nothing : get(lc,fld,nothing))
            val===nothing || (ratings[fld]=val)
        end
        for (side,bus,terms,sgn) in (("from",f,fm,1),("to",t,tm,-1))
            G=_bmatrix(data,"G_$(side)_",n,src;scale=_brat(len)); B=_bmatrix(data,"B_$(side)_",n,src;scale=_brat(len))
            currents=NTuple{3,Int}[]; powers=NTuple{3,Int}[]
            for k in 1:n
                sh=_bcomplex!(e,:line,id,:shunt_current,(side,terms[k]),:A)
                row=Any[(sh,1)]
                append!(row,[(_bvoltage(c,bus,terms[j]),-complex(G[k,j],B[k,j])) for j in 1:n])
                _brow!(e,row,sources)
                i=_bcomplex!(e,:line,id,:current,(side,terms[k]),:A)
                _brow!(e,[(i,1),(series[k],-sgn),(sh,-1)],sources)
                push!(currents,i); _binject!(c,bus,terms[k],i,-1)
                push!(powers,_bpower!(c,:line,id,(side,terms[k]),_bvoltage(c,bus,terms[k]),i,path))
            end
            for (field,items) in (("i_max",currents),("s_max",powers))
                origin=haskey(line,field) ? path : _bpath("linecode",get(line,"linecode",""))
                _blimits!(e,[x[3] for x in items],ratings,field,origin)
            end
        end
        _bcover!(e,path,:supported,"Full matrix KVL, both π-shunts, both terminal current/power limits.")
    end
end

function _bswitches!(c)
    e=c.engine
    for id in _bsorted(_bcollection(c.net,"switch"))
        sw=c.net["switch"][id]; path=_bpath("switch",id)
        f,t,fm,tm=_bports(c,sw,path); currents=NTuple{3,Int}[]; powers=NTuple{3,Int}[]
        opened=get(sw,"open_switch",false)
        opened isa Bool || throw(ArgumentError("$path/open_switch must be boolean"))
        for k in eachindex(fm)
            i=_bcomplex!(e,:switch,id,:current,("from",fm[k]),:A); push!(currents,i)
            _binject!(c,f,fm[k],i,-1); _binject!(c,t,tm[k],i,1)
            if opened
                _bdeclare!(e,i[3],0,0,"$path/open_switch")
            else
                _brow!(e,[(_bvoltage(c,f,fm[k]),1),(_bvoltage(c,t,tm[k]),-1)],[path])
            end
            push!(powers,_bpower!(c,:switch,id,("from",fm[k]),_bvoltage(c,f,fm[k]),i,path))
        end
        _blimits!(e,[i[3] for i in currents],sw,"i_max",path)
        _blimits!(e,[s[3] for s in powers],sw,"s_max",path)
        _bcover!(e,path,:supported,"Fixed open/closed state, mapped voltage equality, conductor currents.")
    end
end

function _bcoils(c,obj,path;source=false)
    bus=string(obj["bus"]); tm=string.(obj["terminal_map"])
    !isempty(tm) && length(unique(tm))==length(tm) || throw(ArgumentError("$path needs distinct terminal_map entries"))
    for t in tm; _bvoltage(c,bus,t); end
    cfg=get(obj,"configuration","WYE"); labels=_neutral_labels(c.net)
    if !source && cfg=="SINGLE_PHASE" && length(tm)==2
        return bus,[(tm[1],tm[2])]
    elseif cfg in ("WYE","SINGLE_PHASE")
        np=_neutral_pos(tm,labels)
        neutral=np===nothing ? (source ? _neutral_terminal(c.net["bus"][bus]) : nothing) : tm[np]
        return bus,[(tm[k],neutral) for k in _phase_positions(tm,labels)]
    elseif !source && cfg=="DELTA"
        return bus,[(tm[k],tm[mod1(k+1,length(tm))]) for k in eachindex(tm)]
    end
    return bus,nothing
end
function _bdevices!(c,family::String)
    e=c.engine
    for id in _bsorted(_bcollection(c.net,family))
        obj=c.net[family][id]; path=_bpath(family,id)
        source=family=="voltage_source"; generator=family=="generator"; load=family=="load"
        bus,coils=_bcoils(c,obj,path;source)
        if coils===nothing
            _bunknown!(c,path,obj,"Unsupported connection configuration."); continue
        end
        # The native load builder treats a two-terminal delta as one coil.
        load && get(obj,"configuration","WYE")=="DELTA" && length(obj["terminal_map"])==2 && (coils=coils[1:1])
        unsupported_load=load && get(obj,"model","constant_power")!="constant_power"
        if source
            tm=string.(obj["terminal_map"]); mags=get(obj,"v_magnitude",[]); angles=get(obj,"v_angle",[])
            length(mags)==length(angles) && length(mags)<=length(tm) || throw(ArgumentError("$path source reference lengths disagree"))
            for k in eachindex(mags)
                mag=_bnumber(mags[k],"$path/v_magnitude"); _bnumber(angles[k],"$path/v_angle")
                mag>=0 || throw(ArgumentError("Negative source magnitude"))
                _bdeclare!(e,_bvoltage(c,bus,tm[k])[3],mag,mag,"$path/v_magnitude/$(k-1)")
                if iszero(angles[k])
                    z=_bvoltage(c,bus,tm[k])
                    refs=["$path/v_magnitude/$(k-1)","$path/v_angle/$(k-1)"]
                    _bupdate!(e,z[1],_bpoint(mag),:source_zero_angle,Int[],refs;force=true)
                    _bupdate!(e,z[2],_bpoint(0),:source_zero_angle,Int[],refs;force=true)
                end
            end
            all(t->(bus,t) in c.grounded,tm[length(mags)+1:end]) || throw(ArgumentError("$path omits an ungrounded voltage reference"))
        end
        if load && !unsupported_load
            for f in ("p_nom","q_nom")
                v=get(obj,f,[])
                v isa AbstractVector && length(v)==length(coils) || throw(ArgumentError("$path/$f must cover every load coil"))
            end
        end
        currents=NTuple{3,Int}[]; powers=NTuple{3,Int}[]
        for (k,(pos,neg)) in enumerate(coils)
            idx=("coil",string(k))
            v=_bconnection!(c,Symbol(family),id,idx,bus,pos,neg,path)
            i=_bcomplex!(e,Symbol(family),id,:current,idx,:A)
            s=_bpower!(c,Symbol(family),id,idx,v,i,path)
            push!(currents,i); push!(powers,s)
            sign=load ? -1 : 1
            _binject!(c,bus,pos,i,sign); _binject!(c,bus,neg,i,-sign)
            if load && !unsupported_load
                for (field,part) in (("p_nom",1),("q_nom",2))
                    v0=_bnumber(obj[field][k],"$path/$field")
                    _bdeclare!(e,s[part],v0,v0,"$path/$field/$(k-1)")
                end
            end
        end
        if !load
            for (field,part,lower) in (("p_min",1,true),("p_max",1,false),("q_min",2,true),("q_max",2,false))
                _blimits!(e,[s[part] for s in powers],obj,field,path;lower)
            end
        end
        if generator
            _blimits!(e,[s[3] for s in powers],obj,"s_max",path)
            limits=get(obj,"i_max",nothing); cfg=get(obj,"configuration","WYE")
            if limits!==nothing
                limits isa AbstractVector || throw(ArgumentError("$path/i_max must be an array"))
                if cfg=="DELTA"
                    _blimits!(e,[i[3] for i in currents],obj,"i_max",path)
                elseif length(coils)==1
                    for (k,v) in enumerate(limits)
                        val=_bnumber(v,"$path/i_max";infinite=true)
                        val>=0 || throw(ArgumentError("Negative current limit"))
                        _bdeclare!(e,currents[1][3],0,val,"$path/i_max/$(k-1)")
                    end
                else
                    n=length(coils)
                    length(limits)<=n+1 || throw(ArgumentError("$path/i_max has too many entries"))
                    _blimits!(e,[i[3] for i in currents],Dict("i_max"=>limits[1:min(n,end)]),"i_max",path)
                    if length(limits)==n+1 && coils[1][2]!==nothing
                        ret=_bcomplex!(e,:generator,id,:current,("neutral",string(coils[1][2])),:A)
                        _brow!(e,vcat(Any[(ret,1)],[(i,1) for i in currents]),[path])
                        val=_bnumber(limits[end],"$path/i_max";infinite=true)
                        val>=0 || throw(ArgumentError("Negative current limit"))
                        _bdeclare!(e,ret[3],0,val,"$path/i_max/$n")
                    end
                end
            end
        end
        reason=source ? "Source magnitude references, exact zero angles and native grounding; nonzero phase angles are relaxed." :
            unsupported_load ? "Voltage-dependent power law omitted; connection currents and power identities retained." :
            "Connection-aware complex power identity and native declared limits."
        _bcover!(e,path,source || unsupported_load ? :relaxed : :supported,reason)
    end
end

function _bshunts!(c)
    e=c.engine
    for id in _bsorted(_bcollection(c.net,"shunt"))
        obj=c.net["shunt"][id]; path=_bpath("shunt",id)
        bus=string(obj["bus"]); tm=string.(obj["terminal_map"]); n=length(tm)
        G=_bmatrix(obj,"G_",n,path); B=_bmatrix(obj,"B_",n,path)
        for (k,t) in enumerate(tm)
            i=_bcomplex!(e,:shunt,id,:current,(t,),:A)
            row=Any[(i,1)]
            append!(row,[(_bvoltage(c,bus,tm[j]),-complex(G[k,j],B[k,j])) for j in 1:n])
            _brow!(e,row,[path]); _binject!(c,bus,t,i,-1)
            _bpower!(c,:shunt,id,(t,),_bvoltage(c,bus,t),i,path)
        end
        _bcover!(e,path,:supported,"Full terminal admittance matrix; nominal ratings are not constraints.")
    end
end

function _bunknown!(c,path,obj,reason)
    # A recognized bus reference narrows the set of unknown injection locations.
    # Unrecognized extension layouts conservatively retain a free injection at
    # every AC terminal, preventing a dropped device from tightening other ones.
    buses=Set{String}()
    function visit(x)
        x isa AbstractDict || return
        for (k,v) in x
            if k in ("bus","bus_from","bus_to") && v isa AbstractString && haskey(c.net["bus"],v)
                push!(buses,String(v))
            elseif v isa AbstractDict
                visit(v)
            elseif v isa AbstractVector
                foreach(visit,v)
            end
        end
    end
    visit(obj)
    isempty(buses) && union!(buses,keys(get(c.net,"bus",Dict())))
    for (bus,t) in sort!(collect(keys(c.voltages)))
        bus in buses || continue
        i=_bcomplex!(c.engine,:unsupported,path,:current,(bus,t),:A)
        _binject!(c,bus,t,i,1)
    end
    _bcover!(c.engine,path,:unsupported,reason*" Free injections retained on connected buses (all AC buses if unresolved).")
end

function _binventory!(c,original)
    e=c.engine
    supported=Set(["bus","line","switch","load","generator","voltage_source","shunt"])
    data=Set(["linecode","wire_data","line_geometry","control_profile","time_series"])
    metadata=Set(["name","meta","extras","terminal_conventions"])
    for family in _bsorted(original)
        startswith(family,"_") && continue
        family in supported && continue
        value=original[family]
        if family in metadata
            _bcover!(e,_bpath(family),:metadata,"Not an operating state.")
        elseif family in data
            for id in _bsorted(_bcollection(original,family))
                _bcover!(e,_bpath(family,id),:data,"Input dependency; no independent operating-state domain.")
            end
        elseif family=="transformer"
            for subtype in _bsorted(value), id in _bsorted(value[subtype])
                _bunknown!(c,_bpath(family,subtype,id),value[subtype][id],"Transformer subtype not yet supported.")
            end
        elseif family in ("dc_bus","dc_branch","dc_grounding","dc_load","dc_source")
            # DC-only objects have no direct AC injection. Any AC/DC converter
            # is independently accounted for by its unsupported IBR port.
            for id in _bsorted(value)
                _bcover!(e,_bpath(family,id),:unsupported,"DC quantities not yet supported; no direct AC port.")
            end
        elseif value isa AbstractDict
            for id in _bsorted(value)
                _bunknown!(c,_bpath(family,id),value[id],"Device family not yet supported.")
            end
        else
            _bunknown!(c,_bpath(family),value,"Unknown network extension.")
        end
    end
end

"""
    propagate_bounds(input; scope=BoundScope(), max_steps=100_000,
                     time_limit=Inf, atol=0.0, rtol=1e-10, validate=true)

Derive outward-valid SI enclosures from a BMOPF dictionary, JSON string or path.
The input is copied and normalized. A referenced time series is resolved through
BMOPFTools for the selected snapshot. Scope is native declared-limit feasibility,
without an objective cutoff. Unsupported physics is reported and relaxed safely.

`atol`/`rtol` only skip small improvements; they never round a bound inward.
`max_steps` and `time_limit` limit propagation work after compilation. Schema
validation runs when a recognized schema is declared; its availability is recorded.
Malformed supported connections throw `ArgumentError`; conflicting valid domains
produce a report with `status == :conflict`.
"""
function propagate_bounds(input;scope::BoundScope=BoundScope(),max_steps::Integer=100_000,
                          time_limit::Real=Inf,atol::Real=0.0,rtol::Real=1e-10,validate::Bool=true)
    max_steps>=0 || throw(ArgumentError("max_steps must be nonnegative"))
    time_limit>=0 && !isnan(time_limit) || throw(ArgumentError("time_limit must be nonnegative"))
    all(x->isfinite(x)&&x>=0,(atol,rtol)) || throw(ArgumentError("Tolerances must be finite and nonnegative"))
    original=_binput(input)
    validation=Dict{String,Any}("requested"=>validate,"adapter"=>_BOUND_ADAPTER_VERSION,
        "arithmetic"=>"outward_binary64_exact_rational_operations","input_semantics"=>"normalized_binary64")
    normalization=deepcopy(get(original,"_meta",Dict{String,Any}()))
    delete!(normalization,"parsed_at")
    validation["normalization"]=normalization
    if validate
        findings=BMOPFTools.Finding[]
        merge!(validation,BMOPFTools.schema_check(original,findings))
        get(validation,"jsonschema_valid",true) || throw(ArgumentError("BMOPF schema validation failed: "*join(string.(findings),"; ")))
    end
    net=BMOPFTools.get_snapshot(original,scope.snapshot)
    e=_BEngine(Float64(atol),Float64(rtol))
    c=_BNetwork(e,net,Dict{Tuple{String,String},NTuple{3,Int}}(),
        Dict{Tuple{String,String},Vector{Tuple{NTuple{3,Int},Int}}}(),
        Set{Tuple{String,String}}(),Set{Tuple{String,String}}())
    _bbuses!(c); _blines!(c); _bswitches!(c)
    for family in ("load","generator","voltage_source"); _bdevices!(c,family); end
    _bshunts!(c); _binventory!(c,original)
    for node in sort!(collect(keys(c.kcl)))
        _brow!(e,c.kcl[node],[_bpath("bus",node...,"kcl")])
    end
    status,steps=_brun!(e;max_steps,time_limit)
    for conflict in e.conflicts
        conflict["quantity"]=_bkeydict(e.keys[conflict["quantity"]])
    end
    domains=Dict(k=>BoundDomain(e.domains[i].lo,e.domains[i].hi,e.units[i],e.lower[i],e.upper[i]) for (i,k) in enumerate(e.keys))
    BoundReport(scope,_bfingerprint(original,scope),net,domains,e.evidence,
        sort!(e.objects;by=x->x.object),e.conflicts,status,steps,length(e.relations),validation)
end
