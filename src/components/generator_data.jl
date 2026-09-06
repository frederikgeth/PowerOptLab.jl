# Local experimental data contract. No upstream record is reinterpreted.
const _GG_SCHEMA_PATH = normpath(joinpath(@__DIR__,"..","..","data","schema","generators.schema.json"))

"""
    GeneratorDataSet

Parsed generator-extension records: `devices` can be passed to
`build_generator_model`; `identifiers` maps each runtime ID to `(family, data_id)`.
Envelope imports qualify every runtime ID as `family:data_id`, preserving IDs
that are repeated across the two data collections. Export retains original IDs.
The native network is supplied separately and is never modified by this codec.
"""
struct GeneratorDataSet
    devices::Vector{_GeneralizedDevice}
    identifiers::Dict{String,Tuple{Symbol,String}}
end

"""
    generator_data_schema()

Return the local experimental JSON Schema as an independent dictionary. Schema
validation checks record syntax; the importer additionally checks dimensions,
topology, finite SI values, passivity, bounds, voltage-law domains and, when a
network is supplied, bus/terminal membership. This is not an upstream release.
"""
generator_data_schema() = JSON3.read(read(_GG_SCHEMA_PATH,String),Dict{String,Any})

Base.include_dependency(_GG_SCHEMA_PATH)
const _GG_DATA_SCHEMA = JSONSchema.Schema(generator_data_schema())

function _gg_validate_data(data)
    problem=JSONSchema.validate(_GG_DATA_SCHEMA,data)
    problem===nothing || throw(ArgumentError("generator data schema: $problem"))
    nothing
end
_gg_family(d::GeneralizedGenerator)=:generalized_generator
_gg_family(d::SourceGenerator)=:source_generator
_gg_data_bound(r,key)=haskey(r,key) ? (r[key] isa AbstractVector ? Float64.(r[key]) : Float64(r[key])) : nothing
function _gg_data_power(r,key)
    a,b=split(key,"_")
    _gg_data_bound(r,haskey(r,key) ? key : "$(a)_total_$(b)")
end
_gg_location(r,key)=Symbol(lowercase(get(r,key,"POC")))
function _gg_data_net(d)
    Dict("bus"=>Dict(d.bus=>Dict("terminal_names"=>_gg_terminals(d))))
end

"""
    generator_from_data(id, record; family=:generalized_generator, net=nothing,
                        voltage_scale=230.0)

Validate and translate one flat specification record to a generator component.
`i_max` in data is the full conductor rating; `i_port_max` is the winding rating.
Generalized configurations are WYE, SINGLE_PHASE, DELTA and oriented PORTS;
sources use WYE with an explicit neutral. Omitted impedance entries are zero.
Pass `net` for bus/terminal validation; absence defers only network checks.
"""
function generator_from_data(id::AbstractString,record::AbstractDict;
        family::Symbol=:generalized_generator,net=nothing,voltage_scale::Real=230.0)
    family in (:generalized_generator,:source_generator) || throw(ArgumentError("unknown generator family $family"))
    r=Dict{String,Any}(String(k)=>v for (k,v) in record)
    _gg_validate_data(Dict(String(family)=>Dict(String(id)=>r)))
    ts=String.(r["terminal_map"]); m=length(ts); config=r["configuration"]
    ports=if config=="WYE"
        [(t,ts[end]) for t in ts[1:end-1]]
    elseif config=="SINGLE_PHASE"
        [(ts[1],ts[2])]
    elseif config=="DELTA"
        [(ts[1],ts[2]),(ts[2],ts[3]),(ts[3],ts[1])]
    else
        all(p->all(k->1<=k<=m,p) && p[1]!=p[2],r["port_map"]) || throw(ArgumentError("invalid port_map indices"))
        pairs=[(ts[Int(p[1])],ts[Int(p[2])]) for p in r["port_map"]]
        Set(vcat(first.(pairs),last.(pairs)))==Set(ts) || throw(ArgumentError("every terminal_map entry must participate"))
        pairs
    end
    n=length(ports); nz=family==:source_generator ? m : n
    # Dimension is determined by topology, never by the largest provided key.
    Z=zeros(ComplexF64,nz,nz)
    for (key,value) in r
        matchkey=match(r"^([RX])_series_([1-9][0-9]*)_([1-9][0-9]*)$",key)
        matchkey===nothing && continue
        row=tryparse(Int,matchkey[2]); col=tryparse(Int,matchkey[3])
        row!==nothing && col!==nothing && row<=nz && col<=nz || throw(ArgumentError("impedance index out of range: $key"))
        Z[row,col]+=(matchkey[1]=="R" ? 1 : im)*Float64(value)
    end
    # Enforce all port vector lengths even when a common law collapses bounds.
    portfields=("v_magnitude","v_angle","angle_offsets","e_min","e_max","v_min","v_max",
        "p_set","q_set","p_min","p_max","q_min","q_max","s_max","i_port_max")
    for key in portfields
        haskey(r,key) && length(r[key])!=n && throw(ArgumentError("$key needs $n entries"))
    end
    haskey(r,"i_max") && length(r["i_max"])!=m && throw(ArgumentError("i_max needs $m conductor entries"))
    mode=Symbol(lowercase(r["voltage_model"]))
    lo,hi=_gg_data_bound(r,"e_min"),_gg_data_bound(r,"e_max")
    mode==:common_magnitude && begin
        lo=lo===nothing ? nothing : maximum(lo)
        hi=hi===nothing ? nothing : minimum(hi)
    end
    law=GeneratorVoltageLaw(mode;phasor=haskey(r,"v_magnitude") ? Float64.(r["v_magnitude"]).*cis.(Float64.(r["v_angle"])) : ComplexF64[],
        angles=haskey(r,"angle_offsets") ? Float64.(r["angle_offsets"]) : Float64[],magnitude_min=lo,magnitude_max=hi)
    ctl=GeneratorControl(p=_gg_data_power(r,"p_set"),q=_gg_data_power(r,"q_set"),
        power_location=_gg_location(r,"power_setpoint_location"),voltage_target=_gg_data_bound(r,"v_target"),
        voltage_metric=get(r,"v_target_measurement","POSITIVE_SEQUENCE")=="PORT" ? :phase : :positive_sequence)
    seqkeys=("v_seq_min","v_seq_max","i_seq_min","i_seq_max")
    seq=any(k->haskey(r,k),seqkeys) ? GeneratorSequenceLimits(location=_gg_location(r,"v_sequence_location"),
        voltage_min=_gg_data_bound(r,"v_seq_min"),voltage_max=_gg_data_bound(r,"v_seq_max"),
        current_min=_gg_data_bound(r,"i_seq_min"),current_max=_gg_data_bound(r,"i_seq_max")) : nothing
    # Runtime terminal order follows outgoing labels then previously unseen returns.
    order=unique(vcat(first.(ports),last.(ports)))
    ratings=_gg_data_bound(r,"i_max")
    ratings===nothing || (ratings=ratings[[findfirst(==(t),ts) for t in order]])
    cap=GeneratorCapability(p_min=_gg_data_power(r,"p_min"),p_max=_gg_data_power(r,"p_max"),
        q_min=_gg_data_power(r,"q_min"),q_max=_gg_data_power(r,"q_max"),s_max=_gg_data_power(r,"s_max"),
        power_location=_gg_location(r,"power_limit_location"),i_max=_gg_data_bound(r,"i_port_max"),
        terminal_i_max=ratings,voltage_min=_gg_data_bound(r,"v_min"),voltage_max=_gg_data_bound(r,"v_max"),
        earth_i_max=_gg_data_bound(r,"ig_max"),sequence=seq)
    args=(id=String(id),bus=String(r["bus"]),impedance=Z,voltage=law,control=ctl,capability=cap,
        voltage_scale=Float64(voltage_scale),cost=Float64(r["cost_total"]))
    d=if family==:source_generator
        ground=r["grounding_model"]
        zg=ground=="OPEN" ? nothing : ground=="IDEAL" ? 0.0 : complex(Float64(get(r,"r_ground",0)),Float64(get(r,"x_ground",0)))
        SourceGenerator(;args...,phase_terminals=ts[1:end-1],neutral=ts[end],grounding=zg)
    else
        GeneralizedGenerator(;args...,connections=ports)
    end
    validate_device(d,(net===nothing ? _gg_data_net(d) : net,))
    d
end

"""
    read_generator_data(input; from_string=false, net=nothing, voltage_scale=230.0)

Read a file, JSON string (`from_string=true`), or parsed dictionary containing
only the `generalized_generator` / `source_generator` collections. Return a
`GeneratorDataSet`. All numeric inputs are SI RMS. Unknown fields and ambiguous
power representations are rejected. No network, ID or record is silently dropped.
"""
function read_generator_data(input;from_string::Bool=false,net=nothing,voltage_scale::Real=230.0)
    data=input isa AbstractDict ? input : JSON3.read(from_string ? input : read(input,String),Dict{String,Any})
    _gg_validate_data(data)
    devices=_GeneralizedDevice[]; identifiers=Dict{String,Tuple{Symbol,String}}()
    for family in (:generalized_generator,:source_generator)
        records=get(data,String(family),Dict())
        for id in sort!(collect(keys(records)))
            runtime_id="$(family):$id"
            d=generator_from_data(runtime_id,records[id];family,net,voltage_scale)
            push!(devices,d); identifiers[runtime_id]=(family,String(id))
        end
    end
    GeneratorDataSet(devices,identifiers)
end

function _gg_put_bound!(r,key,value,n;aggregate=false)
    value===nothing && return
    if aggregate && value isa Real
        a,b=split(key,"_");r["$(a)_total_$(b)"]=Float64(value)
    else
        r[key]=_gg_vector(value,n)
    end
end

"""
    generator_data(device)
    generator_data(dataset)
    generator_data(devices::AbstractVector)

Export canonical flat SI records (one device) or an extension envelope. Scalars
that apply to each entry become explicit vectors, matrix zeros are omitted, and
relative-angle defaults become explicit. Canonicalization preserves the physical
model, not source JSON formatting. Numerical `voltage_scale` is intentionally
not serialized. Dataset exports preserve the original family and data ID.
"""
function generator_data(d::_GeneralizedDevice)
    validate_device(d,(_gg_data_net(d),))
    ts=_gg_terminals(d); pairs=_gg_connections(d); n=_gg_n(d); m=length(ts)
    config=d isa SourceGenerator ? "WYE" : _gg_delta(d) ? "DELTA" : n==1 ? "SINGLE_PHASE" :
        length(unique(last.(pairs)))==1 ? "WYE" : "PORTS"
    r=Dict{String,Any}("bus"=>d.bus,"terminal_map"=>ts,"configuration"=>config,
        "voltage_model"=>uppercase(String(d.voltage.mode)),"cost_total"=>d.cost)
    if config=="PORTS"
        r["port_map"]=[[findfirst(==(a),ts),findfirst(==(b),ts)] for (a,b) in pairs]
    end
    Z=_gg_matrix(d.impedance,d isa SourceGenerator ? m : n)
    for j in axes(Z,2),i in axes(Z,1)
        iszero(real(Z[i,j])) || (r["R_series_$(i)_$(j)"]=real(Z[i,j]))
        iszero(imag(Z[i,j])) || (r["X_series_$(i)_$(j)"]=imag(Z[i,j]))
    end
    law=d.voltage
    if law.mode in (:fixed_phasor,:fixed_magnitudes)
        r["v_magnitude"]=abs.(law.phasor); r["v_angle"]=angle.(law.phasor)
    elseif law.mode in (:common_magnitude,:phase_magnitudes)
        r["angle_offsets"]=_gg_angles(law,n)
    end
    _gg_put_bound!(r,"e_min",law.magnitude_min,n);_gg_put_bound!(r,"e_max",law.magnitude_max,n)
    c=d.control; cap=d.capability
    _gg_put_bound!(r,"p_set",c.p,n;aggregate=true);_gg_put_bound!(r,"q_set",c.q,n;aggregate=true)
    # Preserve even inactive declared locations; they have defined enum semantics.
    r["power_setpoint_location"]=uppercase(String(c.power_location))
    r["power_limit_location"]=uppercase(String(cap.power_location))
    if c.voltage_target!==nothing
        r["v_target"]=c.voltage_target;r["v_target_measurement"]=c.voltage_metric==:phase ? "PORT" : "POSITIVE_SEQUENCE"
    end
    for field in (:p_min,:p_max,:q_min,:q_max,:s_max)
        _gg_put_bound!(r,String(field),getproperty(cap,field),n;aggregate=true)
    end
    for (key,value,count) in (("i_port_max",cap.i_max,n),("i_max",cap.terminal_i_max,m),
            ("v_min",cap.voltage_min,n),("v_max",cap.voltage_max,n))
        _gg_put_bound!(r,key,value,count)
    end
    seq=cap.sequence
    if seq!==nothing
        r["v_sequence_location"]=uppercase(String(seq.location))
        for (key,value) in (("v_seq_min",seq.voltage_min),("v_seq_max",seq.voltage_max),
                ("i_seq_min",seq.current_min),("i_seq_max",seq.current_max))
            _gg_put_bound!(r,key,value,3)
        end
    end
    if d isa SourceGenerator
        r["grounding_model"]=d.grounding===nothing ? "OPEN" : iszero(d.grounding) ? "IDEAL" : "IMPEDANCE"
        if d.grounding!==nothing && !iszero(d.grounding)
            r["r_ground"]=real(d.grounding);r["x_ground"]=imag(d.grounding)
        end
        cap.earth_i_max===nothing || (r["ig_max"]=cap.earth_i_max)
    end
    _gg_validate_data(Dict(String(_gg_family(d))=>Dict(d.id=>r)))
    r
end
function generator_data(data::GeneratorDataSet)
    out=Dict{String,Any}("generalized_generator"=>Dict{String,Any}(),"source_generator"=>Dict{String,Any}())
    length(unique(d.id for d in data.devices))==length(data.devices) || throw(ArgumentError("duplicate runtime ID"))
    for d in data.devices
        haskey(data.identifiers,d.id) || throw(ArgumentError("missing identity for $(d.id)"))
        family,id=data.identifiers[d.id]
        family==_gg_family(d) || throw(ArgumentError("family mismatch for $(d.id)"))
        records=out[String(family)]
        haskey(records,id) && throw(ArgumentError("duplicate data ID $family:$id"))
        records[id]=generator_data(d)
    end
    _gg_validate_data(out)
    out
end
function generator_data(devices::AbstractVector)
    ds=_GeneralizedDevice[devices...]
    # IDs need only be unique within each data family on export.
    out=Dict{String,Any}("generalized_generator"=>Dict{String,Any}(),"source_generator"=>Dict{String,Any}())
    for d in ds
        records=out[String(_gg_family(d))]
        haskey(records,d.id) && throw(ArgumentError("duplicate data ID $(_gg_family(d)):$(d.id)"))
        records[d.id]=generator_data(d)
    end
    out
end

"""
    write_generator_data(path, devices_or_dataset)

Validate and write a generator-extension JSON envelope. Validation completes
before opening the destination. Returns `path`; no native network data is written.
"""
function write_generator_data(path::AbstractString,data)
    envelope=generator_data(data)
    text=JSON3.write(envelope)
    open(path,"w") do io; JSON3.pretty(io,text); println(io); end
    path
end
