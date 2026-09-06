using JSON3, JSONSchema

@testset "Generator data: schema and physical round trips under unbalance" begin
    net=inv_grid3_unbal()
    emf=[235,225,232].*cis.([0.02,-2.13,2.05])
    for source in (false,true),mode in (:none,:fixed_phasor,:fixed_magnitudes,:common_magnitude,:phase_magnitudes)
        law=mode in (:fixed_phasor,:fixed_magnitudes) ? GeneratorVoltageLaw(mode;phasor=emf) :
            mode==:none ? GeneratorVoltageLaw() : GeneratorVoltageLaw(mode;magnitude_min=180.0,magnitude_max=280.0)
        ctl=mode==:none ? GeneratorControl(p=[400.0,700.0,200.0],q=[40.0,-20.0,60.0]) :
            mode==:fixed_phasor ? GeneratorControl() : mode==:fixed_magnitudes ? GeneratorControl(p=1300.0) :
            GeneratorControl(p=1300.0,q=100.0)
        cap=GeneratorCapability(power_location=:internal,p_min=-1e5,p_max=[1e5,1e5,1e5],q_min=-1e5,q_max=1e5,
            s_max=3e5,i_max=100.0,terminal_i_max=[150.0,140.0,130.0,120.0],
            voltage_min=100.0,voltage_max=400.0,earth_i_max=source ? 100.0 : nothing,
            sequence=GeneratorSequenceLimits(location=:internal,voltage_min=[0.0,0.0,0.0],voltage_max=[500.0,500.0,500.0],
                current_min=[0.0,0.0,0.0],current_max=[200.0,200.0,200.0]))
        args=(id="g",bus="poc",voltage=law,control=ctl,capability=cap,cost=0.2)
        d=source ? SourceGenerator(;args...,impedance=fill(0.4+0.8im,4),grounding=0.8+0.1im) :
            GeneralizedGenerator(;args...,connections=_GG_PORTS,impedance=0.4+0.8im)
        envelope=generator_data([d])
        @test JSONSchema.validate(JSONSchema.Schema(generator_data_schema()),envelope)===nothing
        parsed=read_generator_data(JSON3.write(envelope);from_string=true,net)
        @test length(parsed.devices)==1
        restored=only(parsed.devices)
        @test restored.id=="$(PowerOptLab._gg_family(d)):g"
        @test restored.capability.terminal_i_max==[150.0,140.0,130.0,120.0]
        @test restored.capability.i_max==fill(100.0,3)
        a=_gg_solve(net,d;per_unit=true,s_base=1e4)
        b=_gg_solve(net,restored;per_unit=true,s_base=1e4)
        @test a.solve.publishable && b.solve.publishable
        # A free phase-magnitude law can have multiple valid dispatches; compare
        # prescribed aggregate quantities and physical residuals, not optimizer tie-breaking.
        if mode!=:phase_magnitudes
            @test a.devices[d.id].terminal_voltage ≈ b.devices[restored.id].terminal_voltage atol=1e-4
            @test a.devices[d.id].p ≈ b.devices[restored.id].p atol=1e-3
        end
        @test abs(b.devices[restored.id].power_balance_error)<1e-4
        @test haskey(generator_data(parsed)[String(PowerOptLab._gg_family(d))],"g")
    end
end

@testset "Generator data: topology, ID and file preservation" begin
    # A generic forest whose data conductor order differs from runtime order.
    record=Dict{String,Any}("bus"=>"poc","terminal_map"=>["n","c","a","b"],
        "configuration"=>"PORTS","port_map"=>[[3,4],[2,1]],"voltage_model"=>"NONE","cost_total"=>0,
        "i_max"=>[4,3,1,2],"i_port_max"=>[9,8],"p_set"=>[500,300],"q_set"=>[50,20],"power_setpoint_location"=>"POC")
    d=generator_from_data("ports",record;net=inv_grid3_unbal())
    @test PowerOptLab._gg_terminals(d)==["a","c","b","n"]
    @test d.capability.terminal_i_max==[1,3,2,4]
    @test d.capability.i_max==[9,8]
    @test generator_data(d)["port_map"]==[[1,3],[2,4]]
    # Source and ordinary generator may have the same original data ID.
    ds=[GeneralizedGenerator(id="same:id",bus="poc",connections=_GG_DELTA,
            control=GeneratorControl(p=[600.0,300.0,900.0],q=[50.0,100.0,150.0])),
        SourceGenerator(id="same:id",bus="poc",grounding=nothing)]
    mktempdir() do dir
        path=joinpath(dir,"generators.json")
        @test write_generator_data(path,ds)==path
        parsed=read_generator_data(path;net=inv_grid3_unbal())
        @test length(unique(d.id for d in parsed.devices))==2
        @test parsed.identifiers["source_generator:same:id"]==(:source_generator,"same:id")
        @test generator_data(parsed)==generator_data(ds)
        @test generator_data(parsed.devices[1])["configuration"]=="DELTA"
        r=_gg_solve(inv_grid3_unbal(),parsed.devices[1];per_unit=true,s_base=1e4)
        @test r.solve.publishable
        @test r.devices[parsed.devices[1].id].p ≈ 1800 atol=1e-4
        # Duplicate IDs within a family are rejected before opening the file.
        before=read(path,String)
        @test_throws ArgumentError write_generator_data(path,[ds[1],ds[1]])
        @test read(path,String)==before
    end
    for (ports,config,angles) in (([("a","n")],"SINGLE_PHASE",[0.0]),
            ([("a","n"),("b","n")],"WYE",[0.0,pi]),(_GG_DELTA,"DELTA",_GG_PHI))
        d=GeneralizedGenerator(id="g",bus="poc",connections=ports,impedance=0.4+0.8im,
            voltage=GeneratorVoltageLaw(:phase_magnitudes;angles,magnitude_min=100.0,magnitude_max=500.0))
        r=generator_data(d)
        @test r["configuration"]==config
        @test generator_from_data("g",r).connections==ports
    end
    for ground in (nothing,0.0,0.8+0.2im)
        d=SourceGenerator(id="s",bus="poc",grounding=ground)
        @test generator_from_data("s",generator_data(d);family=:source_generator).grounding==ground
    end
    common=Dict{String,Any}("bus"=>"poc","terminal_map"=>["a","b","n"],"configuration"=>"WYE",
        "voltage_model"=>"COMMON_MAGNITUDE","angle_offsets"=>[0,pi],"e_min"=>[100,110],"e_max"=>[130,120],"cost_total"=>0)
    d=generator_from_data("c",common)
    @test d.voltage.magnitude_min==110
    @test d.voltage.magnitude_max==120
end

@testset "Generator data rejects schema and semantic mistakes" begin
    base=generator_data(GeneralizedGenerator(id="g",bus="poc",connections=_GG_DELTA))
    delete!(base,"power_setpoint_location");delete!(base,"power_limit_location")
    changes=[Dict("unexpected"=>1),Dict("cost_total"=>"free"),Dict("p_set"=>[1,2,3]),
        Dict("p_set"=>[1,2,3],"p_total_set"=>6,"power_setpoint_location"=>"POC"),
        Dict("p_min"=>[1,2,3]),Dict("s_total_max"=>-1,"power_limit_location"=>"POC"),
        Dict("v_seq_max"=>[0,400,10]),Dict("v_seq_max"=>[0,400],"v_sequence_location"=>"POC"),
        Dict("v_target"=>400),Dict("v_target_measurement"=>"POSITIVE_SEQUENCE"),
        Dict("voltage_model"=>"FIXED_PHASOR"),Dict("angle_offsets"=>[0,1,2]),
        Dict("configuration"=>"DELTA","terminal_map"=>["a","b"]),
        Dict("i_max"=>[1,2]),Dict("i_port_max"=>[1,2]),Dict("i_port_max"=>Any[1,true,2]),
        Dict("R_series_4_4"=>1),Dict("R_series_1_1"=>-1),Dict("R_series_1_1"=>Inf),
        Dict("R_series_999999999999999999999999_1"=>1),Dict("port_map"=>[[1,2]]),
        Dict("ig_max"=>1),Dict("voltage_model"=>"PHASE_MAGNITUDES","angle_offsets"=>[0,1,2],"e_min"=>[0,1,1]),
        Dict("v_min"=>[400,400,400],"v_max"=>[300,300,300]),Dict("bus"=>"missing")]
    for change in changes
        @test_throws ArgumentError generator_from_data("g",merge(base,change);net=inv_grid3_unbal())
    end
    for key in ("bus","terminal_map","configuration","voltage_model","cost_total")
        r=copy(base);delete!(r,key)
        @test_throws ArgumentError generator_from_data("g",r)
    end
    source=generator_data(SourceGenerator(id="s",bus="poc",grounding=nothing))
    @test_throws ArgumentError generator_from_data("s",merge(source,Dict("r_ground"=>1));family=:source_generator)
    @test_throws ArgumentError generator_from_data("s",merge(source,Dict("configuration"=>"DELTA"));family=:source_generator)
    @test_throws ArgumentError read_generator_data(Dict("generator"=>Dict("g"=>base)))
    @test_throws ArgumentError read_generator_data(Dict("generalized_generator"=>Dict(""=>base)))
    @test_throws ArgumentError generator_from_data("g",base;family=:load)
end

@testset "Fixed phasor polar round trips at exact magnitude bounds" begin
    # Cartesian -> polar -> Cartesian can move a norm by a few ulps. A fixed
    # template at its exact capability boundary must survive that conversion.
    for mode in (:fixed_phasor,:fixed_magnitudes),theta in range(-3,3;length=61)
        e=230cis(theta)
        d=GeneralizedGenerator(id="g",bus="poc",connections=[("a","n")],
            voltage=GeneratorVoltageLaw(mode;phasor=[e],magnitude_min=[abs(e)],magnitude_max=[abs(e)]))
        restored=generator_from_data("g",generator_data(d))
        @test restored.voltage.phasor ≈ [e] atol=1e-12
    end
    d=GeneralizedGenerator(id="g",bus="poc",connections=[("a","n")],
        voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=[230cis(0.5)],magnitude_max=230.0-1e-8))
    @test_throws ArgumentError validate_device(d,(inv_grid3_unbal(),))
    for location in (:poc,:internal)
        d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_DELTA,impedance=0.8+1.2im,
            control=GeneratorControl(p=[600.0,300.0,900.0],q=300.0,power_location=location),
            capability=GeneratorCapability(power_location=location,p_min=[0.0,0.0,0.0],p_max=2000.0,
                q_min=-1000.0,q_max=[1000.0,1000.0,1000.0],s_max=[1500.0,1500.0,1500.0]))
        restored=generator_from_data("g",generator_data(d);net=inv_grid3_unbal())
        @test restored.control.power_location==location
        r=_gg_solve(inv_grid3_unbal(),restored;per_unit=true,s_base=1e4)
        @test r.solve.publishable
        g=r.devices["g"]
        @test (location==:poc ? g.p : g.p_internal) ≈ 1800 atol=1e-4
    end
end
