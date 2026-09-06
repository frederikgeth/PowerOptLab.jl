using BMOPFTools, LinearAlgebra

const _GG_PORTS = [("a","n"),("b","n"),("c","n")]
const _GG_PHI = [0.0,-2pi/3,2pi/3]
const _GG_E = 230 .* cis.(_GG_PHI)
const _GG_OPTIONS = Pair{String,Any}["tol"=>1e-9,"constr_viol_tol"=>1e-9,"bound_relax_factor"=>0.0,"max_iter"=>1200]

function _gg_fixture(;grounded=true)
    parse_bmopf("""{"bus":{"poc":{"terminal_names":["a","b","c","n"],
        "perfectly_grounded_terminals":$(grounded ? "[\"n\"]" : "[]")}}}""";from_string=true)
end

function _gg_imposed_current(currents)
    (ctx,handles) -> begin
        bs=BMOPFTools.opf_coordinate_bases(ctx,"poc")
        for (phase,i) in zip(["a","b","c"],currents)
            BMOPFTools.add_terminal_injection!(ctx,"poc",phase,-real(i)/bs.current,-imag(i)/bs.current)
        end
        ret=sum(currents)
        BMOPFTools.add_terminal_injection!(ctx,"poc","n",real(ret)/bs.current,imag(ret)/bs.current)
    end
end
function _gg_solve(net,d;kwargs...)
    solve_generator_opf(net,[d];solver_options=_GG_OPTIONS,kwargs...)
end

@testset "Every capability field under unbalance" begin
    currents=ComplexF64[10+2im,-3-5im,-2+4im]
    eph=[235,225,232].*cis.([0.02,-2.13,2.05])
    law=GeneratorVoltageLaw(:fixed_phasor;phasor=eph,magnitude_min=210.0,magnitude_max=250.0)
    Z=Matrix(Diagonal(ComplexF64[0.1+0.2im,0.15+0.1im,0.12+0.3im,0.2]))
    Z[1,4]=Z[4,1]=0.02+0.01im
    make(cap)=SourceGenerator(id="s",bus="poc",impedance=Z,grounding=0.8,voltage=law,capability=cap)
    run(cap;pu=true)=_gg_solve(_gg_fixture(),make(cap);per_unit=pu,s_base=1e4,
        generator_hook! = _gg_imposed_current(currents))
    g=run(GeneratorCapability()).devices["s"]
    @test g.solve.publishable
    # Four-conductor mutual coupling with separate neutral/earth paths.
    expected_n=(-0.8sum(currents)-sum(Z[4,k]*currents[k] for k in 1:3))/(Z[4,4]+0.8)
    @test g.terminal_current[end] ≈ expected_n atol=1e-6
    @test g.emf.+g.star_voltage.-g.terminal_voltage[1:3] ≈ (Z*g.terminal_current)[1:3] atol=1e-6
    for location in (:poc,:internal)
        sv=location==:poc ? g.voltage_sequence : g.emf_sequence
        pp=location==:poc ? g.port_voltage.*conj.(g.port_current) : g.emf.*conj.(g.port_current)
        sp=location==:poc ? complex(g.p,g.q) : complex(g.p_internal,g.q_internal)
        seq=GeneratorSequenceLimits(location=location,
            voltage_min=max.(0,abs.(sv).-0.1),voltage_max=abs.(sv).+0.1,
            current_min=max.(0,abs.(g.current_sequence).-0.1),current_max=abs.(g.current_sequence).+0.1)
        for phasepowers in (false,true)
            p=phasepowers ? real.(pp) : real(sp); q=phasepowers ? imag.(pp) : imag(sp)
            cap=GeneratorCapability(power_location=location,p_min=p.-0.1,p_max=p.+0.1,
                q_min=q.-0.1,q_max=q.+0.1,s_max=(phasepowers ? abs.(pp) : abs(sp)).+1,
                i_max=abs.(currents).+0.1,terminal_i_max=abs.(g.terminal_current).+0.1,
                voltage_min=abs.(g.port_voltage).-0.1,voltage_max=abs.(g.port_voltage).+0.1,
                earth_i_max=abs(g.earth_current)+0.1,sequence=seq)
            for pu in (false,true)
                r=run(cap;pu)
                @test r.solve.publishable
                @test r.devices["s"].terminal_current ≈ g.terminal_current atol=1e-5
                @test r.devices["s"].p ≈ g.p atol=1e-4
            end
        end
    end
    # Each limit can independently exclude the otherwise unique physical point.
    badcaps=[GeneratorCapability(p_min=g.p+10),GeneratorCapability(p_max=g.p-10),
        GeneratorCapability(q_min=g.q+10),GeneratorCapability(q_max=g.q-10),
        GeneratorCapability(s_max=0.95abs(complex(g.p,g.q))),
        GeneratorCapability(i_max=0.95abs.(currents)),
        GeneratorCapability(terminal_i_max=0.95abs.(g.terminal_current)),
        GeneratorCapability(earth_i_max=0.95abs(g.earth_current)),
        GeneratorCapability(voltage_min=abs.(g.port_voltage).+1),
        GeneratorCapability(voltage_max=abs.(g.port_voltage).-1),
        GeneratorCapability(sequence=GeneratorSequenceLimits(voltage_min=abs.(g.voltage_sequence).+0.1)),
        GeneratorCapability(sequence=GeneratorSequenceLimits(voltage_max=0.95abs.(g.voltage_sequence))),
        GeneratorCapability(sequence=GeneratorSequenceLimits(current_min=abs.(g.current_sequence).+0.1)),
        GeneratorCapability(sequence=GeneratorSequenceLimits(current_max=0.95abs.(g.current_sequence)))]
    for cap in badcaps
        @test !run(cap).solve.publishable
    end
end

@testset "Power locations, fixed-zero limits and physical units" begin
    for location in (:poc,:internal), pu in (true,false), base in (1e4,1e6)
        d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,impedance=0.4+0.2im,
            control=GeneratorControl(p=[400.0,700.0,200.0],q=[40.0,-20.0,60.0],power_location=location),
            capability=GeneratorCapability(power_location=location,p_min=[390.0,690.0,190.0],p_max=1400.0,
                q_min=-100.0,q_max=[50.0,-10.0,70.0],s_max=[500.0,800.0,300.0]))
        r=_gg_solve(inv_grid3_unbal(),d;per_unit=pu,s_base=base)
        @test r.solve.publishable
        g=r.devices["g"]
        @test (location==:poc ? g.p : g.p_internal) ≈ 1300 atol=1e-4
        @test g.p_internal>g.p
    end
    # One dead port among live unbalanced ports. S=0 gives P=Q=0 rows;
    # I=0 gives fixed current components, without a zero-gradient norm cone.
    for cap in (GeneratorCapability(i_max=[0.0,20.0,20.0]),GeneratorCapability(s_max=[0.0,3000.0,3000.0]))
        d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,capability=cap)
        hook! = (ctx,h) -> begin
            m=BMOPFTools.opf_model(ctx)
            for k in 2:3
                JuMP.fix(h["g"].ir[k],k/h["g"].bases.current;force=true)
                JuMP.fix(h["g"].ii[k],-k/h["g"].bases.current;force=true)
            end
        end
        r=_gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4,generator_hook! = hook!)
        @test r.solve.publishable
        @test abs(r.devices["g"].port_current[1])<1e-6
    end
end

@testset "Connections and stateless multi-period integration" begin
    # Single phase-to-phase component inside an unbalanced three-phase network.
    d=GeneralizedGenerator(id="ll",bus="poc",connections=[("a","b")],impedance=0.3+0.2im,
        voltage=GeneratorVoltageLaw(:none),control=GeneratorControl(p=1000.0,q=100.0))
    r=_gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4)
    @test r.solve.publishable
    @test r.devices["ll"].terminal_current[1] ≈ -r.devices["ll"].terminal_current[2]
    @test r.devices["ll"].voltage_sequence===nothing
    @test abs(r.devices["ll"].port_voltage[1])>300
    # Split phase has explicit 180-degree offsets and unequal leg loading.
    net=inv_grid3_unbal()
    for mode in (:fixed_phasor,:fixed_magnitudes,:common_magnitude,:phase_magnitudes)
        law=mode in (:fixed_phasor,:fixed_magnitudes) ?
            GeneratorVoltageLaw(mode;phasor=ComplexF64[230,-230]) :
            GeneratorVoltageLaw(mode;angles=[0.0,pi],magnitude_min=180.0,magnitude_max=290.0)
        # A split-phase grid is materialized explicitly; the c lateral remains
        # present, demonstrating mixed conductor widths in one native network.
        n=deepcopy(net); n["voltage_source"]["vs"]["v_angle"]=[0.03,pi+0.03,2.0]
        ctl=mode==:fixed_phasor ? GeneratorControl() : mode==:fixed_magnitudes ? GeneratorControl(p=1300.0) :
            mode==:common_magnitude ? GeneratorControl(p=1300.0,q=100.0) : GeneratorControl(p=1300.0,q=[80.0,20.0])
        sd=GeneralizedGenerator(id="split",bus="poc",connections=[("a","n"),("b","n")],
            impedance=[0.4+0.7im,0.5+0.8im],voltage=law,control=ctl)
        rr=_gg_solve(n,sd;per_unit=true,s_base=1e4)
        @test rr.solve.publishable
        g=rr.devices["split"]
        @test abs(abs(g.port_voltage[1])-abs(g.port_voltage[2]))>10
        @test abs(angle(-g.emf[2]/g.emf[1]))<1e-6
        @test abs(sum(g.terminal_current))<1e-6
    end
    one=GeneralizedGenerator(id="one",bus="poc",connections=[("a","n")],impedance=0.3+0.5im,
        voltage=GeneratorVoltageLaw(:common_magnitude;magnitude_min=180.0,magnitude_max=280.0),
        control=GeneratorControl(p=700.0,q=50.0),cost=0.15)
    multi=solve_multiperiod_opf([inv_grid3_unbal(),inv_grid3_unbal()],[one];
        time_grid=TimeGrid([0.5,1.5]),s_base=1e4,solver_options=_GG_OPTIONS)
    @test multi.solve.publishable
    @test multi.objective ≈ 0.21 atol=1e-5
    @test multi.snapshots[1]["custom_injection"]["p"] ≈ 700 atol=1e-5
    @test length(multi.dispatch["one"].snapshots)==2
end

@testset "Native ownership, costs, bounds and non-mutation" begin
    net=inv_grid3_unbal()
    net["generator"]=Dict("native"=>Dict("bus"=>"poc","terminal_map"=>["a","b","c","n"],
        "configuration"=>"WYE","p_min"=>[400.0,200.0,100.0],"p_max"=>[400.0,200.0,100.0],
        "q_min"=>[10.0,20.0,-10.0],"q_max"=>[10.0,20.0,-10.0],"cost"=>[0.1,0.1,0.1]))
    d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,
        control=GeneratorControl(p=[400.0,200.0,100.0],q=[10.0,20.0,-10.0]),cost=0.1)
    original=deepcopy(net)
    native=BMOPFTools.solve_opf(net;per_unit=true,s_base=1e4,solver_options=_GG_OPTIONS)
    replacement=Dict("g"=>(:generator,"native"))
    r=_gg_solve(net,d;replacements=replacement,per_unit=true,s_base=1e4)
    @test r.solve.publishable
    @test net==original
    @test r.objective ≈ native["objective"] atol=1e-6
    @test !haskey(r.network["generator"],"native")
    @test r.build_manifest.component_owners[(:generator,"native")]==:PowerOptLab
    @test r.network["custom_injection"]["p"] ≈ 700 atol=1e-5
    built=build_generator_model(net,[d];replacements=replacement,per_unit=true,s_base=1e4)
    for k in 1:3,component in (:real,:imag)
        v=BMOPFTools.opf_object(built.context,BMOPFTools.opf_generator_current_key("native",k;component))
        @test JuMP.is_fixed(v) && JuMP.fix_value(v)==0
    end
    # Source replacement must not inherit native source-bus bound suppression.
    n=inv_grid3_unbal(); n["bus"]["grid"]["perfectly_grounded_terminals"]=String[]
    n["bus"]["grid"]["v_min"]=[200.0,200.0,200.0]
    n["bus"]["grid"]["v_max"]=[260.0,260.0,260.0]
    src=SourceGenerator(id="s",bus="grid",voltage=GeneratorVoltageLaw(:fixed_phasor;
        phasor=[245,215,230].*cis.([0.05,-2.15,2.0])))
    repl=Dict("s"=>(:voltage_source,"vs"))
    rr=_gg_solve(n,src;replacements=repl,per_unit=true,s_base=1e4)
    @test rr.solve.publishable
    @test !haskey(rr.network["voltage_source"],"vs")
    bb=build_generator_model(n,[src];replacements=repl,per_unit=true,s_base=1e4)
    @test BMOPFTools.OpfModelKey(:constraint,:bus_voltage_upper,("grid",1)) in BMOPFTools.opf_object_keys(bb.context)
    # A real contradiction is caught, not silently dropped as a source bound.
    n["bus"]["grid"]["v_max"]=[240.0,260.0,260.0]
    @test !_gg_solve(n,src;replacements=repl,per_unit=true,s_base=1e4).solve.publishable
end

@testset "Validation rejects ambiguous or unsupported physics before stamping" begin
    net=inv_grid3_unbal()
    invalid=[
        GeneralizedGenerator(id="",bus="poc",connections=_GG_PORTS),
        GeneralizedGenerator(id="g",bus="missing",connections=_GG_PORTS),
        GeneralizedGenerator(id="g",bus="poc",connections=[("a","a")]),
        GeneralizedGenerator(id="g",bus="poc",connections=[("a","b"),("b","c"),("c","a")]),
        GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,impedance=-0.1),
        GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,impedance=ones(2,2)),
        GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,voltage=GeneratorVoltageLaw(:unknown)),
        GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=[0,1,2])),
        GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,voltage=GeneratorVoltageLaw(:phase_magnitudes)),
        GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,voltage=GeneratorVoltageLaw(:none;angles=[0,1,2])),
        GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,capability=GeneratorCapability(i_max=-1.0)),
        GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,capability=GeneratorCapability(earth_i_max=1.0)),
        GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,capability=GeneratorCapability(p_min=5.0,p_max=4.0)),
        GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,control=GeneratorControl(voltage_target=230.0)),
        GeneralizedGenerator(id="g",bus="poc",connections=[("a","n")],capability=GeneratorCapability(sequence=GeneratorSequenceLimits())),
        SourceGenerator(id="g",bus="poc",grounding=-1.0),
        SourceGenerator(id="g",bus="poc",grounding=0.0)]
    for d in invalid
        @test_throws ArgumentError validate_device(d,(net,))
    end
    @test_throws ArgumentError generator_sequence_impedance(-1,1,1)
    Z=generator_sequence_impedance(0.2+0.1im,0.1+0.4im,0.15+0.2im)
    a=cis(2pi/3); T=ComplexF64[1 1 1;1 a a^2;1 a^2 a]/3
    @test T*Z/T ≈ Diagonal(ComplexF64[0.2+0.1im,0.1+0.4im,0.15+0.2im]) atol=1e-12
    d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS)
    @test_throws ArgumentError build_generator_model(net,[d,d])
    @test_throws ArgumentError build_generator_model(net,[d];replacements=Dict("g"=>(:line,"l1")))
    b=build_generator_model(net,[d])
    count=JuMP.num_variables(BMOPFTools.opf_model(b.context))
    @test_throws ArgumentError stamp_device!(b.context,d)
    @test count==JuMP.num_variables(BMOPFTools.opf_model(b.context))
end
@testset "Generalized generator: voltage freedoms under unbalance" begin
    cases = [
        (GeneratorVoltageLaw(:none),GeneratorControl(p=[1200.0,700.0,300.0],q=[90.0,-50.0,130.0])),
        (GeneratorVoltageLaw(:fixed_phasor;phasor=_GG_E),GeneratorControl()),
        (GeneratorVoltageLaw(:fixed_magnitudes;phasor=[235,225,230].*cis.(_GG_PHI)),GeneratorControl(p=1500.0)),
        (GeneratorVoltageLaw(:common_magnitude;magnitude_min=190.0,magnitude_max=270.0),GeneratorControl(p=1500.0,q=300.0)),
        (GeneratorVoltageLaw(:phase_magnitudes;magnitude_min=190.0,magnitude_max=270.0),GeneratorControl(p=1500.0,q=[40.0,30.0,20.0]))]
    for (law,control) in cases
        @testset "$(law.mode)" begin
            Z=ComplexF64[0.4+0.8im 0.03+0.02im 0.01;0.03+0.02im 0.45+0.7im 0.02;0.01 0.02 0.35+0.9im]
            d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,impedance=Z,
                voltage=law,control=control)
            r=_gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4)
            @test r.solve.publishable
            g=r.devices["g"]
            @test abs(g.voltage_sequence[3])>1
            @test g.emf-g.port_voltage ≈ Z*g.port_current atol=1e-7
            @test abs(sum(g.terminal_current))<1e-7
            @test abs(g.power_balance_error)<1e-5
            @test g.series_loss ≈ real(dot(g.port_current,Z*g.port_current)) atol=1e-5
            @test complex(g.p,g.q) ≈ sum(g.port_voltage.*conj.(g.port_current)) atol=1e-5
            if control.p isa Real
                @test g.p ≈ control.p atol=1e-4
            elseif control.p isa Vector
                @test real.(g.port_voltage.*conj.(g.port_current)) ≈ control.p atol=1e-4
            end
            if law.mode==:fixed_phasor
                @test g.emf ≈ law.phasor atol=1e-6
            elseif law.mode==:fixed_magnitudes
                @test abs.(g.emf) ≈ abs.(law.phasor) atol=1e-5
                @test g.emf./law.phasor ≈ fill(g.emf[1]/law.phasor[1],3) atol=1e-6
            elseif law.mode in (:common_magnitude,:phase_magnitudes)
                rotated=g.emf.*cis.(-_GG_PHI)
                @test angle.(rotated./rotated[1]) ≈ zeros(3) atol=1e-6
                @test all(190-1e-5 .<= abs.(g.emf) .<= 270+1e-5)
                law.mode==:common_magnitude && @test maximum(abs.(g.emf))-minimum(abs.(g.emf))<1e-5
                law.mode==:phase_magnitudes && @test maximum(abs.(g.emf))-minimum(abs.(g.emf))>1
            end
        end
    end
end

@testset "Grounded source: independent return paths and physical oracle" begin
    currents=ComplexF64[10+2im,-3-5im,-2+4im]
    eph=[235,225,232].*cis.([0.02,-2.13,2.05])
    Z=Matrix(Diagonal(ComplexF64[0.1+0.2im,0.15+0.1im,0.12+0.3im,0.2]))
    # All phase/neutral currents, earth paths, and voltages are independently
    # computed from the circuit, not by calling the production stamping helper.
    for ground in (nothing,0.0,0.8+0.1im), pu in (false,true)
        @testset "ground=$ground pu=$pu" begin
            d=SourceGenerator(id="s",bus="poc",impedance=Z,grounding=ground,
                voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=eph))
            r=_gg_solve(_gg_fixture(),d;per_unit=pu,s_base=1e4,
                generator_hook! = _gg_imposed_current(currents))
            @test r.solve.publishable
            g=r.devices["s"]
            star=ground===nothing ? -0.2sum(currents) :
                iszero(ground) ? 0.0im : -sum(currents)/(1/0.2+1/ground)
            ineutral=star/0.2
            ig=-sum(currents)-ineutral
            @test g.port_current ≈ currents atol=1e-7
            @test g.star_voltage ≈ star atol=1e-7
            @test g.terminal_current[end] ≈ ineutral atol=1e-7
            @test g.earth_current ≈ ig atol=1e-7
            @test abs(sum(g.terminal_current)+g.earth_current)<1e-7
            @test g.terminal_voltage[1:3] ≈ eph.+star.-diag(Z)[1:3].*currents atol=1e-6
            @test abs(g.power_balance_error)<1e-6
            @test abs(3g.current_sequence[1]+g.terminal_current[end]+g.earth_current)<1e-7
            if ground===nothing
                @test abs(g.earth_current)<1e-7
            else
                @test abs(sum(g.terminal_current))>1
                @test g.ground_loss ≈ real(ground)*abs2(ig) atol=1e-6
            end
        end
    end
    # Nonzero PCC neutral potential exposes the correction missed by a
    # phase-to-neutral-only source power ledger.
    d=SourceGenerator(id="s",bus="poc",impedance=Z,grounding=0.8,
        voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=eph))
    hook! = (ctx,h) -> begin
        _gg_imposed_current(currents)(ctx,h)
        for (component,value) in ((:real,1.0),(:imag,0.3))
            v=BMOPFTools.opf_object(ctx,BMOPFTools.opf_bus_voltage_key("poc","n";component))
            JuMP.fix(v,value/h["s"].bases.voltage;force=true)
        end
        # External ideal neutral reference has its own return current.
        m=BMOPFTools.opf_model(ctx)
        nr=JuMP.@variable(m,start=0.0); ni=JuMP.@variable(m,start=0.0)
        BMOPFTools.add_terminal_injection!(ctx,"poc","n",nr,ni)
    end
    r=_gg_solve(_gg_fixture(grounded=false),d;per_unit=true,s_base=1e4,generator_hook! = hook!)
    @test r.solve.publishable
    g=r.devices["s"]
    @test abs(g.power_balance_error)<1e-5
    phasepower=sum(g.port_voltage.*conj.(g.port_current))
    @test complex(g.p,g.q) ≈ phasepower-g.terminal_voltage[end]*conj(g.earth_current) atol=1e-5
    @test abs(complex(g.p,g.q)-phasepower)>0.1
end

@testset "Ideal/singular series primitives and all source laws" begin
    for z in (0.0,Diagonal(ComplexF64[0,0.2+0.1im,0.3+0.2im]))
        d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,impedance=z,
            control=GeneratorControl(p=[500.0,800.0,200.0],q=[20.0,60.0,-10.0]))
        b=build_generator_model(inv_grid3_unbal(),[d];per_unit=true,s_base=1e4)
        @test b.handles["g"].internal_voltage_variables==0
        @test all(!(x isa JuMP.VariableRef) for x in b.handles["g"].er)
        @test all(JuMP.constraint_object(c).func isa Union{JuMP.AffExpr,JuMP.QuadExpr,JuMP.VariableRef}
            for c in values(b.handles["g"].constraints))
        r=_gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4)
        @test r.solve.publishable
        iszero(z) && @test r.devices["g"].emf ≈ r.devices["g"].port_voltage
    end
    # Exactly ideal source with a free PCC neutral: no ground-current shortcut.
    d=SourceGenerator(id="s",bus="poc",voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=_GG_E))
    r=_gg_solve(_gg_fixture(grounded=false),d;per_unit=false,
        generator_hook! = _gg_imposed_current(ComplexF64[10,-2-4im,-1+2im]))
    @test r.solve.publishable
    @test r.devices["s"].emf ≈ _GG_E atol=1e-6
    @test abs(r.devices["s"].star_voltage)<1e-8
    @test abs(r.devices["s"].power_balance_error)<1e-5
    # Additional source modes operate against a genuinely unbalanced native grid.
    for mode in (:none,:fixed_magnitudes,:common_magnitude,:phase_magnitudes)
        law=mode==:fixed_magnitudes ? GeneratorVoltageLaw(mode;phasor=_GG_E) :
            GeneratorVoltageLaw(mode;magnitude_min=190.0,magnitude_max=270.0)
        control=mode==:none ? GeneratorControl(p=[400.0,200.0,500.0],q=[30.0,-10.0,40.0]) :
            mode==:fixed_magnitudes ? GeneratorControl(p=1100.0) :
            mode==:common_magnitude ? GeneratorControl(p=1100.0,q=100.0) :
            GeneratorControl(p=1100.0,q=[30.0,-10.0,40.0])
        d=SourceGenerator(id="s",bus="poc",impedance=ComplexF64[0.4+0.8im,0.45+0.7im,0.35+0.9im,0.2],
            grounding=0.8,voltage=law,control=control)
        r=_gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4)
        @test r.solve.publishable
        @test abs(r.devices["s"].voltage_sequence[3])>1
        @test abs(r.devices["s"].power_balance_error)<1e-4
    end
end

@testset "Terminal PV differs from fixed internal excitation" begin
    z=0.5+1im
    net=inv_grid3_unbal()
    for k in 1:3
        net["linecode"]["lc"]["X_series_$(k)_$(k)"]=0.25
    end
    net["linecode"]["lc"]["X_series_4_4"]=0.0
    pv=GeneralizedGenerator(id="pv",bus="poc",connections=_GG_PORTS,impedance=z,
        voltage=GeneratorVoltageLaw(:common_magnitude;magnitude_min=180.0,magnitude_max=280.0),
        control=GeneratorControl(p=2400.0,voltage_target=232.0))
    r=_gg_solve(net,pv;per_unit=true,s_base=1e4)
    @test r.solve.publishable
    g=r.devices["pv"]
    @test g.p ≈ 2400 atol=1e-4
    @test abs(g.voltage_sequence[2]) ≈ 232 atol=1e-5
    @test abs(g.voltage_sequence[3])>1
    @test abs(abs(g.emf_sequence[2])-232)>1
    fixed=GeneralizedGenerator(id="fixed",bus="poc",connections=_GG_PORTS,impedance=z,
        voltage=GeneratorVoltageLaw(:fixed_magnitudes;phasor=232 .*cis.(_GG_PHI)),
        control=GeneratorControl(p=2400.0))
    f=_gg_solve(net,fixed;per_unit=true,s_base=1e4).devices["fixed"]
    @test abs(abs(f.voltage_sequence[2])-232)>0.1
end

@testset "Independent sequence capability and normalized limits" begin
    currents=ComplexF64[12,-6,-6] # I1=I2=6, so 8 A sequence caps permit 12 A in a.
    for phasecap in (13.0,10.0), base in (1e4,1e8)
        cap=GeneratorCapability(i_max=phasecap,sequence=GeneratorSequenceLimits(
            current_max=[0.1,8.0,8.0],voltage_max=[20.0,260.0,20.0]))
        d=SourceGenerator(id="s",bus="poc",impedance=[0.1,0.2,0.15,0.3],grounding=0.8,
            voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=_GG_E),capability=cap)
        r=_gg_solve(_gg_fixture(),d;per_unit=true,s_base=base,generator_hook! = _gg_imposed_current(currents))
        @test r.solve.publishable == (phasecap==13.0)
        if phasecap==13.0
            @test abs.(r.devices["s"].current_sequence) ≈ [0,6,6] atol=1e-5
        else
            @test isnan(r.devices["s"].p)
        end
    end
    # Inspect exact-zero sequence lowering: two affine rows, no squared-zero row.
    d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,impedance=0.3+0.6im,
        control=GeneratorControl(p=1800.0,q=300.0),
        capability=GeneratorCapability(sequence=GeneratorSequenceLimits(current_max=[0.0,20.0,10.0])))
    b=build_generator_model(inv_grid3_unbal(),[d];per_unit=true,s_base=1e4)
    rows=b.handles["g"].constraints
    @test haskey(rows,"sequence_i_1_r") && haskey(rows,"sequence_i_1_i")
    @test JuMP.constraint_object(rows["sequence_i_1_r"]).func isa JuMP.AffExpr
    @test !haskey(rows,"sequence_i_1_hi")
    r=_gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4,objective=:loss)
    @test r.solve.publishable
    @test abs(r.devices["g"].current_sequence[1])<1e-7
    @test abs(r.devices["g"].voltage_sequence[3])>1
end

@testset "Internal sequence zeros respect source-law degrees of freedom" begin
    for mode in (:none,:fixed_phasor,:fixed_magnitudes,:common_magnitude,:phase_magnitudes)
        law=mode in (:fixed_phasor,:fixed_magnitudes) ? GeneratorVoltageLaw(mode;phasor=_GG_E) :
            GeneratorVoltageLaw(mode;magnitude_min=190.0,magnitude_max=270.0)
        control=mode==:fixed_phasor ? GeneratorControl() :
            mode==:fixed_magnitudes ? GeneratorControl(p=1800.0) : GeneratorControl(p=1800.0,q=300.0)
        d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,impedance=0.3+0.6im,
            voltage=law,control=control,capability=GeneratorCapability(sequence=GeneratorSequenceLimits(
                location=:internal,voltage_min=[0.0,190.0,0.0],voltage_max=[0.0,270.0,0.0])))
        b=build_generator_model(inv_grid3_unbal(),[d];per_unit=true,s_base=1e4)
        @test any(startswith(k,"sequence_v") for k in keys(b.handles["g"].constraints)) == (mode==:none)
        r=_gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4)
        @test r.solve.publishable
        @test maximum(abs,r.devices["g"].emf_sequence[[1,3]])<1e-6
        @test abs(r.devices["g"].voltage_sequence[3])>1
        @test abs(r.devices["g"].power_balance_error)<1e-4
    end
    for mode in (:fixed_magnitudes,:common_magnitude,:phase_magnitudes)
        law=mode==:fixed_magnitudes ? GeneratorVoltageLaw(mode;phasor=_GG_E) :
            GeneratorVoltageLaw(mode;magnitude_min=190.0,magnitude_max=270.0)
        d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,voltage=law,
            capability=GeneratorCapability(sequence=GeneratorSequenceLimits(
                location=:internal,voltage_max=[0.0,180.0,0.0])))
        @test_throws ArgumentError build_generator_model(inv_grid3_unbal(),[d])
    end
end


@testset "Capability lowering keeps independent physical constraints" begin
    for location in (:poc,:internal), per_port in (false,true)
        p=per_port ? [500.0,200.0,800.0] : 1500.0
        q=per_port ? [40.0,10.0,50.0] : 100.0
        d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,impedance=0.3+0.6im,
            control=GeneratorControl(;p,q,power_location=location),
            capability=GeneratorCapability(power_location=location,p_min=p,p_max=p,q_min=q,q_max=q,s_max=hypot.(p,q)))
        b=build_generator_model(inv_grid3_unbal(),[d];per_unit=true,s_base=1e4)
        @test !any(startswith(k,"cap_") for k in keys(b.handles["g"].constraints))
        r=_gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4,objective=:loss)
        @test r.solve.publishable
        @test abs(r.devices["g"].voltage_sequence[3])>1
    end
    d=GeneralizedGenerator(id="one",bus="poc",connections=[("a","n")],
        control=GeneratorControl(p=700.0,q=50.0),
        capability=GeneratorCapability(i_max=20.0,terminal_i_max=[15.0,10.0]))
    b=build_generator_model(inv_grid3_unbal(),[d];per_unit=true,s_base=1e4)
    rows=b.handles["one"].constraints
    @test count(k -> occursin("_i_",k) && endswith(k,"_hi") && !occursin("_box_",k),keys(rows))==1
    @test JuMP.upper_bound(b.handles["one"].ir[1])*b.handles["one"].bases.current ≈ 10
    @test _gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4).solve.publishable
    # A positive S limit on free dispatch uses quadratic P/Q lifts, never quartic.
    free=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,
        capability=GeneratorCapability(s_max=2000.0))
    b=build_generator_model(inv_grid3_unbal(),[free])
    rows=b.handles["g"].constraints
    @test JuMP.constraint_object(rows["cap_S_p_1"]).func isa JuMP.QuadExpr
    @test JuMP.constraint_object(rows["cap_S_q_1"]).func isa JuMP.QuadExpr
    @test JuMP.constraint_object(rows["cap_S_circle_1"]).func isa JuMP.QuadExpr
    invalid=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,
        control=GeneratorControl(p=1000.0,q=0.0),capability=GeneratorCapability(s_max=900.0))
    @test_throws ArgumentError build_generator_model(inv_grid3_unbal(),[invalid])
end

@testset "Single-phase PV and grounded sources on mixed-width networks" begin
    net=inv_grid3_unbal()
    for k in 1:4
        net["linecode"]["lc"]["X_series_$(k)_$(k)"]=k==4 ? 0.0 : 0.25
    end
    d=GeneralizedGenerator(id="pv",bus="poc",connections=[("a","n")],impedance=0.3+0.5im,
        voltage=GeneratorVoltageLaw(:common_magnitude;magnitude_min=180.0,magnitude_max=300.0),
        control=GeneratorControl(p=700.0,voltage_target=244.0,voltage_metric=:phase))
    r=_gg_solve(net,d;per_unit=true,s_base=1e4)
    @test r.solve.publishable
    @test abs(r.devices["pv"].port_voltage[1]) ≈ 244 atol=1e-5
    @test r.devices["pv"].p ≈ 700 atol=1e-4
    @test r.devices["pv"].voltage_sequence===nothing
    for phases in (["a"],["a","b"])
        n=deepcopy(net); n["voltage_source"]["vs"]["v_angle"]=[0.03,pi+0.03,2.0]
        d=SourceGenerator(id="s",bus="poc",phase_terminals=phases,
            impedance=vcat(fill(0.4+0.7im,length(phases)),[0.2]),grounding=0.8,
            voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=ComplexF64[230,-230][1:length(phases)]))
        r=_gg_solve(n,d;per_unit=true,s_base=1e4)
        @test r.solve.publishable
        g=r.devices["s"]
        @test g.terminals==vcat(phases,["n"])
        @test g.voltage_sequence===nothing
        @test abs(sum(g.terminal_current)+g.earth_current)<1e-7
        @test abs(g.earth_current)>0.1
        @test abs(g.power_balance_error)<1e-4
    end
end

@testset "Exact redundant declarations do not duplicate equality rows" begin
    d=SourceGenerator(id="s",bus="poc",grounding=nothing,
        impedance=ComplexF64[0.4+0.7im,0.4+0.7im,0.4+0.7im,0.2],
        voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=_GG_E),
        capability=GeneratorCapability(earth_i_max=0.0))
    b=build_generator_model(inv_grid3_unbal(),[d];per_unit=true,s_base=1e4)
    @test !haskey(b.handles["s"].constraints,"earth_i_r")
    r=_gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4)
    @test r.solve.publishable
    @test abs(r.devices["s"].earth_current)<1e-7
    @test abs(r.devices["s"].voltage_sequence[3])>1
    d=GeneralizedGenerator(id="g",bus="poc",connections=_GG_PORTS,
        capability=GeneratorCapability(p_min=[500.0,200.0,800.0],p_max=[500.0,200.0,800.0],
            q_min=[40.0,10.0,50.0],q_max=[40.0,10.0,50.0]))
    b=build_generator_model(inv_grid3_unbal(),[d];per_unit=true,s_base=1e4)
    @test count(startswith("cap_"),keys(b.handles["g"].constraints))==6
    @test all(JuMP.constraint_object(row).set isa JuMP.MOI.EqualTo for row in values(b.handles["g"].constraints))
    r=_gg_solve(inv_grid3_unbal(),d;per_unit=true,s_base=1e4)
    @test r.solve.publishable
    @test real.(r.devices["g"].port_voltage.*conj.(r.devices["g"].port_current)) ≈ [500,200,800] atol=1e-5
end
