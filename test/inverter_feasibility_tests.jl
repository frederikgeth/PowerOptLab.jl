@testset "Continuous rail enclosure" begin
    a=cis(2pi/3)
    u=230 .* [1,a^2,a]
    common=(id="audit",bus="poc",phase_terminals=["a","b","c"],neutral="n",
        s_max=20000.,c_dc=.001,v_dc=700.,m_max=.96)
    for topology in (:THREE_LEG,:FOUR_LEG,:SPLIT_DC)
        inv=AdvancedInverter(;common...,topology=topology)
        r=continuous_rail_margin(inv,u)
        exact=.96*700-(topology===:SPLIT_DC ? 2sqrt(2)*230 : sqrt(6)*230)
        topology===:SPLIT_DC && (exact/=2)
        @test r.converged
        @test r.lower_V <= exact <= r.upper_V
        @test r.upper_V-r.lower_V <= 1e-6
    end
    # Two harmonics with arbitrary phases: use the old, independent dense
    # waveform implementation only as an upper bound, not the new certificate.
    for topology in (:THREE_LEG,:FOUR_LEG,:SPLIT_DC), rotation in (0.,.137,.71)
        inv=AdvancedInverter(;common...,topology=topology,
            c_dc_upper=.001,c_dc_lower=.0017)
        v=(u .+ [8+3im,-2+5im,1-4im]) .* cis(rotation)
        d=(13+7im)*cis(2rotation)
        mid=(1.3-.7im)*cis(rotation)
        bus=(2-.4im)*cis(rotation)
        r=continuous_rail_margin(inv,v,d;midpoint=mid,midpoint_mean=2.,bus_fundamental=bus)
        sampled=PowerOptLab._switching_margin(inv,real.(v),imag.(v),real(d),imag(d),
            real(mid),imag(mid),2.;bre=real(bus),bim=imag(bus),n_eval=20000)
        @test r.converged
        @test r.lower_V <= sampled+1e-10
        @test sampled-r.upper_V < 1e-4
        @test r.upper_V-r.lower_V <= 1e-6
    end
    # Conservative enclosures do not invent feasibility at a tangency.
    c=PowerOptLab._trig_minimum(1.,-1.,0.,0.,0.)
    @test c.lower <= 0 <= c.upper
    @test c.upper-c.lower <= 1e-6
    flat=PowerOptLab._trig_minimum(3.,0.,0.,0.,0.)
    @test flat.lower <= 3 <= flat.upper
    @test flat.converged
    @test !PowerOptLab._trig_minimum(NaN,1.,0.,0.,0.).converged
    @test !PowerOptLab._trig_minimum(1.,1.,0.,0.,0.;atol=1e-30).converged
    @test_throws ArgumentError PowerOptLab._trig_minimum(1.,0.,0.,0.,0.;atol=0.)
end

@testset "Independent SI engineering ledger" begin
    a=cis(2pi/3); voltage=230 .* [1,a^2,a]; current=5 .* [1,a^2,a]
    common=(id="audit",bus="poc",phase_terminals=["a","b","c"],neutral="n",
        topology=:THREE_LEG,s_max=6000.,v_dc=750.,c_dc=.0003)
    inv=AdvancedInverter(;common...,i_max=10.,i_positive_max=6.,v_int_max=250.,
        v_int_min=210.,i_cap_max=5.,dv2_max=10.)
    report=inverter_physical_feasibility(inv;converter_voltage=voltage,
        converter_current=current,grid_voltage=voltage,dc_max=4000.)
    @test report.checked && report.passed
    @test all(c.passed for c in report.checks)
    failnames(r)=Set(c.name for c in r.checks if !c.passed)
    for (kw,name) in (((i_max=4.,),"converter_current_1"),
            ((i_grid_max=4.,),"grid_current_1"),
            ((i_positive_max=4.,),"positive_sequence_current"),
            ((s_max=3000.,),"converter_apparent_power"),
            ((v_int_max=220.,),"internal_voltage_max_1"),
            ((v_int_min=240.,),"internal_voltage_min_1"),
            ((v_dc=500.,),"continuous_rail_headroom"))
        plant=AdvancedInverter(;common...,kw...)
        r=inverter_physical_feasibility(plant;converter_voltage=voltage,
            converter_current=current,grid_voltage=voltage)
        @test !r.passed
        @test name in failnames(r)
    end
    # Independently derive ripple for a known negative-sequence current.
    i2=2 .* [1,a,a^2]; mixed=current+i2
    inv=AdvancedInverter(;common...,i_negative_max=.4,dv2_max=.1,i_cap_max=.1)
    r=inverter_physical_feasibility(inv;converter_voltage=voltage,
        converter_current=mixed,grid_voltage=voltage)
    @test Set(["negative_sequence_current","dc_voltage_ripple","capacitor_current"]) ⊆ failnames(r)
    ripplecheck=only(c for c in r.checks if c.name=="dc_voltage_ripple")
    @test ripplecheck.value ≈ 3*230*2/(4pi*50*.0003*750) rtol=1e-12
    # Neither successful solver metadata nor inconsistent lifted quantities
    # can hide a physical violation.
    success=PowerOptLab._result_solve_status("LOCALLY_SOLVED",true)
    rejected=PowerOptLab._physical_solve_status(success,r)
    @test rejected.termination_status==success.termination_status
    @test rejected.optimal && rejected.has_primal
    @test !rejected.publishable && rejected.feasible
    mismatch=inverter_physical_feasibility(inv;converter_voltage=voltage,
        converter_current=current,grid_voltage=voltage,reported=(p_conv=0.,),dc_max=3000.)
    @test "reported_p_conv_residual" in failnames(mismatch)
    @test "dc_source_power_max" in failnames(mismatch)
    broken=inverter_physical_feasibility(inv;converter_voltage=voltage.+1,
        converter_current=current,grid_voltage=voltage)
    @test "converter_filter_kvl_1" in failnames(broken)
    @test_throws ArgumentError inverter_physical_feasibility(inv;converter_voltage=voltage,
        converter_current=current,grid_voltage=voltage,
        tolerances=InverterFeasibilityTolerances(current_A=-1.))
    @test_throws DimensionMismatch inverter_physical_feasibility(inv;converter_voltage=voltage[1:2],
        converter_current=current,grid_voltage=voltage)
end

@testset "Physical acceptance rejects a missed inter-sample rail peak" begin
    inv=AdvancedInverter(id="audit",bus="poc",phase_terminals=["a","b","c"],neutral="n",
        topology=:THREE_LEG,s_max=20000.,v_dc=550.,c_dc=.05,m_max=.96,n_samples=6,c_loss=.01)
    result=solve_advanced_inverter(inv_grid3_bal(),inv;objective=:min_loss,p_set=0.,rail_refinement=false)
    @test result.termination_status in ("LOCALLY_SOLVED","OPTIMAL")
    @test !solve_status(result).publishable
    report=physical_feasibility(result)
    @test report.checked && !report.passed
    @test report.rail.upper_V < -1.
    @test "continuous_rail_headroom" in Set(c.name for c in report.checks if !c.passed)
    # No primal from the solver means no physical feasibility claim.
    failed=PowerOptLab._unchecked_inverter_feasibility()
    @test !failed.checked && !failed.passed
end

@testset "Adaptive rail cuts reject a false sampled optimum" begin
 inv=AdvancedInverter(id="cuts",bus="poc",phase_terminals=["a","b","c"],neutral="n",
  topology=:THREE_LEG,s_max=20000.,v_dc=550.,c_dc=.05,m_max=.96,n_samples=6,c_loss=.01)
 r=solve_advanced_inverter(inv_grid3_bal(),inv;objective=:min_loss,p_set=0.)
 @test !solve_status(r).publishable
 @test r.termination_status ∉ ("OPTIMAL","LOCALLY_SOLVED")
end

@testset "Four-wire, split-bank and single-phase certificates" begin
    a=cis(2pi/3); v=230 .* [1,a^2,a]
    current=5 .* [1,a^2,a] .+ 2 .* [1,a,a^2] .+ 1
    common=(id="fourwire",bus="poc",phase_terminals=["a","b","c"],neutral="n",s_max=6000.)
    failnames(r)=Set(c.name for c in r.checks if !c.passed)
    four=AdvancedInverter(;common...,topology=:FOUR_LEG,v_dc=900.,c_dc=.002,
        In_max=2.,i_zero_max=.4)
    report=inverter_physical_feasibility(four;converter_voltage=v,converter_current=current,grid_voltage=v)
    @test Set(["neutral_current","zero_sequence_current"]) ⊆ failnames(report)
    for source in ((;),(dc_source_r=0.,),(dc_source_r=.2,),(dc_source_r=.2,dc_source_l=.001))
        split=AdvancedInverter(;common...,topology=:SPLIT_DC,v_dc=900.,c_dc=.003,
            c_dc_upper=.002,c_dc_lower=.004,q_mid_balance_max=.5,v_mid_mean_max=1.,
            i_cap_upper_max=10.,i_cap_lower_max=10.,source...)
        # The exact charge required to cancel the natural mean displacement.
        r=inverter_physical_feasibility(split;converter_voltage=v,converter_current=current,
            grid_voltage=v,midpoint_charge=.45)
        @test r.passed
        @test r.rail.converged
    end
    split=AdvancedInverter(;common...,topology=:SPLIT_DC,v_dc=900.,c_dc=.003,
        c_dc_upper=.002,c_dc_lower=.004,q_mid_balance_max=.5,v_mid_mean_max=1.,
        dv_mid_max=.1,i_cap_upper_max=.1,i_cap_lower_max=.1)
    r=inverter_physical_feasibility(split;converter_voltage=v,converter_current=current,
        grid_voltage=v,midpoint_charge=.6)
    @test Set(["midpoint_fundamental_ripple","midpoint_mean","midpoint_balance_charge",
        "upper_capacitor_own_rating","lower_capacitor_own_rating"]) ⊆ failnames(r)
    single=AdvancedInverter(id="single",bus="poc",s_max=6000.,v_dc=600.,
        modulation_max=.5,p_ripple_max=500.)
    r=inverter_physical_feasibility(single;converter_voltage=[230.],converter_current=[5.],grid_voltage=[230.])
    @test Set(["single_phase_power_ripple","single_phase_modulation_1"]) ⊆ failnames(r)
end

@testset "Adaptive rail cuts retain feasible dispatch" begin
    for topology in (:THREE_LEG,:FOUR_LEG,:SPLIT_DC)
        inv=AdvancedInverter(id="refine",bus="poc",phase_terminals=["a","b","c"],neutral="n",
            topology=topology,s_max=20000.,i_max=40.,In_max=40.,
            v_dc=topology===:SPLIT_DC ? 700. : 590.,c_dc=.01,
            m_max=.96,n_samples=6,r_filter=.05,x_filter=.5,c_loss=.01)
        r=solve_advanced_inverter(inv_grid3_bal(),inv)
        @test solve_status(r).publishable
        @test physical_feasibility(r).passed
        @test physical_feasibility(r).rail.lower_V >= -1e-5
        @test r.p_poc > 10000.
    end
end

@testset "Engineering bounds reject fixed-current contradictions across bases" begin
    a=cis(2pi/3)
    for base in (1e4,1e6,1e8), bound in (:converter_current,:grid_current,:converter_va,:poc_va,:neutral)
        kw=bound===:converter_current ? (i_max=.4,) :
           bound===:grid_current ? (i_grid_max=.4,) :
           bound===:converter_va ? (s_max=1000.,) :
           bound===:poc_va ? (s_max=1000.,power_limit_location=:poc) :
           (topology=:FOUR_LEG,In_max=.4,)
        inv=AdvancedInverter(;id="bounds",bus="poc",phase_terminals=["a","b","c"],neutral="n",
            topology=:THREE_LEG,s_max=6000.,v_dc=1000.,c_dc=.001,kw...)
        current=5 .* [1,a^2,a] .+ (bound===:neutral ? 2. : 0.)
        hook! = ctx -> begin
            h=stamp_device!(ctx,inv;retain_lossless_magnitude_epigraph=false)
            m=PowerOptLab._opf_model(ctx)
            for k in 1:3
                @constraint(m,h.cri[k]*h.ib==real(current[k]))
                @constraint(m,h.cii[k]*h.ib==imag(current[k]))
            end
            @objective(m,Min,0.)
        end
        ctx=PowerOptLab.build_opf_model(inv_grid3_bal();per_unit=true,s_base=base,
            add_objective=false,model_hook! = hook!,optimizer=Ipopt.Optimizer)
        PowerOptLab._set_solver_options!(PowerOptLab._opf_model(ctx),
            ("tol"=>1e-8,"bound_relax_factor"=>0.,"max_iter"=>400))
        PowerOptLab.BMOPFTools.optimize_opf!(ctx)
        @test !solve_status(PowerOptLab._opf_model(ctx)).publishable
    end
    for base in (1e4,1e6,1e8)
        inv=AdvancedInverter(id="ripple",bus="poc",s_max=6000.,p_ripple_max=500.)
        r=solve_advanced_inverter(inv_grid(),inv;objective=:min_loss,p_set=3000.,q_set=0.,
            per_unit=true,s_base=base)
        @test !solve_status(r).publishable
    end
end

@testset "Filter initialization obeys KVL, KCL and three-leg topology" begin
    a=cis(2pi/3);v=230 .* [1,a^2,a] .+ 3 .* [1,a,a^2] .+ 2
    net=inv_grid3_src(mags=abs.(v),angs=angle.(v))
    for lcl in (false,true)
        inv=AdvancedInverter(id="start",bus="poc",phase_terminals=["a","b","c"],neutral="n",
            topology=:THREE_LEG,s_max=6000.,i_max=10.,v_dc=1000.,c_dc=.001,
            r_filter=.1,x_filter=.3,r_filter_neutral=.03,x_filter_neutral=.05,
            r_filter_grid=lcl ? .07 : 0.,x_filter_grid=lcl ? .2 : 0.,
            c_filter_mid=lcl ? 30e-6 : 0.)
        ref=Ref{Any}()
        PowerOptLab.build_opf_model(net;per_unit=true,s_base=1e6,
            add_objective=false,optimizer=Ipopt.Optimizer,
            model_hook! = ctx -> (ref[]=stamp_device!(ctx,inv;retain_lossless_magnitude_epigraph=false)))
        h=ref[]
        start(e)=JuMP.value(v->something(JuMP.start_value(v),0.),e)
        ph(re,im,scale)=complex.(start.(re),start.(im)).*scale
        r=inverter_physical_feasibility(inv;converter_voltage=ph(h.vrint,h.viint,h.vb),
            converter_current=ph(h.cri,h.cii,h.ib),grid_voltage=ph(h.vrpoc,h.vipoc,h.vb),
            grid_current=ph(h.gri,h.gii,h.ib),filter_voltage=ph(h.vrmid,h.vimid,h.vb))
        @test all(c.value<1e-10 for c in r.checks if occursin("kvl",c.name)||occursin("kcl",c.name)||c.name=="three_leg_zero_current")
        @test maximum(abs,imag.(ph(h.vrint,h.viint,h.vb)))>100.
    end
end

@testset "Strict LCL grid balance distinguishes source common mode" begin
    a=cis(2pi/3)
    balanced=230 .* [1,a^2,a]
    @test abs(sum(balanced))<1e-12
    request=InverterControlRequest(p_available=5000.,q_scale=0.)
    inv=AdvancedInverter(id="strict",bus="poc",phase_terminals=["a","b","c"],neutral="n",
        topology=:THREE_LEG,s_max=6000.,i_max=10.,v_dc=750.,c_dc=.001,
        r_filter=.1,x_filter=.3,r_filter_grid=.05,x_filter_grid=.15,c_filter_mid=30e-6)
    options=("tol"=>1e-9,"max_iter"=>500,"mu_strategy"=>"adaptive","nlp_scaling_method"=>"none")
    for offset in (0.,1.),target in (:converter,:grid)
        voltage=balanced .+ offset
        net=inv_grid3_src(mags=abs.(voltage),angs=angle.(voltage))
        ctrl=SequenceController(AverageVoltageVoltVarWatt();current_target=target)
        result=solve_controlled_inverter(net,ControlledDevice(inv,ctrl),request;solver_options=options)
        @test solve_status(result).publishable == (offset==0. || target===:converter)
        if offset==1. && target===:converter
            # A physically feasible three-leg plant carries the passive LCL
            # grid I0. Requiring grid I0=0 instead would force Ymid*Vmid,0=0,
            # incompatible with this source and a sequence-decoupled circuit.
            @test abs(sum(result.grid_phase_current)/3)>1e-3
            @test physical_feasibility(result).passed
        end
    end
end

@testset "Zero-PV PCC law has loss import without fictitious active service" begin
    a=cis(2pi/3);voltage=230 .* [1,a^2,a] .+ 2.3 .* [1,a,a^2]
    net=inv_grid3_src(mags=abs.(voltage),angs=angle.(voltage))
    options=("tol"=>1e-9,"max_iter"=>600,"mu_strategy"=>"adaptive","nlp_scaling_method"=>"none")
    for priority in (:common,:power,:balance),target in (:converter,:grid)
        inv=AdvancedInverter(id="night",bus="poc",phase_terminals=["a","b","c"],neutral="n",
            topology=:THREE_LEG,power_limit_location=:poc,s_max=6000.,i_max=.5,
            v_dc=750.,c_dc=.001,r_filter=.1,x_filter=.3,p_loss_fixed=20.,c_loss=.01)
        limiter=priority===:common ? CommonScaleLimiter(pq_priority=:var) :
            SequencePriorityLimiter(priority=priority,pq_priority=:var)
        controller=SequenceController(AverageVoltageVoltVarWatt(
            volt_var=PiecewiseLinearLaw([0.,400.],[.4,.4];smoothing_epsilon=.01)),
            negative_sequence_virtual_impedance(1.,0.),limiter;current_target=target)
        request=InverterControlRequest(p_available=0.,p_rated=6000.,q_scale=6000.)
        r=solve_controlled_inverter(net,ControlledDevice(inv,controller),request;solver_options=options)
        @test solve_status(r).publishable
        if solve_status(r).publishable
            @test abs(r.plant.p_dc)<.01
            @test r.plant.p_poc < -20.
            @test maximum(r.plant.i_mag) <= .50001
            @test r.exact_smooth_current_residual < .005
            @test all(c.passed for c in physical_feasibility(r).checks if startswith(c.name,"command_current"))
        end
    end
end

@testset "PCC allocator starts preserve a feasible priority anchor" begin
    m=Model()
    # Allocate P to a circle with a fixed Q offset, then move Q inward.
    # The second stage must begin inside the same circle, including the
    # saturated-anchor case that an arbitrary 0.9 scale start violates.
    tp=PowerOptLab._pcc_segment_smooth!(m,(0.,.6),(.9,0.),1.,1e-6;
        normalize=true,initialize=true)
    point=(.9tp,.6)
    tq=PowerOptLab._pcc_segment_smooth!(m,point,(0.,-.7),1.,1e-6;
        normalize=true,initialize=true)
    start(e)=JuMP.value(v->something(JuMP.start_value(v),0.),e)
    p=start(tp);q=start(tq)
    @test p ≈ 8/9 atol=1e-7
    @test q ≈ 1. atol=1e-7
    @test hypot(.9p,.6) <= 1.0 + 1e-10
    @test hypot(.9p,.6-.7q)<=1.
    # Initial values satisfy the defining root/output equalities, rather
    # than merely passing through an easy final dispatch solve.
    for (F,S) in JuMP.list_of_constraint_types(m)
        S===JuMP.MOI.EqualTo{Float64} || continue
        for cref in JuMP.all_constraints(m,F,S)
            row=JuMP.constraint_object(cref)
            @test abs(start(row.func)-row.set.value)<1e-9
        end
    end
end

@testset "Unsupported single-phase bounds cannot be silently ignored" begin
    for kw in ((dv2_max=1.,),(In_max=1.,),(dc_source_r=.1,),
            (dc_source_r=0.,dc_source_l=.001,))
        inv=AdvancedInverter(;id="unsupported",bus="poc",s_max=6000.,kw...)
        @test_throws ArgumentError validate_device(inv)
    end
end

@testset "Inverter numerical defaults respect optimizer-specific options" begin
    m=Model(Ipopt.Optimizer)
    PowerOptLab._set_inverter_solver_options!(m,())
    @test get_attribute(m,"bound_relax_factor")==0.
    PowerOptLab._set_inverter_solver_options!(m,("bound_relax_factor"=>1e-6,))
    @test get_attribute(m,"bound_relax_factor")==1e-6
    mock=Model(()->JuMP.MOI.Utilities.MockOptimizer(JuMP.MOI.Utilities.Model{Float64}()))
    PowerOptLab._set_inverter_solver_options!(mock,())
    @test JuMP.MOI.RawOptimizerAttribute("bound_relax_factor") ∉
        JuMP.MOI.get(backend(mock),JuMP.MOI.ListOfOptimizerAttributesSet())
end
