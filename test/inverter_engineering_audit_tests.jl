const _ENGINEERING_OPTIONS=("tol"=>1e-9,"max_iter"=>1000,"mu_strategy"=>"adaptive",
    "nlp_scaling_method"=>"none","bound_relax_factor"=>0.)
_engineering_plant(;kwargs...)=AdvancedInverter(;id="pv",bus="poc",topology=:THREE_LEG,
    phase_terminals=["a","b","c"],neutral="n",s_max=6000.,i_max=nothing,
    v_dc=750.,c_dc=300e-6,kwargs...)
_engineering_balanced(x)=x.*[1,cis(-2pi/3),cis(2pi/3)]

function _engineering_composed_solve(net,inv;per_unit,s_base,allocation=:balanced_converter)
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    device=ConstrainedInverter(inv,InverterOperatingConstraints(reactive=FixedReactivePower(-600.));
        current_allocation=allocation,dc_source=DCSourceBounds(p_min=0.))
    handles=Ref{Any}()
    hook! = ctx -> begin
        handles[]=stamp_device!(ctx,device;request)
        h=handles[].plant
        @objective(PowerOptLab._opf_model(ctx),Max,h.p_poc/(inv.s_max/h.sb))
    end
    ctx=PowerOptLab.build_opf_model(net;per_unit,s_base,add_objective=false,model_hook! = hook!,optimizer=Ipopt.Optimizer)
    PowerOptLab._set_inverter_solver_options!(PowerOptLab._opf_model(ctx),_ENGINEERING_OPTIONS)
    PowerOptLab._optimize_inverter_rails!(ctx,((inv,handles[].plant),))
    raw=solve_status(PowerOptLab._opf_model(ctx))
    extracted=extract_device(device,handles[],raw)
    (raw=raw,plant=extracted.plant,status=PowerOptLab._physical_solve_status(raw,extracted.plant.physical_feasibility))
end

@testset "Engineering audit: DC-link boundary and source line-frequency dissipation" begin
    voltage=_engineering_balanced(230.) .+ (2+im).*[1,cis(2pi/3),cis(-2pi/3)]
    for topology in (:THREE_LEG,:FOUR_LEG,:SPLIT_DC),source in
            ((;),(dc_source_r=0.,),(dc_source_r=.2,),(dc_source_r=.2,dc_source_l=.001),
             (dc_source_r=0.,dc_source_l=.001))
        banks=topology===:SPLIT_DC ? (c_dc_upper=600e-6,c_dc_lower=1200e-6) : (;)
        inv=_engineering_plant(;topology,banks...,source...)
        current=_engineering_balanced(4-im) .+ .5 .* [1,cis(2pi/3),cis(-2pi/3)]
        topology===:THREE_LEG || (current.+=.3+.1im)
        d=dc_source_line_frequency_diagnostics(inv;converter_voltage=voltage,converter_current=current)
        @test d.resolved && d.finite
        @test d.mean_power_boundary===:dc_link && !d.upstream_mean_path_loss_included
        @test d.source_line_frequency_loss_W>=0
        w=2pi*inv.f;excitation=sum(voltage.*current)/inv.v_dc
        if inv.dc_source_r===nothing
            @test d.source_fundamental_current_rms_A==d.source_second_harmonic_current_rms_A==0.
        elseif inv.dc_source_r==0. && inv.dc_source_l==0.
            @test d.source_second_harmonic_peak_current_A ≈ excitation
            if topology===:SPLIT_DC
                @test d.source_fundamental_current_A ≈ (1200e-6-600e-6)*sum(current)/(2*(600e-6+1200e-6))
            end
        else
            # Independently solve capacitor-node KCL, rather than using either
            # transfer function from the implementation.
            z=complex(inv.dc_source_r,2w*inv.dc_source_l);y=1/z
            if topology===:SPLIT_DC
                cu,cl=600e-6,1200e-6
                nodes=ComplexF64[im*2w*cu+y -y; -y im*2w*cl+y]\[-excitation,excitation]
                source2=-(nodes[1]-nodes[2])/z
                z1=complex(inv.dc_source_r,w*inv.dc_source_l);y1=1/z1
                nodes1=ComplexF64[im*w*cu+y1 -y1; -y1 im*w*cl+y1]\fill(-sum(current)/2,2)
                @test d.source_fundamental_current_A ≈ -(nodes1[1]-nodes1[2])/z1 rtol=1e-11
            else
                node=-excitation/(im*2w*inv.c_dc+y)
                source2=-node/z
                @test d.source_fundamental_current_A==0.
            end
            @test d.source_second_harmonic_peak_current_A ≈ source2 rtol=1e-11
        end
        # Direct period integration checks RMS/peak conventions and orthogonal
        # frequency loss addition. The mean source-current path is absent.
        samples=[sqrt(2)*real(d.source_fundamental_current_A*cis(theta))+
            real(d.source_second_harmonic_peak_current_A*cis(2theta)) for theta in range(0,2pi;length=4097)[1:end-1]]
        rms2=sum(abs2,samples)/length(samples)
        @test rms2 ≈ d.source_fundamental_current_rms_A^2+d.source_second_harmonic_current_rms_A^2 atol=1e-12
        @test d.source_line_frequency_loss_W ≈ something(inv.dc_source_r,0.)*rms2 atol=1e-12
    end
    inv=_engineering_plant(dc_source_r=0.,dc_source_l=1/((4pi*50)^2*300e-6))
    d=dc_source_line_frequency_diagnostics(inv;converter_voltage=voltage,converter_current=_engineering_balanced(4.))
    @test !d.resolved && d.network_margin<64eps(Float64)
    @test_throws DimensionMismatch dc_source_line_frequency_diagnostics(inv;converter_voltage=[1im],converter_current=[1im])
end

@testset "Engineering audit: physical current characteristic and absolute loss budget" begin
    inv=_engineering_plant()
    reference=nothing
    for (pu,base) in ((false,1e5),(true,1e4),(true,1e6),(true,1e8))
        ctx=PowerOptLab.build_opf_model(inv_grid3_unbal();per_unit=pu,s_base=base,add_objective=false)
        actual=PowerOptLab._inverter_current_characteristic_at_start(ctx,inv)
        @test 7. < actual < 10.
        reference===nothing ? (reference=actual) : (@test actual ≈ reference rtol=1e-12)
    end
    @test PowerOptLab._inverter_current_characteristic(_engineering_plant(i_max=10.),[0im,0im,0im])==10.
    @test PowerOptLab._inverter_current_characteristic(_engineering_plant(i_grid_max=7.),[0im,0im,0im])==7.
    @test PowerOptLab._inverter_current_characteristic(inv,[0im,0im,0im])>0
    for topology in (:THREE_LEG,:FOUR_LEG),rating in (nothing,10.),coefficient in (2.,100.,1e5)
        plant=_engineering_plant(;topology,i_max=rating,a_loss=coefficient)
        epsilon=PowerOptLab._loss_magnitude_epsilon(plant,rating)
        legs=topology===:FOUR_LEG ? 4 : 3
        @test 0<epsilon<=PowerOptLab._magnitude_epsilon(rating,1.)
        @test coefficient*legs*epsilon <= 1e-4*(1+eps(Float64))
        m=Model();@variable(m,x)
        for ib in (1.,4000.)
            expression=PowerOptLab._push_magnitude(m,nothing,x,0.,rating,ib;
                name="absolute_loss_budget",epsilon_SI=epsilon,scale_SI=10.)
            for radius in (0.,.1,10.)
                deficit=coefficient*legs*(radius-value(_->radius/ib,expression)*ib)
                @test -1e-8<=deficit<=1.001e-4
            end
        end
    end
    # The former unrated 1 mA width produced ~0.3 W bias here, exceeding
    # the unchanged independent 0.01 W acceptance floor.
    old=100*3*(5-(sqrt(5^2+1e-6)-.001))
    @test old>.29
    @test PowerOptLab._loss_magnitude_epsilon(_engineering_plant(i_max=10.,a_loss=2.),10.)==1e-8
end

@testset "Engineering audit: omitted current rating across SI and numerical bases" begin
    for extra in ((a_loss=100.,),(r_filter=.05,x_filter=.15,r_filter_grid=.03,x_filter_grid=.1,
            c_filter_mid=20e-6,r_filter_damping=.5,a_loss=2.))
        inv=_engineering_plant(;extra...)
        reference=nothing
        for (pu,base) in ((false,1e5),(true,1e4),(true,1e6),(true,1e8))
            r=_engineering_composed_solve(inv_grid3_unbal(),inv;per_unit=pu,s_base=base)
            @test r.status.publishable
            @test r.plant.physical_feasibility.passed
            @test r.plant.q_poc ≈ -600. atol=.01
            report=r.plant.physical_feasibility
            @test all(c.passed for c in report.checks if startswith(c.name,"operating_"))
            @test any(c.name=="operating_balanced_negative_sequence" for c in report.checks)
            pdc=only(c.value for c in report.checks if c.name=="dc_source_power_max")
            @test pdc<=5000.00011
            @test 0<=pdc-r.plant.p_dc<=.00011
            if reference!==nothing
                @test r.plant.p_poc ≈ reference.plant.p_poc atol=.01
                @test collect(r.plant.poc_phase_voltage) ≈ collect(reference.plant.poc_phase_voltage) atol=1e-5
            end
            reference=r
        end
    end
end

@testset "Engineering audit: independently reconstructed operating policies" begin
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    inv=_engineering_plant()
    function audit(policy,p,q;allocation=:balanced_converter,negative=0.)
        voltage=_engineering_balanced(230.)
        current=_engineering_balanced(complex(p,-q)/(3*230)) .+
            negative.*[1,cis(2pi/3),cis(-2pi/3)]
        device=ConstrainedInverter(inv,policy;current_allocation=allocation)
        PowerOptLab._constrained_operating_checks(device,request;
            grid_voltage=voltage,grid_current=current,converter_current=current,
            sensed_voltage=voltage,p_dc=p)
    end
    fixed=InverterOperatingConstraints(reactive=FixedReactivePower(-600.),p_max=3500.)
    @test all(c.passed for c in audit(fixed,3000.,-600.))
    @test "operating_reactive_relation" in [c.name for c in audit(fixed,3000.,-500.) if !c.passed]
    @test "operating_active_power" in [c.name for c in audit(fixed,4000.,-600.) if !c.passed]
    @test "operating_balanced_negative_sequence" in [c.name for c in audit(fixed,3000.,-600.;negative=1.) if !c.passed]
    @test all(c.passed for c in audit(fixed,3000.,-600.;allocation=:free,negative=1.))
    pf=InverterOperatingConstraints(reactive=FixedPowerFactor(.95,:absorbing))
    @test all(c.passed for c in audit(pf,3000.,-3000tan(acos(.95))))
    @test "operating_power_factor" in [c.name for c in audit(pf,3000.,0.) if !c.passed]
    for basis in (:rated,:available)
        policy=InverterOperatingConstraints(AverageVoltageVoltVarWatt(
            volt_watt=PiecewiseLinearLaw([200.,260.],[.5,.5];smoothing_epsilon=.01),volt_watt_basis=basis))
        @test all(c.passed for c in audit(policy,2000.,0.))
        @test "operating_volt_watt" in [c.name for c in audit(policy,4000.,0.) if !c.passed]
    end
    curve=VoltVarWattConstraints(AverageVoltageVoltVarWatt(
        volt_var=PiecewiseLinearLaw([200.,260.],[-.1,-.1];smoothing_epsilon=.01)))
    @test all(c.passed for c in audit(curve,3000.,-600.))
    @test "operating_reactive_relation" in [c.name for c in audit(curve,3000.,0.) if !c.passed]
    r=solve_constrained_inverter(inv_grid3_bal(),ConstrainedInverter(inv,fixed),request;
        objective=:feasibility,p_set=2500.,solver_options=_ENGINEERING_OPTIONS)
    @test r.solve.publishable
    @test only(c for c in physical_feasibility(r).checks if c.name=="operating_fixed_dispatch").passed
end

@testset "Engineering audit: raw solver success cannot bypass operating relationships" begin
    inv=_engineering_plant(i_max=10.)
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    coarse=("tol"=>1e3,"constr_viol_tol"=>1e3,"dual_inf_tol"=>1e3,
        "compl_inf_tol"=>1e3,"bound_relax_factor"=>0.)
    r=solve_constrained_inverter(inv_grid3_bal(),ConstrainedInverter(inv,
        InverterOperatingConstraints(reactive=FixedReactivePower(-600.))),request;
        objective=:feasibility,p_set=2500.,solver_options=coarse)
    @test r.solve.feasible && r.solve.optimal && r.solve.has_primal
    @test !r.solve.publishable
    failed=[c.name for c in physical_feasibility(r).checks if !c.passed]
    @test Set(failed)==Set(["operating_reactive_relation","operating_fixed_dispatch"])
    @test abs(r.plant.q_poc+600)>500
    @test abs(r.plant.p_poc-2500)>2000
    zero_watt=AverageVoltageVoltVarWatt(volt_watt=PiecewiseLinearLaw([200.,260.],[0.,0.];smoothing_epsilon=.01))
    for (policy,name) in ((InverterOperatingConstraints(reactive=FixedPowerFactor(.95,:absorbing)),"operating_power_factor"),
            (InverterOperatingConstraints(zero_watt),"operating_volt_watt"))
        r=solve_constrained_inverter(inv_grid3_bal(),ConstrainedInverter(inv,policy),request;solver_options=coarse)
        @test r.solve.feasible && r.solve.optimal
        @test !r.solve.publishable
        @test [c.name for c in physical_feasibility(r).checks if !c.passed]==[name]
    end
end
