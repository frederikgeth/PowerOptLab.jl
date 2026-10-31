using LinearAlgebra

_projection_test_plant(;kwargs...)=AdvancedInverter(;
    id="pv",bus="poc",phase_terminals=["a","b","c"],neutral="n",
    topology=:THREE_LEG,s_max=6000.,i_max=10.,v_dc=900.,c_dc=300e-6,
    m_max=.96,kwargs...)

const _PROJECTION_A=cis(2pi/3)
_projection_balanced(x)=x.*[1,_PROJECTION_A^2,_PROJECTION_A]
const _PROJECTION_OPTIONS=("tol"=>1e-9,"max_iter"=>1000,
    "mu_strategy"=>"adaptive","nlp_scaling_method"=>"none")
const _PROJECTION_REQUEST=InverterControlRequest(p_available=10000.,p_rated=6000.,q_scale=6000.)

@testset "Offline current projection: independent capability oracle and bases" begin
    v=_projection_balanced(230.)
    inv=_projection_test_plant()
    for base in (1e4,1e5,1e6)
        r=project_inverter_current(inv,_PROJECTION_REQUEST;grid_voltage=v,
            requested_current=_projection_balanced(12.),s_base=base,solver_options=_PROJECTION_OPTIONS)
        @test solve_status(r).publishable
        @test physical_feasibility(r).passed
        @test r.request_feasibility===:violated
        @test r.plant_feasibility===:feasible
        # Lossless balanced case: the VA circle has the analytic radial projection.
        expected=6000/(3*230)
        @test collect(r.allocation.converter_terminal.phase_current) ≈ _projection_balanced(expected) atol=2e-5
        @test r.distance_A ≈ sqrt(3)*(12-expected) atol=2e-5
        @test r.allocation.plant.p_poc ≈ 6000. atol=.01
        @test solve_diagnostics(r).request_feasibility===:violated
    end
    desired=_projection_balanced(4+2im)
    r=project_inverter_current(inv,_PROJECTION_REQUEST;grid_voltage=v,requested_current=desired,
        solver_options=_PROJECTION_OPTIONS)
    @test r.request_feasibility===:feasible
    @test r.plant_feasibility===:feasible
    @test solve_status(r).publishable
    @test r.distance_A < 1e-4
    @test collect(r.requested_current) == desired
    # The source voltage is a fixed boundary, not a network variable to exploit.
    @test collect(r.allocation.plant.poc_phase_voltage) ≈ v atol=1e-8
    optional=project_inverter_current(_projection_test_plant(i_max=nothing),_PROJECTION_REQUEST;
        grid_voltage=v,requested_current=_projection_balanced(12.),s_base=1e8,
        solver_options=_PROJECTION_OPTIONS)
    @test solve_status(optional).publishable
    @test optional.allocation.plant.p_poc ≈ 6000. atol=.01
end

@testset "Offline current projection: rail-infeasible request, feasible plant witness" begin
    inv=_projection_test_plant(v_dc=620.,r_filter=.13225,x_filter=1.3225,
        p_loss_fixed=20.,c_loss=.01,dv2_max=7.5,i_cap_max=1.)
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    v=_projection_balanced(253.)
    desired=_projection_balanced(4800/(3*253))
    r=project_inverter_current(inv,request;grid_voltage=v,requested_current=desired,
        target=:grid,solver_options=_PROJECTION_OPTIONS)
    @test r.request_feasibility===:violated
    @test r.requested_physical.rail.upper_V < -26.
    @test [c.name for c in r.requested_physical.checks if !c.passed]==["continuous_rail_headroom"]
    @test r.plant_feasibility===:feasible
    @test solve_status(r).publishable
    @test physical_feasibility(r).passed
    @test physical_feasibility(r).rail.lower_V >= -1e-5
    @test r.allocation.plant.q_poc < -5000.
    @test r.allocation.plant.p_dc <= request.p_available + .01
    @test minimum(r.allocation.operation.converter_current_margin_A) >= -1e-5
    # Same capability with a rotated electrical reference; no angle-specific safety.
    rotation=cis(.37)
    rotated=project_inverter_current(inv,request;grid_voltage=v.*rotation,
        requested_current=desired.*rotation,target=:grid,solver_options=_PROJECTION_OPTIONS)
    @test solve_status(rotated).publishable
    @test physical_feasibility(rotated).passed
    @test rotated.distance_A ≈ r.distance_A atol=2e-4
end

@testset "Offline current projection: LCL arm distinction and neutral primitive" begin
    Rc=Matrix(Diagonal([.08,.08,.08,.1]))+fill(.005,4,4)
    Xc=Matrix(Diagonal([.30,.30,.30,.2]))+fill(.01,4,4)
    Rg=Matrix(Diagonal([.04,.04,.04,.05]))+fill(.003,4,4)
    Xg=Matrix(Diagonal([.15,.15,.15,.12]))+fill(.007,4,4)
    # Unequal mutual terms survive the neutral reduction.
    Rc[1,2]=Rc[2,1]=.007; Xg[2,4]=Xg[4,2]=.009
    inv=_projection_test_plant(topology=:FOUR_LEG,
        r_filter_matrix=Rc,x_filter_matrix=Xc,r_filter_grid_matrix=Rg,x_filter_grid_matrix=Xg,
        c_filter_mid=20e-6,r_filter_damping=1.,b_filter_shunt=1e-4,In_max=8.)
    v=ComplexF64[235cis(.03),220cis(-2pi/3-.02),228cis(2pi/3)]
    desired=_projection_balanced(2.) .+ [.1+.2im,-.05+.1im,.03-.05im]
    r=project_inverter_current(inv,_PROJECTION_REQUEST;grid_voltage=v,
        requested_current=desired,target=:grid,solver_options=_PROJECTION_OPTIONS)
    @test solve_status(r).publishable
    @test r.request_feasibility===:feasible
    @test r.distance_A < 2e-4
    # Independent conductor-domain construction, including the neutral return.
    zc=complex.(Rc,Xc)
    zg=complex.(Rg,Xg)
    drop(Z,i)=let d=Z*[i;-sum(i)]; d[1:3].-d[4] end
    vm=v+drop(zg,desired)
    ic=desired+vm/(1-im/(2pi*50*20e-6))
    u=vm+drop(zc,ic)
    @test r.requested_state.filter_voltage ≈ vm atol=1e-10
    @test r.requested_state.converter_current ≈ ic atol=1e-10
    @test r.requested_state.converter_voltage ≈ u atol=1e-10
    @test collect(r.allocation.converter_terminal.phase_current) ≈ ic atol=2e-4
    @test abs(r.allocation.plant.q_poc-imag(sum(v.*conj.(desired)))) > 10.
    # Converter-arm target must invert the transfer rather than equating arms.
    converter=project_inverter_current(inv,_PROJECTION_REQUEST;grid_voltage=v,
        requested_current=ic,target=:converter,solver_options=_PROJECTION_OPTIONS)
    @test solve_status(converter).publishable
    @test converter.request_feasibility===:feasible
    @test converter.requested_state.grid_current ≈ desired atol=1e-10
    @test collect(converter.allocation.grid_phase_current) ≈ desired atol=2e-4
end

@testset "Offline current projection: midpoint state, source bounds and inconclusive failures" begin
    inv=_projection_test_plant(topology=:SPLIT_DC,v_dc=1500.,q_mid_balance_max=.02,
        v_mid_mean_max=30.,p_loss_fixed=20.)
    v=_projection_balanced(230.)
    r=project_inverter_current(inv,_PROJECTION_REQUEST;grid_voltage=v,
        requested_current=_projection_balanced(2.),midpoint_charge=.001,
        solver_options=_PROJECTION_OPTIONS)
    @test solve_status(r).publishable
    @test r.request_feasibility===:feasible
    @test r.allocation.plant.q_mid_balance ≈ .001 atol=1e-12
    @test r.requested_state.midpoint_charge == .001
    night=project_inverter_current(_projection_test_plant(p_loss_fixed=20.),
        InverterControlRequest(p_available=0.,p_rated=6000.,q_scale=6000.);
        grid_voltage=v,requested_current=zeros(3),solver_options=_PROJECTION_OPTIONS)
    @test solve_status(night).publishable
    @test night.request_feasibility===:violated
    @test night.plant_feasibility===:feasible
    @test night.allocation.plant.p_dc ≈ 0. atol=.01
    @test night.allocation.plant.p_poc ≈ -20. atol=.01

    failed=project_inverter_current(_projection_test_plant(),_PROJECTION_REQUEST;
        grid_voltage=v,requested_current=_projection_balanced(12.),solver_options=("max_iter"=>0,))
    @test !solve_status(failed).publishable
    @test failed.request_feasibility===:violated
    @test failed.plant_feasibility===:unresolved
    @test isnan(failed.distance_A)
    witnessed=project_inverter_current(_projection_test_plant(),_PROJECTION_REQUEST;
        grid_voltage=v,requested_current=_projection_balanced(2.),solver_options=("max_iter"=>0,))
    @test !solve_status(witnessed).publishable
    @test witnessed.request_feasibility===:feasible
    @test witnessed.plant_feasibility===:feasible

    resonant=_projection_test_plant(x_filter_grid=1.,c_filter_mid=1/(2pi*50))
    unknown=project_inverter_current(resonant,_PROJECTION_REQUEST;
        grid_voltage=v,requested_current=zeros(3),solver_options=("max_iter"=>0,))
    @test unknown.requested_state===nothing
    @test unknown.requested_physical===nothing
    @test unknown.request_feasibility===:unresolved
    @test unknown.plant_feasibility===:unresolved
end

@testset "Offline current projection: input contracts and enclosure uncertainty" begin
    inv=_projection_test_plant(); v=_projection_balanced(230.); i=_projection_balanced(2.)
    for bad in (:poc,:neutral)
        @test_throws ArgumentError project_inverter_current(inv,_PROJECTION_REQUEST;
            grid_voltage=v,requested_current=i,target=bad)
    end
    @test_throws DimensionMismatch project_inverter_current(inv,_PROJECTION_REQUEST;
        grid_voltage=[230.],requested_current=i)
    @test_throws ArgumentError project_inverter_current(inv,_PROJECTION_REQUEST;
        grid_voltage=v,requested_current=[NaN,0.,0.])
    @test_throws ArgumentError project_inverter_current(inv,_PROJECTION_REQUEST;
        grid_voltage=v,requested_current=i,midpoint_charge=.01)
    @test_throws ArgumentError project_inverter_current(inv,_PROJECTION_REQUEST;
        grid_voltage=v,requested_current=i,s_base=0.)
    for charge in (.01,Inf)
        @test_throws ArgumentError project_inverter_current(_projection_test_plant(topology=:SPLIT_DC),
            _PROJECTION_REQUEST;grid_voltage=v,requested_current=i,midpoint_charge=charge)
    end
    @test_throws ArgumentError project_inverter_current(_projection_test_plant(pwm_strategy=:SPWM,f_sw=1e4),
        _PROJECTION_REQUEST;grid_voltage=v,requested_current=i)
    # A failed enclosure alone must not become a proof of request infeasibility.
    check=InverterLimitCheck("continuous_rail_headroom",:V,-1.,0.,Inf,1e-5,-1.,false)
    report=InverterFeasibilityReport(true,false,[check],
        (lower_V=-1.,upper_V=1.,angle_rad=0.,branch="test",converged=false,evaluations=1))
    @test PowerOptLab._projection_request_feasibility(report)===:unresolved
    # A physically feasible inverter at a different boundary is not an accepted
    # projection at the requested voltage/actuator state, even after solver success.
    good=inverter_physical_feasibility(inv;converter_voltage=v,converter_current=i,
        grid_voltage=v,grid_current=i,filter_voltage=v)
    @test good.passed
    wrong_voltage=PowerOptLab._projection_boundary_report(good,v.+[.001,0.,0.],v,0.,0.)
    @test !wrong_voltage.passed
    @test !PowerOptLab._physical_solve_status(
        PowerOptLab._result_solve_status("LOCALLY_SOLVED",true),wrong_voltage).publishable
    wrong_charge=PowerOptLab._projection_boundary_report(good,v,v,1e-6,0.)
    @test !wrong_charge.passed
end
