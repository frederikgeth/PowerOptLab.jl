_energy_bal(x)=x.*ComplexF64[1,cis(-2pi/3),cis(2pi/3)]
_energy_inv(;kwargs...)=AdvancedInverter(;id="pv",bus="poc",topology=:THREE_LEG,
    phase_terminals=["a","b","c"],neutral=nothing,s_max=6000.,i_max=10.,
    v_dc=750.,c_dc=300e-6,i_cap_max=1.,dv2_max=7.5,power_limit_location=:poc,
    dc_energy_model=:nonlinear,kwargs...)

@testset "DC energy: mean voltage, capacitor RMS and nonlinear Fourier coefficient" begin
    inv=_energy_inv();u=_energy_bal(230.)
    for ripple in (0.,1.,7.5,30.),angle in (0.,.4,1.7)
        s2=4pi*50*300e-6*750*ripple*cis(angle)
        i=_energy_bal(3.).+s2/690 .* ComplexF64[1,cis(2pi/3),cis(-2pi/3)]
        d=dc_link_energy_diagnostics(inv;converter_voltage=u,converter_current=i)
        @test d.resolved
        theta=2pi.*(0:8191)./8192
        voltage=sqrt.(d.center_V2.+2inv.v_dc.*real.(d.energy_coordinate.*cis.(2theta)))
        cap=-real.(sum(u.*i).*cis.(2theta))./voltage
        @test sum(voltage)/length(voltage) ≈ 750. atol=1e-9
        @test sum(abs2,cap)/length(cap) ≈ d.capacitor_lf_sq atol=1e-10
        @test 2sum(voltage.*cis.(-2theta))/length(voltage) ≈ d.voltage_2omega atol=1e-9
        @test minimum(voltage)>=d.minimum_V-1e-9
        @test maximum(voltage)<=d.maximum_V+1e-9
        ripple==7.5 && @test abs(d.voltage_2omega)>7.50005
    end
end

@testset "DC energy: full-cycle nonlinear rail enclosure" begin
    for four in (false,true),factor in (1.,1.5),angle in (0.,.6)
        inv=_energy_inv(topology=four ? :FOUR_LEG : :THREE_LEG,neutral=four ? "n" : nothing)
        u=_energy_bal(230factor*cis(angle));i=_energy_bal(3.) .+ 0.7 .* [1,cis(2pi/3),cis(-2pi/3)]
        energy=PowerOptLab._dc_energy_state(inv,sum(u.*i))
        r=continuous_rail_margin(inv,u,energy.voltage_2omega;energy)
        gaps=Float64[]
        for theta in 2pi.*(0:8191)./8192
            bus=sqrt(energy.center_V2+2inv.v_dc*real(energy.energy_coordinate*cis(2theta)))
            volt=sqrt(2).*real.(u.*cis(theta))
            needed=maximum(abs(volt[a]-volt[b]) for (a,b) in ((1,2),(2,3),(3,1)))
            four && (needed=max(needed,maximum(abs,volt)))
            push!(gaps,inv.m_max*bus-needed)
        end
        @test r.converged
        @test r.lower_V<=minimum(gaps)+1e-7
        @test r.upper_V>=minimum(gaps)-1e-4
        @test (r.lower_V>=0)==(factor==1.)
    end
end

@testset "DC energy: bounded allocation, extrema and coordinate invariance" begin
    inv=_energy_inv(r_filter=.13225,x_filter=1.3225,p_loss_fixed=15.,a_loss=2.,c_loss=.3)
    req=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    v=_energy_bal(230.);ref=nothing
    for base in (1e4,1e5,1e6)
        r=allocate_inverter_services(inv,req;grid_voltage=v,reactive_request=-1800.,
            active_fraction=.7,negative_sequence_request=2.0+0im,s_base=base)
        @test solve_status(r).publishable
        @test r.allocation.plant.dv2<=7.50001
        @test r.allocation.plant.i_cap_thermal<=1.00001
        @test r.allocation.plant.p_poc ≈ 4200. atol=.001
        @test all(solve_status(s).publishable for (_,s) in r.stages)
        ref===nothing || @test r.allocation.plant.dv2 ≈ ref atol=1e-5
        ref=r.allocation.plant.dv2
    end
    bounded=_energy_inv(v_dc_min=745.,v_dc_max=755.)
    r=allocate_inverter_services(bounded,req;grid_voltage=v,reactive_request=-1800.,
        active_fraction=.7,negative_sequence_request=2.0+0im)
    @test solve_status(r).publishable
    diagnostic=dc_link_energy_diagnostics(bounded;converter_voltage=collect(r.allocation.converter_terminal.phase_voltage),
        converter_current=collect(r.allocation.converter_terminal.phase_current))
    @test diagnostic.minimum_V ≈ 745. atol=.001
    @test diagnostic.maximum_V<=755.00001
    for four in (false,true)
        plant=_energy_inv(topology=four ? :FOUR_LEG : :THREE_LEG,neutral=four ? "n" : nothing,
            pwm_strategy=:CENTERED,f_sw=10000.,r_filter=.13225,x_filter=1.3225,
            p_loss_fixed=15.,a_loss=2.,c_loss=.3,In_max=four ? 10. : nothing)
        result=four ? allocate_phase_inverter_services(plant,req;grid_voltage=v,reactive_request=fill(-600.,3)) :
            allocate_inverter_services(plant,req;grid_voltage=v,reactive_request=-1800.)
        if four
            # Preserve the original zero-neutral-inductance hardware. Its
            # all-duty coupled thermal bound exhausts the 1 A budget.
            @test !solve_status(result).publishable
            @test result.current_command===nothing
            continue
        end
        @test solve_status(result).publishable
        @test result.allocation.plant.i_cap_thermal<=1.00001
        @test all(c.passed for c in physical_feasibility(result).checks)
        night=InverterControlRequest(p_available=0.,p_rated=6000.,q_scale=6000.)
        offsun=four ? allocate_phase_inverter_services(plant,night;grid_voltage=v,reactive_request=fill(-600.,3)) :
            allocate_inverter_services(plant,night;grid_voltage=v,reactive_request=-1800.)
        @test solve_status(offsun).publishable
        @test abs(offsun.allocation.plant.p_dc)<.001
        @test offsun.allocation.plant.p_poc < -15.
        @test all(c.passed for c in physical_feasibility(offsun).checks)
    end
    @test_throws ArgumentError validate_device(_energy_inv(topology=:SPLIT_DC,neutral="n"))
    @test_throws ArgumentError validate_device(_energy_inv(dc_source_r=.1))
    @test_throws ArgumentError validate_device(_energy_inv(v_dc_min=751.))
    @test_throws ArgumentError validate_device(_energy_inv(dc_energy_model=:linearized,v_dc_max=760.))
    @test_throws ArgumentError continuous_rail_margin(inv,v)
end

@testset "DC energy: SI and per-unit engineering voltage bounds" begin
    inv=_energy_inv(v_dc_min=749.4,v_dc_max=751.,r_filter=.13225,x_filter=1.3225)
    v=_energy_bal(230.).+4.6.*[1,cis(2pi/3),cis(-2pi/3)]
    req=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    net=PowerOptLab._projection_network(inv,v)
    device=ConstrainedInverter(inv,InverterOperatingConstraints(reactive=FixedReactivePower(-1800.));
        dc_source=DCSourceBounds(p_min=0.))
    ref=nothing
    for pu in (false,true)
        handles=Ref{Any}()
        hook! = ctx -> begin
            handles[]=stamp_device!(ctx,device;request=req)
            h=handles[].plant
            @objective(PowerOptLab._opf_model(ctx),Max,h.p_poc/(inv.s_max/h.sb))
        end
        ctx=PowerOptLab.build_opf_model(net;per_unit=pu,s_base=1e4,add_objective=false,model_hook! = hook!,optimizer=Ipopt.Optimizer)
        PowerOptLab._set_inverter_solver_options!(PowerOptLab._opf_model(ctx),
            ("tol"=>1e-9,"max_iter"=>1000,"nlp_scaling_method"=>"none","mu_strategy"=>"adaptive"))
        PowerOptLab._optimize_inverter_rails!(ctx,((inv,handles[].plant),))
        status=solve_status(PowerOptLab._opf_model(ctx));result=extract_device(device,handles[],status)
        @test status.publishable && result.plant.physical_feasibility.passed
        @test 3000<result.plant.p_poc<4500
        @test all(c.passed for c in result.plant.physical_feasibility.checks)
        minimum=only(c.value for c in result.plant.physical_feasibility.checks if c.name=="dc_energy_positive_voltage")
        @test minimum ≈ 749.4 atol=.001
        ref===nothing || @test result.plant.p_poc ≈ ref atol=.01
        ref=result.plant.p_poc
    end
end
