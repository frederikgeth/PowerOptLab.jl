using BMOPFTools: parse_bmopf
import BMOPFTools

const _SEQ_LIMIT_OPTIONS=("tol"=>1e-8,"max_iter"=>600,"bound_relax_factor"=>0.,
    "mu_strategy"=>"adaptive","nlp_scaling_method"=>"none")

function _sequence_limit_plant(;kwargs...)
    AdvancedInverter(;id="pv",bus="poc",phase_terminals=["a","b","c"],
        neutral=nothing,topology=:THREE_LEG,s_max=6000.,i_max=10.,
        v_dc=750.,c_dc=300e-6,r_filter=.13225,x_filter=1.3225,kwargs...)
end

@testset "Sequence ratings: equal ratings and fleet margin reporting" begin
    a=cis(2pi/3)
    v=230 .* [1,a^2,a] .+ 2.3 .* [1,a,a^2]
    net=inv_grid3_src(mags=abs.(v),angs=angle.(v))
    net["ibr"]=Dict("pv"=>Dict("bus"=>"poc","terminal_map"=>["a","b","c"],
        "topology"=>"THREE_LEG","prime_mover"=>"PV","s_max"=>fill(2000.,3),
        "p_min"=>zeros(3),"p_max"=>fill(2000.,3),"q_min"=>fill(-2000.,3),"q_max"=>fill(2000.,3)))
    inv=_sequence_limit_plant(power_limit_location=:poc,i_positive_max=4.,i_negative_max=4.)
    control=SequenceController(AverageVoltageVoltVarWatt(),negative_sequence_virtual_impedance(1.,0.),
        SequencePriorityLimiter(priority=:balance,pq_priority=:var);sensing_reference="n")
    req=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    spec=ControlledInverterFleetSpec(Dict("pv"=>ControlledDevice(inv,control)),Dict("pv"=>req))
    result=solve_controlled_inverter_fleet(net,spec;s_base=1e5,solver_options=_SEQ_LIMIT_OPTIONS)
    @test solve_status(result).publishable
    if solve_status(result).publishable
        row=only(controlled_inverter_rows(result))
        @test row.positive_sequence_current_rating_A==4.
        @test row.negative_sequence_current_rating_A==4.
        @test -.00001 <= row.positive_sequence_current_margin_A < .01
        @test row.negative_sequence_current_margin_A > 1.5
        @test row.p_poc_W < 3000.
    end
end

@testset "Sequence ratings: hard plant bounds survive unit scaling" begin
    for s_base in (1e4,1e6,1e8), component in (:zero,:positive,:negative), limit in (.4,6.)
        limits=component===:positive ? (i_positive_max=limit,) :
            component===:negative ? (i_negative_max=limit,) : (i_zero_max=limit,)
        inv=_sequence_limit_plant(;topology=component===:zero ? :FOUR_LEG : :THREE_LEG,
            neutral=component===:zero ? "n" : nothing,limits...)
        handles=Ref{Any}()
        hook! = ctx -> begin
            h=stamp_device!(ctx,inv;retain_lossless_magnitude_epigraph=false)
            handles[]=h
            m=PowerOptLab._opf_model(ctx)
            z,p,n=PowerOptLab._sequence_pair(collect(h.cri),collect(h.cii))
            # Independently force 5 A positive and 2 A negative sequence;
            # include 2 A zero sequence for the four-leg fixture. A 0.4 A
            # rating is contradictory in each channel; 6 A is feasible.
            @constraint(m,p[1]*h.ib==5.)
            @constraint(m,p[2]*h.ib==0.)
            @constraint(m,n[1]*h.ib==2.)
            @constraint(m,n[2]*h.ib==0.)
            if component===:zero
                @constraint(m,z[1]*h.ib==2.)
                @constraint(m,z[2]*h.ib==0.)
            end
            @objective(m,Min,0.)
        end
        ctx=BMOPFTools.build_opf_model(inv_grid3_bal();per_unit=true,s_base=s_base,
            add_objective=false,model_hook! = hook!,optimizer=Ipopt.Optimizer)
        m=PowerOptLab._opf_model(ctx)
        PowerOptLab._set_solver_options!(m,_SEQ_LIMIT_OPTIONS)
        BMOPFTools.optimize_opf!(ctx)
        outcome=PowerOptLab._solve_outcome(m)
        @test outcome.optimal == (limit==6.)
        if outcome.optimal
            r=extract_device(inv,handles[],SolveStatus(outcome))
            @test r.i_negative ≈ 2. atol=1e-5
        end
    end
end

@testset "Sequence ratings: controller saturation and service priority" begin
    a=cis(2pi/3)
    v=230 .* [1,a^2,a] .+ 2.3 .* [1,a,a^2]
    net=inv_grid3_src(mags=abs.(v),angs=angle.(v))
    request=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
    for location in (:converter,:poc), target in (:converter,:grid),
        priority in (:common,:power,:balance), bound in (:positive,:negative)
        @testset "$location $target $priority $bound" begin
            limiter=priority===:common ? CommonScaleLimiter() :
                SequencePriorityLimiter(priority=priority)
            control=SequenceController(AverageVoltageVoltVarWatt(),
                negative_sequence_virtual_impedance(1.,0.),limiter;
                current_target=target,sensing_reference="n")
            limits=bound===:positive ? (i_positive_max=4.,) : (i_negative_max=.4,)
            inv=_sequence_limit_plant(;power_limit_location=location,limits...)
            r=solve_controlled_inverter(net,ControlledDevice(inv,control),request;
                s_base=1e5,solver_options=_SEQ_LIMIT_OPTIONS)
            @test r.solve.publishable
            if r.solve.publishable
                # Independent Fortescue sums from physical converter phasors.
                phase=r.converter_terminal.phase_current
                i1=abs((phase[1]+a*phase[2]+a^2*phase[3])/3)
                i2=abs((phase[1]+a^2*phase[2]+a*phase[3])/3)
                @test i1 ≈ r.plant.i_positive atol=1e-8
                @test i2 ≈ r.plant.i_negative atol=1e-8
                @test r.exact_smooth_current_residual < .005
                @test r.plant.p_dc <= request.p_available+.01
                if bound===:negative
                    @test .395 < i2 <= .40001
                    priority===:common ? (@test r.plant.p_poc<1500.) :
                        (@test r.plant.p_poc>4900.)
                else
                    @test 3.99 < i1 <= 4.00001
                    @test r.plant.p_poc<3000.
                    if priority===:balance || (location===:converter && priority===:power)
                        @test i2>2.
                    end
                end
            end
        end
    end
end

@testset "Sequence ratings: converter boundary with an LCL shunt" begin
    for location in (:converter,:poc), target in (:converter,:grid)
        inv=_sequence_limit_plant(;neutral="n",power_limit_location=location,
            i_positive_max=4.,r_filter_grid=.1,x_filter_grid=.2,c_filter_mid=30e-6)
        controller=SequenceController(AverageVoltageVoltVarWatt(),NoUnbalanceControl(),
            CommonScaleLimiter(pq_priority=:var);current_target=target,sensing_reference="n")
        req=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)
        r=solve_controlled_inverter(inv_grid3_bal(),ControlledDevice(inv,controller),req;
            s_base=1e5,solver_options=_SEQ_LIMIT_OPTIONS)
        @test r.solve.publishable
        if r.solve.publishable
            @test 3.99 < r.plant.i_positive <= 4.00001
            @test r.exact_smooth_current_residual < .005
            grid_i1=abs(PowerOptLab._phasor_sequences(r.grid_phase_current)[2])
            @test abs(grid_i1-r.plant.i_positive)>.1
            location===:poc && @test abs(r.plant.q_poc)<.01
        end
    end
end
