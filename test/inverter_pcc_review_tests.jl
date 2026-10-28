# Reproducible steady-state fixture for PR #74's external review. These are
# model-behaviour checks using the declared Australia A curves, not certification.
using BMOPFTools: parse_bmopf

function _pcc_review_grid(v)
    parse_bmopf("""
    {"bus":{"poc":{"terminal_names":["a","b","c","n"],
      "perfectly_grounded_terminals":["n"],"v_min":[170,170,170],"v_max":[280,280,280]}},
     "voltage_source":{"vs":{"bus":"poc","terminal_map":["a","b","c"],
      "v_magnitude":[$v,$v,$v],"v_angle":[0,-2.0943951023931953,2.0943951023931953]}}}
    """; from_string=true)
end
_pcc_review_plant(; kwargs...) = AdvancedInverter(;
    id="pv", bus="poc", phase_terminals=["a","b","c"], neutral="n",
    topology=:THREE_LEG, power_limit_location=:poc, s_max=10e3,
    v_dc=800., c_dc=2e-3, r_filter=.2, x_filter=.5,
    r_filter_grid=.1, x_filter_grid=.2, c_filter_mid=30e-6,
    p_loss_fixed=20., kwargs...)
function _pcc_review_controller(epsilon=.05; headroom=1e-3)
    SequenceController(AverageVoltageVoltVarWatt(
        volt_var=PiecewiseLinearLaw([207.,220.,240.,258.],[.44,0.,0.,-.6];
            smoothing_epsilon=epsilon),
        volt_watt=PiecewiseLinearLaw([253.,260.],[1.,.2];smoothing_epsilon=epsilon),
        volt_watt_basis=:rated), NoUnbalanceControl(),
        CommonScaleLimiter(pq_priority=:var,priority_headroom_fraction=headroom);
        current_target=:converter,sensing_reference="n")
end
const _PCC_REVIEW_OPTIONS = ("tol"=>1e-10,"bound_relax_factor"=>0.,"max_iter"=>1000,"mu_strategy"=>"adaptive")
# Use one fixed solver configuration, including at the saturated pure-P axis;
# the adaptive barrier retains tight accuracy without per-case retries.
# Independent stable softplus hinge and analytic curve oracles.
_review_hinge(x,e) = max(x,0.) + e*log1p(exp(-abs(x)/e))
_review_q(v,e) = .44-.44/13*_review_hinge(v-207,e)+.44/13*_review_hinge(v-220,e)-
    .6/18*_review_hinge(v-240,e)+.6/18*_review_hinge(v-258,e)
_review_p(v,e) = 1.0 - .8/7*_review_hinge(v-253,e)+.8/7*_review_hinge(v-260,e)
_review_q_exact(v) = v <= 207 ? .44 : v < 220 ? .44*(220-v)/13 :
    v <= 240 ? 0. : v < 258 ? -.6*(v-240)/18 : -.6
_review_p_exact(v) = clamp(1-.8*(v-253)/7,.2,1.)

@testset "PCC review: corner bias and epsilon refinement" begin
    for (which, knots, values, jumps) in (
            (:q,[207.,220.,240.,258.],[.44,0.,0.,-.6],[-.44/13,.44/13,-.6/18,.6/18]),
            (:p,[253.,260.],[1.,.2],[-.8/7,.8/7]))
        previous = nothing
        for epsilon in (.05,.025,.0125)
            law = PiecewiseLinearLaw(knots,values;smoothing_epsilon=epsilon)
            errors = Float64[]
            for (v,value,jump) in zip(knots,values,jumps)
                exact = PowerOptLab._pwl_exact(law,v)
                smooth = PowerOptLab._pwl_numeric_smooth(law,v)
                @test exact ≈ value atol=1e-14
                @test smooth ≈ (which === :q ? _review_q(v,epsilon) : _review_p(v,epsilon)) atol=1e-13
                @test smooth-exact ≈ jump*epsilon*log(2) atol=1e-13
                push!(errors,smooth-exact)
            end
            previous === nothing || @test errors ≈ previous ./ 2 atol=1e-13
            previous = errors
        end
    end
    # The rated volt-watt ceiling is rounded upward at its lower corner.
    @test _review_p(260.,.05) ≈ .203960841031771 atol=1e-13
end

@testset "PCC review: Australia A LCL voltage and availability sweep" begin
    inv = _pcc_review_plant()
    controller = _pcc_review_controller()
    # Include every curve knot, interior slopes, endpoints and both sides of
    # the analytic rating/volt-watt crossover. Keep case labels on failures.
    a=(.8/7)^2+(1/30)^2; b=-2*.8/7+26/900; c=169/900
    crossover = 253+(-b-sqrt(b^2-4a*c))/(2a)
    @test crossover ≈ 254.01324903528123 atol=1e-10
    for available in (1.2,.5,0.), v in (180.,207.,220.,230.,240.,250.,253.,
                                      crossover-.1,crossover+.1,258.,260.)
        @testset "availability=$available voltage=$v" begin
            request = InverterControlRequest(p_available=available*inv.s_max,
                p_rated=inv.s_max,q_scale=inv.s_max)
            r = solve_controlled_inverter(_pcc_review_grid(v),ControlledDevice(inv,controller),request;
                solver_options=_PCC_REVIEW_OPTIONS)
            @test r.solve.publishable
            if r.solve.publishable
                q = inv.s_max*_review_q(v,.05)
                # Independently account for losses from component outputs,
                # without replaying either controller allocator.
                losses = r.plant.p_loss+r.plant.p_cap_loss+r.plant.p_filter_loss
                circle = sqrt(inv.s_max^2-q^2)
                ceiling = request.p_rated*_review_p(v,.05)
                expected_p = min(ceiling,request.p_available-losses,circle)
                @test r.plant.q_poc ≈ q atol=.1
                @test r.plant.p_poc ≈ expected_p atol=.1
                @test hypot(r.plant.p_poc,r.plant.q_poc) <= inv.s_max+.02
                @test r.plant.p_dc <= request.p_available+.02
                exact_q = inv.s_max*_review_q_exact(v)
                exact_p = min(request.p_rated*_review_p_exact(v),request.p_available-losses,
                    sqrt(inv.s_max^2-exact_q^2))
                # At curve knots the expected current error is not near zero:
                # derive it from the signed P/Q smoothing errors, in amperes.
                expected_residual = hypot(expected_p-exact_p,q-exact_q)/(3v)
                @test r.exact_smooth_current_residual ≈ expected_residual atol=5e-4
                @test sum(r.plant.poc_phase_voltage .* conj.(r.grid_phase_current)) ≈
                    complex(r.plant.p_poc,r.plant.q_poc) atol=1e-5
                if available == 1.2 && abs(v-crossover) < .11
                    if v < crossover
                        @test r.plant.p_poc ≈ circle atol=.1
                        @test ceiling > circle+10.
                    else
                        @test r.plant.p_poc ≈ ceiling atol=.1
                        @test circle > ceiling+10.
                    end
                end
                if available == 1.2 && v == 250.
                    @test hypot(r.plant.p_conv,r.plant.q_conv) > inv.s_max+20.
                end
                if available == 0.
                    @test r.plant.p_poc < 0. # loss import, not a charging request
                end
            end
        end
    end
end

@testset "PCC review: independent arm ratings and feasible-Q priority" begin
    v = 207.
    controller = _pcc_review_controller()
    request = InverterControlRequest(p_available=24e3,p_rated=20e3,q_scale=10e3)
    results = Dict()
    for (name,ic,ig) in ((:converter,8.,nothing),(:grid,nothing,7.5),(:both,8.,7.5))
        @testset "$name current rating" begin
            inv = _pcc_review_plant(s_max=20e3,i_max=ic,i_grid_max=ig)
            r = solve_controlled_inverter(_pcc_review_grid(v),ControlledDevice(inv,controller),request;
                solver_options=_PCC_REVIEW_OPTIONS)
            @test r.solve.publishable
            if r.solve.publishable
                results[name] = r
                @test r.plant.q_poc ≈ 10e3*_review_q(v,.05) atol=.1
                @test r.plant.p_poc < 10e3
                @test hypot(r.plant.p_poc,r.plant.q_poc) < inv.s_max-100.
                ic === nothing || @test maximum(r.plant.i_mag) <= ic+1e-5
                ig === nothing || @test maximum(r.plant.i_grid_mag) <= ig+1e-5
                exact_q = 10e3*_review_q_exact(v)
                exact_p = Inf
                if ic !== nothing
                    shunt = r.converter_terminal.phase_current .- r.grid_phase_current
                    for k in 1:3
                        local_shunt = shunt[k]*conj(r.plant.poc_phase_voltage[k]/v)
                        exact_p = min(exact_p,3v*(sqrt(ic^2-(imag(local_shunt)-exact_q/(3v))^2)-real(local_shunt)))
                    end
                end
                if ig !== nothing
                    exact_p = min(exact_p,sqrt((3v*ig)^2-exact_q^2))
                end
                @test r.exact_smooth_current_residual ≈
                    hypot(r.plant.p_poc-exact_p,r.plant.q_poc-exact_q)/(3v) atol=5e-4
            end
        end
    end
    if length(results) == 3
        @test maximum(results[:converter].plant.i_mag) ≈ 8. atol=1e-4
        @test maximum(results[:converter].plant.i_grid_mag) > 8.05
        @test maximum(results[:both].plant.i_grid_mag) ≈ 7.5 atol=1e-4
        @test results[:both].plant.p_poc < results[:converter].plant.p_poc-100.
    end
    # A retained compatibility keyword must not silently change PCC dispatch.
    inv = _pcc_review_plant(s_max=20e3,i_max=8.,i_grid_max=7.5)
    changed = solve_controlled_inverter(_pcc_review_grid(v),
        ControlledDevice(inv,_pcc_review_controller(.05;headroom=.2)),request;
        solver_options=_PCC_REVIEW_OPTIONS)
    @test changed.solve.publishable
    if changed.solve.publishable && haskey(results,:both)
        @test changed.plant.p_poc ≈ results[:both].plant.p_poc atol=.02
        @test changed.plant.q_poc ≈ results[:both].plant.q_poc atol=.02
    end
end
