@testset "BMOPFTools native power API and controlled reactive capability" begin
    B = PowerOptLab.BMOPFTools
    cp = ConnectionPoint(id="pv", bus="poc", ibr_id="pv", export_max=6000.)
    for pf in (0.8, -0.8), per_unit in (false, true)
        net = inv_grid()
        net["ibr"] = Dict("pv" => Dict{String,Any}(
            "bus"=>"poc", "terminal_map"=>["1", "n"],
            "topology"=>"SINGLE_PHASE", "prime_mover"=>"PV",
            "p_min"=>[0.], "p_max"=>[6000.],
            "q_min"=>[-1000.], "q_max"=>[1000.],
            "control_profile"=>"pf"))
        net["control_profile"] = Dict("pf" => Dict(
            "power_factor" => Dict("pf"=>pf)))

        # Power access must work without s_max and use the native object itself.
        ctx = B.build_opf_model(net; per_unit, add_objective=false)
        @test PowerOptLab._ibr_active_power(ctx, cp) ===
              B.opf_object(ctx, B.opf_ibr_power_key("pv", 1))

        result = solve_operating_envelope(net, [cp]; per_unit, fairness=:sum)
        @test solve_status(result).publishable
        dispatch = result.snapshots[1]["ibr"]["pv"]["1"]
        # Q = -sign(pf)*0.75P and |Q| <= 1000 var, hence P <= 4000/3 W.
        # Before upstream #438 the PF law bypassed these Q bounds entirely.
        @test result.envelope["pv"][1] ≈ 4000/3 atol=0.1
        @test dispatch["pg"] ≈ result.envelope["pv"][1] atol=0.1
        @test dispatch["qg"] ≈ -sign(pf)*1000 atol=0.1

        # Upstream minimum PF composes with a compatible prescribed law, and
        # rejects contradictory requests rather than dropping either policy.
        net["extras"] = Dict("bmopftools"=>Dict("minimum_power_factor"=>Dict(
            "schema"=>"bmopftools.minimum_power_factor/v1",
            "ibr"=>Dict("pv"=>Dict("pf_min"=>0.8,
                "direction"=>"generation", "scope"=>"per_coil")))))
        compatible = solve_operating_envelope(net, [cp]; per_unit, fairness=:sum)
        @test solve_status(compatible).publishable
        @test compatible.envelope["pv"][1] ≈ result.envelope["pv"][1] atol=0.1
        net["extras"]["bmopftools"]["minimum_power_factor"]["ibr"]["pv"]["pf_min"] = 0.95
        @test_throws ArgumentError solve_operating_envelope(net, [cp]; per_unit)
    end

    # FOUR_LEG power follows native phase order with a non-trailing neutral.
    net = inv_grid3_bal()
    net["ibr"] = Dict("pv"=>Dict{String,Any}(
        "bus"=>"poc", "terminal_map"=>["c", "n", "a", "b"],
        "topology"=>"FOUR_LEG", "prime_mover"=>"PV",
        "p_min"=>zeros(3), "p_max"=>fill(2000., 3),
        "q_min"=>zeros(3), "q_max"=>zeros(3)))
    ctx = B.build_opf_model(net; add_objective=false)
    cp3 = ConnectionPoint(id="pv", bus="poc", ibr_id="pv",
        phase_terminals=["c", "a", "b"], neutral="n", export_max=6000.)
    power = PowerOptLab._ibr_active_power(ctx, cp3)
    @test isequal(power, sum(B.opf_object(ctx, B.opf_ibr_power_key("pv", k)) for k in 1:3))
end
