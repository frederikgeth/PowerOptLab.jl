import BMOPFTools

function _fleet_native_ibr(; bus="poc", terminal_map=["a", "b", "c", "n"],
                           topology="FOUR_LEG", p=0.0, s_max=1.0)
    nphase = topology == "SINGLE_PHASE" ? 1 : 3
    Dict{String,Any}(
        "bus" => bus,
        "terminal_map" => terminal_map,
        "topology" => topology,
        "prime_mover" => "PV",
        "s_max" => fill(Float64(s_max), nphase),
        "p_min" => fill(Float64(p), nphase),
        "p_max" => fill(Float64(p), nphase),
        "q_min" => zeros(nphase),
        "q_max" => zeros(nphase),
    )
end

function _fleet_controlled(id, positive; neutral="n", c_filter_mid=0.0,
                           r_filter_grid=0.0, b_filter_shunt=0.0)
    inverter = AdvancedInverter(
        id=id, bus="poc", phase_terminals=["a", "b", "c"], neutral=neutral,
        topology=:THREE_LEG, s_max=20e3, i_max=40.0,
        v_dc=700.0, c_dc=1.1e-3, r_filter=0.05, x_filter=0.15,
        c_filter_mid=c_filter_mid, r_filter_grid=r_filter_grid,
        b_filter_shunt=b_filter_shunt, m_max=0.96)
    ControlledDevice(inverter, SequenceController(positive))
end

@testset "Controlled-inverter fleet: specification validation" begin
    controlled = _fleet_controlled("pv", AverageVoltageVoltVarWatt())
    request = InverterControlRequest(p_available=8e3, q_scale=0.0)

    @test_throws ArgumentError ControlledInverterFleetSpec(Dict(), Dict())
    @test_throws ArgumentError ControlledInverterFleetSpec(
        Dict("pv" => controlled), Dict("other" => request))
    @test_throws ArgumentError ControlledInverterFleetSpec(
        Dict("other" => controlled), Dict("other" => request))

    spec = ControlledInverterFleetSpec(
        Dict("pv" => controlled), Dict("pv" => request))
    @test_throws ArgumentError PowerOptLab._validate_controlled_fleet(
        inv_grid3_bal(), spec)

    wrong_bus = inv_grid3_bal()
    wrong_bus["ibr"] = Dict("pv" => _fleet_native_ibr(bus="grid"))
    @test_throws ArgumentError PowerOptLab._validate_controlled_fleet(
        wrong_bus, spec)

    wrong_order = inv_grid3_bal()
    wrong_order["ibr"] = Dict(
        "pv" => _fleet_native_ibr(terminal_map=["b", "a", "c", "n"]))
    @test_throws ArgumentError PowerOptLab._validate_controlled_fleet(
        wrong_order, spec)

    dc_coupled = inv_grid3_bal()
    dc_coupled["ibr"] = Dict("pv" => merge(
        _fleet_native_ibr(), Dict{String,Any}("dc_bus" => "dc")))
    @test_throws ArgumentError PowerOptLab._validate_controlled_fleet(
        dc_coupled, spec)
    three_leg = inv_grid3_bal()
    three_leg["ibr"] = Dict("pv" => _fleet_native_ibr(
        terminal_map=["a", "b", "c"], topology="THREE_LEG"))
    @test_throws ArgumentError PowerOptLab._validate_controlled_fleet(
        three_leg, spec)
    delta_controlled = _fleet_controlled(
        "pv", AverageVoltageVoltVarWatt(); neutral=nothing)
    delta_spec = ControlledInverterFleetSpec(
        Dict("pv" => delta_controlled), Dict("pv" => request))
    @test PowerOptLab._validate_controlled_fleet(
        three_leg, delta_spec)["pv"] == 3
    for shunted in (
        _fleet_controlled("pv", AverageVoltageVoltVarWatt();
                          neutral=nothing, c_filter_mid=1e-6,
                          r_filter_grid=0.01),
        _fleet_controlled("pv", AverageVoltageVoltVarWatt();
                          neutral=nothing, b_filter_shunt=1e-6),
    )
        shunted_spec = ControlledInverterFleetSpec(
            Dict("pv" => shunted), Dict("pv" => request))
        @test_throws ArgumentError PowerOptLab._validate_controlled_fleet(
            three_leg, shunted_spec)
    end
    @test_throws ArgumentError solve_controlled_inverter_fleet(
        merge(inv_grid3_bal(), Dict("ibr" => Dict(
            "pv" => _fleet_native_ibr()))), spec; per_unit=false)
    @test_throws ArgumentError solve_controlled_inverter_fleet(
        merge(inv_grid3_bal(), Dict("ibr" => Dict(
            "pv" => _fleet_native_ibr()))), spec;
        selection_objective=:unsupported)

    # Masking can narrow a homogeneous top-level result dictionary. Publishing
    # the pruned network must widen it explicitly before replacing `"ibr"`.
    narrowed = Dict(
        "bus" => Dict("poc" => 1),
        "ibr" => Dict("pv" => 1, "native" => 2),
    )
    @test narrowed isa Dict{String,Dict{String,Int}}
    pruned = PowerOptLab._network_without_replaced_ibrs(narrowed, ("pv",))
    @test pruned isa Dict{String,Any}
    @test pruned["ibr"] == Dict{String,Any}("native" => 2)
end

@testset "Controlled-inverter fleet: unpublished and neutral-arity regressions" begin
    request = InverterControlRequest(p_available=8e3, q_scale=0.0)
    controlled = _fleet_controlled("pv", AverageVoltageVoltVarWatt())
    spec = ControlledInverterFleetSpec(
        Dict("pv" => controlled), Dict("pv" => request))
    net = inv_grid3_bal()
    net["ibr"] = Dict("pv" => _fleet_native_ibr())

    limited = solve_controlled_inverter_fleet(
        net, spec; solver_options=("max_iter" => 0,))
    @test !solve_status(limited).publishable
    @test isnan(limited.devices["pv"].plant.p_poc)
    @test all(isnan(phase["vm"])
              for phase in values(limited.network["bus"]["poc"]))

    # BMOPFTools resolves neutral labels when declaring native IBR variables.
    # Rename the grounded neutral to an intentionally non-conventional label:
    # BMOPFTools then declares four semantic current pairs for the selected
    # FOUR_LEG placeholder. The custom owner must fix all four, including when
    # the native-cost expression still visits those variables.
    unusual = inv_grid3_bal()
    for bus in values(unusual["bus"])
        bus["terminal_names"] = replace.(bus["terminal_names"], "n" => "0")
        bus["perfectly_grounded_terminals"] = replace.(
            bus["perfectly_grounded_terminals"], "n" => "0")
        haskey(bus, "v_min") && push!(bus["v_min"], last(bus["v_min"]))
        haskey(bus, "v_max") && push!(bus["v_max"], last(bus["v_max"]))
    end
    for line in values(unusual["line"])
        line["terminal_map_from"] = replace.(
            line["terminal_map_from"], "n" => "0")
        line["terminal_map_to"] = replace.(
            line["terminal_map_to"], "n" => "0")
    end
    unusual_native = _fleet_native_ibr(
        terminal_map=["a", "b", "c", "0"])
    unusual_native["cost"] = ones(4)
    unusual["ibr"] = Dict("pv" => unusual_native)
    observed_arity = Ref(0)
    arity_builder = BMOPFTools.OpfDeviceBuilder(
        :PowerOptLab, function (ctx, ids)
            observed_arity[] = PowerOptLab._native_ibr_declared_current_count(
                ctx, only(ids))
            PowerOptLab._disable_native_ibr_port!(ctx, only(ids))
            ctx
        end)
    arity_spec = BMOPFTools.OpfBuildSpec(component_builders=Dict(
        (:ibr, "pv") => arity_builder))
    arity_ctx = BMOPFTools.build_opf_model(
        unusual; build_spec=arity_spec, add_objective=false)
    @test observed_arity[] == 4
    @test all(JuMP.is_fixed(PowerOptLab._opf_ibr_current(
                  arity_ctx, "pv", phase; component=component))
              for phase in 1:4, component in (:real, :imag))

    unusual_controlled = _fleet_controlled(
        "pv", AverageVoltageVoltVarWatt(); neutral="0")
    unusual_spec = ControlledInverterFleetSpec(
        Dict("pv" => unusual_controlled), Dict("pv" => request))
    unusual_result = solve_controlled_inverter_fleet(
        unusual, unusual_spec; selection_objective=:network_cost,
        solver_options=("max_iter" => 500, "tol" => 1e-8))
    @test solve_status(unusual_result).publishable
    @test unusual_result.devices["pv"].plant.p_poc > 1000.0
end

@testset "Controlled-inverter fleet: replacement, coexistence, and extraction" begin
    law = PiecewiseLinearLaw(
        [230.0, 240.0, 250.0], [1.0, 1.0, 0.2];
        smoothing_epsilon=0.05)
    average = _fleet_controlled(
        "pv_average", AverageVoltageVoltVarWatt(volt_watt=law))
    worst = _fleet_controlled(
        "pv_worst", WorstPhaseVoltVarWatt(
            volt_watt=law, extrema_epsilon=0.05))
    request = InverterControlRequest(p_available=8e3, q_scale=0.0)

    net = inv_grid3_unbal()
    # The selected native records deliberately permit essentially no output.
    # A successful multi-kW result therefore discriminates replacement from a
    # second stamp constrained in parallel by the native IBR equations.
    net["ibr"] = Dict(
        "pv_worst" => _fleet_native_ibr(s_max=1.0),
        "pv_native" => _fleet_native_ibr(p=100.0, s_max=500.0),
        "pv_average" => _fleet_native_ibr(s_max=1.0),
    )
    spec = ControlledInverterFleetSpec(
        Dict("pv_worst" => worst, "pv_average" => average),
        Dict("pv_worst" => request, "pv_average" => request))
    result = solve_controlled_inverter_fleet(
        net, spec; solver_options=("max_iter" => 500, "tol" => 1e-8))

    @test solve_status(result).publishable
    @test result.build_manifest.component_owners[(:ibr, "pv_average")] ==
          :PowerOptLab
    @test result.build_manifest.component_owners[(:ibr, "pv_worst")] ==
          :PowerOptLab
    @test result.build_manifest.component_owners[(:ibr, "pv_native")] ==
          :BMOPFTools
    @test Set(keys(result.devices)) == Set(("pv_average", "pv_worst"))
    @test result.spec !== spec
    @test Set(keys(result.network["ibr"])) == Set(("pv_native",))
    @test result.devices["pv_average"].bus === result.network["bus"]
    @test result.devices["pv_worst"].bus === result.network["bus"]
    @test sum(phase["pg"] for phase in
              values(result.network["ibr"]["pv_native"])) ≈ 300.0 atol=1e-4
    @test result.devices["pv_average"].plant.p_poc > 1000.0
    @test result.devices["pv_worst"].plant.p_poc > 1000.0

    # Both laws see the same local voltage state. The average law misses the
    # high-voltage phase, while the worst-phase guard materially curtails.
    @test result.devices["pv_worst"].control.p_request <
          result.devices["pv_average"].control.p_request - 1000.0

    rows = controlled_inverter_rows(result)
    phase_rows = controlled_inverter_phase_rows(result)
    @test getproperty.(rows, :device_id) == ["pv_average", "pv_worst"]
    @test all(getproperty.(rows, :p_available_W) .== 8e3)
    @test all(getproperty.(rows, :dc_capacitance_F) .== 1.1e-3)
    @test all(isapprox(row.dc_stored_energy_J, 0.5*1.1e-3*700.0^2)
              for row in rows)
    @test length(phase_rows) == 6
    @test getproperty.(phase_rows, :device_id) ==
          repeat(["pv_average", "pv_worst"], inner=3)
    @test getproperty.(phase_rows, :phase) == [1, 2, 3, 1, 2, 3]

    diagnostics = solve_diagnostics(result)
    @test diagnostics.controlled_device_count == 2
    @test diagnostics.total_p_poc ≈
          sum(device.plant.p_poc for device in values(result.devices))
    @test diagnostics.maximum_converter_current <= 40.0 + 1e-4

    for device in values(result.devices)
        @test abs(sum(device.control.phase_current)) < 1e-6
        @test maximum(device.plant.i_mag) <= 40.0 + 1e-4
        @test real(device.converter_terminal.total_power) ≈
              device.plant.p_conv atol=1e-5
        @test imag(device.converter_terminal.total_power) ≈
              device.plant.q_conv atol=1e-5
        @test device.plant.p_conv ≈
              device.plant.p_poc + device.plant.p_filter_loss atol=1e-4
        @test device.plant.p_dc ≈
              device.plant.p_conv + device.plant.p_loss +
              device.plant.p_cap_loss atol=1e-4
        @test device.exact_smooth_current_residual < 0.1
    end
end
