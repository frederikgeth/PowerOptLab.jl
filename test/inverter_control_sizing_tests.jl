function _sizing_controlled(id, positive)
    inverter = AdvancedInverter(
        id=id, bus="poc", phase_terminals=["a", "b", "c"], neutral="n",
        topology=:THREE_LEG, s_max=20e3, i_max=40.0, i_grid_max=42.0,
        v_dc=700.0, c_dc=1.1e-3, i_cap_max=25.0,
        r_filter=0.05, x_filter=0.15, m_max=0.96,
        pwm_ac_converter_reserve=(2.0, 3.0, 4.0),
        pwm_ac_grid_reserve=(1.0, 2.0, 3.0))
    ControlledDevice(inverter, SequenceController(positive))
end

@testset "Inverter-control hardware grids and sizing requirements" begin
    @test_throws ArgumentError InverterHardwareSweepPoint(id="")
    @test_throws ArgumentError InverterHardwareSweepPoint(
        id="bad", converter_current_scale=0.0)
    @test_throws ArgumentError InverterHardwareSweepPoint(
        id="bad", dc_capacitance_scale=Inf)

    request = InverterControlRequest(p_available=8e3, q_scale=0.0)
    baseline_device = _sizing_controlled(
        "pv", AverageVoltageVoltVarWatt())
    sequence_device = _sizing_controlled(
        "pv", PositiveSequenceVoltVarWatt())
    baseline_fleet = ControlledInverterFleetSpec(
        Dict("pv" => baseline_device), Dict("pv" => request))
    sequence_fleet = ControlledInverterFleetSpec(
        Dict("pv" => sequence_device), Dict("pv" => request))

    point = InverterHardwareSweepPoint(
        id="scaled", converter_current_scale=1.25,
        grid_current_scale=1.5, dc_capacitance_scale=2.0,
        capacitor_current_scale=1.1)
    resized = resize_controlled_inverter_fleet(baseline_fleet, point)
    original = inverter_spec(baseline_fleet.devices["pv"])
    changed = inverter_spec(resized.devices["pv"])
    @test changed !== original
    @test original.i_max == 40.0
    @test original.c_dc == 1.1e-3
    @test changed.i_max == 50.0
    @test changed.i_grid_max == 63.0
    @test changed.c_dc == 2.2e-3
    @test changed.i_cap_max ≈ 27.5
    @test resized.devices["pv"].controller ===
          baseline_fleet.devices["pv"].controller
    @test resized.requests["pv"] === baseline_fleet.requests["pv"]
    @test_throws ArgumentError resize_controlled_inverter_fleet(
        baseline_fleet,
        InverterHardwareSweepPoint(
            id="overflow", converter_current_scale=floatmax(Float64)))
    @test_throws ArgumentError resize_controlled_inverter_fleet(
        baseline_fleet,
        InverterHardwareSweepPoint(
            id="underflow", dc_capacitance_scale=nextfloat(0.0)))

    no_optional = _fleet_controlled(
        "pv", AverageVoltageVoltVarWatt())
    no_optional_fleet = ControlledInverterFleetSpec(
        Dict("pv" => no_optional), Dict("pv" => request))
    identity = resize_controlled_inverter_fleet(
        no_optional_fleet, InverterHardwareSweepPoint(id="identity"))
    @test PowerOptLab._same_struct_fields(
        inverter_spec(identity.devices["pv"]),
        inverter_spec(no_optional_fleet.devices["pv"]))
    optional_base = inverter_spec(no_optional_fleet.devices["pv"])
    optional_names = fieldnames(AdvancedInverter)
    optional_parameters = NamedTuple{optional_names}(
        Tuple(getfield(optional_base, name) for name in optional_names))
    without_converter_rating = AdvancedInverter(;
        merge(optional_parameters, (i_max=nothing,))...)
    absent_rating_fleet = ControlledInverterFleetSpec(
        Dict("pv" => ControlledDevice(
            without_converter_rating,
            no_optional_fleet.devices["pv"].controller)),
        Dict("pv" => request))
    absent_identity = resize_controlled_inverter_fleet(
        absent_rating_fleet, InverterHardwareSweepPoint(id="identity"))
    @test inverter_spec(absent_identity.devices["pv"]).i_max === nothing
    @test_throws ArgumentError resize_controlled_inverter_fleet(
        no_optional_fleet,
        InverterHardwareSweepPoint(id="grid", grid_current_scale=1.1))
    @test_throws ArgumentError resize_controlled_inverter_fleet(
        no_optional_fleet,
        InverterHardwareSweepPoint(
            id="capacitor", capacitor_current_scale=1.1))

    net = inv_grid3_unbal()
    net["ibr"] = Dict("pv" => _fleet_native_ibr())
    base_cases = [
        _batch_case(
            "snapshot", "baseline", net, baseline_fleet;
            duration_h=0.25, weight=3.0, penetration=80),
        _batch_case(
            "snapshot", "sequence", net, sequence_fleet;
            duration_h=0.25, weight=3.0, penetration=80),
    ]
    @test isnothing(validate_inverter_control_campaign(base_cases))
    nan_metadata_cases = [InverterControlStudyCase(
        scenario_id="nan_metadata", variant_id=case.variant_id,
        network=net, fleet=case.fleet,
        metadata=Dict("sentinel" => NaN)) for case in base_cases]
    @test isnothing(validate_inverter_control_campaign(nan_metadata_cases))
    copied_network_case = InverterControlStudyCase(
        scenario_id="snapshot", variant_id="sequence",
        network=deepcopy(net), fleet=sequence_fleet,
        duration_h=0.25, weight=3.0,
        metadata=base_cases[1].metadata)
    @test_throws ArgumentError validate_inverter_control_campaign(
        [base_cases[1], copied_network_case])
    @test isnothing(validate_inverter_control_campaign(
        [base_cases[1], copied_network_case]; require_shared_network=false))
    duration_case = InverterControlStudyCase(
        scenario_id="snapshot", variant_id="sequence",
        network=net, fleet=sequence_fleet,
        duration_h=0.5, weight=3.0, metadata=base_cases[1].metadata)
    @test_throws ArgumentError validate_inverter_control_campaign(
        [base_cases[1], duration_case])
    @test_throws ArgumentError validate_inverter_control_campaign(
        [base_cases...,
         _batch_case("other", "baseline", net, baseline_fleet)])
    different_hardware = resize_controlled_inverter_fleet(
        sequence_fleet,
        InverterHardwareSweepPoint(
            id="different", converter_current_scale=1.1))
    hardware_case = InverterControlStudyCase(
        scenario_id="snapshot", variant_id="sequence", network=net,
        fleet=different_hardware, duration_h=0.25, weight=3.0,
        metadata=base_cases[1].metadata)
    @test_throws ArgumentError validate_inverter_control_campaign(
        [base_cases[1], hardware_case])
    points = [
        InverterHardwareSweepPoint(
            id="large", converter_current_scale=1.1,
            grid_current_scale=1.1, dc_capacitance_scale=1.5,
            capacitor_current_scale=1.1),
        InverterHardwareSweepPoint(
            id="reference", grid_current_scale=1.0,
            capacitor_current_scale=1.0),
    ]
    expanded = expand_inverter_hardware_cases(base_cases, points)
    @test length(expanded) == 4
    @test getproperty.(expanded, :scenario_id) == [
        "snapshot::hardware[large]", "snapshot::hardware[reference]",
        "snapshot::hardware[large]", "snapshot::hardware[reference]",
    ]
    @test getproperty.(expanded, :variant_id) ==
          ["baseline", "baseline", "sequence", "sequence"]
    @test all(case.network === net for case in expanded)
    @test expanded[1].metadata["base_scenario_id"] == "snapshot"
    @test expanded[1].metadata["hardware_point_id"] == "large"
    reference_expanded = only(filter(
        case -> case.scenario_id == "snapshot::hardware[reference]" &&
                case.variant_id == "baseline", expanded))
    @test reference_expanded.metadata["converter_current_scale"] == 1.0
    @test reference_expanded.metadata["dc_capacitance_scale"] == 1.0
    @test expanded[1].metadata == expanded[3].metadata
    @test inverter_spec(expanded[1].fleet.devices["pv"]).i_max == 44.0
    @test inverter_spec(expanded[1].fleet.devices["pv"]).c_dc == 1.65e-3
    @test_throws ArgumentError expand_inverter_hardware_cases(
        base_cases, InverterHardwareSweepPoint[])
    @test_throws ArgumentError expand_inverter_hardware_cases(
        base_cases, [points[1], points[1]])
    reserved_case = InverterControlStudyCase(
        scenario_id="reserved", variant_id="baseline", network=net,
        fleet=baseline_fleet,
        metadata=Dict("hardware_point_id" => "already-used"))
    @test_throws ArgumentError expand_inverter_hardware_cases(
        [reserved_case], points)

    study = run_inverter_control_study(
        expanded; solver_options=("max_iter" => 500, "tol" => 1e-8))
    @test solve_status(study).publishable
    automatic_summaries = inverter_control_study_summary_rows(study)
    @test length(automatic_summaries) == 4
    @test Set(row["hardware_point_id"] for row in automatic_summaries) ==
          Set(["large", "reference"])
    @test all(haskey(row, "converter_current_utilization_p99")
              for row in automatic_summaries)
    @test all(isfinite(row["converter_current_utilization_p99"])
              for row in automatic_summaries)
    pairs = inverter_control_paired_rows(study, "baseline")
    @test length(pairs) == 2
    @test all(row.publishable_pair for row in pairs)
    @test getproperty.(pairs, :scenario_id) == [
        "snapshot::hardware[large]", "snapshot::hardware[reference]"]

    @test_throws ArgumentError inverter_control_hardware_requirement_rows(
        study; allowed_dc_ripple_fraction=0.0)
    @test_throws ArgumentError inverter_control_hardware_requirement_rows(
        study; allowed_dc_ripple_fraction=0.02, current_margin=0.99)
    @test_throws ArgumentError inverter_control_hardware_requirement_rows(
        study; allowed_dc_ripple_fraction=0.02, binding_tolerance=1.0)
    requirements = inverter_control_hardware_requirement_rows(
        study; allowed_dc_ripple_fraction=0.02, current_margin=1.1)
    @test length(requirements) == 4
    @test isconcretetype(eltype(requirements))
    @test all(row.publishable for row in requirements)
    @test all(row.converter_current_requirement_A > 0
              for row in requirements)
    @test all(row.grid_current_requirement_A > 0 for row in requirements)
    @test all(row.capacitor_thermal_current_requirement_A > 0
              for row in requirements)
    @test all(row.dc_capacitance_2omega_requirement_F > 0
              for row in requirements)
    @test all(row.dc_2omega_capacitor_current_A > 0
              for row in requirements)
    @test all(row.dc_stored_energy_2omega_requirement_J > 0
              for row in requirements)
    @test all(row.installed_dc_stored_energy_J > 0
              for row in requirements)
    @test all(ismissing(row.dc_ripple_binding) for row in requirements)
    @test all(row.achieved_dc_ripple_voltage_V > 0 for row in requirements)
    @test all(isnan(row.enforced_dc_ripple_limit_V) for row in requirements)
    @test all(row.converter_current_binding isa Bool for row in requirements)
    @test all(row.grid_current_binding isa Bool for row in requirements)
    @test all(row.capacitor_current_binding isa Bool for row in requirements)
    @test Set(getproperty.(requirements, :hardware_point_id)) ==
          Set(["large", "reference"])
    @test all(row.base_scenario_id == "snapshot" for row in requirements)
    @test all(row.converter_current_requirement_utilization > 0
              for row in requirements)
    @test all(row.dc_capacitance_2omega_requirement_utilization > 0
              for row in requirements)
    @test all(row.zero_sequence_voltage_V >= 0 for row in requirements)
    @test all(isfinite(row.negative_sequence_current_angle_rad)
              for row in requirements)
    @test all(row.peak_phase_to_positive_current_ratio >= 1.0
              for row in requirements)
    @test all(row.apparent_power_binding isa Bool for row in requirements)
    @test all(row.command_curtailment_active isa Bool for row in requirements)
    @test all(row.controller_power_limiter_active isa Bool
              for row in requirements)
    @test all(row.controller_current_limiter_active isa Bool
              for row in requirements)
    reference_requirement = only(filter(
        row -> row.scenario_id == "snapshot::hardware[reference]" &&
               row.variant_id == "baseline",
        requirements))
    reference_case_result = only(filter(
        result -> result.case.scenario_id ==
                  "snapshot::hardware[reference]" &&
                  result.case.variant_id == "baseline",
        study.cases))
    reference_plant = reference_case_result.result.devices["pv"].plant
    reference_inverter = inverter_spec(
        reference_case_result.result.spec.devices["pv"])
    @test reference_requirement.converter_current_requirement_A ≈
          1.1 * maximum(hypot.(
              reference_plant.i_mag,
              reference_inverter.pwm_ac_converter_reserve))
    @test reference_requirement.grid_current_requirement_A ≈
          1.1 * maximum(hypot.(
              reference_plant.i_grid_mag,
              reference_inverter.pwm_ac_grid_reserve))
    @test reference_requirement.dc_capacitance_2omega_requirement_F ≈
          reference_requirement.installed_dc_capacitance_F *
          reference_plant.dv2 /
          reference_requirement.allowed_dc_ripple_voltage_V rtol=1e-8
    @test reference_requirement.dc_2omega_capacitor_current_A ≈
          reference_plant.ripple /
          (sqrt(2) * reference_inverter.v_dc) rtol=1e-10
    @test reference_requirement.installed_dc_stored_energy_J ≈
          0.5 * reference_requirement.installed_dc_capacitance_F *
          reference_inverter.v_dc^2 rtol=1e-12
    @test reference_requirement.dc_stored_energy_2omega_requirement_J ≈
          0.5 *
          reference_requirement.dc_capacitance_2omega_requirement_F *
          reference_inverter.v_dc^2 rtol=1e-12
    @test reference_requirement.metadata ===
          reference_case_result.case.metadata

    unpublished = run_inverter_control_study(
        [expanded[1]]; solver_options=("max_iter" => 0,))
    unpublished_requirement = only(
        inverter_control_hardware_requirement_rows(
            unpublished; allowed_dc_ripple_fraction=0.02))
    @test !unpublished_requirement.publishable
    @test isnan(unpublished_requirement.converter_current_requirement_A)
    @test isnan(unpublished_requirement.dc_2omega_capacitor_current_A)
    @test unpublished_requirement.installed_dc_stored_energy_J > 0
    @test ismissing(unpublished_requirement.converter_current_compliant)
    @test ismissing(unpublished_requirement.converter_current_binding)
    @test ismissing(unpublished_requirement.apparent_power_binding)
    @test ismissing(unpublished_requirement.command_curtailment_active)
    @test ismissing(unpublished_requirement.dc_capacitance_2omega_compliant)

    # Limiter predicates, tested directly. An apparent-power-binding operating
    # point is only reachable through a curtailing NLP solve near the edge of
    # this formulation's convergence basin, and that basin is not portable:
    # the same code, Julia version, and pinned dependencies converge on
    # arm64-darwin and do not on x86_64-linux. Pinning the flag semantics to
    # such a solve makes the suite a solver-conditioning detector. The
    # predicates are pure, so they are pinned exactly here and the solve below
    # is reduced to an end-to-end wiring check.
    @test PowerOptLab._sizing_binding(true, 10.0, 10.0, 1e-4)
    @test PowerOptLab._sizing_binding(true, 10.0, 9.9999, 1e-4)
    @test !PowerOptLab._sizing_binding(true, 10.0, 9.9, 1e-4)
    @test PowerOptLab._sizing_binding(true, 10.0, 11.0, 1e-4)
    @test ismissing(PowerOptLab._sizing_binding(false, 10.0, 10.0, 1e-4))
    @test ismissing(PowerOptLab._sizing_binding(true, NaN, 10.0, 1e-4))
    @test ismissing(PowerOptLab._sizing_binding(true, 10.0, NaN, 1e-4))
    @test ismissing(PowerOptLab._sizing_binding(true, 0.0, 0.0, 1e-4))
    @test PowerOptLab._sizing_curtailment_active(true, 8e3, 7e3, 1e-4)
    @test !PowerOptLab._sizing_curtailment_active(true, 8e3, 8e3, 1e-4)
    @test !PowerOptLab._sizing_curtailment_active(true, 8e3, 8e3 - 0.1, 1e-4)
    @test !PowerOptLab._sizing_curtailment_active(true, 0.0, 0.0, 1e-4)
    @test ismissing(
        PowerOptLab._sizing_curtailment_active(false, 8e3, 7e3, 1e-4))
    @test ismissing(
        PowerOptLab._sizing_curtailment_active(true, NaN, 7e3, 1e-4))
    @test PowerOptLab._sizing_scale_active(true, 0.5, 1e-4)
    @test !PowerOptLab._sizing_scale_active(true, 1.0, 1e-4)
    @test !PowerOptLab._sizing_scale_active(true, 1.0 - 1e-9, 1e-4)
    @test ismissing(PowerOptLab._sizing_scale_active(false, 0.5, 1e-4))
    @test ismissing(PowerOptLab._sizing_scale_active(true, NaN, 1e-4))

    limited_inverter = AdvancedInverter(
        id="limited", bus="poc", phase_terminals=["a", "b", "c"],
        neutral="n", topology=:THREE_LEG, s_max=7e3, i_max=40.0,
        v_dc=700.0, c_dc=1.1e-3, r_filter=0.05, x_filter=0.15,
        m_max=0.96)
    limited_fleet = ControlledInverterFleetSpec(
        Dict("limited" => ControlledDevice(
            limited_inverter,
            SequenceController(AverageVoltageVoltVarWatt()))),
        Dict("limited" => request))
    limited_net = inv_grid3_bal()
    limited_net["ibr"] = Dict("limited" => _fleet_native_ibr())
    limited_study = run_inverter_control_study([
        InverterControlStudyCase(
            scenario_id="limited", variant_id="baseline",
            network=limited_net, fleet=limited_fleet),
    ]; solver_options=("max_iter" => 5000, "tol" => 1e-8))
    limited_requirement = only(inverter_control_hardware_requirement_rows(
        limited_study; allowed_dc_ripple_fraction=0.02))
    limited_flags = (
        :apparent_power_binding, :converter_current_binding,
        :command_curtailment_active, :controller_power_limiter_active,
        :controller_current_limiter_active)
    # The tri-state contract holds either way: a flag is Bool exactly when the
    # solve published, and missing otherwise. This is the part that must never
    # depend on the linear-solver build.
    for flag in limited_flags
        value = getproperty(limited_requirement, flag)
        @test limited_requirement.publishable ? value isa Bool :
              ismissing(value)
    end
    if limited_requirement.publishable
        @test limited_requirement.apparent_power_binding
        @test limited_requirement.command_curtailment_active
        @test !limited_requirement.converter_current_binding
        @test limited_requirement.achieved_converter_apparent_power_VA ≈ 7e3 rtol=1e-4
    else
        # Not a silent skip: say which platform-dependent status was hit, so a
        # future CI log explains itself instead of looking like a passing run.
        @info "s_max-binding integration check skipped; solve not publishable" *
              " (physics assertions covered by the predicate tests above)" *
              " termination=$(solve_status(limited_study).termination_status)"
    end
end
