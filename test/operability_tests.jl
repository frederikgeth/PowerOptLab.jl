using BMOPFTools
using LinearAlgebra
using SparseArrays

@testset "post-OPF operability checker" begin
    net = single_bus_net(pload=100.0)
    pf = solve_pf(net; per_unit=false)

    scope_audit = operability_scope_audit(net)
    @test scope_audit["status"] == :supported
    @test scope_audit["closure"] == :frozen_dispatch
    @test scope_audit["control_closure"] == "frozen_dispatch_native_static"
    @test scope_audit["generator_count"] == 0
    @test scope_audit["ibr_count"] == 0
    @test scope_audit["topology"]["has_voltage_source"] === true
    @test scope_audit["topology"]["voltage_source_count"] == 1
    @test scope_audit["topology"]["source_buses"] == ["sourcebus"]
    @test isempty(scope_audit["topology"]["missing_source_buses"])

    no_source_net = deepcopy(net)
    delete!(no_source_net, "voltage_source")
    no_source_audit = operability_scope_audit(no_source_net)
    @test no_source_audit["status"] == :not_applicable
    @test no_source_audit["topology"]["has_voltage_source"] === false
    @test no_source_audit["topology"]["voltage_source_count"] == 0
    @test any(occursin("no voltage source", reason)
              for reason in no_source_audit["unsupported_reasons"])

    dangling_source_net = deepcopy(net)
    dangling_source_net["voltage_source"]["vs"]["bus"] = "missing_bus"
    dangling_source_audit = operability_scope_audit(dangling_source_net)
    @test dangling_source_audit["status"] == :not_applicable
    @test dangling_source_audit["topology"]["has_voltage_source"] === false
    @test dangling_source_audit["topology"]["missing_source_buses"] == ["missing_bus"]
    @test any(occursin("missing bus", reason)
              for reason in dangling_source_audit["unsupported_reasons"])

    @test_throws ArgumentError check_opf_operability(net, pf)

    spec = OperabilitySpec(
        scaling_policy=SIUnitsScaling(),
        voltage_min=800.0,
        voltage_max=1100.0,
    )
    report = check_opf_operability(net, pf; spec)
    @test report.status == :pass
    @test report.checks["endpoint"].status == :pass
    @test report.checks["endpoint"].value == report.endpoint_residual_normalized
    @test report.endpoint_residual < 1e-6
    @test report.checks["jacobian_regular"].status == :pass
    @test report.checks["load_scale_sensitivity"].status == :pass
    connection = report.load_connections["ld1/1"]
    @test connection["requested_power"] ≈ connection["realized_power"] atol=1e-8
    load_scale_connection = report.sensitivities["load_scale"]["load_connections"]["ld1/1"]
    @test load_scale_connection["magnitude_derivative"] < 0.0
    @test isfinite(load_scale_connection["path_dP_dV"])
    @test load_scale_connection["path_dP_dV"] < 0.0
    @test report.branch_evidence["dP_dV"]["connections"]["ld1/1"]["classification"] ==
          "negative_high_side_indicator"
    @test haskey(report.sensitivities["directions"], "P")
    @test haskey(report.sensitivities["directions"], "Q")
    p_direction = report.sensitivities["directions"]["P"]["ld1/1"]
    @test p_direction["units"] == "W"
    @test p_direction["load_connections"]["ld1/1"]["magnitude_derivative"] < 0.0
    @test report.provenance["operability"]["scope"] == "static_ybus_linearized"
    @test report.provenance["operability"]["closure"] == "frozen_dispatch"
    @test report.provenance["operability"]["model_inventory"]["load_models"] ==
          ["constant_power"]
    snapshot_row = operability_snapshot_row(report; snapshot_id="base")
    @test snapshot_row.snapshot_id == "base"
    @test snapshot_row.status == :pass
    @test snapshot_row.scope == "single_snapshot_static_ybus"
    @test snapshot_row.minimum_terminal_voltage ≈ report.load_connections["ld1/1"]["magnitude"]
    @test snapshot_row.high_side_indicator_count == 1
    @test snapshot_row.fixed_point_certificate_status == :not_applicable
    @test snapshot_row.jacobian_spectrum_mode == "full"
    @test snapshot_row.jacobian_storage_mode == "dense"
    @test snapshot_row.jacobian_nonzero_count > 0
    @test snapshot_row.jacobian_storage_bytes_estimate > 0
    @test snapshot_row.zbus_storage_mode == "not_requested"
    # Scaling changes coordinates, not the physical snapshot conclusion.  The
    # checker receives the same SI-valued solution under both policies; the
    # normalized residual is expected to change with the declared current
    # base, while physical voltages, branch indicators, and claims must not.
    classic_spec = OperabilitySpec(
        scaling_policy=ClassicPerUnitScaling(1.0e6),
        scaling_bases=Dict(
            "sourcebus" => (voltage=1000.0, current=1000.0),
            "bus1" => (voltage=1000.0, current=1000.0)),
        voltage_min=800.0,
        voltage_max=1100.0,
        compute_fixed_point_certificate=true)
    classic_report = check_opf_operability(net, pf; spec=classic_spec)
    @test classic_report.status == report.status
    @test classic_report.checks["endpoint"].status == report.checks["endpoint"].status
    @test classic_report.checks["terminal_voltage_bounds"].status ==
          report.checks["terminal_voltage_bounds"].status
    @test classic_report.checks["fixed_point_certificate"].status ==
          :pass
    @test classic_report.load_connections["ld1/1"]["magnitude"] ≈
          report.load_connections["ld1/1"]["magnitude"] atol=1e-10
    @test classic_report.branch_evidence["dP_dV"]["connections"]["ld1/1"]["classification"] ==
          report.branch_evidence["dP_dV"]["connections"]["ld1/1"]["classification"]
    @test classic_report.branch_evidence["dP_dV"]["connections"]["ld1/1"]["path_dP_dV"] ≈
          report.branch_evidence["dP_dV"]["connections"]["ld1/1"]["path_dP_dV"] atol=1e-8
    @test classic_report.endpoint_residual_normalized ≈
          report.endpoint_residual_normalized / 1000.0 atol=1e-12
    @test_throws ArgumentError check_opf_operability(net, pf;
        spec=OperabilitySpec(scaling_policy=ClassicPerUnitScaling(1.0e6),
            scaling_bases=Dict("bus1" => (voltage=1000.0, current=1000.0))))

    for policy in (
            ConsistentPerUnitScaling(s_base=1.0e6,
                voltage_bases=Dict("sourcebus" => 1000.0, "bus1" => 1000.0)),
            ZonePerUnitScaling(
                voltage_bases=Dict("sourcebus" => 1000.0, "bus1" => 1000.0),
                power_bases=Dict("sourcebus" => 1.0e6, "bus1" => 1.0e6)))
        scaled = check_opf_operability(net, pf; spec=OperabilitySpec(
            scaling_policy=policy,
            scaling_bases=Dict(
                "sourcebus" => (voltage=1000.0, current=1000.0),
                "bus1" => (voltage=1000.0, current=1000.0)),
            voltage_min=800.0, voltage_max=1100.0,
            compute_fixed_point_certificate=true))
        @test scaled.status == report.status
        @test scaled.checks["terminal_voltage_bounds"].status ==
              report.checks["terminal_voltage_bounds"].status
        @test scaled.checks["fixed_point_certificate"].status == :pass
        @test scaled.load_connections["ld1/1"]["magnitude"] ≈
              report.load_connections["ld1/1"]["magnitude"] atol=1e-10
        @test scaled.branch_evidence["dP_dV"]["connections"]["ld1/1"]["classification"] ==
              report.branch_evidence["dP_dV"]["connections"]["ld1/1"]["classification"]
    end

    certificate_report = check_opf_operability(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_fixed_point_certificate=true))
    snapshot_rows = operability_snapshot_rows(Dict(
        "z_certificate" => certificate_report,
        "a_base" => report))
    @test [row.snapshot_id for row in snapshot_rows] == ["a_base", "z_certificate"]
    @test snapshot_rows[1].status == report.status
    @test snapshot_rows[2].fixed_point_certificate_status == :pass
    @test_throws ArgumentError operability_snapshot_rows([report])
    @test certificate_report.checks["fixed_point_certificate"].status == :pass
    @test certificate_report.checks["fixed_point_local_region"].status == :pass
    @test certificate_report.checks["fixed_point_euclidean_region"].status == :pass
    certificate = certificate_report.branch_evidence["fixed_point_certificate"]
    @test certificate["method"] == "bernstein_style_zbus_contraction"
    @test certificate["candidate_inside_region"] === true
    @test certificate["contraction_factor"] < 1.0
    @test certificate["candidate_distance"] <= certificate["radius"]
    @test certificate["selected_region"] == "componentwise"
    @test length(certificate["region_radii"]) == length(report.state_nodes)
    @test haskey(certificate, "uniform_region")
    @test haskey(certificate, "componentwise_region")
    @test certificate["law_bound_validation"]["status"] == :pass
    @test maximum(certificate["law_bound_validation"]["connections"]["ld1/1"]["ratio"]) <= 1.001
    @test certificate["local_candidate_region"]["status"] == :pass
    @test certificate["local_candidate_region"]["condition_margin"] > 0.0
    @test certificate["euclidean_region"]["status"] == :pass
    @test certificate["euclidean_region"]["condition_margin"] > 0.0
    complexity = certificate_report.branch_evidence["complexity"]
    @test complexity["free_node_count"] == length(certificate_report.state_nodes)
    @test complexity["real_state_dimension"] == length(certificate_report.state)
    @test complexity["load_connection_count"] == 1
    @test complexity["sensitivity_direction_count"] == 3
    @test complexity["linear_solver_factorization_reused"] === true
    @test complexity["jacobian_storage_bytes_dense"] ==
          sizeof(Float64) * length(certificate_report.state)^2
    @test complexity["zbus_storage_bytes_dense"] > 0
    @test complexity["zbus_storage_mode"] in
          ("implicit_dense_factorization", "implicit_sparse_factorization")
    @test complexity["fixed_point_scan_points_per_geometry"] == 33
    @test complexity["fixed_point_refinement_iterations"] == 12
    @test complexity["fixed_point_scan_strategy"] == "coarse_grid_bisection"
    @test certificate_report.status == :pass

    reduced_spectrum_report = check_opf_operability(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            voltage_min=800.0, voltage_max=1100.0,
            jacobian_spectrum=:extremes))
    @test reduced_spectrum_report.status == :pass
    @test length(reduced_spectrum_report.singular_values) == 2
    @test reduced_spectrum_report.singular_values[1] >=
          reduced_spectrum_report.singular_values[2] > 0.0
    @test reduced_spectrum_report.branch_evidence["critical_mode"]["status"] ==
          :not_applicable
    @test reduced_spectrum_report.branch_evidence["complexity"]["jacobian_spectrum_mode"] ==
          "extremes"
    sparse_spectrum_report = check_opf_operability(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            voltage_min=800.0, voltage_max=1100.0,
            compute_sensitivity=false,
            jacobian_spectrum=:extremes, jacobian_storage=:sparse))
    @test sparse_spectrum_report.status == :pass
    @test issparse(sparse_spectrum_report.jacobian)
    @test sparse_spectrum_report.singular_values ≈
          reduced_spectrum_report.singular_values rtol=1e-6
    @test sparse_spectrum_report.branch_evidence["complexity"]["jacobian_storage_mode"] ==
          "sparse"
    @test sparse_spectrum_report.branch_evidence["complexity"]["jacobian_nonzero_count"] > 0
    sparse_pattern = sparse_spectrum_report.branch_evidence["complexity"]["jacobian_pattern"]
    @test length(sparse_pattern["rows"]) ==
          sparse_spectrum_report.branch_evidence["complexity"]["jacobian_nonzero_count"]
    @test length(sparse_pattern["columns"]) == length(sparse_pattern["rows"])
    sparse_unrecorded = check_opf_operability(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            voltage_min=800.0, voltage_max=1100.0,
            compute_sensitivity=false, jacobian_spectrum=:extremes,
            jacobian_storage=:sparse, record_jacobian_pattern=false))
    @test sparse_unrecorded.status == :pass
    @test sparse_unrecorded.branch_evidence["complexity"]["jacobian_pattern"] === nothing
    @test sparse_unrecorded.branch_evidence["complexity"]["jacobian_pattern_recorded"] === false
    sparse_certificate_report = check_opf_operability(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            voltage_min=800.0, voltage_max=1100.0,
            compute_sensitivity=false,
            compute_fixed_point_certificate=true,
            jacobian_spectrum=:extremes, jacobian_storage=:sparse))
    @test sparse_certificate_report.status == :pass
    @test issparse(sparse_certificate_report.jacobian)
    @test sparse_certificate_report.checks["fixed_point_certificate"].status == :pass
    @test sparse_certificate_report.branch_evidence["fixed_point_certificate"]["zbus_storage_mode"] in
          ("implicit_sparse_factorization", "implicit_dense_factorization")
    certificate_row = operability_snapshot_row(certificate_report)
    @test certificate_row.fixed_point_certificate_status == :pass
    @test certificate_row.fixed_point_condition_margin > 0.0
    @test certificate_row.fixed_point_local_region_status == :pass
    @test certificate_row.fixed_point_local_condition_margin > 0.0
    @test certificate_row.fixed_point_euclidean_region_status == :pass
    @test certificate_row.fixed_point_euclidean_condition_margin > 0.0
    @test certificate_row.jacobian_spectrum_mode == "full"
    @test certificate_row.jacobian_storage_mode == "dense"
    @test certificate_row.jacobian_storage_bytes_estimate > 0
    @test certificate_row.zbus_storage_mode in
          ("implicit_dense_factorization", "implicit_sparse_factorization")
    direction = OperabilityStressDirection(:phase_selective;
        p_scale=1.0, q_scale=0.0, connection_weights=Dict("ld1" => [0.0]))
    stressed = operability_stress_network(net, 0.5, direction)
    @test stressed["load"]["ld1"]["p_nom"] == [0.0]
    @test stressed["load"]["ld1"]["q_nom"] == [0.0]
    stress_rows = operability_stress_rows(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_fixed_point_certificate=true),
        directions=[OperabilityStressDirection(:uniform),
            OperabilityStressDirection(:reactive_removed; p_scale=1.0, q_scale=0.0)],
        lambdas=[1.0, 0.0], solve=n -> solve_pf(n; per_unit=false))
    @test length(stress_rows) == 4
    @test [row.lambda for row in stress_rows] == [0.0, 1.0, 0.0, 1.0]
    @test all(row.status == :pass for row in stress_rows)
    @test all(row.endpoint_status == :pass for row in stress_rows)
    vararg_stress_rows = operability_stress_rows(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_fixed_point_certificate=true),
        directions=[OperabilityStressDirection(:vararg)], lambdas=[1.0],
        solve=(args...) -> solve_pf(args[1]; per_unit=false), solve_arity=1)
    @test only(vararg_stress_rows).status == :pass
    three_arg_stress_rows = operability_stress_rows(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_fixed_point_certificate=true),
        directions=[OperabilityStressDirection(:three_arg)], lambdas=[1.0],
        solve=(network, previous, direction) -> solve_pf(network; per_unit=false))
    @test only(three_arg_stress_rows).status == :pass
    int32_arity_rows = operability_stress_rows(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_fixed_point_certificate=true),
        directions=[OperabilityStressDirection(:int32_arity)], lambdas=[1.0],
        solve=(args...) -> solve_pf(args[1]; per_unit=false), solve_arity=Int32(1))
    @test only(int32_arity_rows).status == :pass
    bad_solution = deepcopy(pf)
    bad_solution["bus"]["bus1"]["1"]["vr"] = NaN
    model_error_rows = operability_stress_rows(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()),
        directions=[OperabilityStressDirection(:bad_data)], lambdas=[1.0],
        solve=network -> bad_solution)
    @test only(model_error_rows).status == :inconclusive
    @test occursin("not finite", only(model_error_rows).error)
    @test_throws UndefVarError operability_stress_rows(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()),
        directions=[OperabilityStressDirection(:broken)], lambdas=[1.0],
        solve=network -> definitely_missing_solver_name(network))
    @test_throws ArgumentError operability_stress_rows(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()),
        directions=[OperabilityStressDirection(:bad_arity)], lambdas=[1.0],
        solve=network -> solve_pf(network; per_unit=false), solve_arity=5)
    stress_summary = operability_stress_summary(stress_rows)
    @test [row.direction for row in stress_summary] == [:reactive_removed, :uniform]
    @test all(row.status == :pass for row in stress_summary)
    @test all(row.boundary_status == :not_observed for row in stress_summary)
    @test all(isfinite(row.minimum_condition_margin) for row in stress_summary)
    failed_rows = [row.direction == :uniform && row.lambda == 1.0 ?
        merge(row, (status=:fail, endpoint_status=:fail,
                    message="declared voltage-limit boundary")) : row
        for row in stress_rows]
    failed_summary = operability_stress_summary(failed_rows)
    uniform_failed = only(filter(row -> row.direction == :uniform, failed_summary))
    @test uniform_failed.status == :fail
    @test uniform_failed.boundary_status == :fail
    @test uniform_failed.boundary_lambda == 1.0
    ensemble = operability_stress_ensemble_rows(
        Dict("base" => stress_rows, "failed" => failed_rows))
    @test length(ensemble) == 4
    @test ensemble[1].model == "base"
    @test any(row.model == "failed" && row.status == :fail for row in ensemble)
    @test isempty(operability_stress_summary(NamedTuple[]))
    @test_throws ArgumentError operability_stress_rows(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()),
        directions=[OperabilityStressDirection(:duplicate),
            OperabilityStressDirection(:duplicate)], solve=n -> solve_pf(n; per_unit=false))
    @test_throws ArgumentError OperabilitySpec(
        scaling_policy=SIUnitsScaling(), jacobian_spectrum=:unknown)
    @test_throws ArgumentError OperabilitySpec(
        scaling_policy=SIUnitsScaling(), jacobian_storage=:sparse)
    critical = report.branch_evidence["critical_mode"]
    @test critical["status"] == :pass
    @test length(critical["left_vector"]) == 2 * length(report.state_nodes)
    @test length(critical["right_node_participation"]) == length(report.state_nodes)
    @test_throws ArgumentError OperabilitySpec(
        scaling_policy=SIUnitsScaling(), closure=:operational)

    validation_report = check_opf_operability(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_sensitivity_validation=true))
    @test validation_report.checks["load_scale_sensitivity_validation"].status == :pass
    validation = validation_report.sensitivities["validation"]["load_scale"]
    @test validation["absolute_error"] < validation["tolerance"]
    @test validation["plus_residual"] < 1e-6
    @test validation["minus_residual"] < 1e-6
    @test validation_report.checks["directional_sensitivity_validation"].status == :pass
    direction_validation = validation_report.sensitivities["validation"]["directions"]
    @test all(v["status"] == :pass for family in values(direction_validation)
              for v in values(family))

    bad = deepcopy(pf)
    bad["bus"]["bus1"]["1"]["vr"] += 10.0
    bad_report = check_opf_operability(net, bad; spec)
    @test bad_report.status == :fail
    @test bad_report.checks["endpoint"].status == :fail

    # The balanced fixture uses the phase ordering/convention shared with the
    # existing inverter symmetrical-component helper.
    net3 = inv_grid3_bal()
    pf3 = solve_pf(net3; per_unit=false)
    report3 = check_opf_operability(net3, pf3;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(), vuf_max=0.01))
    @test report3.status == :pass
    @test report3.checks["sequence_unbalance"].status == :pass
    @test report3.sequences["poc"]["vuf"] < 1e-4
    @test abs(report3.sequences["poc"]["v1"]) > 200.0
    sequence_sensitivity = report3.branch_evidence["sequence_sensitivity"]["buses"]["poc"]
    @test isfinite(sequence_sensitivity["positive_sequence_magnitude_derivative"])
    @test isfinite(sequence_sensitivity["vuf_derivative"])

    delta_net = inv_grid3_bal()
    delta_net["load"] = Dict("dΔ" => Dict{String,Any}(
        "bus" => "poc", "terminal_map" => ["a", "b", "c"],
        "configuration" => "DELTA", "p_nom" => [100.0, 100.0, 100.0],
        "q_nom" => [0.0, 0.0, 0.0]))
    delta_pf = solve_pf(delta_net; per_unit=false)
    delta_report = check_opf_operability(delta_net, delta_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()))
    @test delta_report.status == :not_applicable
    @test length(delta_report.load_connections) == 3
    @test all(r["positive"] !== r["negative"] for r in values(delta_report.load_connections))
    @test delta_report.provenance["operability"]["model_inventory"]["load_configurations"] ==
          ["DELTA"]
    delta_certificate = check_opf_operability(delta_net, delta_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_fixed_point_certificate=true))
    @test delta_certificate.checks["fixed_point_certificate"].status == :pass
    @test delta_certificate.branch_evidence["fixed_point_certificate"]["edge_count"] == 3
    @test length(delta_certificate.branch_evidence["fixed_point_certificate"]["region_radii"]) ==
          length(delta_certificate.state_nodes)
    unbalanced_delta_net = deepcopy(delta_net)
    unbalanced_delta_net["load"]["dΔ"]["p_nom"] = [100.0, 200.0, 300.0]
    unbalanced_delta_pf = solve_pf(unbalanced_delta_net; per_unit=false)
    unbalanced_delta_certificate = check_opf_operability(unbalanced_delta_net,
        unbalanced_delta_pf; spec=OperabilitySpec(
            scaling_policy=SIUnitsScaling(), compute_fixed_point_certificate=true))
    unbalanced_region = unbalanced_delta_certificate.branch_evidence[
        "fixed_point_certificate"]["region_radii"]
    @test unbalanced_delta_certificate.checks["fixed_point_certificate"].status == :pass
    @test length(unbalanced_region) == 3
    @test maximum(unbalanced_region) > minimum(unbalanced_region)

    # A floating-neutral WYE connection has one nonlinear edge mapped into
    # four free node states; exercise the Euclidean incidence bound without
    # assuming one edge per state.
    floating_net = inv_grid3_floating_neutral()
    floating_pf = solve_pf(floating_net; per_unit=false)
    floating_certificate = check_opf_operability(floating_net, floating_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_fixed_point_certificate=true))
    @test floating_certificate.checks["fixed_point_certificate"].status == :pass
    # The floating-neutral incidence has ||a||₂=√2; the assembly-level probe
    # below guards the corrected Euclidean bound for this multi-state edge.
    @test floating_certificate.checks["fixed_point_euclidean_region"].status == :pass
    @test length(floating_certificate.state_nodes) == 4

    four_wire_net = operability_four_wire_radial()
    four_wire_pf = solve_pf(four_wire_net; per_unit=false)
    four_wire_report = check_opf_operability(four_wire_net, four_wire_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            voltage_min=180.0, voltage_max=280.0,
            compute_fixed_point_certificate=true))
    @test four_wire_report.status == :pass
    @test four_wire_report.checks["fixed_point_certificate"].status == :pass
    @test length(four_wire_report.state_nodes) == 8
    @test four_wire_report.branch_evidence["fixed_point_certificate"]["edge_count"] == 3
    four_wire_extremes = check_opf_operability(four_wire_net, four_wire_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            voltage_min=180.0, voltage_max=280.0,
            compute_sensitivity=false, jacobian_spectrum=:extremes))
    @test four_wire_extremes.checks["jacobian_regular"].status == :pass
    @test four_wire_extremes.singular_values ≈
          [first(four_wire_report.singular_values), last(four_wire_report.singular_values)] rtol=1e-5
    @test four_wire_extremes.branch_evidence["complexity"]["jacobian_spectrum_effective_mode"] ==
          "extremes"

    # Assembly-level regression: on a lightly loaded copy of the multi-bus
    # four-wire fixture, probe several phase-rotated incidence directions near
    # the lower connection-voltage boundary and ensure the reported q bounds
    # the fixed-point map's real Jacobian. The non-zero source angles matter.
    four_wire_light_net = operability_stress_network(four_wire_net, 0.05,
        OperabilityStressDirection(:light))
    four_wire_light_pf = solve_pf(four_wire_light_net; per_unit=false)
    four_wire_light = check_opf_operability(four_wire_light_net, four_wire_light_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_fixed_point_certificate=true))
    euclidean = four_wire_light.branch_evidence["fixed_point_certificate"]["euclidean_region"]
    @test euclidean["status"] == :pass
    lin_light = BMOPFTools.ybus_linearized(four_wire_light_net; fold=:constant_z)
    vmap_light = PowerOptLab._operability_solution_map(four_wire_light_pf, lin_light.nodes)
    meta_light = PowerOptLab._operability_state_meta(lin_light, vmap_light,
        PowerOptLab._operability_source_buses(four_wire_light_net))
    nf_light = length(meta_light.free)
    yll_light = ComplexF64.(lin_light.Y[meta_light.free, meta_light.free])
    yls_light = ComplexF64.(lin_light.Y[meta_light.free, meta_light.fixed])
    yll_factor_light = lu(yll_light)
    w_light = -(yll_factor_light \ (yls_light * meta_light.fixed_values))
    map_light(vfree) = begin
        vfull = zeros(ComplexF64, length(lin_light.nodes))
        vfull[meta_light.fixed] = meta_light.fixed_values
        vfull[meta_light.free] = vfree
        w_light + yll_factor_light \ ComplexF64.(lin_light.i_comp(vfull)[meta_light.free])
    end
    sl_light = only(PowerOptLab._load_subloads(four_wire_light_net["load"]["la"], four_wire_light_net))
    a_light = zeros(Float64, nf_light)
    free_index_light = Dict(gi => li for (li, gi) in enumerate(meta_light.free))
    pi_light = lin_light.index[sl_light.pos]
    ni_light = lin_light.index[sl_light.neg]
    haskey(free_index_light, pi_light) && (a_light[free_index_light[pi_light]] += 1.0)
    haskey(free_index_light, ni_light) && (a_light[free_index_light[ni_light]] -= 1.0)
    dv0_light = (haskey(free_index_light, pi_light) ? w_light[free_index_light[pi_light]] : 0.0im) -
        (haskey(free_index_light, ni_light) ? w_light[free_index_light[ni_light]] : 0.0im)
    phase_light = dv0_light / abs(dv0_light)
    pack_light(v) = vcat(real.(v), imag.(v))
    unpack_light(z) = ComplexF64.(z[1:nf_light] .+ im .* z[nf_light+1:end])
    probe_factors = (phase_light, -phase_light, im * phase_light, -im * phase_light)
    probe_norms = Float64[]
    for factor in probe_factors
        probe = w_light .+ (0.999 * euclidean["radius"]) .*
            ComplexF64.(a_light ./ norm(a_light)) .* factor
        map_jacobian = finite_difference_jacobian(
            z -> pack_light(map_light(unpack_light(z))), pack_light(probe); step=1e-6)
        push!(probe_norms, opnorm(map_jacobian, 2))
    end
    @test maximum(probe_norms) <= euclidean["contraction_factor"] * (1.0 + 1e-3)

    helm_report = check_opf_operability(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(), compute_helm=true))
    @test helm_report.checks["helm_reachability"].status == :pass
    @test helm_report.branch_evidence["reachability"]["base"] == "energized_no_load_germ"
    @test helm_report.branch_evidence["reachability"]["endpoint_mismatch"] < 1e-4

    current_net = single_bus_net(pload=100.0)
    current_net["load"]["ld1"]["model"] = "constant_current"
    current_net["load"]["ld1"]["v_nom"] = [1000.0]
    current_pf = solve_pf(current_net; per_unit=false)
    current_report = check_opf_operability(current_net, current_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(), compute_helm=true))
    @test current_report.checks["helm_reachability"].status == :not_applicable
    @test any(occursin("constant_current", reason) for reason in
        current_report.branch_evidence["reachability"]["reasons"])
    current_certificate = check_opf_operability(current_net, current_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_fixed_point_certificate=true))
    @test current_certificate.checks["fixed_point_certificate"].status == :pass
    @test current_certificate.branch_evidence["fixed_point_certificate"]["contraction_factor"] < 1.0
    @test current_certificate.branch_evidence["fixed_point_certificate"][
        "law_bound_validation"]["status"] == :pass

    for model in ("zip", "exponential")
        model_net = single_bus_net(pload=100.0)
        load = model_net["load"]["ld1"]
        load["model"] = model
        load["v_nom"] = [1000.0]
        if model == "zip"
            load["alpha_z"] = [0.0]; load["alpha_i"] = [0.0]
            load["alpha_p"] = [1.0]; load["beta_z"] = [0.0]
            load["beta_i"] = [0.0]; load["beta_p"] = [1.0]
        else
            load["gamma_p"] = [1.5]; load["gamma_q"] = [1.5]
        end
        model_pf = solve_pf(model_net; per_unit=false)
        model_certificate = check_opf_operability(model_net, model_pf;
            spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
                compute_fixed_point_certificate=true))
        @test model_certificate.checks["fixed_point_certificate"].status == :pass
        @test model_certificate.branch_evidence["fixed_point_certificate"]["edge_count"] == 1
        @test model_certificate.branch_evidence["fixed_point_certificate"][
            "law_bound_validation"]["status"] == :pass
    end

    for model in ("constant_current", "constant_impedance")
        model_net = single_bus_net(pload=100.0)
        model_net["load"]["ld1"]["model"] = model
        model_net["load"]["ld1"]["v_nom"] = [1000.0]
        model_pf = solve_pf(model_net; per_unit=false)
        model_report = check_opf_operability(model_net, model_pf;
            spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
                compute_sensitivity_validation=true))
        @test model_report.status == :not_applicable
        @test model_report.provenance["operability"]["model_inventory"]["load_models"] == [model]
        @test model_report.checks["load_scale_sensitivity_validation"].status == :pass
        @test model_report.checks["directional_sensitivity_validation"].status == :pass
        model_trace = continue_opf_operability(model_net, model_pf;
            spec=OperabilitySpec(scaling_policy=SIUnitsScaling()),
            continuation=OperabilityContinuationSpec(initial_step=0.2))
        @test model_trace.status == :pass
        @test last(model_trace.lambdas) == 1.0
    end

    q_only_net = single_bus_net(pload=0.0)
    q_only_net["load"]["ld1"]["p_nom"] = Float64[]
    q_only_net["load"]["ld1"]["q_nom"] = [100.0]
    q_only_pf = solve_pf(q_only_net; per_unit=false)
    q_only_report = check_opf_operability(q_only_net, q_only_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            voltage_min=800.0, voltage_max=1100.0))
    @test haskey(q_only_report.sensitivities["directions"], "Q")
    @test haskey(q_only_report.sensitivities["directions"]["Q"], "ld1/1")

    trace = continue_opf_operability(net, pf;
        spec=spec,
        continuation=OperabilityContinuationSpec(initial_step=0.2, max_step=0.25))
    @test trace.status == :pass
    @test first(trace.lambdas) == 0.0
    @test last(trace.lambdas) == 1.0
    @test trace.endpoint_match === true
    @test trace.endpoint_distance < 1e-4
    @test maximum(trace.residuals) < 1e-5
    @test trace.provenance["continuation"]["pseudo_arclength"] === false
    @test_throws ArgumentError OperabilityContinuationSpec(initial_step=0.01, min_step=0.1)

    pseudo_trace = continue_opf_operability_pseudo_arclength(net, pf;
        spec=spec,
        continuation=OperabilityPseudoArclengthSpec(
            initial_step=0.05, max_step=0.1, target_lambda_tol=0.01))
    @test pseudo_trace.status == :pass
    @test first(pseudo_trace.lambdas) == 0.0
    @test abs(last(pseudo_trace.lambdas) - 1.0) <= 0.01
    @test pseudo_trace.endpoint_match === true
    @test pseudo_trace.endpoint_distance < 1e-4
    @test maximum(pseudo_trace.residuals) < 1e-5
    @test pseudo_trace.provenance["continuation"]["pseudo_arclength"] === true
    @test all(>(0.0), pseudo_trace.provenance["continuation"]["arclength_state_scale"])
    @test !isempty(pseudo_trace.provenance["continuation"]["curvature_history"])
    @test all(isfinite, pseudo_trace.provenance["continuation"]["curvature_history"])
    @test length(pseudo_trace.provenance["continuation"]["curvature_history"]) ==
          length(pseudo_trace.provenance["continuation"]["arclength_steps"])
    @test pseudo_trace.provenance["continuation"]["curvature_control"] ==
          Dict("low" => 2.0, "high" => 5.0)
    pseudo_rows = operability_continuation_rows(pseudo_trace)
    @test length(pseudo_rows) == length(pseudo_trace.lambdas)
    @test pseudo_rows[1].index == 1
    @test pseudo_rows[1].lambda == 0.0
    @test isnan(pseudo_rows[1].curvature)
    @test [row.lambda for row in pseudo_rows] == pseudo_trace.lambdas
    @test pseudo_trace.provenance["continuation"]["margin"]["status"] == :not_observed
    @test_throws ArgumentError OperabilityPseudoArclengthSpec(initial_step=0.01,
                                                               min_step=0.1)
    @test_throws ArgumentError OperabilityPseudoArclengthSpec(
        curvature_low=2.0, curvature_high=1.0)
    @test_throws ArgumentError OperabilityPseudoArclengthSpec(lambda_min=1.0)

    # Analytic two-bus resistive feeder: P = V(E - V) / R has a high branch
    # at 600 V and a low branch at 400 V for P = 2.4 MW.  Both are valid
    # equilibria, but only the high branch is connected to the energized germ
    # before the λ=1 target crossing.
    nose_net = single_bus_net(pload=2.4e6)
    nose_pf = solve_pf(nose_net; per_unit=false)
    @test nose_pf["bus"]["bus1"]["1"]["vr"] ≈ 600.0 atol=1e-6
    nose_trace = continue_opf_operability_pseudo_arclength(nose_net, nose_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()),
        continuation=OperabilityPseudoArclengthSpec(
            initial_step=0.05, max_step=0.1, max_steps=40, target_lambda_tol=0.01))
    @test nose_trace.status == :pass
    @test nose_trace.endpoint_match === true
    @test any(get(event, "kind", "") == "target_crossing" for event in nose_trace.events)
    @test any(get(event, "kind", "") == "target_refinement" &&
              get(event, "status", nothing) == :pass for event in nose_trace.events)

    low_solution = deepcopy(nose_pf)
    low_solution["bus"]["bus1"]["1"]["vr"] = 400.0
    low_report = check_opf_operability(nose_net, low_solution;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()))
    @test low_report.status == :not_applicable
    @test low_report.branch_evidence["dP_dV"]["connections"]["ld1/1"]["classification"] ==
          "positive_low_side_indicator"
    low_certificate_report = check_opf_operability(nose_net, low_solution;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(),
            compute_fixed_point_certificate=true))
    @test low_certificate_report.checks["fixed_point_certificate"].status == :inconclusive
    @test low_certificate_report.checks["fixed_point_local_region"].status == :inconclusive
    @test low_certificate_report.checks["fixed_point_euclidean_region"].status == :inconclusive
    @test low_certificate_report.branch_evidence["fixed_point_certificate"][
        "candidate_inside_region"] === false
    low_trace = continue_opf_operability_pseudo_arclength(nose_net, low_solution;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()),
        continuation=OperabilityPseudoArclengthSpec(
            initial_step=0.05, max_step=0.1, max_steps=40, target_lambda_tol=0.01))
    @test low_trace.status == :fail
    @test low_trace.endpoint_match === false
    @test low_trace.endpoint_distance ≈ 200.0 atol=1e-6
    @test any(get(event, "kind", "") == "target_refinement" &&
              get(event, "status", nothing) == :pass for event in low_trace.events)

    stress_trace = continue_opf_operability_pseudo_arclength(nose_net, nose_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()),
        continuation=OperabilityPseudoArclengthSpec(
            initial_step=0.05, max_step=0.1, max_steps=40),
        stop_at_target=false)
    @test stress_trace.status == :pass
    @test stress_trace.provenance["continuation"]["stop_at_target"] === false
    @test stress_trace.provenance["continuation"]["lambda_min"] == 0.0
    @test count(event -> get(event, "kind", "") == "target_crossing",
                stress_trace.events) == 2
    @test count(event -> get(event, "kind", "") == "target_refinement" &&
                         get(event, "status", nothing) == :pass,
                stress_trace.events) == 2
    @test any(get(event, "kind", "") == "model_domain_failure" &&
              get(event, "status", nothing) == :pass for event in stress_trace.events)
    @test last(stress_trace.lambdas) == 0.0
    fold_events = filter(event -> get(event, "kind", "") == "fold_candidate",
                         stress_trace.events)
    @test !isempty(fold_events)
    localized_fold = fold_events[1]["fold_localization"]
    @test localized_fold["status"] == :pass
    @test localized_fold["lambda"] ≈ 25 / 24 atol=1e-6
    @test localized_fold["sigma_min"] < 1e-7
    @test stress_trace.provenance["continuation"]["margin"]["mechanism"] == :fold_candidate
    @test stress_trace.provenance["continuation"]["margin"]["lambda"] ≈ 25 / 24 atol=1e-6
    @test stress_trace.provenance["continuation"]["margin"]["parameter_margin"] ≈ 1 / 24 atol=1e-6
    @test operability_continuation_margin(stress_trace)["lambda"] ≈ 25 / 24 atol=1e-6
    @test solve_diagnostics(stress_trace).margin["mechanism"] == :fold_candidate
    stress_rows = operability_continuation_rows(stress_trace)
    @test any(:fold_candidate in row.event_kinds for row in stress_rows)
    @test count(:target_crossing in row.event_kinds for row in stress_rows) == 2
    @test_throws ArgumentError operability_continuation_margin(stress_trace;
                                                               reference_lambda=0.0)

    limit_trace = continue_opf_operability_pseudo_arclength(nose_net, nose_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(), voltage_min=800.0),
        continuation=OperabilityPseudoArclengthSpec(
            initial_step=0.05, max_step=0.1, max_steps=40),
        stop_on_voltage_limit=true)
    @test limit_trace.status == :fail
    @test limit_trace.provenance["continuation"]["stop_on_voltage_limit"] === true
    @test any(get(event, "kind", "") == "voltage_limit" for event in limit_trace.events)
    @test last(limit_trace.lambdas) < 1.0
    @test limit_trace.provenance["continuation"]["margin"]["mechanism"] == :voltage_limit
    @test limit_trace.provenance["continuation"]["margin"]["pre_reference"] === true

    fold_guess = deepcopy(nose_pf)
    fold_guess["bus"]["bus1"]["1"]["vr"] = 500.0
    fold = locate_opf_operability_fold(nose_net, fold_guess;
        lambda=1.04, spec=OperabilitySpec(scaling_policy=SIUnitsScaling()))
    @test fold.status == :pass
    @test fold.lambda ≈ 25 / 24 atol=1e-7
    @test fold.node_voltages[("bus1", "1")] ≈ 500.0 + 0.0im atol=1e-5
    @test fold.residual_norm < 1e-8
    @test fold.sigma_min < 1e-8
    @test fold.critical_mode["coordinate_order"] == "[real(state_nodes); imag(state_nodes)]"
    @test fold.provenance["fold_localization"]["equations"] == "F=0,Jv=0,norm(v)=1"

    # Scope exclusions are explicit evidence rather than a silent partial pass.
    ibr_net = doe_ibr_feeder()
    ibr_pf = solve_pf(ibr_net; per_unit=false)
    ibr_scope_audit = operability_scope_audit(ibr_net)
    @test ibr_scope_audit["status"] == :not_applicable
    @test ibr_scope_audit["ibr_count"] == 1
    @test ibr_scope_audit["control_closure"] == "outside_native_static_seam"
    @test any(occursin("IBR", reason) for reason in ibr_scope_audit["unsupported_reasons"])
    ibr_report = check_opf_operability(ibr_net, ibr_pf; spec)
    @test ibr_report.status == :not_applicable
    @test !isempty(ibr_report.unsupported)
    @test ibr_report.checks["scope"].status == :not_applicable
    @test ibr_report.provenance["operability"]["status"] == :not_applicable
    @test any(occursin("IBR", reason) for reason in ibr_report.unsupported)
    ibr_trace = continue_opf_operability(ibr_net, ibr_pf; spec=spec)
    @test ibr_trace.status == :not_applicable
    @test ibr_trace.events[1]["kind"] == "unsupported_physics"
    ibr_pseudo_trace = continue_opf_operability_pseudo_arclength(
        ibr_net, ibr_pf; spec=spec)
    @test ibr_pseudo_trace.status == :not_applicable
    @test ibr_pseudo_trace.events[1]["kind"] == "unsupported_physics"
    ibr_fold = locate_opf_operability_fold(ibr_net, ibr_pf;
        spec=spec, lambda=1.0)
    @test ibr_fold.status == :not_applicable
    @test_throws ArgumentError locate_opf_operability_fold(
        net, pf; spec=spec, lambda=0.0)

    unsupported_config_net = deepcopy(net)
    unsupported_config_net["load"]["ld1"]["configuration"] = "BOGUS"
    unsupported_config_report = check_opf_operability(unsupported_config_net, pf; spec)
    @test unsupported_config_report.status == :not_applicable
    @test unsupported_config_report.checks["scope"].status == :not_applicable
    @test any(occursin("unsupported configuration", reason)
              for reason in unsupported_config_report.unsupported)
    @test unsupported_config_report.provenance["operability"]["model_inventory"][
        "load_configurations"] == ["BOGUS"]
end
