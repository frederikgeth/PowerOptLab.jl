# Export quantities from the executed tutorial module; do not repeat its solves.
using JSON3

function study_figure_exporter(directory)
    mkpath(directory)
    function export_page(filename, tutorial)
        value(name) = getfield(tutorial, name)
        data = if filename == "battery_storage_models.md"
            ideal, resistive = value(:ivq0), value(:ivqr)
            ns, np, resistance = value(:ns), value(:np), value(:resistance)
            soc_grid = collect(range(0.05, 0.95; length=181))
            # Independent analytic reference derived in the published study.
            ocv = 3.0 .+ 0.6 .* soc_grid
            current_bound = min.(5.0, (ocv .- 3.0) ./ resistance)
            cell_bound = ns*np .* (ocv .- resistance .* current_bound) .* current_bound
            Dict(
                "time_h" => [0.0, 1.0, 2.0],
                "ideal_power_w" => ideal.p_poc,
                "resistive_power_w" => resistive.p_poc,
                "ideal_pe_power_w" => value(:pe0).p_net,
                "calibrated_pe_power_w" => value(:per).p_net,
                "ideal_soc" => ideal.soc, "resistive_soc" => resistive.soc,
                "soc_grid" => soc_grid, "cell_export_bound_w" => cell_bound,
                "export_reference_w" => min.(value(:p_rating), cell_bound),
                "inverter_rating_w" => value(:p_rating),
                "voltage_inverter_crossover_soc" => value(:p_rating)*resistance/(ns*np*3.0*0.6),
                "solved_soc" => [0.1, 0.5, 0.9],
                "linear_export_w" => value(:export_limits)[1:2:end],
                "tabulated_export_w" => value(:export_limits)[2:2:end],
            )
        elseif filename == "constrained_nlls_state_estimation.md"
            function neutral(s, x)
                k = s.free_state_map[("load", "n")]
                [x[k], x[length(s.free_state_map) + k]]
            end
            cov = value(:neutral_covariance)
            meters = value(:clean_readings)
            Dict(
                "true_neutral_v" => [real(value(:truth)[("load", "n")]),
                                      imag(value(:truth)[("load", "n")])],
                "fitted_neutral_v" => neutral(value(:structure), value(:sparse_fit).state),
                "voltage_only_neutral_v" => neutral(value(:weak_structure), value(:weak_fit).state),
                "neutral_component_sigma_v" => sqrt.([cov[1,1], cov[2,2]]),
                "weighted_residuals" => value(:sparse_fit).evaluation.residual,
                "voltage_only_residuals" => value(:weak_fit).evaluation.residual,
                "voltage_only_positions" => [k for (k,m) in enumerate(meters) if
                    m in value(:voltage_only)],
                "meter_labels" => ["Vr1", "Vi1", "P1", "Q1", "Vr2", "Vi2", "P2", "Q2",
                                   "Vr3", "Vi3", "P3", "Q3", "Vrn", "Vin", "|I1|"],
                "unobserved_directions" => value(:weak_obs).unobservable_dimension,
            )
        else
            return
        end
        path = joinpath(directory, replace(filename, ".md" => ".json"))
        open(path, "w") do io
            JSON3.pretty(io, data)
            println(io)
        end
        println("Exported solved figure data: ", path)
    end
end
