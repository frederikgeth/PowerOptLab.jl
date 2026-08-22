#!/usr/bin/env julia

include(joinpath(@__DIR__, "common.jl"))
using .ENWLAdvancedIBR
using PowerOptLab

function usage(io=stdout)
    println(io, """
Usage:
  julia --project=. scripts/enwl_advanced_ibr/run_campaign.jl [options]

Options:
  --root=PATH          ENWLsnapshots directory (defaults to sibling BMOPFDraftData)
  --feeder=NAME        Feeder directory (default: 30bus_LN)
  --snapshots=LIST     Comma-separated full IDs or suffixes (default: t13_1400)
  --variants=LIST      baseline,worst_phase,sequence_droop (default: all three)
  --devices=N          Evenly spaced PVs per case; 0 means every PV (default: 3)
  --output=DIR         Output directory (default: results/enwl_advanced_ibr)
  --s-base=VA          OPF power base (default: 100000)
  --max-iter=N         Ipopt iteration limit (default: 2000)
  --tol=VALUE          Ipopt tolerance (default: 1e-6)
  --dry-run=true       Build and validate the campaign without solving
  --help               Show this help

The ENWL inputs contain only single-phase PVs. Each selected record is treated
as a synthetic balanced three-leg retrofit at the same bus. Its original total
VA rating and available W are preserved; unselected PVs remain native.
""")
end

_list(value) = String.(filter(!isempty, strip.(split(value, ','))))

function main(args=ARGS)
    defaults = Dict(
        "root" => default_data_root(), "feeder" => "30bus_LN",
        "snapshots" => "t13_1400",
        "variants" => "baseline,worst_phase,sequence_droop",
        "devices" => "3", "output" => "results/enwl_advanced_ibr",
        "s-base" => "100000", "max-iter" => "2000", "tol" => "1e-6",
        "dry-run" => "false")
    options = parse_options(args; defaults=defaults)
    get(options, "help", "false") == "true" && (usage(); return)
    allowed = Set([keys(defaults)..., "help"])
    unknown = setdiff(Set(keys(options)), allowed)
    isempty(unknown) || throw(ArgumentError(
        "unknown options: $(join(sort!(collect(unknown)), ", "))"))

    snapshot_filter = options["snapshots"] == "all" ? nothing :
        _list(options["snapshots"])
    paths = snapshot_paths(options["root"]; feeder=options["feeder"],
                            snapshots=snapshot_filter)
    variants = _list(options["variants"])
    device_count = parse(Int, options["devices"])
    cases = build_study_cases(paths; variants=variants,
                              device_count=device_count)
    print_campaign_plan(cases)
    lowercase(options["dry-run"]) in ("true", "1", "yes") && return

    result = run_inverter_control_study(cases;
        s_base=parse(Float64, options["s-base"]),
        continue_on_error=:all,
        solver_options=(
            "max_iter" => parse(Int, options["max-iter"]),
            "tol" => parse(Float64, options["tol"]),
            "bound_relax_factor" => 0.0))

    output = abspath(options["output"])
    case_rows = inverter_control_study_case_rows(result)
    device_rows = inverter_control_study_device_rows(result)
    phase_rows = inverter_control_study_phase_rows(result)
    summary_rows = inverter_control_study_summary_rows(
        result; group_by=["feeder", "controlled_device_count"])
    compare_variants = "baseline" in variants && length(unique(variants)) > 1
    paired_rows = compare_variants ?
        inverter_control_paired_rows(result, "baseline") : NamedTuple[]
    paired_summary_rows = compare_variants ?
        inverter_control_paired_summary_rows(
            result, "baseline";
            group_by=["feeder", "controlled_device_count"]) : Dict[]

    outputs = [
        write_csv(joinpath(output, "cases.csv"), case_rows),
        write_csv(joinpath(output, "devices.csv"), device_rows),
        write_csv(joinpath(output, "phases.csv"), phase_rows),
        write_csv(joinpath(output, "summary.csv"), summary_rows),
        write_csv(joinpath(output, "paired_devices.csv"), paired_rows),
        write_csv(joinpath(output, "paired_summary.csv"), paired_summary_rows),
    ]
    diagnostics = solve_diagnostics(result)
    println("  publishable cases: ", diagnostics.publishable_case_count,
            "/", diagnostics.case_count)
    println("  validation errors: ", diagnostics.validation_error_case_count)
    println("  unexpected errors: ", diagnostics.unexpected_error_case_count)
    println("  elapsed solver time [s]: ",
            round(diagnostics.elapsed_seconds; digits=2))
    for path in filter(!isnothing, outputs)
        println("  wrote ", path)
    end
    solve_status(result).publishable || exit(2)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
