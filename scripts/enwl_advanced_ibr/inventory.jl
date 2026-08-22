#!/usr/bin/env julia

include(joinpath(@__DIR__, "common.jl"))
using .ENWLAdvancedIBR

function usage(io=stdout)
    println(io, """
Usage:
  julia --project=. scripts/enwl_advanced_ibr/inventory.jl [options]

Options:
  --root=PATH       ENWLsnapshots directory (defaults to sibling BMOPFDraftData)
  --feeder=NAME     Restrict to one feeder, e.g. 30bus_LN
  --output=PATH     Also write one-row-per-snapshot CSV
  --help            Show this help
""")
end

function main(args=ARGS)
    options = parse_options(args; defaults=Dict("root" => default_data_root()))
    get(options, "help", "false") == "true" && (usage(); return)
    allowed = Set(["root", "feeder", "output", "help"])
    unknown = setdiff(Set(keys(options)), allowed)
    isempty(unknown) || throw(ArgumentError(
        "unknown options: $(join(sort!(collect(unknown)), ", "))"))
    paths = snapshot_paths(options["root"];
        feeder=get(options, "feeder", nothing))
    rows = inventory_rows(paths)
    print_inventory(rows)
    if haskey(options, "output")
        write_csv(options["output"], rows)
        println("wrote ", abspath(options["output"]))
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
