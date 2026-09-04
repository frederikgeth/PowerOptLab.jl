#!/usr/bin/env julia

using PowerOptLab

function _tsv_cell(value)
    text = value isa AbstractString ? String(value) : repr(value)
    return replace(text, '\t' => ' ', '\n' => ' ', '\r' => ' ')
end

function _write_tsv(path, rows)
    isempty(rows) && error("benchmark produced no rows")
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(String.(columns), '\t'))
        for row in rows
            println(io, join((_tsv_cell(getproperty(row, column))
                              for column in columns), '\t'))
        end
    end
end

function _usage()
    println(stderr,
        "usage: julia --project=. scripts/run_doe_benchmark.jl CASE.jl OUTPUT.tsv")
    exit(2)
end

length(ARGS) == 2 || _usage()
case_path = abspath(ARGS[1])
output_path = abspath(ARGS[2])
isfile(case_path) || error("case file not found: $case_path")

include(case_path)
isdefined(Main, :doe_benchmark_case) || error(
    "case file must define doe_benchmark_case()")
case = Main.doe_benchmark_case()
case isa NamedTuple || error("doe_benchmark_case() must return a NamedTuple")

nets = case.nets
connection_points = case.connection_points
methods = get(case, :methods, ["default" => NamedTuple()])
metadata = get(case, :metadata, Dict{String,Any}())
seeds = get(case, :seeds, Dict{String,Int}())

rows = NamedTuple[]
for entry in methods
    entry isa Pair || error("methods must contain label => NamedTuple pairs")
    label = string(first(entry))
    keywords = last(entry)
    keywords isa NamedTuple || error("method '$label' options must be a NamedTuple")
    control_policy = get(keywords, :control_policy, PerfectRecourse())
    fairness = get(keywords, :fairness, :equal)
    direction = get(keywords, :direction, :export)
    security = get(keywords, :security, :bound_point)
    utilizations = get(keywords, :utilizations, nothing)
    solver_options = get(keywords, :solver_options, NamedTuple())
    spec = DOEStudySpec(nets, connection_points;
        direction, security, utilizations, control_policy, fairness,
        solver_options, seeds,
        metadata=merge(Dict{String,Any}(string(key) => value
                       for (key, value) in pairs(metadata)),
                       Dict("method" => label)))
    result = solve_operating_envelope(nets, connection_points; keywords...)
    append!(rows, doe_benchmark_rows(spec, result; method_label=label))
end

mkpath(dirname(output_path))
_write_tsv(output_path, rows)
println("wrote $(length(rows)) DOE benchmark rows to $output_path")
