#!/usr/bin/env julia

# Diagnostic size-scaling benchmark for the single-snapshot checker.
#
# Usage:
#   julia --project=. scripts/benchmark_post_opf_operability_scaling.jl [n ...]
#
# The generated radial feeders are deliberately simple. Results are timing and
# allocation observations for the local machine, not performance guarantees or
# a replacement for the scientific operability claims.

using PowerOptLab
using BMOPFTools: parse_bmopf, solve_pf, SIUnitsScaling

function radial_net(n::Int)
    n >= 1 || throw(ArgumentError("radial feeder size must be >= 1"))
    io = IOBuffer()
    print(io, "{\"bus\":{")
    print(io, "\"source\":{\"terminal_names\":[\"1\",\"n\"],\"perfectly_grounded_terminals\":[\"n\"]}")
    for i in 1:n
        print(io, ",\"bus", i,
            "\":{\"terminal_names\":[\"1\",\"n\"],\"perfectly_grounded_terminals\":[\"n\"],",
            "\"v_min\":[180.0],\"v_max\":[280.0]}")
    end
    print(io, "},\"voltage_source\":{\"vs\":{\"bus\":\"source\",",
        "\"terminal_map\":[\"1\"],\"v_magnitude\":[230.0],\"v_angle\":[0.0]}},")
    print(io, "\"linecode\":{\"lc\":{\"R_series_1_1\":0.05,\"R_series_2_2\":0.05}},\"line\":{")
    for i in 1:n
        i > 1 && print(io, ",")
        from = i == 1 ? "source" : "bus$(i - 1)"
        print(io, "\"l", i, "\":{\"bus_from\":\"", from,
            "\",\"bus_to\":\"bus", i,
            "\",\"terminal_map_from\":[\"1\",\"n\"],",
            "\"terminal_map_to\":[\"1\",\"n\"],\"linecode\":\"lc\",\"length\":1.0}")
    end
    print(io, "},\"load\":{")
    for i in 1:n
        i > 1 && print(io, ",")
        print(io, "\"ld", i, "\":{\"bus\":\"bus", i,
            "\",\"terminal_map\":[\"1\",\"n\"],\"configuration\":\"WYE\",",
            "\"model\":\"constant_power\",\"p_nom\":[100.0],\"q_nom\":[20.0]}")
    end
    print(io, "}}")
    parse_bmopf(String(take!(io)); from_string=true)
end

function benchmark_case(n::Int, spectrum::Symbol)
    net = radial_net(n)
    solution = solve_pf(net; per_unit=false)
    spec = OperabilitySpec(
        scaling_policy=SIUnitsScaling(),
        voltage_min=180.0,
        voltage_max=280.0,
        compute_sensitivity=false,
        compute_fixed_point_certificate=false,
        jacobian_spectrum=spectrum)
    timed = @timed check_opf_operability(net, solution; spec=spec)
    report = timed.value
    complexity = report.branch_evidence["complexity"]
    (nodes=n, spectrum=String(spectrum), status=String(report.status),
     seconds=timed.time, allocated_bytes=timed.bytes,
     real_state=complexity["real_state_dimension"],
     jacobian_kib=complexity["jacobian_storage_bytes_dense"] / 1024)
end

sizes = isempty(ARGS) ? [4, 8, 16, 32, 64] : parse.(Int, ARGS)
all(>(0), sizes) || throw(ArgumentError("benchmark sizes must be positive integers"))
println("n,spectrum,status,seconds,allocated_bytes,real_state,jacobian_kib")
for n in sizes
    for spectrum in (:full, :extremes)
        row = benchmark_case(n, spectrum)
        println(row.nodes, ",", row.spectrum, ",", row.status, ",",
            row.seconds, ",", row.allocated_bytes, ",", row.real_state, ",",
            row.jacobian_kib)
    end
end
