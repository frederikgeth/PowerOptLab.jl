# Run: julia --project=. scripts/benchmark_compiled_state_estimation.jl
# Synthetic sparse chains isolate overhead. These are not operational feeder
# benchmarks: no nonlinear devices, noise, or difficult initialisation.
using PowerOptLab, LinearAlgebra, SparseArrays

function benchmark_se_chain(n)
    net = Dict{String,Any}(
        "bus" => Dict{String,Any}("src" => Dict("terminal_names" => ["1"])),
        "voltage_source" => Dict("s" => Dict("bus" => "src", "terminal_map" => ["1"],
            "v_magnitude" => [230.], "v_angle" => [0.])),
        "linecode" => Dict("lc" => Dict("R_series_1_1" => .1)),
        "line" => Dict{String,Any}())
    ms = Measurement[]
    for i in 1:n
        bus = "b$i"
        net["bus"][bus] = Dict("terminal_names" => ["1"])
        net["line"]["l$i"] = Dict("bus_from" => i == 1 ? "src" : "b$(i-1)",
            "bus_to" => bus, "terminal_map_from" => ["1"], "terminal_map_to" => ["1"],
            "linecode" => "lc", "length" => 1.)
        push!(ms, Measurement(kind=:vr, bus=bus, reference=nothing, value=230., sigma=1.),
                  Measurement(kind=:vi, bus=bus, reference=nothing, value=0., sigma=1.))
    end
    s = compile_state_estimator(net, ms)
    s, SEParameters(s, ms), vcat(fill(230., n), zeros(n))
end

function benchmark_compiled_se(sizes=(50, 200, 800))
    println("Julia $(VERSION); BLAS threads $(BLAS.get_num_threads())")
    println("states,operation,min_seconds,min_allocated_bytes")
    for n in sizes
        s, p, x = benchmark_se_chain(n)
        for (label, f) in (("evaluate", evaluate_state_estimator),
                           ("jacobian", residual_jacobian),
                           ("constraints", constraint_jacobian),
                           ("exact_warm_start", solve_sparse_state_estimator))
            f(s, p, x) # exclude first-call compilation
            runs = [@timed f(s, p, x) for _ in 1:3]
            println("$(2n),$label,$(minimum(r.time for r in runs)),$(minimum(r.bytes for r in runs))")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_compiled_se()
end
