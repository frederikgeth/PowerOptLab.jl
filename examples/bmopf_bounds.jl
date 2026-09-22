using PowerOptLab, BMOPFTools, JSON3

net = parse_bmopf("""
{
  "bus": {
    "grid": {"terminal_names": ["a", "n"], "perfectly_grounded_terminals": ["n"]},
    "customer": {"terminal_names": ["a", "n"], "perfectly_grounded_terminals": ["n"],
                 "vpn_min": [210.0], "vpn_max": [250.0]}
  },
  "voltage_source": {
    "grid": {"bus": "grid", "terminal_map": ["a"],
             "v_magnitude": [240.0], "v_angle": [0.0]}
  },
  "line": {
    "service": {"bus_from": "grid", "bus_to": "customer",
                "terminal_map_from": ["a"], "terminal_map_to": ["a"],
                "R_series_1_1": 0.1}
  },
  "load": {
    "customer": {"bus": "customer", "terminal_map": ["a", "n"],
                 "configuration": "SINGLE_PHASE", "p_nom": [3000.0], "q_nom": [4000.0]}
  }
}
"""; from_string=true)

report = propagate_bounds(net)
key = QuantityKey(:load, "customer", :current, ("coil", "1"))
current = bound(report, key)
println("Stop: ", report.status)
println("Load current magnitude: [", current.lower, ", ", current.upper, "] A")
println("Endpoint explanation: ", length(explain_bound(report, key)), " evidence records")

filled, patch = materialize_bounds(net, report)
println("Filled line current limit: ", filled["line"]["service"]["i_max"], " A")
println("Changed fields: ", length(patch))

# Optional output directory; the example is read-only when no argument is given.
if !isempty(ARGS)
    directory = only(ARGS)
    mkpath(directory)
    write_bmopf(filled, joinpath(directory, "network.bmopf.json"))
    open(joinpath(directory, "bounds.json"), "w") do io
        JSON3.pretty(io, bounds_sidecar(report))
    end
    open(joinpath(directory, "patch.json"), "w") do io
        JSON3.pretty(io, patch)
    end
end
