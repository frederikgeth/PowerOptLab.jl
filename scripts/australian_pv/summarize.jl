# Build results/summary.csv and regenerate the README results tables from the
# committed result CSVs, so that no number in the README is copied by hand.
using Printf
const ROOT = normpath(joinpath(@__DIR__, "..", "..", "studies", "australian_pv"))
const LABELS = [
    "unity" => "Unity PF",
    "mean_vv" => "Mean VV",
    "mean_vv_vw" => "Mean VV+VW",
    "worst_vv_vw" => "Worst-phase VV+VW",
    "sequence_vv_vw" => "Sequence VV+VW",
    "sequence_droop" => "Sequence + I₂ droop (common scale)",
    "sequence_droop_power" => "Sequence + I₂ droop (power priority)",
    "gfm_virtual_impedance" => "Grid-forming virtual impedance",
    "single_phase_native" => "Single-phase PV, native AS/NZS 4777.2",
]

function read_csv(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    [Dict(String(h) => String(v) for (h, v) in zip(header, split(l, ',')))
     for l in lines[2:end] if !isempty(l)]
end
num(r, k) = parse(Float64, r[k])
select(rows, mode; level=1.0, source=1.03) = filter(r -> r["mode"] == mode &&
    num(r, "availability") == level && num(r, "source_pu") == source, rows)
finite_max(x) = (y = filter(isfinite, x); isempty(y) ? NaN : maximum(y))
finite_sum(x) = (y = filter(isfinite, x); isempty(y) ? NaN : sum(y))

function metrics(customers, devices, mode)
    c = select(customers, mode); d = select(devices, mode)
    isempty(c) && return nothing
    (max_V_PN=maximum(num.(c, "vmax_pn_V")),
     max_VUF_percent=maximum(num.(c, "vuf_percent")),
     max_neutral_V=maximum(num.(c, "neutral_V")),
     delivered_P_kW=sum(num.(d, "p_poc_W"))/1e3,
     absorbed_Q_kvar=-sum(num.(d, "q_poc_var"))/1e3,
     converter_loss_kW=finite_sum(num.(d, "converter_loss_W"))/1e3,
     max_current_A=maximum(num.(d, "maximum_converter_current_A")),
     max_I2_A=finite_max(num.(d, "negative_sequence_current_A")),
     max_dv2_V=finite_max(num.(d, "dc_ripple_voltage_V")))
end

fmt(x; digits=2) = isfinite(x) ? @sprintf("%.*f", digits, x) : "—"

function table(dir)
    customers = read_csv(joinpath(dir, "customers.csv"))
    devices = read_csv(joinpath(dir, "devices.csv"))
    rows = [(mode, label, metrics(customers, devices, mode)) for (mode, label) in LABELS]
    filter(r -> r[3] !== nothing, rows)
end

function markdown(rows)
    io = IOBuffer()
    println(io, "| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Peak converter current (A) | Peak I₂ (A) | Peak 2ω ripple (V) |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for (_, label, m) in rows
        println(io, "| ", label, " | ", fmt(m.max_V_PN), " | ", fmt(m.max_VUF_percent; digits=3),
                " | ", fmt(m.max_neutral_V; digits=3), " | ", fmt(m.delivered_P_kW),
                " | ", fmt(m.absorbed_Q_kvar), " | ", fmt(m.converter_loss_kW; digits=3),
                " | ", fmt(m.max_current_A), " | ", fmt(m.max_I2_A), " | ", fmt(m.max_dv2_V), " |")
    end
    String(take!(io))
end

function write_summary(path, rows)
    open(path, "w") do io
        println(io, "mode,max_V_PN,max_VUF_percent,max_neutral_V,delivered_P_kW,absorbed_Q_kvar,converter_loss_kW,max_current_A,max_I2_A,max_dv2_V")
        for (mode, _, m) in rows
            println(io, join([mode; [string(getproperty(m, k)) for k in propertynames(m)]], ','))
        end
    end
end

function replace_block(text, name, content)
    start = "<!-- $name:start -->"; stop = "<!-- $name:end -->"
    i = findfirst(start, text); j = findfirst(stop, text)
    (i === nothing || j === nothing) && error("README is missing the $name markers")
    text[1:last(i)] * "\n" * content * text[first(j):end]
end

function summarize()
    main = table(joinpath(ROOT, "results"))
    write_summary(joinpath(ROOT, "results", "summary.csv"), main)
    readme = read(joinpath(ROOT, "README.md"), String)
    readme = replace_block(readme, "results", markdown(main))
    sensitivity = IOBuffer()
    for dir in sort(filter(d -> startswith(d, "sensitivity_"),
                           readdir(joinpath(ROOT, "results"))))
        rows = table(joinpath(ROOT, "results", dir))
        write_summary(joinpath(ROOT, "results", dir, "summary.csv"), rows)
        println(sensitivity, "**", replace(dir, "sensitivity_" => ""), "**\n")
        println(sensitivity, markdown(rows))
    end
    readme = replace_block(readme, "sensitivity", String(take!(sensitivity)))
    write(joinpath(ROOT, "README.md"), readme)
    println("Updated summary.csv and README tables")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    summarize()
end
