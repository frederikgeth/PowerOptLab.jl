using Documenter
using PowerOptLab

makedocs(
    sitename = "PowerOptLab.jl",
    modules  = [PowerOptLab],
    repo     = Documenter.Remotes.GitHub("frederikgeth", "PowerOptLab.jl"),
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link  = "main",
    ),
    pages = [
        "Home"               => "index.md",
        "Storage & EVs"      => "devices.md",
        "Multi-period OPF"   => "multiperiod.md",
        "State estimation"   => "state_estimation.md",
        "Parameter estimation" => "parameter_estimation.md",
        "Operating envelopes" => "operating_envelope.md",
        "Advanced inverter"  => "advanced_inverter.md",
        "API reference"      => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/frederikgeth/PowerOptLab.jl.git",
    devbranch = "main",
    push_preview = false,
)
