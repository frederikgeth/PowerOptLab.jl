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
        "Home"            => "index.md",
        "Concepts"        => "concepts.md",
        "Component models" => [
            "Storage & EVs"     => "components/devices.md",
            "Advanced inverter" => "components/advanced_inverter.md",
            "IVQ battery"       => "components/ivq_battery.md",
        ],
        "Problem specifications" => [
            "Multi-period OPF"     => "problems/multiperiod.md",
            "Legacy WLS state estimation" => "problems/state_estimation.md",
            "Constrained NLLS state estimation" => "problems/constrained_state_estimation.md",
            "Parameter estimation" => "problems/parameter_estimation.md",
            "Inverse Carson"       => "problems/inverse_carson.md",
            "Operating envelopes"  => "problems/operating_envelope.md",
        ],
        "Bespoke algorithms" => [
            "Overview"        => "algorithms/index.md",
            "HELM power flow" => "algorithms/helm.md",
        ],
        "Tutorials" => [
            "Dynamic operating envelopes" => "tutorials/dynamic_operating_envelopes.md",
            "Constrained NLLS state estimation" => "tutorials/constrained_nlls_state_estimation.md",
        ],
        "Contributing"       => "contributing.md",
        "API reference"      => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/frederikgeth/PowerOptLab.jl.git",
    devbranch = "main",
    push_preview = false,
)
