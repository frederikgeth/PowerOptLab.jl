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
        "Research program" => "research_program.md",
        "Inverter-based resources" => [
            "Overview"                 => "ibr/index.md",
            "Scientific foundations"   => "ibr/foundations.md",
            "Phase-aware controls"      => "ibr/phase_aware_control_laws.md",
            "Control-study methodology" => "ibr/control_study_methodology.md",
            "Advanced inverter API"    => "components/advanced_inverter.md",
            "Modelling tutorial"       => "tutorials/advanced_inverter_modelling.md",
            "Topology under unbalance" => "tutorials/ibr_topology_under_unbalance.md",
            "DC source & split link"   => "tutorials/ibr_dc_source_and_split_link.md",
            "Carrier harmonics & LCL"  => "tutorials/ibr_ac_harmonics_lcl.md",
            "Verification & benchmarks" => "ibr/verification.md",
            "References"               => "ibr/references.md",
        ],
        "Component models" => [
            "Storage & EVs"     => "components/devices.md",
            "IVQ battery"       => "components/ivq_battery.md",
        ],
        "Problem specifications" => [
            "Multi-period OPF"     => "problems/multiperiod.md",
            "Legacy WLS state estimation" => "problems/state_estimation.md",
            "Constrained NLLS state estimation" => "problems/constrained_state_estimation.md",
            "Parameter estimation" => "problems/parameter_estimation.md",
            "Inverse Carson"       => "problems/inverse_carson.md",
            "Operating envelopes"  => "problems/operating_envelope.md",
            "Bilevel PV/tap POC" => "problems/bilevel.md",
        ],
        "Bespoke algorithms" => [
            "Overview"        => "algorithms/index.md",
            "HELM power flow" => "algorithms/helm.md",
        ],
        "Tutorials" => [
            "Dynamic operating envelopes" => "tutorials/dynamic_operating_envelopes.md",
            "Constrained NLLS state estimation" => "tutorials/constrained_nlls_state_estimation.md",
            "Inverse Carson reconstruction" => "tutorials/inverse_carson_reconstruction.md",
            "Battery storage models" => "tutorials/battery_storage_models.md",
            "HELM versus nonlinear power flow" => "tutorials/helm_vs_nonlinear_power_flow.md",
            "Learning smart-inverter controls" => "tutorials/learning_smart_inverter_controls.md",
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
