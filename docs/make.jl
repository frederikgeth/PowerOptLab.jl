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
        "Research program" => [
            "Overview" => "research_program.md",
            "DOE quantification review" => "problems/doe_quantification_review.md",
            "Post-OPF operability roadmap" => "problems/post_opf_operability_roadmap.md",
        ],
        "Inverter-based resources" => [
            "Overview"                 => "ibr/index.md",
            "Scientific foundations"   => "ibr/foundations.md",
            "Phase-aware controls"      => "ibr/phase_aware_control_laws.md",
            "Control-study methodology" => "ibr/control_study_methodology.md",
            "Network control studies"   => "ibr/network_control_studies.md",
            "Closed-loop evidence"      => "ibr/closed_loop_evidence.md",
            "Phase-aware control API"   => "components/inverter_controls.md",
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
        "State estimation" => [
            "Overview"                     => "estimation/index.md",
            "Choosing a formulation"       => "estimation/comparison.md",
            "Likelihood, loss, and priors" => "estimation/maximum_likelihood.md",
            "Legacy WLS estimator"         => "problems/state_estimation.md",
            "Constrained NLLS estimator"   => "problems/constrained_state_estimation.md",
            "Modelling tutorial"           => "tutorials/constrained_nlls_state_estimation.md",
            "Observability & under-observed" => "estimation/observability.md",
            "Current-magnitude measurements" => "estimation/current_magnitude.md",
            "Background & roadmap"         => "estimation/state_of_the_art.md",
        ],
        "Problem specifications" => [
            "Multi-period OPF"     => "problems/multiperiod.md",
            "Post-OPF operability" => "problems/operability.md",
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
            "Single-snapshot post-OPF operability study" => "tutorials/post_opf_operability_study.md",
            "Dynamic operating envelopes" => "tutorials/dynamic_operating_envelopes.md",
            "Inverse Carson reconstruction" => "tutorials/inverse_carson_reconstruction.md",
            "Battery storage models" => "tutorials/battery_storage_models.md",
            "HELM versus nonlinear power flow" => "tutorials/helm_vs_nonlinear_power_flow.md",
            "Learning smart-inverter controls" => "tutorials/learning_smart_inverter_controls.md",
        ],
        "Contributing"       => "contributing.md",
        "API reference"      => "api.md",
        "Symbol index"       => "symbol_index.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/frederikgeth/PowerOptLab.jl.git",
    devbranch = "main",
    push_preview = false,
)
