# Documentation tutorial audit

Date: 2026-10-11. Branch: `codex/reorganize-component-docs`.
This is a developer record, outside the published documentation.

## Coverage of the 23 pages in docs/src/tutorials

| Coverage | Pages | Check |
|:--|:--|:--|
| Study runner | battery_storage_models; constrained_nlls_state_estimation; ibr_topology_under_unbalance; ibr_dc_source_and_split_link; ibr_ac_harmonics_lcl; helm_vs_nonlinear_power_flow; post_opf_operability_study | Exact published Julia fences, fresh module per page, explicit expected results and scientific assertions; CI floor and stable |
| DOE runner | doe_analytic_reference; doe_faithful_replay; doe_getting_started; doe_range_guarantees; doe_scenario_design; doe_statistical_validation; doe_uncertainty_coverage; doe_uncertainty_models | Existing exact-fence DOE runner, CI floor and stable |
| Documenter | generalized_generator_models; generalized_generator_tradeoffs; generalized_generator_validation | Published @example blocks execute the analytic, trade-off and validation companion scripts in documentation CI |
| Schematic guides | advanced_inverter_modelling; inverse_carson_reconstruction; learning_smart_inverter_controls; dynamic_operating_envelopes; doe_research_workflow | Explicitly labelled model-selection or research workflows with caller-supplied inputs; links to complete worked studies |

The five schematic guides are deliberately not advertised as sequential scripts.
Adding an executable marker to them requires complete study inputs and scientific
assertions; syntax alone would not be an adequate check. The inverse Carson
reference contains a complete single-candidate fit, but the tutorial's catalog
comparison and profile fragments still require study-specific choices.

## Corrections in this pass

- Added three scientific SVGs from the battery and state-estimation studies.
  Julia exports quantities from the executed modules. Matplotlib renders the
  figures, and stable CI checks fresh data against the committed snapshot,
  renderer/input fingerprints and a rendering smoke run. Physical tolerances
  permit small cross-platform solver differences; no pixel equality is claimed.
- The inverter studies already ran, but CI had not executed their plain fences.
  Added coverage and checks for neutral/zero-sequence identities, sequence limits,
  DC current sharing, reserve coverage and AC/grid total RMS ratings. The
  DC-source comparison now holds AC P/Q fixed with explicit power scaling and a
  declared quadratic converter loss: maximum-export cases had both a Julia 1.10
  convergence failure and changing operating points that confounded source
  current-sharing comparisons. Strict status checks remain in place.
- Completed HELM and operability setup. Added analytic load voltage/neutral rise,
  explicit nominal voltage for voltage-dependent loads, finite stress assertions
  and an intentionally unsupported network case. Scope diagnostics are conditional
  on rejection, not a key present in every successful report.
- Corrected the native IBR descriptions: Volt-VAr is a reactive-power equality;
  Volt-Watt is an active-power upper bound, allowing further curtailment. The
  pinned solver supports fractional droop power units; absolute VAR/W variants
  warn and skip. Volt-VAr and Volt-Watt can coexist; fixed PF takes precedence.
  Evidence: pinned BMOPFTools BMOPFOpfExt/ibr.jl, curve resolvers and
  `_ibr_pq_bounds!` stamping.
- Restricted the HELM fallback recommendation: the static Ybus fixed-point map
  supports additional native load laws, not generator/IBR-control equations.
- Distinguished observability rank from conditioning and removed misleading
  meter/forecast and topology-error language in the estimation overview.
- Brought the complete battery and estimation studies ahead of detailed model
  references in their groups. Kept existing URLs and single page registration.

## Limits of this audit

This pass checked tutorial execution coverage, explicit setup, selected capability
claims against implementation, terminology and navigation. It is not a complete
line-by-line review of all 108 published pages or a new physical-validation study.
Existing DOE and generator validation remain evidence under their declared scope.
No roadmap or other workflow record was removed in this pass.

Useful subsequent work is a complete inverse-line reconstruction comparison with
measured-data provenance, and a worked identification/replay example for a fitted
local inverter law. Those need explicit data and inference/evaluation targets;
turning schematic fragments into arbitrary numerical demonstrations would not
supply those scientific assumptions.
