# Dynamic operating envelopes API

The public types and functions for defining, solving, verifying, and recording
dynamic operating-envelope studies. See
[Reproducible DOE recourse, verification, and search](@ref) for
an end-to-end example and guidance on interpreting finite numerical evidence.

## Connections, fairness, and control recourse

```@docs
ConnectionPoint
FairnessPolicy
DOEControlRegistration
DOEControlRule
DOEControlPolicy
PerfectRecourse
IssueFixedControls
IssuePlusLocalLaws
```

## Study specification and evidence records

```@docs
DOEStudySpec
doe_study_manifest
doe_benchmark_rows
doe_context_benchmark_rows
OperatingEnvelopeResult
OperatingEnvelopeVerification
OperatingEnvelopeContextResult
OperatingEnvelopeSearchResult
DOEAdversarialSearchResult
SearchStableOperatingEnvelopeResult
OperatingEnvelopeMultistartResult
```

## Solvers and verification

```@docs
solve_operating_envelope
solve_operating_envelope_multistart
compare_operating_envelope_policies
verify_operating_envelope
search_operating_envelope_utilizations
search_operating_envelope_adversarial
solve_search_stable_operating_envelope
```
