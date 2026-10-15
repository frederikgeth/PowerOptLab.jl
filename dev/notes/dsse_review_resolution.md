# DSSE review: assessment and implementation

The blocking findings were valid. An ordinary loaded four-wire feeder exposed
both Ipopt stopping/acceptance problems and inadequate L-BFGS curvature. The
legacy current-variable model had also inherited restrictions that belong only
to voltage elimination. Agreement on small component cases did not establish
robustness on a feeder or equivalence of the two input contracts.

The implementation followed four steps: restore injection capabilities; repair
the Ipopt backend and reproduce feeder failures; align statistical/API semantics;
then broaden independent and statistical validation. The table records the
result and the boundaries that remain.

| Concern | Assessment and action |
|---|---|
| Ipopt rejects an optimum because of its termination label | Valid. Every available primal now receives independent feasibility, stationarity, rank and curvature checks. The optimizer's original status remains separate. Failed physical checks still prevent publication. |
| L-BFGS stalls on a larger feeder | Valid. The default supplies Gauss–Newton curvature, sparse symbolic derivative patterns and per-iterate caches. Independent 5/20/40-bus four-wire circuits compare states, J and dof against dense, sparse and native WLS. |
| Generic compiled solve selects the least robust backend | Valid default-selection concern. The overload now selects the sparse solver; explicit names remain preferable for comparisons. Native and compiled overloads return their documented formulation-specific result types. |
| Legacy overlapping paths and different merged returns regress | Valid. Only conflicting returns for the same original customer are rejected. The legacy current-variable formulation retains its broader path support. |
| Compiled mixed delta/wye customers and partially exact shared neutrals are excluded | Valid model gap. Default compilation adds explicit unknown currents for those cases and for merged customers. Exact and unknown device currents coexist in conductor KCL, including an additional explicitly declared unknown on the same phase/neutral endpoints. The shared-bus ZIP derivative regression is restored at its original bus. |
| Merged profiling imposes an artificial low-voltage domain | Valid for a general API. Default merged-current equations use finite currents and `S=U*conj(I)` at zero voltage. Optional `formulation=:voltage` still requires nonzero voltage for analytical power profiling; its domain guard is retained and tested as an explicit narrower contract. Physical nonzero constant-power device laws remain undefined at zero connection voltage. |
| Undeclared conductors have different meanings | Valid. Native KCL continues to impose zero injection and warns about undeclared conductors. Compiled preflight reports `:undeclared_injection` and leaves those currents free. Comparisons declare all intended zero injections explicitly. |
| Public objective and degrees of freedom differ | Valid. All solver diagnostics report J as the full squared residual sum, its measurement/prior components, residual count, nuisance dimension and dof. Explicit current directions enter the measurement rank; profiled nuisance powers are subtracted separately. The tutorial uses the same scale. |
| Legacy finite differences defeat a roundoff-level projection floor | Valid. Ordinary legacy diagnostics use ForwardDiff. Overlapping current paths use analytic Jacobians and eliminate current/nuisance tangent directions for the voltage diagnostic. |
| The winding identity proof can fail open | Valid. Unknown element groups prevent the structural reduction. Finite electrodes retain their equations regardless of strength. |
| Two reductions change `constraint_rank` semantics | The two proofs have different purposes, but reported rank is computed from the same complete physical Jacobian in every backend. Structural identities are removed for all solvers; Ipopt additionally selects optimizer rows. Feeder comparisons assert equal reported ranks, and optimizer row mappings remain separate diagnostics. |
| Open-delta connectivity belongs upstream | Valid. BMOPFTools PR [#436](https://github.com/frederikgeth/BMOPFTools.jl/pull/436) owns the shared-phase copper bond in passive node assembly and the native model. PowerOptLab removes the synthetic-switch/deep-copy patch and pins that upstream fix. |
| Repeated dense equation reduction and evaluator work | Valid. The coefficient signature and QR are sparse. A bounded synchronized cache reuses a verified proof for reading/uncertainty/prior updates; physical parameter changes invalidate it. Residual and constraint Jacobians are cached separately at an iterate. Final SVD rank/covariance remains dense and is not a large-network scaling solution. |
| Current-device angle convention is ambiguous | Valid gap. `ConstantCurrentDevice` documents the fixed-phasor default and supports `angle_reference=:voltage` for constant magnitude and voltage-relative power factor. An independent two-conductor equilibrium checks the latter. ZIP's existing fixed-phasor current term is explicit; exact models can be composed. |
| N-winding excitation is absent | The implementation claim was incorrect: `g_no_load`/`b_no_load` already stamp a winding-2 coil branch. The coverage concern was valid. The independent star circuit and OpenDSS terminal equations now exercise nonzero excitation across mixed connections and DELTA rotations. Nonlinear saturation is not implemented. |
| N-winding tap fields and ideal equipment are unsupported | Correct limitations. Fixed winding ratios use `v_nom`; ignored tap fields are rejected. Ideal transformers, ideal-bond reaction currents and floating phase-neutral source laws need further equipment equations/current states. Unknown customer-current states do not manufacture those laws. |
| Transformer subtype leaks into meter construction | Valid API concern. IDs are looked up across subtypes; `variant` is optional and only needed for duplicated IDs. |
| Dead result compatibility constructor | Removed. No compatibility claim is made for a type newly introduced in the same PR. |
| Time-series kwargs target the wrong backend | Valid API concern. `solver_kwargs` holds per-backend settings. Trust-region-only options passed to Ipopt produce an explanatory `ArgumentError`, rather than an internal method error or silent omission. |
| Sign conventions invite mistakes | The conventions are retained explicitly: positive customer power injects into the network, positive exact-device power/current is consumption. Tests put both on one electrical node/shared neutral and check independent KCL and voltage expectations. |
| Ipopt coverage omits older injection/quality/regression suites | Fixed. Those solver loops now include shared Ipopt. Trust-region-specific assertions remain specific to the algorithms that implement them. |
| Private tests are brittle | Public states, predictions, selective covariance and diagnostics cover end-to-end behavior. A few private coefficient-identity unit tests remain deliberately: they test proof construction that a successful solve alone cannot establish. |
| No statistical calibration | A seeded 500-draw nonlinear merged-customer case checks E[J] against dof and empirical voltage covariance against the local calculation. It also compares explicit-current and profiled likelihoods. This is a finite low-noise regression, not universal coverage or bad-data validation. |
| Style and assertion counts obscure quality | Reworked numerical files expand dense statements. Validation descriptions lead with circuit diversity, independent equations, uncertainty calibration and feeder size; assertion totals are only execution records. |
| Local manifest disagrees with the declared dependency | The dependency pin includes the upstream fix while preserving the declared API base. Reproducible runs resolve/instantiate that pin; the ignored local manifest is refreshed to this declared commit. Runs on a stale manifest are not evidence for compatibility with a different upstream API. |

Remaining work includes ideal-device equipment equations, nonlinear core models,
correlated errors, calibrated robust/bad-data processing, realistic sparse
telemetry/dropout studies, multi-start behavior and sparse rank/covariance
algorithms. These are recorded feature boundaries, rather than capabilities
inferred from solver agreement.

## Validation record

The 16 October 2026 run uses BMOPFTools compatibility commit
`8ddec1e557e6bb405567b26b810b5b14cbab5274`. Both source declarations and ignored
local manifests are refreshed with `scripts/instantiate_pinned.jl`. Companion
PR #436 also integrates the same fix with current upstream main; PowerOptLab
keeps the existing upstream API base rather than adopting that separate API
migration.

On Julia 1.12.6 the full package suite passes 15,718 checks plus 11 Aqua checks.
The complete documentation build/doctests, seven executable study tutorials
and three reproducible figure checks pass. The final review regressions also
pass directly against the root manifest's immutable dependency source.

The independent nonlinear calibration has mean J=2.90376 for 3 dof and relative
voltage-covariance error 0.16066 over 500 draws. Its finite regression budgets are
recorded in the test; it does not establish universal coverage.
