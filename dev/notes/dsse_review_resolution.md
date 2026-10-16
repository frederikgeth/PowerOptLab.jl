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
| Public objective and degrees of freedom differ | Valid. All solver diagnostics report J as the full squared residual sum, its measurement/prior components, residual count, nuisance dimension and dof. Explicit current directions enter the measurement rank; profiled nuisance powers are subtracted separately. Fully fixed states retain all residual degrees of freedom, checked against J=13/dof=2 in every solver. The tutorial uses the same scale. |
| Legacy finite differences defeat a roundoff-level projection floor | Valid. Ordinary legacy diagnostics use ForwardDiff. Overlapping current paths use analytic Jacobians and eliminate current/nuisance tangent directions for the voltage diagnostic. |
| The winding identity proof can fail open | Valid. Unknown element groups prevent the structural reduction. Finite electrodes retain their equations regardless of strength. |
| Two reductions change `constraint_rank` semantics | The two proofs have different purposes, but reported rank is computed from the same complete physical Jacobian in every backend. Structural identities are removed for all solvers; Ipopt additionally selects optimizer rows. Feeder comparisons assert equal reported ranks, and optimizer row mappings remain separate diagnostics. |
| Open-delta connectivity belongs upstream | Valid. BMOPFTools PR [#436](https://github.com/frederikgeth/BMOPFTools.jl/pull/436) owns the shared-phase copper bond in passive node assembly and the native model. PowerOptLab removes the synthetic-switch/deep-copy patch and pins that upstream fix. |
| Repeated dense equation reduction and evaluator work | Valid. The coefficient signature and QR are sparse. A bounded synchronized cache reuses a verified proof for reading/uncertainty/prior updates; physical parameter changes invalidate it. Residual and constraint Jacobians are cached separately at an iterate. Final SVD rank/covariance remains dense and is not a large-network scaling solution. Native WLS current-observability diagnostics also compile the model twice and use dense Jacobians, null spaces and `pinv(A)`; that diagnostic is a separate runtime/memory limit. |
| Current-device angle convention is ambiguous | Valid gap. `ConstantCurrentDevice` documents the fixed-phasor default and supports `angle_reference=:voltage` for constant magnitude and voltage-relative power factor. An independent two-conductor equilibrium checks the latter. ZIP's existing fixed-phasor current term is explicit; exact models can be composed. |
| N-winding excitation is absent | The implementation claim was incorrect: `g_no_load`/`b_no_load` already stamp a winding-2 coil branch. The coverage concern was valid. The independent star circuit and OpenDSS terminal equations now exercise nonzero excitation across mixed connections and DELTA rotations. Nonlinear saturation is not implemented. |
| N-winding tap fields and ideal equipment are unsupported | Correct limitations. Fixed winding ratios use `v_nom`; ignored tap fields are rejected. Ideal transformers, ideal-bond reaction currents and floating phase-neutral source laws need further equipment equations/current states. Unknown customer-current states do not manufacture those laws. |
| Transformer subtype leaks into meter construction | Valid API concern. IDs are looked up across subtypes; `variant` is optional and only needed for duplicated IDs. |
| Dead result compatibility constructor | Removed. No compatibility claim is made for a type newly introduced in the same PR. |
| Time-series kwargs target the wrong backend | Valid API concern. `solver_kwargs` holds per-backend settings. Every backend validates the selected options and produces an explanatory `ArgumentError` for unsupported keys, rather than an internal method error or silent omission. |
| Sign conventions invite mistakes | The conventions are retained explicitly: positive customer power injects into the network, positive exact-device power/current is consumption. Tests put both on one electrical node/shared neutral and check independent KCL and voltage expectations. |
| Ipopt coverage omits older injection/quality/regression suites | Fixed. Those solver loops now include shared Ipopt. Trust-region-specific assertions remain specific to the algorithms that implement them. |
| Private tests are brittle | Public states, predictions, selective covariance and diagnostics cover end-to-end behavior. A few private coefficient-identity unit tests remain deliberately: they test proof construction that a successful solve alone cannot establish. |
| No statistical calibration | A StableRNG 500-draw nonlinear merged-customer case checks E[J] against dof and mean whitened state error against its local chi-square prediction. Independent 5/20/40-bus feeders assert noiseless voltage truth and noisy state-error calibration with the constrained covariance pseudoinverse. It also compares explicit-current and profiled likelihoods. This is a finite low-noise regression, not universal coverage or bad-data validation. |
| Style and assertion counts obscure quality | Reworked numerical files expand dense statements. Validation descriptions lead with circuit diversity, independent equations, uncertainty calibration and feeder size; assertion totals are only execution records. |
| Local manifest disagrees with the declared dependency | The dependency pin includes the upstream fix while preserving the declared API base. Reproducible runs resolve/instantiate that pin; the ignored local manifest is refreshed to this declared commit. Runs on a stale manifest are not evidence for compatibility with a different upstream API. |

Remaining work includes ideal-device equipment equations, nonlinear core models,
correlated errors, calibrated robust/bad-data processing, realistic sparse
telemetry/dropout studies, multi-start behavior and sparse rank/covariance
algorithms. These are recorded feature boundaries, rather than capabilities
inferred from solver agreement.

## Validation record

The earlier 16 October 2026 runs used BMOPFTools compatibility commit
`8ddec1e557e6bb405567b26b810b5b14cbab5274`. Both source declarations and ignored
local manifests are refreshed with `scripts/instantiate_pinned.jl`. Companion
PR #436 subsequently integrated current upstream main and was merged. The
latest-pin follow-up below supersedes this historical compatibility source.

On Julia 1.12.6 the full package suite passes 15,718 checks plus 11 Aqua checks.
Julia 1.10.11 passes the 15,720-check suite plus its 10 enabled Aqua checks;
the two explicit-zero warning checks were added between those runs. After the
full runs, all 2,527 review checks pass on both versions, including 17 new
fully fixed-state checks. The Julia 1.12 run uses the actual root manifest's
immutable dependency source. Both versions pass the complete documentation
build/doctests, seven executable study tutorials and three reproducible figure
checks.

Before the StableRNG follow-up, Julia 1.12.6 nonlinear calibration had mean J=2.90376 for
3 dof and relative voltage-covariance error 0.16066 over 500 draws; Julia 1.10.11
gave 3.04283 and 0.04348 respectively. Those MersenneTwister streams differed
across Julia versions, so these historical values are not comparable draws. The
follow-up uses StableRNGs and a mean whitened-error test instead of a relative
Frobenius-error threshold. The chi-square approximation is local for this
nonlinear estimator; it does not establish universal coverage.

The upstream compatibility tree passes its full suite and required scientific,
Finding-registry, metadata and execution-interface gates. The conflict-free
integration with upstream main passes those gates and its full suite as well.

## Final critical pass

The final PR review found two additional estimator boundaries and one CI setup
failure:

- A fully fixed voltage model could retain declared unknown current states while
  `initial_state_estimator` returned an empty vector. Initialization now preserves
  those states. Regressions check their unobservability and identification by a
  current prior across all three shared solvers.
- A zero-injection declaration touching an unknown customer path changed meaning
  between voltage elimination and explicit currents: a reproduction gave J=10000
  versus essentially zero. Preflight now rejects the conflicting declarations,
  including electrical aliases and persistent paths without telemetry. Known
  exact devices remain supported; net meters can observe a known zero without
  declaring an unknown customer.
- Current-Julia optional CI failed during incremental Pkg precompilation with
  `Unknown package: MathOptIIS`. Optional setup now resolves/develops its complete
  environment with automatic precompilation deferred, then explicitly precompiles
  that environment once. The production dependency pin is unchanged. Setup has
  been verified on Julia 1.13.2, the version from the failing CI job.

The expanded DSSE/parameter numerical suite passes 6,454 checks on Julia 1.12.6.
Final validation also exercises the optional solver integrations and reruns the
full package suite, documentation and executable tutorials; the PR body and
head-specific CI checks record those results.

Linux Julia 1.10 CI subsequently exposed a platform-dependent oracle issue:
an impedance-grounded delta phase was exactly zero, whereas local SVD roundoff
left it tiny and nonzero. The mixed-winding test had requested an unsmoothed
voltage-magnitude derivative at that point. Electrode phase telemetry now uses
rectangular voltage channels, retaining displaced-neutral magnitudes. A new
deterministic exact-zero regression checks the physical solution, rank and all
shared solvers, and explicitly tests rejection/smoothing of the nonsmooth
magnitude channel. The production derivative guard and numerical tolerances
are unchanged.

## Independent-truth and calibration follow-up

The additional review correctly identified that the multi-bus fixture returned
independent conductor voltages without asserting against them. The 5/20/40-bus
suite now verifies all four backends against noiseless truth at 1e-8 V and checks
each shared backend's noisy error against the constrained covariance
pseudoinverse. The covariance rank comes from the reported tangent dimension, checked against
6n independent voltage directions for n downstream buses; full observability
and null-direction error are also checked explicitly.
A five-standard-deviation budget is an approximate nonlinear regression bound,
not a claim of exact finite-sample coverage.

The 500-draw merged-customer test uses a test-only StableRNGs dependency. Its
mean squared whitened voltage error is compared with 2 using a four-standard-error
budget (2/sqrt(500) per standard error under the linear Gaussian model). Relative
Frobenius covariance error remains logged for diagnosis without an arbitrary
pass/fail threshold. Both fixtures now use the same random streams across Julia
versions.

Preflight emits `:info` / `:reinterpreted_customer_measurement` for each original
P/Q row whose known-device meaning changes when an additional unknown customer
is explicitly declared on the same endpoints. This preserves the documented
semantics while exposing the change. Time-series option validation is symmetric
across dense, sparse and Ipopt; unsupported options, misspelled backend keys and
invalid settings containers produce `ArgumentError`. Remaining multiple-statement
lines in the identified device/preflight paths were expanded, and the native
observability diagnostic's double compilation and dense pseudoinverse are now
listed among scaling limits.

Validation of this follow-up: the full package suite passes 15,903 assertions on
both Julia 1.10.11 and 1.12.6, plus 10 and 11 Aqua checks respectively. The same
StableRNG draws give mean J=2.8476751924 and mean q=1.8738165957 on both runtimes
(to roundoff), with q's reference standard error 0.0894427191. The noiseless
shared-solver probe reaches about 1e-11 V maximum coordinate error or better;
the regression retains the requested 1e-8 V bound on complex conductor voltage.
Full Documenter/doctests, seven executable study tutorials and three reproducible
figure checks also pass. This supplements the earlier CI record; new CI must
be checked on the follow-up commit itself.


## Merged upstream dependency pin

After upstream PR #436 merged, both root and documentation source declarations
were advanced to BMOPFTools main `df6bc42bfd64fc5475658c438adcdb28a3469777`.
The merged source retains the shared-phase open-delta bond and also includes
PR #435's canonical `open_switch` semantics. PowerOptLab preflight now uses
that Boolean field (absent means closed), reports the corresponding topology,
and rejects non-Boolean values or the unsupported legacy `status` field.
Fixtures now use `open_switch=false`; a regression checks default/closed/open
node aliasing, distinct recovered voltages across open terminals, residual dof,
and invalid declarations. The production version bound remains 0.1.0 and the
source is an immutable merged main commit, not an unmerged compatibility branch.

Validation at this merged source: full package suites pass 15,963 checks on both
Julia 1.10.11 and current stable 1.13.2, plus 10 and 11 Aqua checks respectively.
The added switch regression contributes 60 checks. Full documentation builds
and doctests pass on both versions, with deployment disabled. Julia 1.13.2
also passes seven executable study tutorials and eight executable DOE tutorials;
all three reproducible figure checks pass. Root, docs and both test manifests
resolve the immutable merged source rather than a path override.
