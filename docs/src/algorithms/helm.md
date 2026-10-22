# HELM power flow

> **Kind:** Bespoke algorithm · **Maturity:** prototype · **Direction:** forward · **Temporal:** single-snapshot

PowerOptLab adds a deterministic, non-iterative multiconductor (4-wire and
beyond) power-flow solver based
on the **Holomorphic Embedding Load-flow Method** (HELM): [`solve_pf_helm`](@ref)
(standard result dictionary) and [`helm_series`](@ref) (programmatic
[`HelmResult`](@ref)). It is a *bespoke algorithm* — a new solution method built
on the BMOPFTools engine's public admittance-matrix API — and complements the
engine's Ipopt-based `solve_pf`: no initial guess and no convergence tuning.
It reports physical mismatch, Padé spread, and coefficient-tail diagnostics;
finite-order divergence is deliberately not presented as a proof that no
power-flow solution exists.

```julia
res = solve_pf_helm(net)
res["termination_status"]        # "HELM_CONVERGED" | "HELM_SERIES_DIVERGED" | "HELM_MAX_ORDER"
res["bus"]["lb"]["a"]["vm"]      # voltage magnitude (V)
res["helm"]["pade_spread"]       # last-two Padé differences, one per series row
res["helm"]["singularity_estimate"] # heuristic Domb–Sykes diagnostic
```

## The augmented (bordered) admittance matrix

Ideal elements — a closed switch, a zero-leakage transformer — have **no finite
admittance**: stamping a huge `1/z` corrupts conditioning, and a singular
zero-leakage block loses the coupling entirely. The engine's `ybus_augmented`
extends the passive Ybus with one **constraint row and auxiliary current
unknown per ideal coupling**:

```math
\mathbf{K} = \begin{bmatrix} \mathbf{Y} & \mathbf{A}^{\mathsf T} \\
\mathbf{A} & \mathbf{0} \end{bmatrix},
\qquad
\mathbf{K} \begin{bmatrix} \mathbf{V} \\ \mathbf{w} \end{bmatrix}
= \begin{bmatrix} \mathbf{I}_{\text{inj}} \\ \mathbf{0} \end{bmatrix}
```

Each row of ``\mathbf{A}`` is one `IdealCoupling` — a linear identity
``\mathbf{a}^{\mathsf T}\mathbf{V} = 0`` that holds *exactly*:

- **closed switch**, per conductor: ``V_{\text{from}} - V_{\text{to}} = 0``
  (with `switches = :constrain`; the default `:alias` keeps the node-fusion
  model of `ybus_passive` — both give identical voltages, `:constrain`
  additionally yields the switch current as a solution unknown);
- **ideal transformer**, per winding core: ``u_{w1} - r\,u_{w2} = 0`` over the
  coil voltages, with the ratio (tap included) from
  `_xfmr_winding_incidence` — the *same* seam the `transformer_yprim`
  builders use, so the admittance and constraint models cannot drift apart.
  Covered subtypes: `single_phase`, `wye_delta`, `delta_wye`; unity-ratio
  zero-leakage units stay node-aliased (an exact identity).

Because the constraint coefficients are real, the same vector is stamped as
row **and** column, preserving the module-wide reciprocity convention
``\mathbf{K} = \mathbf{K}^{\mathsf T}`` (plain transpose — never the adjoint).
Constraint rows are scaled by a nominal admittance for conditioning; the
physical coupling current is `coupling.scale * w`.

!!! note "Promoted transformers keep their shunt branches"
    `ybus_augmented` drops the whole Yprim of a transformer it promotes to an
    ideal coupling (or node-aliases at unity ratio), including its finite
    no-load (magnetising) shunt and its neutral-grounding impedance
    (`r/x_neutral_*`). HELM stamps that shunt-only part back on top of the
    augmented matrix, so an ideal Dyn or single-phase unit whose secondary is
    earthed only through its neutral resistor keeps that earth reference.

## The holomorphic embedding

The power-flow equations are embedded in a complex parameter ``s``. Voltage
sources stay fixed for all ``s``; loads enter through one of two embeddings:

| `embedding` | germ (``s = 0``) | load at ``s`` |
|---|---|---|
| `:load_scaling` (default) | energized **no-load** network | every load part scaled by ``s`` — ``s`` is a load multiplier |
| `:impedance_germ` | network loaded by **equivalent admittances** | constant-Z parts fixed in the matrix; other parts ``(1-s)\,y_g\Delta V + s\,I_{\text{law}}(\Delta V)`` |

In both cases the germ comes from one linear solve and ``s = 1`` is the
actual problem. The impedance germ determines nodes that the no-load germ
leaves undetermined (a star point held only by its loads); its ``s`` is not
a load multiplier, so the singularity estimate is then only a series
diagnostic.

**Load laws.** A sub-load across ``(p, q)`` draws
``I = \sum_\gamma \kappa_\gamma |\Delta V|^\gamma / \overline{\Delta V}``
(``\gamma = 0`` constant power, ``1`` constant current, ``2`` constant
impedance, any other exponent exponential). With the reflected series
``\tilde{\Delta V}(s) = \overline{\Delta V(\bar s)}`` — holomorphic, and equal to
``\overline{\Delta V}`` on the real axis — the product
``m = \Delta V\,\tilde{\Delta V}`` is a holomorphic stand-in for ``|\Delta V|^2``,
and ``|\Delta V|^\gamma = m^{\gamma/2}`` follows from J. C. P. Miller's
power-series recursion. Every ZIP and exponential model is therefore inside
the embedding; nothing is linearised.

**Source reference.** `source_reference = :earth` (default) fixes source
phase terminals to their phasors and the source-bus neutral to 0 V, as the
OPF does — a Dirichlet boundary eliminated from the unknowns.
`source_reference = :neutral` imposes each phasor between the phase and the
source neutral as a bordered constraint row (with the source current as its
unknown), so the source neutral floats with the network's actual earthing —
a grounding resistor at the source is honoured instead of being shorted.

Every voltage and coupling/source current is a holomorphic function of ``s``, expanded
as a power series whose order-``n`` coefficients solve

```math
\mathbf{K}_{UU}\, \mathbf{x}[n] = \mathbf{rhs}\big(\mathbf{x}[0..n{-}1]\big)
```

— **one LU factorization** serves the germ and every order. The nonlinear
terms use the conjugate-reflection trick with incremental inverse- and
power-series convolutions; phase-to-phase and delta-connected
loads come out of the same machinery (the load split is shared with
`ybus_linearized`, so the load model is byte-equivalent to the OPF's).
Constraint rows have identically zero right-hand side at every order: the
ideal-element identities hold term-by-term, hence exactly in the summed
solution.

The series is evaluated at ``s = 1`` by **Padé analytic continuation** (Wynn's
epsilon algorithm). Convergence is judged *physically*: the returned voltages
are plugged back into the full nonlinear current-mismatch equations.

## Convergence and singularity diagnostics

The implementation makes one positive claim: `HELM_CONVERGED` means the
Padé-evaluated state satisfies the full nonlinear current mismatch tolerance on
the no-load-connected branch. Its other statuses describe the finite series
that was computed; they do not certify non-existence:

| `termination_status` | meaning |
|---|---|
| `HELM_CONVERGED` | the operational solution (the branch continuously connected to no-load) |
| `HELM_SERIES_DIVERGED` | nonlinear mismatch failed and every exposed coefficient-tail norm ratio exceeds one |
| `HELM_MAX_ORDER` | nonlinear mismatch failed without the growing-tail classification; retrying with a larger `max_order` may help |

The raw diagnostics are available in both APIs:

- `pade_spread`: absolute difference between the final two Wynn-epsilon Padé
  estimates for every coefficient row;
- `coefficient_tail_norms`: max-norm of each final coefficient order;
- `coefficient_tail_ratios`: adjacent ratios used by `HELM_SERIES_DIVERGED`;
- `singularity_estimate`: Domb–Sykes extrapolation of the nearest
  coefficient-dominating singularity in ``s``.

On a simple radial constant-power saddle-node, `singularity_estimate` tracks
the analytic collapse multiplier. On a general network, another real or
complex singularity, finite-order effects, or cancellation can dominate the
coefficients. Treat it as a heuristic and validate any collapse interpretation
with continuation power flow or an analytically known system. `load_margin`
remains a compatibility alias for this value, but its old name must not be read
as a guarantee.

## Four-wire and n-wire modelling

Nothing is Kron-reduced: every declared `(bus, terminal)` is a node unless it
is perfectly grounded or fused by an ideal unity coupling. Neutral, PE or any
further conductor is an ordinary node with its own KCL row, coupled through
the full ``n \times n`` line impedance matrices. The recursion is generic over
node indices and terminal pairs, so the solver has no conductor-count
assumption — the test suite solves a 5-conductor feeder with separate N and PE
and a phase-to-PE leakage load. The remaining 4-wire-specific conventions live
at the data boundary, not in the solver:

- WYE loads and sources identify **one** neutral per bus by its label
  (`_neutral_terminal` / `_neutral_labels`); a bus with two neutral-like
  conductors needs `SINGLE_PHASE` per-pair loads to choose between them.
- Voltage sources are WYE phase-to-reference phasors (the BMOPF schema has no
  source configuration or internal impedance).
- Connection configurations (`WYE`, `DELTA`, `SINGLE_PHASE`) are expanded to
  terminal pairs by the shared `_load_subloads` decomposition. The solver
  consumes only those pairs, so an explicit per-pair connection list in the
  data model would make HELM fully n-wire with no solver change.

## Requirements and limitations

- **At least one voltage source** (WYE / SINGLE_PHASE; DELTA references raise
  an error).
- **Generators and IBRs raise an `ArgumentError`.** The BMOPF model gives them
  bounds and control profiles, not fixed power-flow injections, so HELM
  refuses them rather than solving a network that silently omits them.
- **Zero-leakage transformers outside the augmented model raise an
  `ArgumentError` naming the unit.** `ybus_augmented` has exact constraint
  rows only for `single_phase`, `wye_delta` and `delta_wye`. A
  `center_tap` with a zero HV or LV star arm, a zero-impedance
  `single_phase_autotransformer` or `open_delta_regulator`, or an `n_winding`
  unit with a singular short-circuit matrix would otherwise fall back to a
  shunt-only stamp (or be skipped) and act as an open circuit, which surfaces
  as a misleading "floating nodes" error. The same applies to non-unity ideal
  two-winding units under `ideal_xfmrs = :stamp`. Leakage at or below
  `xfmr_z_min_ohm` (default 1e-6 Ω) counts as zero. Just above that threshold
  the stamped admittance (~1e5 S) can leave a round-off floor of ~1e-6 A in
  the absolute current mismatch, so near-ideal units may report
  `HELM_MAX_ORDER` with accurate voltages; from about 1e-4 Ω they converge
  normally.
- **Every node needs a linear path to a reference at the germ.** Islands and
  floating ideal deltas are diagnosed with an error naming the floating
  `(bus, terminal)` nodes. Nodes held only by loads need
  `embedding = :impedance_germ`.
- **Unearthed star points with constant-power or constant-current loads are
  often ill-posed.** Their neutral can be non-unique (balanced constant-power
  stars admit a continuum of displaced neutrals) or non-existent (constant
  current magnitudes that violate the triangle inequality cannot sum to zero).
  HELM then reports `HELM_SERIES_DIVERGED` with a singularity estimate below
  one instead of returning one branch arbitrarily.
- Time series are not applied; HELM solves the snapshot values in the network
  dictionary.

## Validation

The test suite (`test/helm_tests.jl`) checks, always:

- Padé continuation against series with known limits — including series that
  diverge at ``s = 1`` but continue analytically, and one with *no* limit
  (the spread indicator must refuse to lock in);
- a closed-form two-node feeder: voltages to `rtol 1e-9`, the neutral-rise
  detail, and `singularity_estimate` against the analytic saddle-node
  ``P^* = E^2/4R`` (`≈ 6.6`, `≈ 2`, and `< 1` past collapse);
- a warm-started natural-parameter continuation trace whose nonlinear Ipopt
  correctors bracket that analytic saddle-node, checked against HELM convergence
  and series-divergence classifications;
- switch `:alias ≡ :constrain` voltage identity with exact conductor currents;
- delta and line-to-line loads against a rotated closed form;
- constant-current, exponential and full-ZIP laws against closed forms on the
  two-node feeder (`rtol 1e-11`, both embeddings);
- `source_reference = :neutral` against an independent MNA solve with source
  and load neutral earthing resistors;
- an ideal single-phase transformer whose secondary is earthed only through
  its neutral resistor, against a near-ideal stamped unit;
- an unearthed star point (impedance germ) and a 5-conductor N + PE feeder;
- zero-leakage `center_tap`, autotransformer, open-delta and n-winding units
  (and an ideal unit under `ideal_xfmrs = :stamp`) are refused by name rather
  than reported as floating nodes;
- the `ybus_linearized` residual oracle
  ``\lVert \mathbf{Y}\mathbf{V} - \mathbf{i}_{\text{comp}}(\mathbf{V})\rVert_\infty \approx 0``
  at every returned solution (two independent code paths).

Gated on OpenDSSDirect: node-to-earth voltage parity on the `pf_comparison`
decks (sub-mV on line/delta/capacitor/ZIP decks; a few mV through Yd/Dy
transformers), and a 3-way HELM ≈ Ipopt `solve_pf` ≈ OpenDSS oracle.

## API reference

The [HELM API](../api/helm.md) lists the public types, functions and result fields
for this guide.
