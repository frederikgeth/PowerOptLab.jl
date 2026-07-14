# HELM power flow

PowerOptLab adds a deterministic, non-iterative 4-wire power-flow solver based
on the **Holomorphic Embedding Load-flow Method** (HELM): [`solve_pf_helm`](@ref)
(standard result dictionary) and [`helm_series`](@ref) (programmatic
[`HelmResult`](@ref)). It is a *bespoke algorithm* — a new solution method built
on the BMOPFTools engine's public admittance-matrix API — and complements the
engine's Ipopt-based `solve_pf`: no initial guess, no convergence tuning, and —
uniquely — a **certificate** when no power-flow solution exists, plus an estimate
of the distance to voltage collapse from every solve.

```julia
res = solve_pf_helm(net)
res["termination_status"]        # "HELM_CONVERGED" | "HELM_NO_SOLUTION" | "HELM_MAX_ORDER"
res["bus"]["lb"]["a"]["vm"]      # voltage magnitude (V)
res["helm"]["load_margin"]       # collapse loading multiplier (see below)
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

!!! note "Promoted transformers are purely ideal"
    A transformer promoted to an ideal coupling contributes *only* the ratio
    identity: its no-load (magnetising) shunt and neutral-grounding branch are
    dropped along with the singular series stamp.

## The holomorphic embedding

The power-flow equations are embedded in a complex parameter ``s``: voltage
sources stay fixed for all ``s`` (a Dirichlet boundary, eliminated from the
unknowns), and every load scales with ``s`` — constant-power draws ``s\,S``,
constant-impedance admittances become ``s\,\mathbf{Y}_Z``.

- ``s = 0`` is the **germ**: the energized no-load network — genuinely
  unbalanced, floating neutrals at their equilibrium — obtained from one
  linear solve.
- ``s = 1`` is the actual problem.

Every voltage and coupling current is a holomorphic function of ``s``, expanded
as a power series whose order-``n`` coefficients solve

```math
\mathbf{K}_{UU}\, \mathbf{x}[n] = \mathbf{rhs}\big(\mathbf{x}[0..n{-}1]\big)
```

— **one LU factorization** serves the germ and every order. The nonlinear
``S^*/\Delta V^*`` term uses the classic conjugate-reflection trick with an
incremental inverse-series convolution; phase-to-phase and delta-connected
loads come out of the same machinery (the load split is shared with
`ybus_linearized`, so the load model is byte-equivalent to the OPF's).
Constraint rows have identically zero right-hand side at every order: the
ideal-element identities hold term-by-term, hence exactly in the summed
solution.

The series is evaluated at ``s = 1`` by **Padé analytic continuation** (Wynn's
epsilon algorithm). Convergence is judged *physically*: the returned voltages
are plugged back into the full nonlinear current-mismatch equations.

## Certified non-existence and the loading margin

By Stahl's theorem the diagonal Padé sequence converges wherever the power-flow
solution exists and provably fails to converge where it does not. HELM
therefore has a three-way outcome, not a "diverged, try another start":

| `termination_status` | meaning |
|---|---|
| `HELM_CONVERGED` | the operational solution (the branch continuously connected to no-load) |
| `HELM_NO_SOLUTION` | **certified voltage collapse** — no power-flow solution exists at this loading |
| `HELM_MAX_ORDER` | series order exhausted; retry with a larger `max_order` |

Because the embedding parameter scales the loads, the series' **radius of
convergence is the collapse loading multiplier** ``\lambda^*``. It is estimated
from the coefficient-ratio tail (Domb–Sykes extrapolation) and reported as
`load_margin` on every solve:

- ``\lambda^* = 2`` — the present operating point is at half its collapse
  loading;
- ``\lambda^* < 1`` — explains a `HELM_NO_SOLUTION` (the solvable range ends
  before ``s = 1``);
- `NaN` — the series is too short or featureless to extrapolate (typically
  very light or constant-impedance-only loading, where the margin is far away).

So every routine HELM power flow doubles as a per-feeder distance-to-collapse
measurement — no continuation power flow required.

## Requirements and v1 limitations

- **At least one voltage source** (WYE / SINGLE_PHASE; DELTA references raise
  an error). Source phase terminals are fixed to their reference phasors and
  the source-bus neutral to 0 V, mirroring the OPF's source semantics.
- **Every node needs a linear path to a reference.** Islands, floating ideal
  deltas, and nodes held only by loads (whose no-load germ voltage is
  genuinely indeterminate) are diagnosed with an error naming the floating
  `(bus, terminal)` nodes.
- **Loads: constant-power + constant-impedance parts only** (including those
  ZIP fractions). Constant-current fractions and non-integer exponential
  models involve ``|\Delta V|``, which is not holomorphic in this embedding —
  they raise a validation error naming the loads. Use the
  `ybus_linearized` fixed-point map or the OPF power flow for those;
  an outer-loop treatment is a planned follow-up.
- Generators and IBRs are not yet injected (as in `ybus_linearized`).

## Validation

The test suite (`test/helm_tests.jl`) checks, always:

- Padé continuation against series with known limits — including series that
  diverge at ``s = 1`` but continue analytically, and one with *no* limit
  (the spread indicator must refuse to lock in);
- a closed-form two-node feeder: voltages to `rtol 1e-9`, the neutral-rise
  detail, `load_margin` against the analytic collapse point (`≈ 6.6`, `≈ 2`,
  and `< 1` past collapse), certified no-solution beyond ``P^* = E^2/4R``;
- switch `:alias ≡ :constrain` voltage identity with exact conductor currents;
- delta and line-to-line loads against a rotated closed form;
- the `ybus_linearized` residual oracle
  ``\lVert \mathbf{Y}\mathbf{V} - \mathbf{i}_{\text{comp}}(\mathbf{V})\rVert_\infty \approx 0``
  at every returned solution (two independent code paths).

Gated on OpenDSSDirect: node-to-earth voltage parity on the `pf_comparison`
decks (sub-mV on line/delta/capacitor decks; a few mV through Yd/Dy
transformers), and a 3-way HELM ≈ Ipopt `solve_pf` ≈ OpenDSS oracle.
