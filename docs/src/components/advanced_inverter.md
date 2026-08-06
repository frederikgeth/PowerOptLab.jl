# Advanced inverter (prototype)

> **Kind:** Component model · **Maturity:** prototype · **Direction:** forward · **Temporal:** single-snapshot

[`AdvancedInverter`](@ref) is a more detailed inverter-based-resource (IBR) than
the BMOPFTools engine's built-in current-injection IBR. The engine models an IBR
as a bounded current source at the point of connection (POC); this prototype adds
the core structural idea from the BMOPFTools
[IBR model extensions design doc](https://github.com/frederikgeth/BMOPFTools.jl/blob/main/docs/ibr_model_extensions.md)
— an explicit **internal AC node** behind the converter — and the five feature
phases layered on it.

```
 POC bus ──[filter r+jx (+ grid shunt b)]── internal node ──[converter]── DC
   network sets V here                       EMF lives here          losses/ripple here
```

On top of that structure it carries exact three-phase **feasible-region models**
from ongoing research — 3-leg (3-wire), 4-leg, and split-DC-link converters —
whose DC-utilisation limit is the exact time-sampled switching-polytope condition
with an endogenous 2ω bus-ripple derating and neutral-current limits (see
[Three-phase topology models](@ref three-phase-topology-models) below).

It is built entirely on the BMOPFTools staged API through a `model_hook!`; it
does **not** modify the engine. Device parameters are SI; the solve runs in SI
(`per_unit=false`) or per-unit (`per_unit=true`), scaling every parameter to model
units via `ctx.bases` — the DC-side quantities (`v_dc`, `c_dc`, `In_max`) stay SI
and the AC↔DC coupling scales through the POC bus's `v_base`/`i_base`/`s_base`.
Results are returned in SI in both modes.

## The five phases

| Phase | Feature | Model |
|---|---|---|
| 0 | Output filter (L/LC) | series `r+jx` from POC to the internal node, optional grid-side shunt `b` |
| 1 | Internal EMF / DC utilisation | `|V_int|` box; single-phase DC-link modulation `|V_int| ≤ modulation_max·v_dc/√3`; or a three-phase `topology`'s exact switching polytope (below) |
| 2 | Grid-forming | balanced 120° internal EMF with a bounded magnitude decision variable |
| 3 | Converter losses | non-branching `P_dc = P_ac + P_loss`, `P_loss = p_loss_fixed + a_loss·|I| + c_loss·|I|²` |
| 4 | Double-frequency ripple | single-phase cap `|Σ_k V_int_k·I_k|² ≤ p_ripple_max²`; three-phase topologies form the 2ω bus-ripple phasor that derates the DC rails (below) |

`p_poc` and `q_poc` are the total grid-side exchange. In particular, `q_poc`
and `q_set` include the optional grid-side shunt; `q_conv` is the converter-side
quantity before the filter and shunt.

Every feature is opt-in: with only `id`, `bus`, and `s_max` the device is a plain
grid-following converter, and the internal node collapses onto the POC when the
filter is zero.

## Key modelling choices (from the design doc)

- **Limits on the converter side.** The apparent-power circle `s_max` and the
  current limit `i_max` are applied on the converter quantities (internal-node
  voltage × current), matching real nameplate — so an output filter reduces the
  power actually delivered to the grid below the converter rating.
- **Non-branching losses.** With AC power positive = injected to grid and DC power
  positive = drawn from the DC source, the single equation `P_dc = P_ac + P_loss`
  (`P_loss ≥ 0`) holds for both discharge and charge — no direction `if`-branch.
- **Grid-forming ≠ slack.** A grid-forming inverter holds a balanced, bounded
  internal EMF *behind the filter*, but does not replace the network's reference;
  the surrounding grid still needs a slack source.

## Worked example

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf

# A stiff grid: slack at "grid", short line to the inverter POC.
net = parse_bmopf("""
{"bus":{
    "grid":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "poc": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],"v_min":[200.0],"v_max":[250.0]}},
 "voltage_source":{"vs":{"bus":"grid","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.05}},
 "line":{"l1":{"bus_from":"grid","bus_to":"poc","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
"""; from_string=true)

# A converter with an output filter and a three-term loss curve; minimise loss
# while delivering 3 kW to the grid.
inv = AdvancedInverter(id="inv", bus="poc", s_max=5000.0,
                       r_filter=0.2, x_filter=0.5,
                       p_loss_fixed=20.0, a_loss=0.3, c_loss=0.02)

r = solve_advanced_inverter(net, inv; objective=:min_loss, p_set=3000.0)

r.p_poc     # ≈ 3000 W delivered at the POC
r.q_poc     # total POC reactive exchange, including any grid-side shunt
r.p_conv    # converter-side active power (> p_poc: filter losses)
r.p_loss    # 20 + 0.3·|I| + 0.02·|I|²
r.p_dc      # = p_conv + p_loss  (the non-branching DC-link balance)
r.v_int_mag # internal EMF magnitude per phase (V)
```

Switch `objective=:max_export` to maximise POC active power and watch the
converter rating, filter, EMF/modulation, or ripple limits bind. For a
three-phase `grid_forming=true` inverter the solved internal EMF magnitudes are
equal across phases (balanced 120°) and the 2ω ripple is ≈ 0.

### Choosing a three-phase topology

Set `topology` to one of `:THREE_LEG`, `:FOUR_LEG`, or `:SPLIT_DC` (with `v_dc`,
`c_dc`, and `In_max` for the 4-wire ones) to use the exact switching-polytope
model:

```julia
inv = AdvancedInverter(id="inv", bus="poc", phase_terminals=["a","b","c"], neutral="n",
                       topology=:FOUR_LEG, s_max=20e3, i_max=40.0,
                       v_dc=700.0, c_dc=1.1e-3, In_max=40.0, m_max=0.96,
                       r_filter=0.05, x_filter=0.15)
r = solve_advanced_inverter(net3, inv)   # net3 = a three-phase grid
r.i_neutral   # neutral current (A) — non-zero only under unbalance
r.dv2         # 2ω bus-ripple amplitude (V) that derated the DC rails
r.i_cap       # capacitor RMS ripple current (A) — compare against i_cap_max
```

On a balanced grid all three topologies coincide (no neutral current, no ripple).
Under unbalance the 4-leg and split-DC draw neutral current (bounded by `In_max`),
and the split-DC needs a **higher `v_dc`** for the same per-phase voltage — the
half-bus utilisation penalty of the split-capacitor structure.

See the API reference for [`AdvancedInverter`](@ref),
[`solve_advanced_inverter`](@ref), and [`InverterResult`](@ref).

## Scope

This is a **prototype** for experimentation, not a validated engine feature. It
implements the design doc's Phases 0–4 plus the three-phase topology models as a
hook-stamped device; the reactive/active priority modes, sequence-current limits,
and grid-forming-as-reference capabilities listed in the design doc's backlog are
not included. If a piece of this matures, it can be folded back into the engine.

## [Three-phase topology models](@id three-phase-topology-models)

For the three-phase topologies (`:THREE_LEG`, `:FOUR_LEG`, `:SPLIT_DC`) the
crude scalar modulation cap is replaced by the **exact time-sampled
switching-polytope** feasibility from ongoing research on converter feasible
regions (fundamental-frequency RMS phasors; the classical pulse-width-modulation
DC-utilisation and DC-link-ripple results of the power-electronics literature).
All constraints apply to the **converter output** `U_x = V_int_x` (the internal
node), so the filter, losses, grid-forming, and `s_max` circle all compose. The
equations below are what the code stamps (shown in SI; per-unit is the same after
base scaling).

**Oscillating (2ω) power.** The unconjugated phase sum

```math
\tilde S = \sum_{x\in\{a,b,c\}} U_x\, I_x,\qquad
\tilde S_{re} = \textstyle\sum_x (U^{re}_x I^{re}_x - U^{im}_x I^{im}_x),\;\;
\tilde S_{im} = \textstyle\sum_x (U^{re}_x I^{im}_x + U^{im}_x I^{re}_x)
```

is the amplitude of the double-frequency power pulsation the DC-link capacitance
absorbs (the two bilinear equalities are the only nonconvexity in the model).

**Bus-ripple phasor.** With DC capacitance `C_eq` (single cap `C_eq = C_dc`;
split link's series caps give `C_eq = C/2`), the 2ω bus voltage ripple is a
phasor `D = j\,\tilde S/(2\omega C_{eq} V_{dc})`, i.e.

```math
D_{re} = -\tilde S_{im}/(2\omega C_{eq} V_{dc}),\qquad
D_{im} =  \tilde S_{re}/(2\omega C_{eq} V_{dc}),
```

giving the instantaneous DC rail `v_{dc}(\theta) = V_{dc} + D_{re}\cos2\theta -
D_{im}\sin2\theta`. An optional `dv2_max` caps `\sqrt{D_{re}^2+D_{im}^2}`.

**Sampled voltage feasibility.** Over a uniform grid `θ_k = 2π(k-1)/N`
(`N = n_samples`, default 36), both signs, with `m = m_max`:

- **3-leg (3-wire)** — pairwise line-to-line references must fit the bus, and no
  zero-sequence current flows:

```math
\pm\sqrt2\big[(U^{re}_x-U^{re}_y)\cos\theta_k - (U^{im}_x-U^{im}_y)\sin\theta_k\big]
   \le m\,v_{dc}(\theta_k),\quad (x,y)\in\{ab,bc,ca\};\qquad \textstyle\sum_x I_x = 0.
```

- **4-leg** — the fourth leg is a movable reference, so the pairwise conditions
  hold **plus** each phase against the neutral leg, and the neutral current is
  limited by the 4th-leg rating:

```math
\pm\sqrt2\big[U^{re}_x\cos\theta_k - U^{im}_x\sin\theta_k\big] \le m\,v_{dc}(\theta_k);
\qquad |I_n| = \Big|\textstyle\sum_x I_x\Big| \le I_{n,\max}.
```

- **split-DC (4-wire)** — each phase is an independent half-bridge against the
  capacitor midpoint (half the bus), and the fundamental midpoint ripple
  `N = I_{ret}/(j2\omega C)` (with `I_{ret} = \sum_x I_x`) merges into the phase
  reference `W_x = U_x + N`:

```math
\pm\sqrt2\big[(U^{re}_x+N_{re})\cos\theta_k - (U^{im}_x+N_{im})\sin\theta_k\big]
   \le \tfrac{m}{2}\,v_{dc}(\theta_k);\qquad
   \underbrace{|I_n| \le I_{n,\max}}_{\text{optional — see below}},
```

Here the standalone `I_{n,\max}` bound is **optional**: the split link's neutral
current flows through the capacitors, so supplying `i_cap_max` bounds it through
the thermal budget instead (and then `|I_n| \le 2\,i_{cap,max}`). One of the two
is required. (`:FOUR_LEG` always needs `I_{n,\max}` — its neutral goes through
the fourth leg, which the capacitor budget says nothing about.)

with `N_{re} = -I^{im}_n/(2\omega C)`, `N_{im} = I^{re}_n/(2\omega C)`. The factor
of two on the rail is the split link's utilisation penalty: it needs roughly
twice the DC voltage of the 4-leg for the same per-phase output.

!!! warning "`In_max` for the split link is net of the 2ω allocation"
    The split link's half-banks carry the fundamental neutral current **and** the
    2ω bus current. They are at different frequencies, so they combine in RMS and
    must be allocated *simultaneously* — the bank rating is not all available for
    neutral current:

    ```math
    \Big(\tfrac{|I_n|}{2}\Big)^2 + I_{2\omega,rms}^2 + I_{sw,rms}^2 \le I_{half,rated}^2,
    \qquad I_{2\omega,rms} = \frac{|\tilde S|}{\sqrt2\,V_{dc}}
    ```

    so `In_max = 2·√(I_half_rated² − I_2ω,rms² − I_sw,rms²)` (the factor 2 because
    `|I_n|` divides equally between the two halves). Passing the raw bank rating
    double-counts the capacitors. Example: a 12 A/half bank at `|S̃|` = 6.6 kVA and
    `V_dc` = 800 V gives `I_2ω,rms` = 5.8 A, hence `In_max ≈ 21 A`, not 24 A.

    The exact statement is a frequency-weighted thermal limit,
    `Σ_k ESR(f_k)·I_k² ≤ P_diss,max`. The equal-weight sum used here is an
    *unweighted* approximation of it, and is conservative **only** when
    `i_cap_max` is referred to the lowest frequency appearing in the sum —
    the 50/60 Hz neutral term — because electrolytic ESR *falls* with
    frequency, so that is the worst case. In practice: take the 100/120 Hz
    datasheet rating and apply the manufacturer's frequency multiplier
    (≈0.7–0.85 at 50 Hz) before passing it in. Passing a raw 100/120 Hz
    rating is **non-conservative** for the 50/60 Hz neutral current.

### Capacitor ripple current: endogenous allocation (`i_cap_max`)

Rather than pre-computing `In_max` by hand, supply the bank's RMS ripple-current
rating as `i_cap_max` (**per half-bank** for `:SPLIT_DC`, the whole DC-link bank
otherwise) and let the solve make the allocation at the operating point:

```math
\underbrace{\Big(\tfrac{|I_n|}{2}\Big)^2}_{\texttt{:SPLIT\_DC}\ \text{only}}
 + \; I_{2\omega,rms}^2 \; + \; i_{sw}^2 \;\le\; i_{cap,max}^2 ,
\qquad I_{2\omega,rms} = \frac{|\tilde S|}{\sqrt2\,V_{dc}} = k\,|D|,
\quad k = \frac{2\omega C_{eq}}{\sqrt2}
```

Note `I_{2ω,rms} = |\tilde S|/(\sqrt2 V_{dc})` is independent of `C`: capacitance
sets the ripple *voltage*, not the ripple *current*. `result.i_cap` reports the
whole left-hand side (square-rooted), `i_sw` included, so it compares directly
against `i_cap_max`.

The neutral term appears **only** for the split link, whose half-banks sit in the
neutral path; the 4-leg's neutral current flows through its fourth leg, so its
bank carries the 2ω component alone. Writing `I_{2ω,rms}` through the ripple
phasor `D` (rather than the bilinear `\tilde S`) keeps the constraint quadratic
in variables the model already has. `i_sw` optionally reserves a constant
switching-frequency allowance out of the budget, for designs where the
electrolytics — not a parallel film capacitor — carry the `f_sw` component.

`i_cap_max` composes with `In_max`; whichever binds, binds. For `:SPLIT_DC` it
may also be supplied *instead of* `In_max` — the bank rating bounds `|I_n|` on
its own. (`:FOUR_LEG` always needs `In_max`: its neutral current flows through
the fourth leg, whose device rating is unrelated to the capacitors.) The solved
bank current is reported as `result.i_cap`, so you can see which limit governed:

```julia
r = solve_advanced_inverter(net, AdvancedInverter(; id="inv", bus="poc",
        phase_terminals=["a","b","c"], neutral="n", s_max=20e3,
        topology=:SPLIT_DC, v_dc=800.0, c_dc=2.8e-3, i_cap_max=12.0))
r.i_cap        # capacitor RMS ripple current (A) vs the 12 A rating
r.i_neutral    # what was left for neutral current after the 2ω share
```

This is the model-side statement of the sizing rule above: tightening
`i_cap_max` on a split link visibly trades away neutral capability, while on a
4-leg it only reduces the bus ripple.

**Current limits** are per-phase `|I_x| \le i_{max}` and the neutral limits above.
The sampling makes these **outer** approximations, exact as `N→∞` with relative
gap `~(π/N)²` (≈0.8 % at N=36). Every per-sample constraint is linear in the
phasor variables and the ripple aux `D`; only the two `\tilde S` equalities are
nonlinear, so the model solves as a smooth NLP (Ipopt).
