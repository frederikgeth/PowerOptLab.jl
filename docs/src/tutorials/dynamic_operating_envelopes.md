# Dynamic operating envelopes for LV networks: a modelling tutorial

> **Audience:** power-system researchers · **Scope:** nonlinear AC DOE studies
> with DSSE snapshots, mandatory DER controls, and network support devices.

This tutorial is about *modelling choices*, not just calling an optimizer. A
dynamic operating envelope (DOE) is an operational promise: a participant may
operate within an advertised active-power range while the LV network remains
within its declared limits. The difficult part is deciding exactly what that
promise means and ensuring the model represents the physical controls that will
be present when it is used.

The runnable end-to-end example is
`scripts/cases/doe_dsse_validation_demo.jl`:

```sh
julia --project=. scripts/validate_doe_from_dsse.jl \
    scripts/cases/doe_dsse_validation_demo.jl
```

It performs DSSE reconstruction, DOE issuance, fixed-capacity verification, an
independent AC power-flow replay, and a with/without-STATCOM comparison. This
tutorial explains how to turn that pattern into a research study.

## 1. Begin with the operational question

State the quantity before choosing a formulation. For export, participant `i`
receives a non-negative active-power limit ``\bar p_i`` such that its injection
is at most ``\bar p_i``. For import, the same positive quantity bounds
withdrawal. These are two different operating studies:

```julia
export = solve_operating_envelope(net, cps; direction=:export)
import = solve_operating_envelope(net, cps; direction=:import)
```

Do not infer an import envelope by negating an export result. In an LV feeder,
the baseline P/Q, voltage profile, converter limits, and controls can make the
two highly asymmetric.

The input `net` must be a *snapshot*: topology, source conditions, DSSE-derived
loads, DER availability, and network limits appropriate to one issuance
interval. A `Vector{Dict}` gives several intervals; each is solved separately
unless rolling fairness is selected.

## 2. Preserve the distinction between estimated injections and controls

The ordinary DSSE load is known P/Q data. It belongs in the network's `load`
block and stays fixed during DOE optimization. A PV inverter or battery is not
a free `(P,Q)` source merely because it can exchange reactive power. Bind a
connection point to its existing BMOPFTools IBR instead:

```julia
pv = ConnectionPoint(id="pv_17", bus="lv17", ibr_id="pv17",
                     export_max=10e3, import_max=5e3)
r = solve_operating_envelope(net, [pv]; direction=:export)
```

The DOE constrains only that IBR's aggregate active power. The network model
continues to enforce its topology, current and apparent-power limits, DC-side
coupling, and mandatory Volt-VAr/Volt-Watt or fixed-power-factor control law.

### Pitfall: independent P and Q envelopes

Giving every participant a freely dispatchable `Q` range is usually physically
wrong for LV operation. It can create a generous-looking P envelope by using
reactive support that the installed PV/battery controller will never provide.
It also obscures who is responsible for voltage support. Use a separate
STATCOM, D-STATCOM, or another explicitly modelled network support device if
reactive dispatch is genuinely available to the operator.

The lightweight `ConnectionPoint` without `ibr_id` is useful for teaching or an
aggregate unity-power-factor port. It is not the preferred representation for a
real controlled inverter.

## 3. Choose the security set deliberately

An allocation does not define its own guarantee. The `security` keyword states
which participant dispatches are embedded in the nonlinear AC optimization:

| Choice | Included dispatches | Suitable interpretation | Principal limitation |
|---|---|---|---|
| `:bound_point` | all participants at their bound | a simultaneous full-export/import operating point | does not establish interior range safety |
| `:corners` | every zero/full corner of the advertised box | a small-N box-range test | `2^N` contexts; still not a global AC robust certificate |

```julia
r_bound = solve_operating_envelope(net, cps; security=:bound_point)
r_box   = solve_operating_envelope(net, cps; security=:corners)
```

For ``N`` participants, corners require ``2^N`` AC contexts per scenario, so
they are intentionally capped (`max_exact_corners=10` by default). Report the
mode together with the capacity. The result diagnostics make this explicit:

```julia
r_box.diagnostics[1]["security_scope"]
r_box.diagnostics[1]["binding_constraints"]
```

### Pitfall: calling an upper-bound point an operating range

`security=:bound_point` is often useful, but it means exactly one joint point
was tested. If an operator tells each participant they can independently choose
anything in `[0, envelope[i]]`, use corners for small studies or develop a
screened/adaptive security set for larger ones. Even corners do not prove that a
non-convex AC feasible set contains every interior point.

## 4. Treat uncertainty as a shared-allocation problem

Forecast and model uncertainty should normally produce one conservative DOE per
interval, not a different DOE chosen after uncertainty resolves. Supply one or
more scenarios per interval:

```julia
scenarios = [[net_central, net_high_load, net_low_source],
             [net_next_central, net_next_high_load, net_next_low_source]]

r = solve_operating_envelope(scenarios, cps; security=:corners)
```

Each interval has a single allocation shared by all of its scenarios. Scenario
sets can represent P/Q forecast error, source-voltage uncertainty, switch
status, or credible feeder-parameter alternatives. If line construction is
uncertain, materialized candidates from [`solve_inverse_carson`](@ref) are a
natural source of network scenarios.

### Pitfall: mixing probability and feasibility claims

This is a scenario-feasibility method, not a chance-constrained guarantee. A
three-scenario set has no implied confidence level. State how scenarios were
generated, whether they are stress cases or samples, and which unmodelled
events remain outside the DOE promise.

## 5. Make fairness units explicit

The objective determines who receives scarce network headroom. Equal kW,
equal fraction of nameplate, and equal fraction of a request are different
policies. Use `FairnessPolicy` rather than relying on an informal description:

```julia
policy = FairnessPolicy(kind=:max_min,
    normalization=:capacity,
    weights=Dict("pv_17" => 1.0, "battery_42" => 1.25))
r = solve_operating_envelope(net, cps; fairness=policy)
```

Useful choices include `:max_total` (efficiency), `:equal` (equal normalized
allocation), `:proportional`, `:alpha`, `:max_min`, and
`:equal_curtailment`. Normalization is `:none`, `:capacity`, `:request`, or
`:custom`. Compare policies rather than claiming one is intrinsically fair:

```julia
frontier = compare_operating_envelope_policies(net, cps,
    ["equal kW" => :equal,
     "efficient" => :sum,
     "proportional capacity" => FairnessPolicy(kind=:proportional,
                                                  normalization=:capacity)])

frontier["efficient"].fairness_metrics[1]
```

The metrics include total capacity, normalized allocations, curtailment
fractions, and Jain's index. They describe the published allocation; they do
not select a social-welfare objective on their own.

### Pitfall: unnormalised fairness silently favours one population

Equal W favours neither participant in an electrical sense, but can be unfair
when systems have very different nameplates or requests. Conversely,
capacity-normalized fairness can allocate substantially more kW to a large
system. Publish the reference and weights alongside every fairness result.

## 6. Separate customer DER behaviour from network support

A STATCOM is represented in the network, not as a customer Q envelope. Compare
otherwise-identical cases:

```julia
using BMOPFTools: add_statcom!, augment_case

base = deepcopy(net)
supported = deepcopy(net)
add_statcom!(supported, "lv17"; s_max=50e3)
supported, _ = augment_case(supported)

r_base = solve_operating_envelope(base, cps)
r_supported = solve_operating_envelope(supported, cps)
gain_W = r_supported.total_capacity[1] - r_base.total_capacity[1]
```

Inspect the STATCOM dispatch in `r_supported.snapshots` and the worst margins
in its diagnostics. A capacity gain is meaningful only with the support-device
rating, location, control policy, source conditions, and binding constraint
reported.

### Pitfall: comparing different cases

Do not change baseline load, source voltage, topology, or DER availability when
adding the STATCOM. Otherwise the difference is not attributable to support.
Also distinguish an *available rating* from a control policy: a free OPF
STATCOM can look better than a deployed device with a fixed or local controller.

## 7. Issue, verify, and audit

The result contains more than a number. Attach issuance metadata and choose an
explicit publication behaviour for failed intervals:

```julia
using Dates

r = solve_operating_envelope(nets, cps;
    issued_at=DateTime(2026, 7, 16, 9),
    interval_seconds=300,
    validity_seconds=600,
    fallback=:missing)

r.schedule[1]
```

`fallback=:zero` and `:last_feasible` are publication policies, not fresh
network-security assertions. A last feasible DOE can be unsafe after topology,
load, or source conditions change.

Before an operational study is trusted, fix the published allocation and solve
again at the intended utilization points:

```julia
check = verify_operating_envelope(nets, cps, r; utilizations=:corners)
all(check.feasible)
```

For a DSSE-to-DOE study, use the standalone validation runner. It compares the
DSSE voltage reconstruction against a truth power flow, checks the issued DOE,
and replays fixed active-power setpoints in a separate AC power flow. The
independent replay currently targets single-phase IBR-bound points; multi-phase
studies must preserve the DOE's phase allocation when replaying it.

### Pitfall: silently publishing infeasible solver iterates

An NLP termination status alone is not an envelope. This implementation returns
`NaN` capacities for intervals without a feasible primal point. Treat that as
an issuance failure requiring an explicit fallback and audit trail, never as a
zero-confidence numerical nuisance to be filled with the last optimizer values.

## 8. A minimum publication checklist

For each research result, record:

1. network model: four-wire/neutral representation, limits, topology, source,
   and parameter assumptions;
2. snapshot provenance: DSSE data, timestamp, load/DER P/Q assumptions, and
   scenario construction;
3. DER and support controls: mandatory Q–V/Q–P laws, STATCOM rating/location,
   and whether its dispatch is free, fixed, or controller-replayed;
4. security semantics: bound point, corners, or another tested set;
5. fairness objective, normalization, weights, and reported fairness metrics;
6. nonlinear solver status, binding margins, and independent replay/validation;
7. issuance interval, validity window, and fallback policy.

With those choices stated, a DOE becomes reproducible and interpretable: not
just an attractive number from an OPF, but a precisely scoped operational claim.

For the full API and field definitions, see [Operating envelopes](../problems/operating_envelope.md).
