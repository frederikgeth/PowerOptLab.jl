# Build your own formulation experiment

Start here to use the toolkit. Read [the mathematical models](theory.md) next,
then [physical error and conditioning](error_budgets.md). The
[controller tutorial](controllers.md) connects these ideas to the existing
three-phase inverter. [Extension methods](extensions.md) let you add a family
without changing PowerOptLab or BMOPFTools.

The package supplies representations, construction hooks, diagnostics and export.
You choose the scientific question, cases, solver settings and acceptance criteria.
An experiment run is one model build and one call to an external optimizer.
External solvers may perform their own homotopy; PowerOptLab adds no retry loop.

## A first experiment, in physical units

Consider an inverter exporting real current through a resistor. Define the
canonical current law in amperes versus voltage in volts. The source voltage and
resistance are independent case parameters, not smoothing parameters.

```@example experiment
using PowerOptLab, JuMP, Ipopt, TOML
curve = PWLFunction([220.,240.,250.,270.], [20.,20.,0.,0.];
    input_unit=:V,output_unit=:A)
case = resistive_control_case(curve;source_voltage=230.,resistance=1.)
reference = resistive_equilibria(curve,230.,1.)
@assert isapprox(only(reference.points).current,40/3)
reference
```

Each affine segment is intersected with `V = Vs + R*I` analytically. The returned
record contains isolated `points` and any equilibrium `intervals`. A unique root
is a property of this example, not an assumption imposed on your case. Endpoint
inclusion uses the caller-adjustable absolute voltage tolerance `atol`. Nearly
parallel intersections can be ill-conditioned; the enumerator does not perform
interval-arithmetic certification.

Choose methods at the same **physical output-error budget**. A representation
callback can read configuration fields and return a different width on each run.
The optimizer factory and its options remain under your control.

```@example experiment
methods = [FormulationMethod(string(family),
    (f,c) -> smoothing_for_error(f,family,c.error_A),
    Ipopt.Optimizer;options=(tol=1e-9,max_iter=300),configure! = set_silent,
    metadata=(family=string(family),optimizer="Ipopt"))
    for family in (SoftplusFormulation,LocalC2Formulation)]
configs = [(error_A=e,input_scale=230.,output_scale=20.,source_voltage=vs,
           start_input=245.) for e in (.01,.1) for vs in (230.,240.)]
rows = run_formulation_experiment([case],methods;configurations=configs,
    on_error=:throw)
@assert length(rows)==8
[(r["method"],r["configuration"].error_A,r["termination_status"],
  only(r["observations"]).exact_graph_error) for r in rows]
```

Here the model coordinate `v` represents `230*v` volts; `h.output` represents
`20*h.output` amperes. Widths and recorded residuals remain physical. Vary
`input_scale`, `output_scale`, `start_input`, `resistance`, or `objective=:zero`
when those changes answer your question. Coordinate changes preserve the intended
mathematical equations but can affect floating-point solver behavior.

## Graph and relation observations

Return graph handles, relation handles, or both in `observations`. Every handle
must belong to the returned model. The runner records canonical curve data and
units, physical coordinate scales, the domain and candidate residuals.
`observation_kind` distinguishes graph and relation records. Relation records
separate `requested_formulation` from `built_formulation`: a requested hull may
build only `:linear_inequalities`. Use `reason_code` for filtering and `reason`
for explanatory text.
`conservative` and `output_shift` identify any contract-based one-sided correction.

Graph records retain `formulation_type`, `exact_graph_error` and their formulation
`contract`.
Relation records instead include `relation`, `strategy`, `semantics`, a signed
`canonical_residual` and a nonnegative `canonical_violation`. Negative residual
for an upper limit is valid slack. A smooth relation also records its
`approximation_contract`, including the domain-specific C2 bounds or global
bounds used for a conservative shift (`error_scope`); an auxiliary `graph_audit` is separate from the relation
check. In particular, an exact projected hull limit must not be labeled a relaxed
relation merely because its auxiliary can lie off the curve. See the
[bounded volt-watt experiment](bounds_and_relations.md) for executable examples.
Use custom `metrics` for circuit, plant or study-specific quantities.

## Assess evidence without rewriting solver status

`strict_solver_success` requires `OPTIMAL` or `LOCALLY_SOLVED` and a feasible
primal status. `candidate_available` can be true after a time/iteration limit or
an almost-solved termination. Candidate audits are useful in either case.
Primal infeasibility/reduction certificates are flagged as certificates and are
not evaluated as operating points.
`run_status="finished"` means that the call and diagnostics finished, not that
the solver succeeded or the physical equations were validated.
`solver_detail_errors` records failures of optional raw-status/dual-status getters
in external backends. Such failures leave the primary termination/primal outcome
and candidate audits intact; unavailable details are exported explicitly.

A researcher may explicitly accept a candidate under a different policy. Keep
that policy and its result separate from the original solver report:

```@example experiment
function assess_current_example(row)
    get(row,"candidate_available",false) || return (accepted=false,reason="no candidate")
    a = only(row["observations"])
    numerical = abs(a.surrogate_equation_error)<=1e-6 &&
        abs(row["metrics"].electrical_residual_V)<=1e-6 && a.domain_violation<=1e-6
    # The graph error includes approximation error as well as numerical error.
    physical_budget = row["configuration"].error_A
    return (accepted=row["strict_solver_success"] && numerical &&
        abs(a.exact_graph_error)<=physical_budget+1e-6,
        policy="example current-law budget plus electrical/domain residuals")
end
assessed = run_formulation_experiment([case],methods;configurations=configs,
    assess=assess_current_example,on_error=:throw)
@assert all(r -> haskey(r,"assessment"),assessed)
path = write_formulation_results(tempname()*".toml",assessed;
    metadata=(purpose="tutorial illustration",acceptance="assess_current_example"))
@assert TOML.parsefile(path)["schema_version"]==1
nothing # hide
```

For a publishable study, supply `sources=[...]` with your configuration script,
case files and resolved manifest. The exporter stores hashes of those files,
core package/runtime versions, options, configurations and metrics. It does not
archive those files, discover arbitrary callback dependencies, or reconstruct a
model from TOML. Preserve the files and environment alongside the bundle. Custom
metrics must be plain data; model and solver objects are rejected. A method's
`metadata` should identify its solver and any custom representation parameters
that are not already explicit in the configuration or observation contract.

## Supply a new case

A case can be a scalar relation, inverter, feeder or DOE. Its build callback
returns a fresh JuMP model, optional PWL handles to audit, and a zero-argument
metric callback evaluated only when a candidate exists:

```@example experiment
custom = FormulationCase("fixed-voltage",(rep,c) -> begin
    m = Model()
    @variable(m,v)
    representation = rep(curve,c)
    h = formulate_pwl!(m,curve,v,representation)
    @constraint(m,v==c.voltage)
    @objective(m,Min,h.output^2)
    (model=m,observations=[h],metrics=()->(current_A=value(h.output),))
end;metadata=(purpose="inspect one breakpoint",))
inspection = run_formulation_experiment([custom],methods;
    configurations=[(voltage=240.,error_A=.01)],solve=false,on_error=:throw)
@assert all(r -> r["run_status"]=="built",inspection)
```

Use `UnsupportedFormulation("reason")` for combinations your case intentionally
does not implement. Other errors record their stage and message under the default
`on_error=:record`; use `:throw` while developing. `on_result(row)` permits live
logging or incremental export. Configurations should be treated as read-only.
Build and solve times include first-use compilation: warm-up, repeat allocation,
randomization and performance statistics are choices for a study, not automatic
claims made by this runner.

For exact graph or MPCC methods, use the optional environment described in
[the backend comparison](theory.md). `FormulationMethod(...; configure! =
MathOptComplements.Bridges.add_all_bridges)` can install the CCOpt bridges before
attaching `CCOpt.Optimizer`. Use `ComplementarityGraph(scale=...)` to choose hinge
normalization independently of voltage coordinates. Preserve physical graph and
complementarity residuals; an MOI status alone is not an MPCC stationarity proof.
