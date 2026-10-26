# Inverter operating constraints: final PR review

Branch: `codex/ibr-constraint-control`, based on `main` at `0955b61`.
The remote was refreshed during review; no newer main commits were present.

## Review findings resolved

1. Staged validation passed per-unit bus voltage bounds into a diagnostic that
   expects SI. A physically significant LCL midpoint shunt (about 2.17 A) was
   therefore reported as negligible. Validation now restores SI bounds in a
   shallow diagnostic copy. The electrical model and its working network remain
   unchanged. Seven regression checks cover the standalone solve, direct staged
   composition, preservation of per-unit data and a truly negligible shunt.
2. The full-network documentation example directly imports `JSON3`, but the docs
   project did not declare it. An isolated documentation build exposed the
   missing dependency. `docs/Project.toml` now declares it with the same minimum
   compatibility as the package.

No unresolved blocking findings remain in the reviewed formulation, public
composition/extraction API, failure masking, examples or documentation.

## Verification

- Before these review fixes, the complete test harness passed on Julia 1.10.12
  and 1.12.6: 18,192 tests each, plus 10 and 11 Aqua package-quality checks
  respectively. Both runs included the declared OpenDSS oracle dependency.
- With the validation fix, the complete Julia 1.12.6 harness passed 18,199 tests
  and all 11 package-quality checks. The affected inverter test file passed all
  192 checks on both Julia versions, including the seven new regression checks.
- The full Documenter build passed on Julia 1.10.12 and 1.12.6 with an isolated
  `JULIA_LOAD_PATH=@:@stdlib`, including the executable network composition,
  independent replay and controller-to-equilibrium examples. Existing size
  warnings for the DOE API page and search index remain.
- All 844 local content links/assets across the six new documentation pages
  were present in the rendered output. Documenter's deployment-only
  `siteinfo.js` / `versions.js` references were excluded from this local check.
- Additional probes checked zero phase voltages and impossible Q requests at
  widely separated power bases; failed solves retained non-publishable status
  and masked numerical results.
- `git diff --check` passed. Linux CI remains to run with the PR.

Tests used the full `test/runtests.jl` harness and local environments with the
reviewed BMOPFTools source revision from `Project.toml`; this was not a fresh
network dependency resolution or a Linux run.

## Scope of readiness

This is ready for review as a nonconvex, fundamental-frequency operating-point
formulation, alongside the retained forward controller. The documentation
distinguishes controller graphs, physical capabilities and objective-selected
dispatch, and limits the agreement claims to matching operating assumptions.
The seven-bus validation includes three inverters, unbalanced four-wire loads,
floating neutrals, mutually coupled lines, losses, binding service constraints
and an independent SI radial power-flow replay.

The documented exclusions remain deliberate: grid-forming operating laws,
single-phase constrained wrappers, automatic PWM reserve closure, exact PWL
inverter command graphs, dedicated constrained-fleet/time/DOE wrappers and
dynamic or standards-compliance certification. No global optimum, global
infeasibility, uniqueness or dynamic-stability claim follows from the tests.
