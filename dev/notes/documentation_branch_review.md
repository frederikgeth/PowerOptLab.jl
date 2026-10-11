# Documentation branch review

Date: 2026-10-11. Branch: `codex/reorganize-component-docs`.
This validation record belongs in developer notes, outside the published docs.

## Base and review scope

Merged main commit `5f6c396` before validation, retaining its staged
`optimize_opf!` integration and BMOPFTools revision
`56e3b2f0261488acd30f1c98ad812c4eff4611f1`. Reviewed the complete branch diff:
component-family and study-topic navigation, API moves, capability corrections,
workflow-note relocation, tutorials, generated figures and CI changes.

The 330 canonical `@docs` entries match main exactly, without duplicate owners.
All 108 published pages are registered once. Guide paths remain stable; canonical
API entries now live in topic references. The relocated workflow records and
extracted sections remain in `dev/notes/`; the published operability/DOE
references retain model equations, scaling, continuation and validation guidance.
Changes under `src/` relative to main are comments/docstrings only; runtime
statements outside documentation are unchanged.

## Findings fixed in this pass

- Reconciled the current upstream finalization API and source pin with the
  reorganized documentation after merging main.
- Corrected tabulated-chemistry bounds: the constructor permits and defaults to
  the outer knots. Strictly interior bounds are a numerical recommendation, not
  enforced constructor behavior. Clamping can introduce an endpoint derivative
  discontinuity; the interior PCHIP C² limitation remains.
- Fixed a repeated word in the moved neutral-reduction reference and identified
  README examples as interface snippets with caller-supplied networks.
- Fixed Julia 1.10 Markdown parsing of bracketed reactive-power units followed by
  a parenthetical phrase, which had become an invalid local link.
- Fixed fresh documentation installation on Julia 1.10. Its Pkg ignores
  `[sources]`, so a temporary environment now seeds the pinned upstream package
  and the local PowerOptLab path before resolving the docs dependencies. Existing
  manifests retain the normal installation path.
- Added documentation CI on Julia 1.10 and current stable. A deployment flag
  ensures only current stable deploys; both versions run Documenter and examples.
- Marked generated study SVGs and the plotting-data snapshot for collapsed
  GitHub review display.

## Clean-environment validation

Used separate full Git snapshots with no initial manifests. Each environment
resolved dependencies from the declared source pin, using a separate writable
Julia depot with the normal download cache as a fallback. This is a fresh
resolution/build test, not a claim of a network-free or cold-cache installation.
A separate fresh Python 3.12 virtual environment installed the plotting pin.
The user's unrelated working-tree results file was not copied into the snapshots.

| Check | Julia 1.10.11 | Julia 1.12.6 |
|:--|:--|:--|
| Package quality | 10 passed | 11 passed, including Aqua's persistent-task probe |
| Full scientific suite | 11,025 passed | 11,025 passed |
| Exact published study tutorials | 7 passed | 7 passed |
| Exact published DOE tutorials | 8 passed | 8 passed |
| Fresh documentation installation | passed after bootstrap fix | passed |
| Documenter build, examples and exported API coverage | passed | passed |

The fresh plotting environment passed the three-figure input/fingerprint and
rendering check. Navigation, local file/image targets and rendered section
anchors were checked separately. The first-run and generator examples execute
through Documenter; tutorial runners do not substitute copied example code.

## Remaining diagnostics and scope

The build retains the existing DOE API page-size and search-index warnings.
They are below Documenter's fatal limits. Stable testing also reports an upstream
OpenDSSDirect precompilation/name-shadowing issue and then loads it without that
precompile cache; the full test suite completes successfully. This docs PR does
not patch or silently replace the external test dependency.

The five schematic tutorial guides remain explicitly labelled rather than being
presented as runnable studies. See `documentation_tutorial_audit.md` for their
coverage and the scientific inputs needed for subsequent worked examples.

## Ubuntu CI follow-up

The first Ubuntu runs on Julia 1.10.12 and 1.13.1 passed the full package and
DOE suites but failed the AC-harmonics tutorial's raw-SI maximum-export solves.
One failed during the L/LCL comparison and the other during the neutral-inductor
comparison. These examples had succeeded on the local macOS environments;
that was insufficient evidence of platform robustness.

The tutorial now uses per-unit coordinates on a 20 kVA base and holds AC P/Q at
15 kW / 0 var. A declared synthetic quadratic converter-loss coefficient selects
a low-current point. Every filter and sampling refinement uses the same setup.
Publishable status, P/Q accuracy, total-RMS ratings, reserve coverage and the
zero-capacitor equivalence checks are retained. This fixes the numerical setup
and the comparison's operating-point assumptions rather than accepting failed
candidates. Julia 1.13's SHA 1.0 stdlib is also admitted by package compatibility.

The next current-stable run stopped earlier in the battery tutorial with an
`ALMOST_LOCALLY_SOLVED` exit. Its 2 kW devices were still using a 1 MVA AC base.
The battery study now sets its AC base to the 2 kVA rating and disables
acceptable-iteration exits. All existing power, energy, SoC, cell-law and
publishability assertions remain. Fresh local figure inputs match the committed
snapshot within its existing physical budgets. Figure CI now depends on a
successful study export so missing inputs do not obscure a primary solve failure.
