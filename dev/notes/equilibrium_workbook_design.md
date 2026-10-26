# Controller-program to equilibrium workbook

Design recorded before implementation of this teaching slice.

Teach the meaning of an equilibrium formulation before solver syntax. Keep
package APIs unchanged. Use a small dimensionless voltage–power feedback loop
whose exact fixed point and branch conditions can be checked analytically.

1. Translate read/compute/update statements into network, controller, hardware,
   state and numerical-method roles. Distinguish requested and delivered power.
2. Solve the same controller equilibrium algebraically and with simultaneous
   JuMP constraints. Show divergent direct iteration, convergent relaxation and
   explicitly declared continuous/discrete physical dynamics. Equilibrium alone
   makes no stability claim.
3. Compare exact clipping graphs with capability inequalities under different
   objectives. Enumerate guarded affine branches for a transparent exact oracle;
   include boundary cases and examples where a missing guard gives a false root.
4. Preserve physical time with a battery schedule and a network equilibrium in
   each period. Replay the optimised schedule sequentially. Demonstrate hysteresis
   history dependence and state the boundary/tie convention explicitly.
5. Explain solved points through equality residuals, signed margins, binding
   constraints and requested/commanded/delivered quantities. Teach the limits of
   these diagnostics and connect back to existing inverter reproductions.

Deliver a workbook with prediction/run/explanation exercises, executable Julia
helpers, independent analytic and replay tests, reproducible scientific figures,
and a small offline browser feedback explorer if it fits the documentation build.
Verify current figures and executable examples in a full Documenter build.

## Completed

Delivered six prediction/run/explanation exercises in
`docs/src/tutorials/controller_program_to_equilibrium.md`, with a line-by-line
program translation, simultaneous JuMP solve, guarded clip-branch oracle,
dispatch counterexamples, iteration/stability distinctions, battery replay,
hysteresis history and an interpretable constraint report. The workbook is
linked before the existing inverter geometry tutorial in the navigation.

The same clipped feedback loop has a unique analytical equilibrium but a
cycling direct numerical iteration. Its capability formulation reproduces
the controller only under the demonstrated selection assumptions; an extra
voltage ceiling provides an explicit non-equivalence example. Declared sampled
and continuous physical dynamics are discussed separately from solver iteration.

Added three reproducible scientific SVG/PDF figures and an offline feedback
explorer. The explorer was browser-tested at default, zero and maximum feedback
gain, and with full numerical relaxation; displayed roots and local multipliers
matched the formulas and no browser console errors were observed. All figures
were visually inspected, with legends moved clear of the plotted relationships.

Validation:

- Final workbook suite: 213 passing checks, including branch knots, zero
  availability/gain/resistance, deduplication, malformed inputs, objective
  counterexamples, missing equilibria, sequential energy conservation, failed
  battery-schedule masking and hysteresis thresholds/history.
- Broader inverter/workbook regression run: 3,676 passing checks and one optional
  OpenDSS skip. It included the first 209 workbook checks; the final four added
  infeasible-battery checks were verified in the focused suite above.
- Full Documenter build passed after the last text edits. All examples executed;
  rendered diagnostic labels, local links and byte-identical copies of the
  explorer plus six SVG/PDF exports were verified.
- Source/generated-SVG whitespace checks passed. Existing documentation API-page
  and search-index size warnings remain.

These are transparent teaching utilities, not new package APIs. The scalar
network, ideal battery and branch enumeration are explicitly scoped; production
inverter, battery and network models remain linked as the next learning step.
