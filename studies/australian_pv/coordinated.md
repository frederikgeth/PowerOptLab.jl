# Coordinated three-leg capability benchmark

> **Historical campaign.** These results describe a pinned earlier model revision. See the [current benchmark](engineering_benchmark.md) and [evidence/replay index](evidence/README.md).

This benchmark maximizes total net POC active injection using the existing
Australian plants and DC-link availability. It releases the local VV+VW laws
and optimizes Q jointly with P. It is a feasible local nonconvex OPF reference,
not a certified global optimum or upper bound. Hardware sizing, converter-arm
current limits, POC apparent-power limits, continuous rails, DC-ripple and
capacitor thermal budgets remain unchanged.

Two current-allocation variants separate different capabilities:

- `balanced_converter` retains zero and negative-sequence converter current at
  zero, matching the current subspace of the three autonomous sensing laws.
- `free` permits complex I₂ within the same three-leg plant limits. I₀ remains
  zero. This releases another control degree of freedom, so its extra export
  must not be attributed solely to coordination or to a changed sensing law.

Each variant is solved with and without an additional 253 V upper
phase-to-neutral constraint at every customer. This is an illustrative planning
criterion, not regional connection conformity. The unrestricted variant removes
only this additional criterion; hardware and native network constraints remain.
No lower phase-to-neutral planning bound is added. Autonomous controller points
must be checked against the additional criterion before comparing their feasible
export with a constrained benchmark. A larger export violating 253 V is not a
feasible improvement under that criterion.

The 24 predeclared cases span the same 0/50/100% available PV and 1/1.03 pu
positive-sequence source conditions as the matched study, with 1% aligned
upstream VUF and the original 60/25/15% phase-load shares. All plants remain
energized, including the zero-PV cases. Q is an optimization variable rather
than a reactive-priority service request. The benchmark does not enforce VV+VW
grid-code behavior or establish admissibility of deliberate unbalanced injection.

```sh
julia --project=. scripts/australian_pv/test_coordinated.jl
julia --project=. scripts/australian_pv/coordinated.jl
julia --project=. scripts/australian_pv/audit_coordinated.jl
julia --project=. scripts/australian_pv/test_coordinated_records.jl
python3 scripts/australian_pv/test_coordinated_csv.py
python3 scripts/australian_pv/plot_coordinated.py
```

An optional output-directory argument selects an empty destination. One common
100 kVA base and 1,000-iteration NLP profile apply to every case, with adaptive
continuous-rail cuts and no per-case retries. Every completed failure retains
its raw status and acceptance reason. The independent inverter-policy, native
SI equation and additional customer-voltage gates must pass, together with a
fixed-PQ cross-formulation replay. That replay uses the same network engine.
The saved-point audit reconstructs inverter physical checks, sequence allocation
and customer voltages from serialized phasors, then repeats the fixed-PQ replay.
Source hashes, settings, all physical/network check rows, points and completion
counts are retained in the output directory.

Saved-record regressions reject modified currents, mismatched source hashes,
missing cases/devices, changed CSV metrics and false inequality-success flags.
An audit begins by clearing its prior success flag. Failed validation cannot
retain an accepted audit report. The Python command (Matplotlib/NumPy required
for plotting) generates the comparison table and standalone PNG/SVG figure only
after the CSV records match the independent phasor audit. It also verifies common
source, hardware and terminal settings against the autonomous snapshot campaign.

Standalone regression tests check the representative constrained/unconstrained
points, three numerical bases, input preservation and rejection of a raw-success
point with violated native network equations. The benchmarks do not supply
dynamic stability, switching, communication or EMT evidence.
