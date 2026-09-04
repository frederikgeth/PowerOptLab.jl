# Current-magnitude measurements

> **Audience:** anyone who has been told that ampere readings are numerically
> dangerous in state estimation.

Current-magnitude (`:imag`) measurements have a bad reputation, usually
summarised as "they make the problem ill-conditioned". That summary is too
coarse to act on. The defensible statement is:

> **Current magnitude does not inevitably cause ill-conditioning. Conditioning
> is operating-point and measurement-set dependent — and there are operating
> points, including very ordinary ones, where it fails badly.**

This page separates the three distinct things that get conflated: a benign
path where conditioning stays bounded, a genuinely degenerate one that arises
at unity power factor, and a non-uniqueness that is not a conditioning
question at all.

Everything below runs against one two-bus feeder, and every table on this page
is pinned by a test in `test/constrained_state_estimation_tests.jl`:

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf

net = parse_bmopf("""
{"bus":{"src":{"terminal_names":["1"]},"b":{"terminal_names":["1"]}},
 "voltage_source":{"s":{"bus":"src","terminal_map":["1"],
     "v_magnitude":[230.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.5}},
 "line":{"l":{"bus_from":"src","bus_to":"b","terminal_map_from":["1"],
     "terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
"""; from_string=true)

vtrue = 220.0 * cis(-0.05)          # the state we will try to recover
itrue = (230.0 - vtrue) / 0.5       # |I| = 30.0981 A
```

The state is `[Vre_b, Vim_b]`. Unless stated otherwise every measurement below
carries `sigma = 0.01` in its own SI unit (volts for `:vmag`, amperes for
`:imag`), so the whitened Jacobian rows are directly comparable.

## Shrinking current alone is benign

Take ``|V_b|`` and ``|I_\ell|`` and let the line current fall by six orders of
magnitude, holding its source-referenced direction fixed so the state stays
genuinely complex:

| ``|I|`` (A) | ``\sigma_{\min}(H)`` | ``\mathrm{cond}(H)`` |
|---:|---:|---:|
| 60 | 4.50e+01 | 4.87 |
| 20 | 4.13e+01 | 5.33 |
| 5 | 4.00e+01 | 5.50 |
| 1 | 3.97e+01 | 5.54 |
| 0.01 | 3.96e+01 | 5.56 |
| 0.0001 | 3.96e+01 | 5.56 |

Along this path the condition number does not blow up — it *converges*, to
about 5.6. With respect to rectangular current components, the gradient of
``\sqrt{I_{re}^2+I_{im}^2}`` is a unit vector along the current direction. The
voltage-state row also contains the fixed branch-current sensitivity; neither
factor shrinks merely because ``|I|`` does along this path.

**This is one path only.** It fixes the current's direction and varies its
magnitude. The next section varies the direction instead, and the outcome is
completely different — which is exactly why the blanket claim in either
direction is unhelpful.

## Real problem 1: sensitivity rows can collapse onto each other

``|V_b|`` senses the voltage state along ``\hat V``. A current-magnitude row
senses it along ``J_I^\top\hat I``, where ``J_I`` maps voltage-state changes to
branch-current changes. On this page's scalar, purely resistive line,
``J_I`` is just a signed scalar: when current and receiving-end voltage are in
phase, the two whitened rows become **parallel**, and the pair loses rank at any
nonzero current magnitude.

For a scalar complex impedance ``Z``, the alignment condition instead involves
the impedance angle: ``\angle V-\angle I = \angle Z`` modulo ``\pi``. With
mutual coupling it depends on the full current-sensitivity row. Unity power
factor is therefore the degeneracy of this resistive example, not a universal
threshold for ampere measurements.

Hold the current at 40 A and sweep its angle relative to the 0-rad source toward
zero. The next column reports the actual receiving-end voltage/current angle:

| source-referenced ``\angle I`` (rad) | ``|\angle V-\angle I|`` | ``\sigma_{\min}(H)`` | ``\mathrm{cond}(H)`` | sd of ``\angle V_b`` |
|---:|---:|---:|---:|---:|
| −0.6 | 0.6529 | 5.61e+01 | 3.86 | 7.2e−05 |
| −0.3 | 0.3280 | 2.91e+01 | 7.63 | 1.6e−04 |
| −0.1 | 0.1095 | 9.78e+00 | 22.8 | 4.8e−04 |
| −0.03 | 0.0329 | 2.94e+00 | 76.1 | 1.6e−03 |
| −0.01 | 0.0110 | 9.80e−01 | 228 | 4.9e−03 |
| −0.003 | 0.0033 | 2.94e−01 | 761 | 1.6e−02 |
| **0** | **0** | **0** | ``\infty`` | undefined |

At exactly unity power factor *on this resistive feeder* the set is rank 1 of 2 while carrying 40 A, and
`selected_state_covariance` refuses to return a covariance. Unity and
near-unity power factor is not exotic — resistive load, and PV inverters
operating at unity — so this degeneracy is reachable in ordinary operation.

The contrast with a rectangular pair is total. Replacing ``|V_b|`` with
``V_{re}``, ``V_{im}`` gives ``\mathrm{cond}(H) = 1.0`` at **every** power-factor
angle in the sweep, including zero.

So: conditioning with ampere measurements depends on the network sensitivity,
measurement placement, whitening, and operating point. It is fine along some
paths and singular along others, and only the diagnostic for the actual model
and operating point can tell you which one you are on.

## Real problem 2: the derivative is undefined at exactly zero current

``\partial |I| / \partial x`` does not exist at ``I = 0``. Not "is large" —
does not exist. And the difficulty is a **single point**, not a
neighbourhood: at ``|I| = 10^{-4}\,\mathrm{A}`` the Jacobian is perfectly
healthy, as the table above shows.

The trouble is that the single bad point is a common start. On an unloaded
radial feeder without shunts, a flat voltage profile puts the line currents at
exactly zero.

```julia
imeas = BranchMeasurement(kind=:imag, line="l", side=:from, terminal="1",
                          value=1.0, sigma=0.01)
s = compile_state_estimator(net, [imeas])
xflat = [230.0, 0.0]          # flat: no current in the line at all

evaluate_state_estimator(s, SEParameters(s, [imeas]; current_epsilon=1e-3), xflat)
```

| `current_epsilon` (A) | predicted ``\|I\|`` | ``\|H\|`` | solver status |
|---|---:|---:|---|
| `0.0` | 0.0000 | *DomainError* | `:undefined_derivative` |
| `1e-3` | 0.0010 | 0.0 | `:converged_underobserved` |
| `1.0` | 1.0000 | 0.0 | `:converged_underobserved` |

Three things to read off this:

* With no smoothing the derivative genuinely does not exist. The estimator
  reports `:undefined_derivative` and returns. It used to throw a `DomainError`
  out of the caller's solve; a property of the measurement set at a point is a
  diagnosis, not a crash.
* `SEParameters(...; current_epsilon=ε)` replaces ``|I|`` with
  ``\sqrt{I_{re}^2+I_{im}^2+\varepsilon^2}``, which is differentiable
  everywhere. The solve now runs. `current_epsilon` is in **amperes**;
  `magnitude_epsilon` is the separate volts-valued knob for ``|V|``, because one
  scalar cannot carry both units.
* **But the smoothed row is information-free at zero current**: its gradient
  there is exactly zero, so ``\|H\| = 0`` and the honest verdict is
  `:converged_underobserved`. Smoothing buys differentiability, not
  information. It does not invent a measurement that the operating point
  cannot support, and the status says so rather than returning a confident
  wrong answer.

The practical remedy is to start somewhere with current in the lines — a
power-flow solution, the previous snapshot, or a nominal-load point — and to
reach for `current_epsilon` knowingly. Sizing it below the ammeter's meaningful
resolution limits the perturbation of the predicted magnitude, but the resulting
sensitivity and estimate should still be checked rather than assumed unchanged.

## Real problem 3: a magnitude carries no direction

This is the difficulty that does not go away, and it is not numerical at all.
It is a property of the measurement function.

Give the estimator ``|V_b|`` and ``|I_\ell|``, both exact, and start it on
either side of the true angle:

```julia
ambiguous = Any[
    Measurement(kind=:vmag, bus="b", reference=nothing,
                value=abs(vtrue), sigma=0.01),
    BranchMeasurement(kind=:imag, line="l", side=:from, terminal="1",
                      value=abs(itrue), sigma=0.01)]
s = compile_state_estimator(net, ambiguous)
p = SEParameters(s, ambiguous)

solve_compiled_state_estimator(s, p, [219.0, -11.0]; initial_radius=1.0)
solve_compiled_state_estimator(s, p, [219.0, +11.0]; initial_radius=1.0)
```

| start | status | ``\|V_b\|`` | ``\angle V_b`` | ``\|r\|`` | observable |
|---|---|---:|---:|---:|---|
| lagging | `:converged_unique` | 220.0000 | **−0.05000** | 2.8e−11 | 2/2 |
| leading | `:converged_unique` | 220.0000 | **+0.05000** | 2.8e−11 | 2/2 |

Two different states. Both fit the data **exactly** on this real-impedance
feeder. Both are reported
`:converged_unique` with a full-rank reduced Jacobian.

Nothing here is a bug, and nothing is ill-conditioned. The two states are a
conjugate pair, and ``|V|`` and ``|I|`` cannot tell them apart because neither
function can see a sign. The local rank test is answering the question it was
asked — *is the estimate locally identified?* — and the answer is genuinely
yes, at each of two points. This is the phenomenon analysed in
Abur & Expósito (1997).

!!! warning "`:converged_unique` is a local statement"
    It reports that the reduced Jacobian ``HZ`` has full rank **at the
    returned point**. It is not a global uniqueness certificate, and no local
    method can supply one. With ampere measurements in the set, solve from
    several starts and compare.

### The fix, when the meter allows it

If the instrument reports a phasor rather than a magnitude, say so. The same
current as a rectangular pair fixes the direction:

```julia
resolved = Any[
    Measurement(kind=:vmag, bus="b", reference=nothing,
                value=abs(vtrue), sigma=0.01),
    BranchMeasurement(kind=:ire, line="l", side=:from, terminal="1",
                      value=real(itrue), sigma=0.01),
    BranchMeasurement(kind=:iim, line="l", side=:from, terminal="1",
                      value=imag(itrue), sigma=0.01)]
```

| start | status | ``\|V_b\|`` | ``\angle V_b`` |
|---|---|---:|---:|
| lagging | `:converged_unique` | 220.0000 | −0.05000 |
| leading | `:converged_unique` | 220.0000 | −0.05000 |

Both starts land on the true state. The ambiguity was in the measurement
model, so that is where it had to be repaired.

## What to actually do

1. **Do not discard ampere measurements on a blanket conditioning argument** —
   but do not assume they are safe either. Check ``\sigma_{\min}`` at your
   operating point. The failure mode is alignment of the ``|V|`` and ``|I|``
   *state-sensitivity rows*. Near-unity power factor is the warning sign for the
   resistive example above; general impedances and coupled lines shift it.
2. **Avoid a zero-current start when `:imag` rows are present.** Warm start from
   a power flow, the previous snapshot, or a nominal-load point. If other rows
   or smoothing move the iterate away from zero, document that mechanism rather
   than relying on an undefined current-magnitude derivative.
3. **Use `current_epsilon` deliberately**, sized below the ammeter's
   resolution, and read the returned status: an `ε`-smoothed row at zero
   current carries no information, and `:converged_underobserved` is telling
   you so.
4. **Prefer `:ire`/`:iim` over `:imag`** whenever the instrument supplies a
   phasor. The magnitude discards exactly the information that resolves the
   ambiguity.
5. **Solve from several starts** and compare, whenever the set is
   magnitude-only. `:converged_unique` at one point says nothing about the
   other.
6. **Anchor the angle.** The mirror pair exists because nothing else in the
   set fixes the sign of the angle. One rectangular voltage or current
   measurement removes it.

The regression tests for all of the above live in
`test/constrained_state_estimation_tests.jl` under
*"|I| admits mirror solutions that |I_re|,|I_im| do not"* and
*"an undefined magnitude derivative is a status, not a throw"*.

## References

Abur, A. and Expósito, A. G., "Detecting multiple solutions in state estimation
in the presence of current magnitude measurements", *IEEE Transactions on Power
Systems*, 12(1):370–375, February 1997. This establishes the **multiple-solution**
result reproduced above. It is not a statement about conditioning, and should
not be cited as one.

Abur, A. and Expósito, A. G., *Power System State Estimation: Theory and
Implementation*, Marcel Dekker, 2004. Discusses the practical difficulties of
ampere measurements, including ill-conditioning near low loading when they are
load-bearing for observability. That low-loading mechanism and the sensitivity
alignment demonstrated above are related practical warnings, not the same
mathematical claim.
