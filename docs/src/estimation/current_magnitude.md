# Current-magnitude measurements

> **Audience:** anyone who has been told that ampere readings are numerically
> dangerous in state estimation.

Current-magnitude (`:imag`) measurements have a bad reputation. The folklore
says they make the estimation problem ill-conditioned. **That is not what
happens**, and believing it leads people to discard perfectly good telemetry
for the wrong reason.

Two things really do go wrong, and neither is conditioning. This page
demonstrates both, and shows what the estimator does about them.

Everything below runs against a two-bus feeder:

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

## The myth: ill-conditioning

It is easy to check. Take ``|V_b|`` and ``|I_\ell|``, and watch the
conditioning of the whitened Jacobian as the line current falls by six orders
of magnitude (holding a 0.9 lagging angle so the state stays genuinely
complex):

| ``|I|`` (A) | ``\sigma_{\min}(H)`` | ``\mathrm{cond}(H)`` |
|---:|---:|---:|
| 60 | 4.50e+01 | 4.87 |
| 20 | 4.13e+01 | 5.33 |
| 5 | 4.00e+01 | 5.50 |
| 1 | 3.97e+01 | 5.54 |
| 0.01 | 3.96e+01 | 5.56 |
| 0.0001 | 3.96e+01 | 5.56 |

The condition number does not blow up. It *converges*, to about 5.6. The
gradient of ``\sqrt{I_{re}^2+I_{im}^2}`` is a unit vector in the direction of
the current — its direction depends on the operating point, but its magnitude
does not shrink with ``|I|``. A measurement set dominated by ampere readings
conditions just as well as any other:

| measurement set | rows | states | observable | cond |
|---|---:|---:|---:|---:|
| ``\|I\|`` on both lines only | 2 | 4 | 2 | ``\infty`` |
| + ``\|V_b\|`` | 3 | 4 | 3 | ``\infty`` |
| + ``\|V_b\|``, ``\|V_c\|`` | 4 | 4 | 4 | 9.71 |

Rank grows one-for-one with rows, and the fully determined set is well
conditioned. The under-observed rows are simply *too few equations* — the
ordinary situation, not a pathology of ampere measurements.

## Real problem 1: the derivative is undefined at exactly zero current

``\partial |I| / \partial x`` does not exist at ``I = 0``. Not "is large" —
does not exist. And the difficulty is a **single point**, not a
neighbourhood: at ``|I| = 10^{-4}\,\mathrm{A}`` the Jacobian is perfectly
healthy, as the table above shows.

The trouble is that the single bad point is the one everybody starts from. A
flat start — every voltage at nominal, no load anywhere — puts *every* line
current at exactly zero.

```julia
s = compile_state_estimator(net, [imeas])
xflat = [230.0, 0.0]          # flat: no current in the line at all
```

| `magnitude_epsilon` | predicted ``\|I\|`` | ``\|H\|`` | solver status |
|---|---:|---:|---|
| `0.0` | 0.0000 | *DomainError* | `:undefined_derivative` |
| `1e-3` | 0.0010 | 0.0 | `:converged_underobserved` |
| `1.0` | 1.0000 | 0.0 | `:converged_underobserved` |

Three things to read off this:

* With no smoothing the derivative genuinely does not exist. The estimator
  reports `:undefined_derivative` and returns. It used to throw a `DomainError`
  out of the caller's solve; a property of the measurement set at a point is a
  diagnosis, not a crash.
* `SEParameters(...; magnitude_epsilon=ε)` replaces ``|I|`` with
  ``\sqrt{I_{re}^2+I_{im}^2+\varepsilon^2}``, which is differentiable
  everywhere. The solve now runs.
* **But the smoothed row is information-free at zero current**: its gradient
  there is exactly zero, so ``\|H\| = 0`` and the honest verdict is
  `:converged_underobserved`. Smoothing buys differentiability, not
  information. It does not invent a measurement that the operating point
  cannot support, and the status says so rather than returning a confident
  wrong answer.

The practical remedy is to start somewhere with current in the lines — a
power-flow solution, the previous snapshot, or a nominal-load point — and to
reach for `magnitude_epsilon` knowingly, sized below the instrument's
meaningful resolution.

## Real problem 2: a magnitude carries no direction

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

Two different states. Both fit the data **exactly**. Both are reported
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

1. **Do not discard ampere measurements for conditioning.** That reason is
   not supported by the numbers.
2. **Never start flat when `:imag` rows are present.** Warm start from a power
   flow, the previous snapshot, or a nominal-load point. A flat start is the
   one operating point where the derivative does not exist.
3. **Use `magnitude_epsilon` deliberately**, sized below the meter's
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

## Reference

Abur, A. and Expósito, A. G., "Detecting multiple solutions in state
estimation in the presence of current magnitude measurements", *IEEE
Transactions on Power Systems*, 12(1):370–375, February 1997.
