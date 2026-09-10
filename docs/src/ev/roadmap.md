# Scope, extensions and research validation

## Implemented foundation

The first implementation provides disconnected conductor isolation, current and
apparent-power ratings, equal-power or explicit independent phase dispatch,
normalized mode relaxation and exact MPCC encoding, fixed AC outlet assignments,
vehicle/equipment/session capability intersections, energy diagnostics, and
optional monotone acceptance envelopes with exact or conservative smooth bounds.

It retains `EVDevice` and the PE stored-energy convention while making the
numerical relaxation explicit. All multi-phase PE ports now require a conductor
rating. These changes should be treated as model changes when reproducing old
studies, not as a numerical performance update.

## Next layers and their intended uses

| Extension | Enables | Required design/validation |
|:--|:--|:--|
| Persistent vehicle trajectories and trip withdrawals | Home–work–home, fleets and depot returns | Boundary event ordering, measured arrival states, trip-energy accounting; never reset energy at each visit |
| DC EVSE and shared conversion groups | Multiple dispensers on one power cabinet | Stamp AC conversion once; couple per-port DC V/I/P and cabinet limits; distinguish cable and converter ratings |
| Measured loss and pilot-response maps | Low-current efficiency, deliverable control schedules | AC/DC boundaries, standby/on losses, pilot versus accepted current, response and sampling times |
| IVQ/thermal/aging composition | Fast-charge feasibility, temperature and lifetime studies | Shared battery/converter interfaces, calibrated data and independent replay; avoid duplicate losses |
| Scenario/MPC driver | Early departures, outages, forecast errors | State transfer, nonanticipativity, explicit shortfalls, out-of-sample evaluation |
| Aggregate service contracts | DOE/ancillary services | Disaggregation, activation duration, recovery, uncertainty and feeder feasibility |
| Minimum/discrete pilots and assignment MPCCs | Operational scheduling and allocation | Exact discrete semantics, multiple starts, numerical residuals, implementation replay |

These are prospective features. In particular, the current model does not
implement trip withdrawals, shared DC cabinets, stochastic control or optimized
assignments. A session's occupancy mask is not a mobile battery.

## NLP and MPCC design principles

1. Keep fixed assignments, equipment schedules and known trips as input data
   unless their optimization answers the study question.
2. Impose multiple capability bounds separately instead of introducing a
   nonsmooth minimum unnecessarily.
3. Reuse exact affine/concave upper-bound formulations before smoothing.
4. Normalize complementarity pairs with physical equipment scales and report
   physical as well as dimensionless residuals.
5. Treat a finite product relaxation, fractional assignment or averaged pilot as
   an approximation with a quantified implementation gap.
6. Keep sampling, model fit, smoothing and solver tolerances separate.
7. Compare economic cases with matched initial/final inventory and distinguish
   physical service from a model's controllable feasible set.

## Research acceptance tests for future extensions

A feature should include a counterexample that demonstrates why it is needed,
an analytical or independently simulated reference, a calibrated/illustrative
data declaration, and a domain where its conclusion is valid. Examples include
trip-energy conservation over repeated visits; shared-cabinet conservation with
two outlets; held-out pilot tracking; temperature trajectories against measured
charging; and explicit disaggregation of every sampled aggregate offer.

Smooth NLP or MPCC compatibility is a mathematical modeling choice. It does not
supply global optimality, field controllability, standards compliance, or a
battery safety certificate on its own.
