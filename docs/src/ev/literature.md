# Literature evidence and scientific boundaries

These are targeted annotated research notes, not a systematic review or a claim
of literature completeness. Dates refer to publications or identified preprint
versions, not repository download dates. Numerical findings from a particular
vehicle, chemistry or feeder should not be promoted to universal constants.

## Charging systems, equipment and data

| Source | Evidence and assumptions | Consequence for this library |
|:--|:--|:--|
| Z. J. Lee, S. Sharma, D. Johansson and S. H. Low, **ACN-Sim: An Open-Source Simulator for Data-Driven Electric Vehicle Charging Research**, [2021 revision](https://arxiv.org/abs/2012.02809) | A modular simulation environment with EV, EVSE, battery, event and infrastructure models; includes unbalanced infrastructure and nonideal charging. Its infrastructure approximation is not our full nonlinear AC formulation. | Separate session, vehicle and equipment; validate scheduling against a response model before claiming controllability. |
| **ACN-Data**, [official schema](https://ev.caltech.edu/dataset) | Session data distinguishes connection/disconnection, charging completion, delivered energy and user requests. A session record does not by itself identify electrochemical state. | Preserve request-versus-delivery provenance. Do not substitute metered AC kWh for battery-energy gain or fabricate initial SoC. |
| K. Sevdari, L. Calearo, B. H. Bakken, P. B. Andersen and M. Marinelli, **Experimental validation of onboard electric vehicle chargers to improve the efficiency of smart charging operation**, SETA 60, 103512 (2023), [DOI](https://doi.org/10.1016/j.seta.2023.103512) | Measurements of 38 vehicle models show charging efficiency and reactive behavior depend on vehicle and operating current. The paper reports increased energy demand under some current-modulated smart-charging strategies. | Fixed efficiencies and Q=0 are declared approximations. Future loss/response maps need vehicle, test conditions and measurement-boundary metadata. Capability is not reactive-power controllability. |
| K. Bopp, J. Bennett and N. Lee, **Electric Vehicle Supply Equipment: An Overview of Technical Standards**, NREL/PR-7A40-78085 (2020), [report](https://docs.nrel.gov/docs/fy21osti/78085.pdf) | Technical overview distinguishing AC and DC charging and equipment/vehicle limits; not a current normative standard text. | AC EVSE is an outlet and permission interface; DC charging needs offboard conversion ownership. Do not claim standards compliance from a scheduling envelope. |

## Batteries, economics and physical fidelity

| Source | Evidence and assumptions | Consequence for this library |
|:--|:--|:--|
| P. Aaslid, F. Geth, M. Korpås, M. M. Belsnes and O. B. Fosso, **Non-linear charge-based battery storage optimization model with bi-variate cubic spline constraints**, Journal of Energy Storage 32, 101979 (2020), [DOI](https://doi.org/10.1016/j.est.2020.101979) | Optimization in voltage/current/charge coordinates with fitted nonlinear battery relations. | Reuse IVQ where voltage/current limits matter. The existing Rint implementation is a simpler model, not a reproduction of all bivariate spline physics. |
| K. Schwenk et al., **Integrating Battery Aging in the Optimization for Bidirectional Charging of Electric Vehicles**, [2021 preprint revision](https://arxiv.org/abs/2009.12201) | Validated battery modeling integrated with charging optimization; examines thermal and aging effects on costs. Quantitative thresholds and cost changes are case-dependent. | Equal-terminal-energy arbitrage examples establish energy economics only. Lifetime claims require thermal/aging calibration and an explicit cost boundary. |
| H. Movahedi et al., **Extra Throughput versus Days Lost in load-shifting V2G services: Influence of dominant degradation mechanism**, [2024 preprint](https://arxiv.org/abs/2408.02139) | Lifetime simulation across chemistries and environments examines calendar/cycle aging trade-offs. | Do not impose a universal damage-per-kWh coefficient as a validated lifetime model. Use sensitivity bands until suitable cell/pack evidence exists. |

A richer model with guessed parameters is not stronger scientific evidence.
For acceptance curves, retain measured versus illustrative status, AC/DC
boundary, energy/charge coordinate, temperature, SoH, operating range, fitting
method, residuals and held-out validation. The executable taper example is
explicitly illustrative and provides no measured charger parameter preset.

## Optimization and flexibility

| Source | Evidence and assumptions | Consequence for this library |
|:--|:--|:--|
| Z. Li et al., **Further Discussions on Sufficient Conditions for Exact Relaxation of Complementarity Constraints for Storage-Concerned Economic Dispatch**, [2015 preprint](https://arxiv.org/abs/1505.02493) | Sufficient conditions for storage-complementarity relaxation in particular economic dispatch formulations. | Efficiency alone is not an exclusion theorem; establish all assumptions before claiming exactness in another model. |
| A. Joshi, H. Kebriaei, V. Mariani and L. Glielmo, **A Sufficient Condition to Guarantee Non-Simultaneous Charging and Discharging of Household Battery Energy Storage**, [2021 preprint](https://arxiv.org/abs/2104.06267) | Derives a criterion for a household demand-response formulation using duality/KKT arguments. | An objective-based guarantee is formulation-dependent. Nonconvex AC network constraints require their own analysis. |
| N. K. Panda and S. H. Tindemans, **Efficient Quantification and Representation of Aggregate Flexibility in Electric Vehicles**, [preprint](https://arxiv.org/abs/2310.02729) | Exact fleet-flexibility representations within specified request windows and linear charging assumptions. | Preserve temporal/individual feasibility and test disaggregation; a summed kWh/kW battery is not generally exact. |
| K. Mukhi, C. Qu, P. You and A. Abate, **Robust Aggregation of Electric Vehicle Flexiblity**, [2025 revision](https://arxiv.org/abs/2405.08232), [conference DOI](https://doi.org/10.1145/3716863.3718054) | Robust aggregate sets with probabilistic tracking guarantees under sampled charging requirements and stated distributional assumptions. | Scenario data and probability guarantees need explicit assumptions; deterministic session tests cannot establish them. |
| A. Nurkanović, A. Pozharskiy and M. Diehl, **Solving mathematical programs with complementarity constraints arising in nonsmooth optimal control**, [2024 revision](https://arxiv.org/abs/2312.11022) | MPCC constraint qualification failures, stationarity concepts and algorithms. | Separate numerical feasibility, complementarity residuals, stationarity and global optimality. |
| A. Pozharskiy, F. Pacaud, M. Diehl and A. Nurkanović, **CCOpt: an Open-Source Solver for Large-Scale Mathematical Programs with Complementarity Constraints**, [2026 preprint](https://arxiv.org/abs/2604.18726) | Specialized relaxation/penalty and active-set algorithms built on MadNLP, evaluated on several benchmark classes. | Delegate continuation to CCOpt; characterize the EV/AC-network composition rather than inheriting unrelated solver performance claims. |

For smoothing foundations see [the formulation bibliography](../formulations/references.md):
Chen–Mangasarian (1996), Nesterov (2005), and Chen (2012). These motivate smooth
approximations and distinctions between function accuracy, derivative behavior
and limiting stationarity. They do not prove convergence or physical validity
for every composed charging/AC-network problem.

## Reproducibility and claims

The audit of the earlier implementation found that an unplugged three-phase
port could carry nonzero cancelling phase currents and that negative prices
could reward simultaneous inefficient charging/discharging. The new regression
cases test those mechanisms, not merely the presence of particular constraints.

`test/evse_tests.jl` checks physical currents, energy conservation, rating
intersections, occupancy conflicts, normalized leakage, endpoint acceptance and
refinement against a continuous analytical trajectory. `scripts/ev/ccopt.jl`
exercises the optional MPCC path. Record Julia/package versions, operating mode,
normalization, solver settings, statuses and physical residuals when publishing
results; the tests are not calibration evidence or a scalability benchmark.
