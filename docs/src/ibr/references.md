# [IBR references](@id ibr-references)

This bibliography is organised by the layer of physics it supports. A reference
appearing here does not imply that every mechanism in it is implemented.

!!! warning "Verification status of individual entries"
    Entries marked **†** were added during a design review to close attribution
    gaps identified below. Their authors, titles, venues, and years are stated
    from the reviewer's knowledge; **volumes, pages, and DOIs have not been
    checked against live sources** and no DOI is given for them. Verify each one
    before it appears in a submitted manuscript, exactly as the unmarked entries
    already were. Where an entry is a standard whose designation itself is
    uncertain, that is stated in the entry.

## Algebraic and optimisation models

1. **R. Heidari and F. Geth**, “Improved algebraic inverter modeling for
   four-wire power flow optimization,” *Electric Power Systems Research*, vol.
   234, 110825, 2024. [doi:10.1016/j.epsr.2024.110825](https://doi.org/10.1016/j.epsr.2024.110825).
   Primary basis for internal/external nodes, four-conductor KVL/KCL, topology,
   sequence limits, losses, and steady-state GFL/GFM constraints.
2. **M. Deakin, R. Heidari, and X. Deng**, “Power converter DC link ripple and
   network unbalance as active constraints in distribution system optimal power
   flow,” arXiv:2512.18293, 2025.
   [doi:10.48550/arXiv.2512.18293](https://doi.org/10.48550/arXiv.2512.18293).
   Primary basis for ``\widetilde S``, the small-ripple DC-link relation, and the
   six equation/PLECS benchmark cases.
3. **M. Deakin, R. Heidari, and X. Deng**, “DC link capacitor ripple constraints
   limit benefits of utility-owned four-wire power converters,” arXiv:2606.21934,
   2026. [arXiv:2606.21934](https://arxiv.org/abs/2606.21934).
   Primary basis for simultaneous capacitor ripple allocation and the distinction
   among 3-leg, 4-leg, and reconfigurable split-link return paths.
4. **E. O. Badmus and A. Pandey**, “Two-stage bidirectional inverter equivalent
   circuit model for distribution grid steady-state analysis and optimization,”
   *IEEE Transactions on Power Systems*, in press, 2026. Follow this work for
   two-stage DC/DC-plus-DC/AC coupling when the public final article is available.

Items 2, 3, and 21 are cited from preprints. Replace each with its published
version, and update the table/figure numbers used by the regression tests, before
any of them supports a manuscript claim.

## Split link, ripple, and switching physics

5. **I. Ziyat, J. Wang, and P. R. Palmer**, “Voltage ripple model and capacitor
   sizing for the three-phase four-wire converter used for power redistribution,”
   *IEEE Journal of Emerging and Selected Topics in Power Electronics*, vol. 12,
   no. 2, pp. 1437-1445, 2024.
   [doi:10.1109/JESTPE.2023.3289485](https://doi.org/10.1109/JESTPE.2023.3289485).
   Basis for treating the split midpoint as a charge-balance problem and for
   retaining the two half-banks as distinct physical components. PowerOptLab's
   unequal-capacitance and bounded-charge equations are quasi-static reductions,
   not a reproduction of the paper's switching/control dynamics.
6. **J. Liang, T. C. Green, C. Feng, and G. Weiss**, “Increasing voltage
   utilization in split-link, four-wire inverters,” *IEEE Transactions on Power
   Electronics*, vol. 24, no. 6, pp. 1562-1569, 2009.
   [doi:10.1109/TPEL.2009.2013351](https://doi.org/10.1109/TPEL.2009.2013351).
7. **A. Viatkin et al.**, “Analysis of the output current ripple of the
   three-phase four-leg inverter with a neutral inductor,” *Energies*, vol. 14,
   1430, 2021. [doi:10.3390/en14051430](https://doi.org/10.3390/en14051430).
   Basis for phase/neutral ripple with an arbitrary neutral-to-phase inductance
   ratio and for the corresponding topology-ordering tests.
8. **R. Mandrioli et al.**, “Prediction of DC-link voltage switching ripple in
   three-phase four-leg PWM inverters,” *Energies*, vol. 14, 1434, 2021.
   [doi:10.3390/en14051434](https://doi.org/10.3390/en14051434).
   Basis for the ideal shared-carrier DC-link audit and the balanced SPWM and
   centered-PWM closed-form regression tests.
9. **M. Vujacic, M. Hammami, M. Srndovic, and G. Grandi**, “Analysis of
   dc-link voltage switching ripple in three-phase PWM inverters,” *Energies*,
   vol. 11, 471, 2018.
   [doi:10.3390/en11020471](https://doi.org/10.3390/en11020471).
   Basis for the finite series R–L source branch, DC-node harmonic current
   sharing, and high-source-impedance limit. The implementation extends that
   circuit beyond the paper's balanced centered-PWM operating cases.
10. **A. Hammami et al.**, “Analysis of input voltage switching ripple in
    three-phase four-wire split capacitor PWM inverters,” *Energies*, vol. 13,
    5076, 2020. [doi:10.3390/en13195076](https://doi.org/10.3390/en13195076).
    Split-capacitor switching-ripple basis and an independent source for future
    modulation/load-angle map regressions.
11. **R. Mandrioli, M. Hammami, A. Viatkin, R. Barbone, D. Pontara, and M.
    Ricco**, “Phase and neutral current ripple analysis in three-phase four-wire
    split-capacitor grid converter for EV chargers,” *Electronics*, vol. 10,
    1016, 2021. [doi:10.3390/electronics10091016](https://doi.org/10.3390/electronics10091016).
    Source of the implemented split-link SPWM phase and neutral RMS regression
    formulas, scaling laws, and experimental parameter set.
12. **R. Mandrioli, F. Lo Franco, M. Ricco, and G. Grandi**, “A generalized
    approach for determining the current ripple RMS in four-leg inverters with
    the neutral inductor,” *Energies*, vol. 16, 1710, 2023.
    [doi:10.3390/en16041710](https://doi.org/10.3390/en16041710).
    General arbitrary-common-mode/neutral-inductance formulation, validated for
    SPWM, SVPWM, DPWM, and THIPWM in simulation and experiment.
13. **† R. Zhang, V. H. Prasad, D. Boroyevich, and F. C. Lee**, “Three-dimensional
    space vector modulation for four-leg voltage-source converters,” *IEEE
    Transactions on Power Electronics*, 2002. The canonical four-leg modulation
    reference. PowerOptLab's four-leg switching hull is a sampled two-level rail
    feasibility region, not a modulation strategy; this is the reference against
    which any future four-leg modulation claim should be positioned.

Items 7-13 support the implemented DC-link and AC-side carrier audits; `i_sw`
remains a separate allowance for residual unmodelled content. The package's
harmonic-network solve generalizes the papers' independent-inductor circuits to
the supplied primitive L/LCL model, but not their full experimental validation
matrix. Items 5-6 remain validation anchors for unequal half-bank dynamics,
active balancing, and improved voltage utilisation beyond the present
steady-state charge abstraction.

## Filters, impedance, and frequency coupling

14. **M. Liserre, F. Blaabjerg, and S. Hansen**, “Design and control of an LCL-
    filter-based three-phase active rectifier,” *IEEE Transactions on Industry
    Applications*, vol. 41, no. 5, pp. 1281-1291, 2005.
    [doi:10.1109/TIA.2005.853373](https://doi.org/10.1109/TIA.2005.853373).
    Basis for separating converter- and grid-side inductive arms, the midpoint
    capacitor/damping branch, and the scalar undamped resonance screening
    formula. PowerOptLab optimises the passive circuit at the fundamental and
    optionally audits its constant R/L/C primitives at carrier harmonics; it does
    not implement the paper's control-design or stability layer.
15. **J. Sun**, “Impedance-based stability criterion for grid-connected
    inverters,” *IEEE Transactions on Power Electronics*, vol. 26, no. 11,
    pp. 3075-3078, 2011.
    [doi:10.1109/TPEL.2011.2136439](https://doi.org/10.1109/TPEL.2011.2136439).
16. **A. Rygg, M. Molinas, C. Zhang, and X. Cai**, “A modified sequence-domain
    impedance definition and its equivalence to the dq-domain impedance
    definition for the stability analysis of AC power electronic systems,”
    *IEEE Journal of Emerging and Selected Topics in Power Electronics*, vol. 4,
    no. 4, pp. 1383-1396, 2016.
    [doi:10.1109/JESTPE.2016.2588733](https://doi.org/10.1109/JESTPE.2016.2588733).

Items 15-16 define a separate future model layer: frequency-dependent
controller, PLL, grid, and sequence-coupling impedances. They should not be
collapsed into the present fundamental-frequency passive LCL algebraic model.

## Numerical methods and smoothing

These support the package's own claims about smoothing bias, conditioning, and
locality of solutions. They were previously recorded only in source comments.

17. **Yu. Nesterov**, “Smooth minimization of non-smooth functions,”
    *Mathematical Programming*, vol. 103, pp. 127-152, 2005.
    [doi:10.1007/s10107-004-0552-5](https://doi.org/10.1007/s10107-004-0552-5).
    The classical accuracy-versus-smoothness trade-off: ``O(\mu)`` approximation
    error against an ``O(1/\mu)`` gradient Lipschitz constant. Its worst-case
    bound is over a function class and is *not* what the package's own
    measurements show, because here the smoothed norm is one term in a problem
    whose conditioning is set by the network equations.
18. **C. Chen and O. L. Mangasarian**, “A class of smoothing functions for
    nonlinear and mixed complementarity problems,” *Computational Optimization
    and Applications*, vol. 5, pp. 97-138, 1996.
    [doi:10.1007/BF00249052](https://doi.org/10.1007/BF00249052).
    The smoothing class that ``\sqrt{x^2+\epsilon^2}`` belongs to, and therefore
    the reference for every smooth min/max selector in the control laws.
19. **P. Charbonnier, L. Blanc-Féraud, G. Aubert, and M. Barlaud**,
    “Deterministic edge-preserving regularization in computed imaging,” *IEEE
    Transactions on Image Processing*, vol. 6, no. 2, pp. 298-311, 1997.
    [doi:10.1109/83.551699](https://doi.org/10.1109/83.551699).
    The same function under its pseudo-Huber/Charbonnier name.
20. **A. Wächter and L. T. Biegler**, “On the implementation of an interior-point
    filter line-search algorithm for large-scale nonlinear programming,”
    *Mathematical Programming*, vol. 106, pp. 25-57, 2006.
    [doi:10.1007/s10107-004-0559-y](https://doi.org/10.1007/s10107-004-0559-y).
    Ipopt. Its tolerance design sets the accuracy floor of every implicit
    magnitude equality here, and its Eqn 35 bound relaxation is the mechanism
    behind the per-unit rating-violation note in the [`AdvancedInverter`](@ref)
    docstring and the publication gate in
    [Verification and benchmark cases](@ref ibr-verification).

## Converter-level sequence control

21. **C. L. Fortescue**, “Method of symmetrical co-ordinates applied to the
    solution of polyphase networks,” *Transactions of the AIEE*, vol. 37, pp.
    1027-1140, 1918.
    [doi:10.1109/T-AIEE.1918.4765570](https://doi.org/10.1109/T-AIEE.1918.4765570).
22. **C. J. O'Rourke, M. M. Qasim, M. R. Overlin, and J. L. Kirtley**, “A
    geometric interpretation of reference frames and transformations: dq0,
    Clarke, and Park,” *IEEE Transactions on Energy Conversion*, vol. 34, no. 4,
    pp. 2070-2083, 2019.
    [doi:10.1109/TEC.2019.2941175](https://doi.org/10.1109/TEC.2019.2941175).
23. **† H.-S. Song and K. Nam**, “Dual current control scheme for PWM converter
    under unbalanced input voltage conditions,” *IEEE Transactions on Industrial
    Electronics*, 1999. The classical dual-sequence current-control scheme, and
    the origin of treating positive- and negative-sequence current references as
    two independently commanded channels.
24. **† H.-S. Suh and T. A. Lipo**, work on instantaneous active and reactive
    power of a PWM AC/DC converter under generalized unbalanced network
    conditions, *IEEE Transactions*, 2006. The derivation of the mean plus
    twice-fundamental power algebra under unbalance, and of the four real
    degrees of freedom available to a three-wire converter. This is the
    provenance of the ``S`` and ``\widetilde S`` decomposition used in
    [Phase-aware local control laws](@ref ibr-phase-aware-control-laws), which
    currently attributes it only to item 30. Confirm which of the two
    Suh–Lipo 2006 papers to cite.
25. **† P. Rodríguez, A. V. Timbus, R. Teodorescu, M. Liserre, and F.
    Blaabjerg**, “Flexible active power control of distributed power generation
    systems during grid faults,” *IEEE Transactions on Industrial Electronics*,
    2007. Flexible positive/negative-sequence current-reference generation with
    a single scalar weighting between the two sequences. This is the direct
    antecedent of the implemented `ripple_blend` ``\lambda`` and should be cited
    wherever that blend is introduced.
26. **† R. Teodorescu, M. Liserre, and P. Rodríguez**, *Grid Converters for
    Photovoltaic and Wind Power Systems*, Wiley, 2011. Textbook treatment of
    sequence extraction, dual-sequence reference generation, and the associated
    power-oscillation algebra.
27. **F. Nejabatkhah, Y. W. Li, and B. Wu**, “Control strategies of three-phase
    distributed generation inverters for grid unbalanced voltage compensation,”
    *IEEE Transactions on Power Electronics*, vol. 31, no. 7, pp. 5228-5241,
    2016. [doi:10.1109/TPEL.2015.2479601](https://doi.org/10.1109/TPEL.2015.2479601).
    Basis for treating positive- and negative-sequence current references as
    separate control degrees of freedom under unbalanced voltage.
28. **Y. Guo, B. C. Pal, and R. A. Jabr**, “On the optimality of voltage
    unbalance attenuation by inverters,” arXiv:2109.10974, 2021.
    [doi:10.48550/arXiv.2109.10974](https://doi.org/10.48550/arXiv.2109.10974).
    Demonstrates why current, active-power, feasibility, and synchronization
    constraints belong in negative-sequence voltage-attenuation design.
29. **N. Helaly**, “A predictive negative sequence current control algorithm
    for voltage imbalance compensation and power oscillation minimization under
    unbalanced conditions,” *Applied Science and Engineering Progress*, vol. 16,
    no. 4, 6562, 2023.
    [doi:10.14416/j.asep.2023.01.003](https://doi.org/10.14416/j.asep.2023.01.003).
    Provides a control-oriented example of explicitly trading voltage-unbalance
    compensation against oscillating power.
30. **† M. Savaghebi, A. Jalilian, J. C. Vasquez, and J. M. Guerrero**,
    “Secondary control scheme for voltage unbalance compensation in an islanded
    droop-controlled microgrid,” *IEEE Transactions on Smart Grid*, 2012.
    Representative of the negative-sequence virtual-admittance/virtual-impedance
    family that the implemented ``I_2^v=-\kappa e^{-j\phi_2}U_2`` law belongs
    to. Cite it so the admittance form is not read as novel.

Items 23-30 motivate the candidate local laws in [Phase-aware local control
laws](@ref ibr-phase-aware-control-laws) and are the converter-control prior art
for the dual-sequence reference generation used here. They do not by themselves
validate the PowerOptLab allocator or its interaction with the implemented
switching and capacitor constraints. The contribution under study is the
fixed-structure algebraic surrogate and its network embedding, not the sequence
reference itself, and the bibliography should make that ordering obvious.

## Grid-forming and validation guardrails

31. **Y. Lin et al.**, *Research Roadmap on Grid-Forming Inverters*,
    NREL/TP-5D00-73476, 2020.
    [NREL report](https://www.nrel.gov/docs/fy21osti/73476.pdf).
32. **S. Shah and D. Ramasubramanian**, *Testing the Performance of Grid-Forming
    Resources: Test Methods and Performance Metrics for Evaluating the Voltage
    Source Behavior of Grid-Forming Resources*, NREL/TP-5D00-94378, 2025.
    [report record](https://research-hub.nrel.gov/en/publications/testing-the-performance-of-grid-forming-resources-test-methods-an/).
33. **† J. Rocabert, A. Luna, F. Blaabjerg, and P. Rodríguez**, “Control of power
    converters in AC microgrids,” *IEEE Transactions on Power Electronics*,
    2012. The standard grid-following/grid-forming taxonomy; useful for keeping
    the steady-state `grid_forming=true` constraint clearly separated from the
    control class it is named after.

Items 31-33 are guardrails on terminology and validation: a steady-state balanced
internal EMF is only one property of a GFM resource. Dynamic voltage-source
behaviour, current limiting, fault recovery, interaction stability, and hardware
testing need their own models and performance tests.

## Local voltage control and unbalance in distribution networks

This layer motivates the *distribution-side* question the control study asks. It
was absent from earlier revisions of this bibliography, which cited only
converter-side control.

34. **† K. Turitsyn, P. Šulc, S. Backhaus, and M. Chertkov**, “Options for
    control of reactive power by distributed photovoltaic generators,”
    *Proceedings of the IEEE*, 2011. Foundational treatment of local reactive
    control by distributed PV.
35. **† M. Jahangiri and D. C. Aliprantis**, “Distributed Volt/VAr control by PV
    inverters,” *IEEE Transactions on Power Systems*, 2013.
36. **† M. Farivar, X. Chen, and S. H. Low**, “Equilibrium and dynamics of local
    voltage control in distribution systems,” *IEEE Conference on Decision and
    Control*, 2013. Existence, uniqueness, and convergence of the equilibrium of
    a local Volt-var droop. This is the theory behind the repeated caveat that a
    solved equilibrium is not a stability result, and behind the argument that a
    discontinuous conflict rule may admit no equilibrium near the tie surface.
37. **† H. Zhu and H. J. Liu**, “Fast local voltage control under limited
    reactive power: optimality and stability analysis,” *IEEE Transactions on
    Power Systems*, 2016. Local droop as a surrogate-gradient method, with
    convergence conditions on the droop slope. Directly relevant to bounding
    the curve slopes stamped here.
38. **† K. Baker, A. Bernstein, E. Dall'Anese, and C. Zhao**,
    “Network-cognizant voltage droop control for distribution grids,” *IEEE
    Transactions on Power Systems*, 2018. Droop gains designed from network
    sensitivities; the closest published relative of the ``\widehat H_2^{-1}``
    sensitivity idea in the negative-sequence section.
39. **† E. Dall'Anese and A. Simonetto**, “Optimal power flow pursuit,” *IEEE
    Transactions on Smart Grid*, 2018. Local controllers that track the solution
    of an optimization problem online. The candidate closest-feasible per-phase
    projection law is a member of this family and should be positioned against
    it rather than presented as a bespoke oracle.
40. **† N. Weckx and J. Driesen**, “Load balancing with EV chargers and PV
    inverters in unbalanced distribution grids,” *IEEE Transactions on
    Sustainable Energy*, 2015. Unbalance mitigation by distributed inverters at
    feeder scale — the distribution-side counterpart of the converter-side
    references above. (Confirm the first author's initials.)
41. **† F. Shahnia, R. Majumder, A. Ghosh, G. Ledwich, and F. Zare**, “Voltage
    imbalance analysis in residential low voltage distribution networks with
    rooftop PV,” *Electric Power Systems Research*, 2011. Establishes the
    magnitude and mechanism of the LV unbalance this work is trying to control.
42. **† P. Pillay and M. Manyage**, “Definitions of voltage unbalance,” *IEEE
    Power Engineering Review*, 2001. The distinction between the true
    ``|U_2|/|U_1|`` factor and the NEMA/IEEE approximations. Cite it wherever
    `voltage_unbalance_factor` is defined, so the reported metric is
    unambiguous.
43. **† V. Girigoudar and L. A. Roald**, work on voltage-unbalance metrics and
    unbalance constraints in distribution-system optimization, *Electric Power
    Systems Research* and related venues, from 2020. The closest prior work on
    putting unbalance metrics inside an optimization model, and therefore the
    natural comparison for the VUF reported by the study layer.

## Standards and reference implementations

44. **IEEE Std 1547-2018**, *IEEE Standard for Interconnection and
    Interoperability of Distributed Energy Resources with Associated Electric
    Power Systems Interfaces*.
    [IEEE record](https://standards.ieee.org/ieee/1547/5915/).
45. **IEEE Std 1547.1-2020**, *IEEE Standard Conformance Test Procedures for
    Equipment Interconnecting Distributed Energy Resources with Electric Power
    Systems and Associated Interfaces*.
    [IEEE record](https://standards.ieee.org/ieee/1547.1/6039/).
46. **† IEEE Std 2800-2022**, *IEEE Standard for Interconnection and
    Interoperability of Inverter-Based Resources (IBRs) Interconnecting with
    Associated Transmission Electric Power Systems*.
    [IEEE record](https://standards.ieee.org/ieee/2800/10453/).
    The only standard cited here that addresses negative-sequence current
    behaviour normatively. It is a bulk-system document and is not an LV
    controller specification, but a discussion of unbalance control that omits
    it is incomplete.
47. **AS/NZS 4777.2**, *Grid connection of energy systems via inverters,
    Part 2: Inverter requirements*.
    [AEMO standards overview](https://www.aemo.com.au/initiatives/major-programs/nem-distributed-energy-resources-der-program/standards-and-connections/as-nzs-4777-2-inverter-requirements-standard).
    The licensed clause specifying which voltage a multi-phase inverter observes
    for Volt-var and Volt-watt has **not** been read; see the open novelty check
    in [Phase-aware local control laws](@ref ibr-phase-aware-control-laws).
48. **† EN 50549-1 and EN 50549-2**, *Requirements for generating plants to be
    connected in parallel with distribution networks*. European connection-code
    comparator. Designation and edition unverified here.
49. **† VDE-AR-N 4105**, *Generators connected to the low-voltage distribution
    network*. German LV connection rule. Designation and edition unverified here.
50. **IEC 61000-4-30:2025**, *Electromagnetic compatibility — Testing and
    measurement techniques — Power quality measurement methods*.
    [IEC record](https://webstore.iec.ch/en/publication/71611).
51. **IEC TR 61000-3-13:2008**, *Assessment of emission limits for the
    connection of unbalanced installations to MV, HV and EHV power systems*.
    [IEC record](https://webstore.iec.ch/en/publication/4145).
52. **† EN 50160**, *Voltage characteristics of electricity supplied by public
    distribution networks*. Cited in the standards-mapping table of
    [Phase-aware local control laws](@ref ibr-phase-aware-control-laws) but
    previously absent from this bibliography. **The designation, title, and
    edition are unverified**; pin them, or remove the mapping-table row, before
    the supply-voltage-characteristics comparison is used in a manuscript.
53. **EPRI OpenDSS documentation**, `InvControl` monitored-voltage,
    Volt-watt-axis, and inverter-priority properties.
    [monitored voltage](https://opendss.epri.com/Commonproperties.html),
    [Volt-watt bases](https://opendss.epri.com/Propertiesofsmartinvertervolt-wa.html),
    [capability priority](https://opendss.epri.com/ViolationofInverterCapabilityCur.html).

Items 44-53 define terminology, test scope, and independent reference-tool
behaviour. PowerOptLab does not reproduce their complete profiles or claim
conformance; any such study must use the current licensed standard and its
prescribed measurement and response-time procedures.

## Literature watch list

The highest-value additions are:

- the licensed AS/NZS 4777.2 and IEEE 1547-2018 clauses that specify the
  monitored voltage for a multi-phase DER, which decide whether the worst-phase
  envelope is a contribution or prior art;
- parameter-identical EMT validation datasets for all three implemented
  topologies;
- dynamic midpoint balancing with leakage, actuator power/loss, and per-half
  overvoltage protection;
- frequency- and temperature-dependent capacitor ESR/ESL electrical and thermal
  models, including lifetime/ageing feedback;
- zero-sequence/common-mode grounding and parasitic-capacitance paths, which
  also determine what a three-wire inverter's voltage sensor actually measures;
- frequency-dependent phase/neutral magnetic models, plus dead time,
  discontinuous PWM, interleaved carriers, and overmodulation;
- measured broadband battery/DC-stage impedance, source-control interaction,
  and validation beyond the implemented constant series R–L approximation;
- control-aware harmonic linearisation with positive/negative-sequence coupling;
- current-limited GFM equilibria linked to dynamic fault-ride-through models;
- multilevel and two-stage DC/DC-plus-DC/AC converters; and
- an explicit reconfigurable four-leg-plus-split-link topology before any soft
  open point work is attempted.
