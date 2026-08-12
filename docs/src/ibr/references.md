# [IBR references](@id ibr-references)

This bibliography is organised by the layer of physics it supports. A reference
appearing here does not imply that every mechanism in it is implemented.

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

Items 7-12 support the implemented DC-link and AC-side carrier audits; `i_sw`
remains a separate allowance for residual unmodelled content. The package's
harmonic-network solve generalizes the papers' independent-inductor circuits to
the supplied primitive L/LCL model, but not their full experimental validation
matrix. Items 5-6 remain validation anchors for unequal half-bank dynamics,
active balancing, and improved voltage utilisation beyond the present
steady-state charge abstraction.

## Filters, impedance, and frequency coupling

13. **M. Liserre, F. Blaabjerg, and S. Hansen**, “Design and control of an LCL-
   filter-based three-phase active rectifier,” *IEEE Transactions on Industry
   Applications*, vol. 41, no. 5, pp. 1281-1291, 2005.
   [doi:10.1109/TIA.2005.853373](https://doi.org/10.1109/TIA.2005.853373).
   Basis for separating converter- and grid-side inductive arms, the midpoint
   capacitor/damping branch, and the scalar undamped resonance screening
   formula. PowerOptLab optimises the passive circuit at the fundamental and
   optionally audits its constant R/L/C primitives at carrier harmonics; it does
   not implement the paper's control-design or stability layer.
14. **J. Sun**, “Impedance-based stability criterion for grid-connected
    inverters,” *IEEE Transactions on Power Electronics*, vol. 26, no. 11,
    pp. 3075-3078, 2011.
    [doi:10.1109/TPEL.2011.2136439](https://doi.org/10.1109/TPEL.2011.2136439).
15. **A. Rygg, M. Molinas, C. Zhang, and X. Cai**, “A modified sequence-domain
    impedance definition and its equivalence to the dq-domain impedance
    definition for the stability analysis of AC power electronic systems,”
    *IEEE Journal of Emerging and Selected Topics in Power Electronics*, vol. 4,
    no. 4, pp. 1383-1396, 2016.
    [doi:10.1109/JESTPE.2016.2588733](https://doi.org/10.1109/JESTPE.2016.2588733).

Items 14-15 define a separate future model layer: frequency-dependent
controller, PLL, grid, and sequence-coupling impedances. They should not be
collapsed into the present fundamental-frequency passive LCL algebraic model.

## Foundations, controls, and validation

16. **C. L. Fortescue**, “Method of symmetrical co-ordinates applied to the
    solution of polyphase networks,” *Transactions of the AIEE*, vol. 37, pp.
    1027-1140, 1918.
    [doi:10.1109/T-AIEE.1918.4765570](https://doi.org/10.1109/T-AIEE.1918.4765570).
17. **C. J. O'Rourke, M. M. Qasim, M. R. Overlin, and J. L. Kirtley**, “A
    geometric interpretation of reference frames and transformations: dq0,
    Clarke, and Park,” *IEEE Transactions on Energy Conversion*, vol. 34, no. 4,
    pp. 2070-2083, 2019.
    [doi:10.1109/TEC.2019.2941175](https://doi.org/10.1109/TEC.2019.2941175).
18. **Y. Lin et al.**, *Research Roadmap on Grid-Forming Inverters*,
    NREL/TP-5D00-73476, 2020.
    [NREL report](https://www.nrel.gov/docs/fy21osti/73476.pdf).
19. **S. Shah and D. Ramasubramanian**, *Testing the Performance of Grid-Forming
    Resources: Test Methods and Performance Metrics for Evaluating the Voltage
    Source Behavior of Grid-Forming Resources*, NREL/TP-5D00-94378, 2025.
    [report record](https://research-hub.nrel.gov/en/publications/testing-the-performance-of-grid-forming-resources-test-methods-an/).
20. **F. Nejabatkhah, Y. W. Li, and B. Wu**, “Control strategies of three-phase
    distributed generation inverters for grid unbalanced voltage compensation,”
    *IEEE Transactions on Power Electronics*, vol. 31, no. 7, pp. 5228-5241,
    2016. [doi:10.1109/TPEL.2015.2479601](https://doi.org/10.1109/TPEL.2015.2479601).
    Basis for treating positive- and negative-sequence current references as
    separate control degrees of freedom under unbalanced voltage.
21. **Y. Guo, B. C. Pal, and R. A. Jabr**, “On the optimality of voltage
    unbalance attenuation by inverters,” arXiv:2109.10974, 2021.
    [doi:10.48550/arXiv.2109.10974](https://doi.org/10.48550/arXiv.2109.10974).
    Demonstrates why current, active-power, feasibility, and synchronization
    constraints belong in negative-sequence voltage-attenuation design.
22. **N. Helaly**, “A predictive negative sequence current control algorithm
    for voltage imbalance compensation and power oscillation minimization under
    unbalanced conditions,” *Applied Science and Engineering Progress*, vol. 16,
    no. 4, 6562, 2023.
    [doi:10.14416/j.asep.2023.01.003](https://doi.org/10.14416/j.asep.2023.01.003).
    Provides a control-oriented example of explicitly trading voltage-unbalance
    compensation against oscillating power.

Items 18-19 are guardrails on terminology and validation: a steady-state balanced
internal EMF is only one property of a GFM resource. Dynamic voltage-source
behaviour, current limiting, fault recovery, interaction stability, and hardware
testing need their own models and performance tests.

Items 20-22 motivate the candidate local laws in [Phase-aware local control
laws](@ref ibr-phase-aware-control-laws). They do not by themselves validate the
proposed PowerOptLab allocator or its interaction with the implemented switching
and capacitor constraints.

## Literature watch list

The highest-value additions are:

- parameter-identical EMT validation datasets for all three implemented
  topologies;
- dynamic midpoint balancing with leakage, actuator power/loss, and per-half
  overvoltage protection;
- frequency- and temperature-dependent capacitor ESR/ESL electrical and thermal
  models, including lifetime/ageing feedback;
- zero-sequence/common-mode grounding and parasitic-capacitance paths;
- frequency-dependent phase/neutral magnetic models, plus dead time,
  discontinuous PWM, interleaved carriers, and overmodulation;
- measured broadband battery/DC-stage impedance, source-control interaction,
  and validation beyond the implemented constant series R–L approximation;
- control-aware harmonic linearisation with positive/negative-sequence coupling;
- current-limited GFM equilibria linked to dynamic fault-ride-through models;
- multilevel and two-stage DC/DC-plus-DC/AC converters; and
- an explicit reconfigurable four-leg-plus-split-link topology before any soft
  open point work is attempted.
