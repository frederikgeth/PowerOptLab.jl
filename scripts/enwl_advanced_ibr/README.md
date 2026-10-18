# ENWL advanced-IBR showcase

For new Australian engineering studies, use the [Australian LV study](../../studies/australian_pv/README.md). This ENWL campaign remains available as the earlier synthetic-retrofit experiment.

These scripts use the 150 BMOPF snapshots under the sibling
`BMOPFDraftData/benchmarks/ENWLsnapshots` checkout as inputs to PowerOptLab's
network-scale, phase-aware advanced-inverter study runner.

Start by checking the available cases:

```sh
julia --project=. scripts/enwl_advanced_ibr/inventory.jl
```

Build a small campaign without solving it:

```sh
julia --project=. scripts/enwl_advanced_ibr/run_campaign.jl --dry-run=true
```

Run the default matched experiment (one midday 30-bus snapshot, three selected
PVs, and three controller variants):

```sh
julia --project=. scripts/enwl_advanced_ibr/run_campaign.jl
```

Scale along either axis. `--devices=0` selects every PV and
`--snapshots=all` selects all 25 half-hour snapshots:

```sh
julia --project=. scripts/enwl_advanced_ibr/run_campaign.jl \
  --feeder=538bus_LN --snapshots=t09_1200,t13_1400,t17_1600 --devices=24
```

The runner is deliberately serial and deterministic. For a large unattended
campaign, partition snapshot lists across Julia processes and give each process
a distinct `--output` directory.

## Interpretation boundary

The source data contains only `SINGLE_PHASE` PV records, whereas
`solve_controlled_inverter_fleet` is a three-phase phase-aware controller study.
The adapter therefore creates a **synthetic retrofit** for each selected PV:

- same bus;
- balanced `THREE_LEG` connection to phases a-b-c;
- original total VA rating and available watts preserved;
- hardware scaled in per unit of each device's rating: a 5% (X/R = 10) series
  filter, 50 µF/kVA DC link with a 1% 2ω ripple limit, 20% current headroom,
  and a loss curve of about 2.2% at rated output (all assumptions; see
  `docs/src/ibr/parameter_selection.md`);
- output-filter, DC-link, current, modulation, switching-polytope, and
  double-frequency ripple physics enabled;
- every unselected single-phase PV remains in the native BMOPFTools model.

This is useful for comparative controller/hardware stress studies, but it is
not a reproduction of the original single-phase PV operating point. The CSV
outputs retain the source path, original phase/terminal map, selected IDs, and
retrofit label in case metadata.

The AS/NZS 4777 curves are read from the snapshot's `control_profile`, and the
declared ordinate units are enforced rather than assumed: `p_ref: S_MAX` selects
`volt_watt_basis=:rated`, so a retrofitted plant curtails against the same base
as the native PVs solved beside it. The three-wire plants sense phase-to-neutral
through a neutral sensing lead (`sensing_reference="n"`), matching the
`PN_PER_PHASE` profiles of the `*_LN` feeders. Both references are recorded in
case metadata (`profile_voltage_reference`, `sensed_voltage_reference`).

The variants are:

- `baseline`: mean-phase Volt-var/Volt-watt, balanced current;
- `worst_phase`: worst-phase voltage guards, balanced current;
- `sequence_droop`: positive-sequence Volt-var/Volt-watt plus a
  negative-sequence admittance droop of 43 pu (the Australian study's measured
  design for 20% attenuation) at an assumed 20° angle, with power-priority
  limiting so unbalance support never curtails active power. The ENWL feeder
  impedances are not measured here, so treat the gain and angle as
  assumptions.

All variants use identical smoothing widths.

All numerical output is written in SI units to `results/enwl_advanced_ibr` by
default. A non-publishable campaign exits with status 2 after preserving its
case, device, phase, summary, and matched-pair evidence.
