# ENWL advanced-IBR showcase

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
- output-filter, DC-link, current, modulation, switching-polytope, and
  double-frequency ripple physics enabled;
- every unselected single-phase PV remains in the native BMOPFTools model.

This is useful for comparative controller/hardware stress studies, but it is
not a reproduction of the original single-phase PV operating point. The CSV
outputs retain the source path, original phase/terminal map, selected IDs, and
retrofit label in case metadata.

The variants are:

- `baseline`: mean-phase Volt-var/Volt-watt, balanced current;
- `worst_phase`: worst-phase voltage guards, balanced current;
- `sequence_droop`: positive-sequence Volt-var/Volt-watt plus negative-sequence
  admittance droop with partial 2ω-ripple cancellation.

All numerical output is written in SI units to `results/enwl_advanced_ibr` by
default. A non-publishable campaign exits with status 2 after preserving its
case, device, phase, summary, and matched-pair evidence.
