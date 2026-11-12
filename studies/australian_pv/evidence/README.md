# Australian PV evidence and provenance

Full raw records (saved phasors, solver stages, network audits, EMT waveforms)
are stored as immutable, content-addressed bundles in `bundles/`. Each manifest
lists every file's SHA256 and size, the bundle SHA256, and the Git commit the
bytes were packed from. The Git tree keeps only compact review files.

| Manifest | Status | Contents |
|---|---|---|
| `current.json` (= `current-25327e0.json`) | Current | Matched benchmark with the declared 7 A / 23 A fixture (15 points, PF/OPF agreement), legacy-fixture sensitivity, 60 matched EMT references, regenerated carrier-domain and DC-boundary validations |

## Restore current records

From the repository root, before running saved-record tests or the report
generators:

```sh
python scripts/australian_pv/evidence.py restore studies/australian_pv/evidence/current.json
julia --project=. scripts/australian_pv/test_engineering_benchmark.jl
python scripts/australian_pv/test_engineering_benchmark_emt.py
python scripts/australian_pv/test_plot_engineering_benchmark.py
python scripts/australian_pv/summarize_evidence.py --check
```

## Earlier campaigns

Records of earlier campaigns are kept by the Git tag `evidence/historical-d2c9a97`
rather than a bundle; see [study history](../HISTORY.md) for how to check them
out. The four night-time thermal-prediction misses are extracted in
`results/historical_night_thermal_misses.csv`.

## How a campaign is sealed

1. Generate from a clean commit (`run_engineering_benchmark.jl`, audit, EMT,
   domain, boundary).
2. Commit the full records once and pack them:
   `python scripts/australian_pv/evidence.py pack <commit> studies/australian_pv/evidence/current-<short>.json <roots...>`.
3. Copy the new manifest to `current.json`, then untrack the raw files so only
   compact files remain in the tree. A squash merge keeps only the bundle.

Saved-record audits compare source hashes: editing any hashed script makes the
records stale until the campaign is regenerated. That is deliberate (no
re-signing of old results), but it means any change to the shared study
scripts requires a rerun (about 15 minutes on a laptop for everything).

Restoration verifies the bundle and destination before writing, and refuses
conflicting files, unsafe paths and checksum failures. If the distributed bundle
is missing, restoration falls back to the local cache, then to the pinned Git
commit.

**Open decision (roadmap D1):** the one remaining bundle is 15 MB of binary in
the repository. Moving it to a GitHub Release or Zenodo, keeping only the
manifest and a download step, would keep the repository small.
