# Australian PV evidence and provenance

Full raw records (saved phasors, solver stages, network audits, EMT waveforms)
are stored as immutable, content-addressed bundles in `bundles/`. Each manifest
lists every file's SHA256 and size, the bundle SHA256, and the Git commit the
bytes were packed from. The Git tree keeps only compact review files.

| Manifest | Status | Contents |
|---|---|---|
| `current.json` (= `current-05db2db.json`) | Current | Matched benchmark with the declared 7 A / 23 A fixture (15 points, PF/OPF agreement), legacy-fixture sensitivity, 60 matched EMT references, regenerated carrier-domain and DC-boundary validations |
| `current-35d5cac.json` | Superseded | First matched benchmark (1 A capacitor budget, 8 points, 32 EMT references) |
| `historical-d2c9a97.json` | Historical | 546 files from all earlier campaigns, including failed predictions and solver controls (see [study history](../HISTORY.md)) |

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

## Restore superseded or historical records

Restore into a separate directory; they deliberately differ from the current
source and must not overwrite it:

```sh
python scripts/australian_pv/evidence.py restore studies/australian_pv/evidence/historical-d2c9a97.json --output /tmp/australian-pv-history
```

Reproducing an old campaign needs the source bytes matching its recorded hashes
(each campaign's `environment.json`), not just the archive commit.
`dependency_manifest-b61f2440.toml` is the dependency lock for the historical
records; current records carry their own `dependency_manifest.toml`.

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
conflicting files, unsafe paths and checksum failures. If a distributed bundle
is missing, restoration falls back to the local cache, then to the pinned Git
history.

**Open decision (roadmap D1):** the bundles are about 34 MB of binary in the
repository. Moving them to a GitHub Release or Zenodo, keeping only manifests
and a download step, would keep the repository small.
