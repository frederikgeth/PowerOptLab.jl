# Australian PV evidence and provenance

The **coupled inductive thermal campaign** is the current model-validation
campaign. It compares local reactive-priority mean/worst/sequence/four-leg
control with coordinated export with Q free, using identical fixed hardware.
Older tables in this study are historical experiments, not current-source
validation of the updated thermal allocator.

Full evidence is stored as immutable, content-addressed bundles. Each manifest
lists every original file's SHA256 and size, the bundle SHA256, and a full Git
commit identifying the original bytes. Compressed bundles in `bundles/` retain
those bytes through squash merges and branch deletion. An ignored local cache
and reconstruction from the pinned Git tree provide additional restore paths;
the distributed bundles do not require that intermediate Git history.
No source manifest is re-signed to make an old result appear current.

| Manifest | Status | Contents |
|---|---|---|
| `current.json` | Current | Complete benchmark, domain/boundary records and matched EMT, including failures and raw statuses |
| `historical-d2c9a97.json` | Pinned historical archive | All 546 pre-update result files, including four night thermal-prediction misses, the unsuccessful solver control and earlier failed hardware controls |

## Restore current records

From a checkout, restore current records into their canonical
paths before running saved-evidence tests:

```sh
python scripts/australian_pv/evidence.py restore studies/australian_pv/evidence/current.json
julia --project=. scripts/australian_pv/test_engineering_benchmark.jl
python scripts/australian_pv/test_engineering_benchmark_emt.py
python scripts/australian_pv/test_plot_engineering_benchmark.py
python scripts/australian_pv/summarize_evidence.py --check
```

## Restore historical records

Restore historical results into a separate directory; they deliberately differ
from the current source and must not overwrite it:

```sh
python scripts/australian_pv/evidence.py restore studies/australian_pv/evidence/historical-d2c9a97.json --output /tmp/australian-pv-history
```

The normal restore uses a distributed bundle or a verified local cache. If both
are unavailable, fetching the pinned Git history enables the reconstruction
fallback. Rebuilding compressed bytes also requires a compatible compression
implementation; distributing the original bundle avoids that dependency.
Restoration verifies the entire bundle and destination before writing any files.
It refuses conflicting files, unsafe paths and checksum failures. Archive tests
exercise cached and Git-only restoration and corruption rejection.

For historical replay, use each campaign's own generation `git_revision`,
`git_status`, Julia version and source hashes. The bundle's `source_revision`
identifies the archived evidence tree; it need not be the source revision that
generated every experiment inside it. Some historical preparations used dirty
working trees: reproducing them requires the source bytes matching their hashes,
not merely checking out the archive commit. Auditing old points against current
equations answers a different question and should fail when sources differ.

`dependency_manifest-b61f2440.toml` preserves the original dependency lock and
matches all 30 historical environment fingerprints. Use it with the recorded
Julia version. Current records also retain their generation
`dependency_manifest.toml` inside the full bundle. Supported-version validation
resolves separate Julia-compatible environments; replay audits verify the
recorded lock and report their separate runtime environment without replacing it.

Normal PR CI runs solver, independent SI record, corruption and circuit tests
on Julia 1.10 and current stable. The `Australian PV evidence` scheduled/manual
workflow regenerates eight coupled points and 32 representative EMT references
at three meshes in averaged and switched modes. It requires every thermal
prediction to enclose its switched reference. Hardware, raw solver, physical/
replay, mesh and prediction acceptance remain separate gates.

`current-35d5cac.json` is the immutable manifest for this campaign;
`current.json` is its active copy. For a future campaign, create a new versioned
manifest and retain prior bundles/manifests before updating the active copy.
Regenerate data when source hashes change; do not edit old provenance to make
it appear current. `summarize_evidence.py --check` independently replays the
audits in temporary directories and checks the compact fixtures and summary
without changing original receipts.
