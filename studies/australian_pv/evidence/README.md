# Australian PV evidence and provenance

The **coupled inductive thermal campaign** is the current model-validation
campaign. It compares local reactive-priority mean/worst/sequence/four-leg
control with coordinated export with Q free, using identical fixed hardware.
Older tables in this study are historical experiments, not current-source
validation of the updated thermal allocator.

Full evidence is stored as immutable, content-addressed bundles. Each manifest
lists every original file's SHA256 and size, the bundle SHA256, and a full Git
commit identifying the original bytes. The bundle cache is local and ignored;
the same bundle can be reconstructed from the pinned Git tree without a remote
artifact service. Commits remain in this branch's ancestry after compaction.
No source manifest is re-signed to make an old result appear current.

| Manifest | Status | Contents |
|---|---|---|
| `current.json` | Current, added after regeneration | Complete benchmark, domain/boundary records and matched EMT, including failures and raw statuses |
| `historical-d2c9a97.json` | Pinned historical archive | All 546 pre-update result files, including four night thermal-prediction misses, the unsuccessful solver control and earlier failed hardware controls |

From a full-history checkout, restore current records into their canonical
paths before running saved-evidence tests:

```sh
python scripts/australian_pv/evidence.py restore studies/australian_pv/evidence/current.json
julia --project=. scripts/australian_pv/test_engineering_benchmark.jl
python scripts/australian_pv/test_engineering_benchmark_emt.py
python scripts/australian_pv/test_plot_engineering_benchmark.py
```

Restore historical results into a separate directory; they deliberately differ
from the current source and must not overwrite it:

```sh
python scripts/australian_pv/evidence.py restore studies/australian_pv/evidence/historical-d2c9a97.json --output /tmp/australian-pv-history
```

If a shallow checkout lacks the pinned commit, fetch the branch history first.
Restoration verifies the entire bundle and destination before writing any files.
It refuses conflicting files, unsafe paths and checksum failures. Seven tests
exercise cached and Git-only restoration and corruption rejection.

For exact historical replay, check out the manifest's `source_revision` in a
separate worktree and use the original campaign's recorded environment and
source hashes. Auditing those points against current equations answers a
different question and should fail when sources differ. Current records retain
their generation `dependency_manifest.toml`; cross-version audits verify that
recorded file and report their separate runtime environment rather than silently
replacing the generation manifest.

Normal PR CI runs solver, independent SI record, corruption and circuit tests
on Julia 1.10 and current stable. The `Australian PV evidence` scheduled/manual
workflow regenerates eight coupled points and 32 representative EMT references
at three meshes in averaged and switched modes. It requires every thermal
prediction to enclose its switched reference. Hardware, raw solver, physical/
replay, mesh and prediction acceptance remain separate gates.
