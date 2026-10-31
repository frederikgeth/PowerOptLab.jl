# Review data and full evidence

The current coupled thermal campaign is summarized in
[`engineering_benchmark/cases.csv`](engineering_benchmark/cases.csv),
[`devices.csv`](engineering_benchmark/devices.csv), and the
[compact evidence summary](../evidence/current_summary.json). Its
[report and teaching figures](../engineering_benchmark.md) describe the service
tradeoff and model scope. Hardware sizing is unchanged.

Other retained case/summary files and metadata belong to the historical
experiments identified in the [evidence index](../evidence/README.md).
They have not been relabeled as current-source validation.

Full raw records, failed outcomes, phasors and waveforms are distributed in
checksum-verified compressed bundles. Restore current records from the repository
root before running saved-record tests or publication scripts:

```sh
python scripts/australian_pv/evidence.py restore studies/australian_pv/evidence/current.json
```

Generated payloads remain ignored locally after restoration. The Git review tree
retains compact tables and fixtures; the bundles preserve the complete result
bytes. Historical records should be restored into a separate directory using
the evidence index instructions.
