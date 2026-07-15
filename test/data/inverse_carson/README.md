# Inverse Carson validation data

These files are versioned validation artifacts, not generic example linecodes.
Every benchmark records the conventions that materially affect its meaning:
frequency, soil resistivity, earth-return model, transposition, neutral
elimination, symmetrical-component ordering, units, and uncertainty status.

- `paper_table_iv.toml` contains the five overhead cases printed in Table IV of
  Tam, Geth & Mithulananthan, arXiv:2404.08210v2. Conductor inputs are reconstructed
  from Tables II and V. The paper does not explicitly report frequency or soil
  resistivity; their inferred status is retained in every case.
- `opendss_mars_triangle.dss` is the independent executable input.
  `opendss_mars_triangle.toml` freezes the resulting primitive and sequence data
  from DSS C-API 0.14.3/OpenDSS SVN 3723.

Do not replace a missing measurement uncertainty with the numerical
reproducibility tolerance of a deterministic forward calculation. Paper
rounding, catalog uncertainty, and earth-model discrepancy are separate error
sources and are labelled separately in the artifacts.

The tests rerun the DSS input when OpenDSSDirect is available. A changed DSS
result should be reviewed as a benchmark-version change, not accepted by merely
loosening the tolerance.
