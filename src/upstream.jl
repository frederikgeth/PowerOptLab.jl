# BMOPFTools 0.1.0 compatibility adapter.
#
# HELM uses the public augmented-admittance API, but also needs the parsed
# constant-power and constant-impedance sub-loads used by BMOPFTools' own
# linearized Y-bus. The pinned upstream release has no public read-only seam for
# that decomposition. Keep every private import in this file so an upstream
# change fails in one obvious place and the adapter can be deleted once the
# public load API lands.
using BMOPFTools: _Node, _SubLoad, _load_subloads, _subload_S, _subload_yz,
    _stamp_pair!, _neutral_terminal, _neutral_labels, _DEFAULT_CONFIG
