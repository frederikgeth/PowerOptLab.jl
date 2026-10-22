# BMOPFTools 0.1.x compatibility adapter.
#
# HELM uses the public augmented-admittance API, but also needs the parsed
# constant-power and constant-impedance sub-loads used by BMOPFTools' own
# linearized Y-bus. The pinned upstream release has no public read-only seam for
# that decomposition. Keep every private import in this file so an upstream
# change fails in one obvious place and the adapter can be deleted once the
# public load API lands.
using BMOPFTools: _Node, _SubLoad, _load_subloads, _subload_S, _subload_S_nz,
    _subload_yz,
    _stamp_pair!, _neutral_terminal, _neutral_labels, _neutral_pos,
    _phase_positions, _DEFAULT_CONFIG

# HELM additionally needs the ideal-transformer classification used by
# `ybus_augmented`, so it can restore the shunt-only part (magnetising branch and
# neutral grounding) of the transformers that `ybus_augmented` promotes to
# constraint rows or node aliases, and refuse zero-leakage units it cannot
# promote (`_nw_zb_matrix` is the n-winding leakage matrix `nwinding_yprim`
# inverts). Kept separate from the operability adapter fingerprint below: these
# imports do not widen the operability closure.
using BMOPFTools: _xfmr_is_ideal, _xfmr_turns_ratio, _domain_thresholds,
    _ybus_xfmr_z_min, _nw_zb_matrix

const _HELM_UPSTREAM_PRIVATE_IMPORTS = (
    "_xfmr_is_ideal", "_xfmr_turns_ratio", "_domain_thresholds",
    "_ybus_xfmr_z_min", "_nw_zb_matrix")

const _OPERABILITY_UPSTREAM_ADAPTER_VERSION =
    "bmopftools-0.1.0-private-load-decomposition/v1"
const _OPERABILITY_UPSTREAM_PRIVATE_IMPORTS = (
    "_Node", "_SubLoad", "_load_subloads", "_subload_S", "_subload_S_nz",
    "_subload_yz", "_stamp_pair!", "_neutral_terminal", "_neutral_labels",
    "_neutral_pos", "_phase_positions", "_DEFAULT_CONFIG")

"""
    operability_upstream_audit()

Return the compatibility boundary used by the native operability and HELM
implementations. The record is intentionally explicit: the pinned upstream
release provides the public `ybus_linearized` equilibrium seam, while the
connection-aware load decomposition remains private and is isolated here.
The record is provenance, not a claim that generator, IBR, or controller
equations are included in the native static closure.
"""
function operability_upstream_audit()
    upstream_version = string(Base.pkgversion(BMOPFTools))
    private_imports = collect(_OPERABILITY_UPSTREAM_PRIVATE_IMPORTS)
    fingerprint = join(("BMOPFTools", upstream_version,
                        _OPERABILITY_UPSTREAM_ADAPTER_VERSION,
                        join(private_imports, ",")), "|")
    public_equilibrium_seam_available = isdefined(BMOPFTools, :ybus_linearized)
    private_imports_available = all(isdefined(BMOPFTools, Symbol(name))
                                    for name in private_imports)
    private_imports_bound = all(isdefined(@__MODULE__, Symbol(name))
                                for name in private_imports)
    private_imports_match_upstream = private_imports_bound && all(
        getfield(@__MODULE__, Symbol(name)) ===
        getfield(BMOPFTools, Symbol(name)) for name in private_imports)
    Dict{String,Any}(
        "upstream_package" => "BMOPFTools",
        "upstream_version" => upstream_version,
        "adapter_version" => _OPERABILITY_UPSTREAM_ADAPTER_VERSION,
        "adapter_fingerprint" => fingerprint,
        "public_equilibrium_seam" => "BMOPFTools.ybus_linearized",
        "public_equilibrium_seam_available" => public_equilibrium_seam_available,
        "private_imports" => private_imports,
        "private_imports_available" => private_imports_available,
        "private_imports_bound" => private_imports_bound,
        "private_imports_match_upstream" => private_imports_match_upstream,
        "compatibility_status" => public_equilibrium_seam_available &&
            private_imports_available && private_imports_match_upstream ?
            :supported : :not_applicable,
        "private_load_decomposition" => true,
        "generator_ibr_controller_residual_seam" => false,
        "replacement_plan" =>
            "replace the single private boundary with a public equilibrium/device residual seam before extending closure scope",
    )
end
