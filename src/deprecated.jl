# Deprecated names kept for one release.

"""
    kron_reduce_bmopf(input; neutral_terminals=nothing, as_json=false,
                      include_coordinates=false)

Deprecated alias of `kron_reduce_neutrals`, the explicit-neutral Kron
reduction that now lives in BMOPFTools (same arguments and behavior). See the
BMOPFTools [network reduction guide](https://frederikgeth.github.io/BMOPFTools.jl/dev/network_reduction/).
"""
function kron_reduce_bmopf(input; kwargs...)
    Base.depwarn("`kron_reduce_bmopf` is deprecated; use `kron_reduce_neutrals` " *
                 "(provided by BMOPFTools and re-exported by PowerOptLab).", :kron_reduce_bmopf)
    kron_reduce_neutrals(input; kwargs...)
end

"""
    reduce_bmopf_neutrals(input; kwargs...)

Deprecated alias of `kron_reduce_neutrals` (BMOPFTools).
"""
function reduce_bmopf_neutrals(input; kwargs...)
    Base.depwarn("`reduce_bmopf_neutrals` is deprecated; use `kron_reduce_neutrals` " *
                 "(provided by BMOPFTools and re-exported by PowerOptLab).", :reduce_bmopf_neutrals)
    kron_reduce_neutrals(input; kwargs...)
end
