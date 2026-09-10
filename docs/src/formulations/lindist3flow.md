# LinDist3Flow

LinDist3Flow has moved to [FormulationLab.jl](https://github.com/frederikgeth/FormulationLab.jl).
Its formulation, component documentation, unit tests, and IEEE/OpenDSS oracles
are maintained there. PowerOptLab retains thin forwarding entry points for
existing studies and supplies BMOPFTools nonlinear replay through an explicit
callback. FormulationLab has no dependency on PowerOptLab or BMOPFTools.

Run `julia --project=. scripts/instantiate_pinned.jl` to install the immutable
FormulationLab and BMOPFTools commits declared in `Project.toml`. For local
FormulationLab development, use `Pkg.develop(path="../FormulationLab.jl")` afterward.

```julia
using PowerOptLab, Clarabel
result = solve_l3f_opf(network, Clarabel.Optimizer;
    options=L3FOptions(validate_nonlinear=true))
```

The result is a `FormulationLab.L3FResult`; `PowerOptLab.solve_status` and
`solve_diagnostics` continue to work. It is an approximation, not an AC lower
bound or a physical-feasibility certificate.
