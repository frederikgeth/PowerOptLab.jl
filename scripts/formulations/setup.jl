# Optional backends live in a separate environment; preserve the production pin.
using Pkg, TOML
length(ARGS) == 1 || error("Usage: julia setup.jl EMPTY_ENV_DIRECTORY")
repo = normpath(joinpath(@__DIR__,"..",".."))
destination = abspath(only(ARGS))
(!isdir(destination) || isempty(readdir(destination))) || error("Destination must be empty")
manifest = joinpath(repo,"Manifest.toml")
isfile(manifest) || error("First run julia --project=. scripts/instantiate_pinned.jl")
project = TOML.parsefile(joinpath(repo,"Project.toml"))
expected_revision = project["sources"]["BMOPFTools"]["rev"]
for key in ("name","uuid","authors","version","extensions")
    pop!(project,key,nothing)
end
mkpath(destination)
open(joinpath(destination,"Project.toml"),"w") do io
    TOML.print(io,project; sorted=true)
end
cp(manifest,joinpath(destination,"Manifest.toml"))
Pkg.activate(destination)
# Optional bridges need newer transitive MOI than the production manifest.
# DiffOpt 0.6.1 caps MOI below that requirement; use 0.6.2 only here.
# Preserve other direct dependencies and the immutable BMOPFTools source.
Pkg.add([
    PackageSpec(name="DiffOpt",version="0.6.2"),
    PackageSpec(name="CCOpt",version="0.1.0"),
    PackageSpec(name="MathOptComplements",version="0.1.1"),
    PackageSpec(name="NLPModelsJuMP",version="0.13.6"),
    PackageSpec(name="PiecewiseLinearOpt",version="0.4.2"),
    PackageSpec(name="HiGHS",version="1.20.1"),
]; preserve=Pkg.PRESERVE_DIRECT)
Pkg.develop(path=repo; preserve=Pkg.PRESERVE_ALL)
info = only(i for i in values(Pkg.dependencies()) if i.name == "BMOPFTools")
info.is_tracking_repo && info.git_revision == expected_revision || error("BMOPFTools pin changed")
println("Optional formulation environment: ",destination)
