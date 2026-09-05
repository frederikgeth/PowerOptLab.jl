# Optional diagnostic tools only: leave PowerOptLab's dependency environment intact.
using Pkg, TOML

repo = normpath(joinpath(@__DIR__, "..", ".."))
length(ARGS) == 1 || error("Usage: julia setup_controller_convergence.jl EMPTY_ENV_DIRECTORY")
destination = abspath(only(ARGS))
(!isdir(destination) || isempty(readdir(destination))) || error("Destination must be empty")
manifest = joinpath(repo, "Manifest.toml")
isfile(manifest) || error("First run julia --project=. scripts/instantiate_pinned.jl")
project = TOML.parsefile(joinpath(repo, "Project.toml"))
expected_revision = project["sources"]["BMOPFTools"]["rev"]
for key in ("name", "uuid", "authors", "version")
    pop!(project, key, nothing)
end
mkpath(destination)
open(joinpath(destination, "Project.toml"), "w") do io
    TOML.print(io, project; sorted=true)
end
cp(manifest, joinpath(destination, "Manifest.toml"))
Pkg.activate(destination)
Pkg.add(PackageSpec(name="MadNLP", version="0.10.1"); preserve=Pkg.PRESERVE_ALL)
Pkg.develop(path=repo; preserve=Pkg.PRESERVE_ALL)
bmopf = only(info for info in values(Pkg.dependencies()) if info.name == "BMOPFTools")
bmopf.git_revision == expected_revision || error("BMOPFTools revision does not match the repository pin")
println("Diagnostic environment: ", destination)
