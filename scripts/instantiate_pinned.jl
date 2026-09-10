#!/usr/bin/env julia

# Use the same immutable upstream source locally, on Julia 1.10 CI (which does
# not consume [sources]), and in the documentation environment. Pkg.add also
# replaces an existing Pkg.develop path in an untracked local Manifest.toml.
using Pkg
using TOML

root = dirname(@__DIR__)
project_paths = [joinpath(root, "Project.toml"), joinpath(root, "docs", "Project.toml")]
names = ["BMOPFTools", "FormulationLab"]
projects = TOML.parsefile.(project_paths)
Base.active_project() in project_paths || error("Run with --project=. or --project=docs")
specs = Pkg.PackageSpec[]
for name in names
    sources = [project["sources"][name] for project in projects]
    sources[1] == sources[2] || error("Root and docs $name source pins differ")
    source = first(sources)
    occursin(r"^[0-9a-f]{40}$", source["rev"]) || error("$name must use a full commit SHA")
    push!(specs, Pkg.PackageSpec(name=name, url=source["url"], rev=source["rev"]))
end

# Resolve both unregistered packages together, including on Julia 1.10.
Pkg.add(specs)
Pkg.instantiate()
for name in names
    source = projects[1]["sources"][name]
    info = only(info for info in values(Pkg.dependencies()) if info.name == name)
    info.is_tracking_repo && info.git_revision == source["rev"] ||
        error("Resolved $name does not match the requested source commit")
    println(name, " ", info.git_revision, " at ", info.source)
end
