#!/usr/bin/env julia

# Use the same immutable upstream source locally, on Julia 1.10 CI (which does
# not consume [sources]), and in the documentation environment. Pkg.add also
# replaces an existing Pkg.develop path in an untracked local Manifest.toml.
using Pkg
using TOML

root = dirname(@__DIR__)
project_paths = [joinpath(root, "Project.toml"), joinpath(root, "docs", "Project.toml")]
sources = [TOML.parsefile(path)["sources"]["BMOPFTools"] for path in project_paths]
sources[1] == sources[2] || error("Root and docs BMOPFTools source pins differ")
source = first(sources)
occursin(r"^[0-9a-f]{40}$", source["rev"]) || error("BMOPFTools must use a full commit SHA")
Base.active_project() in project_paths || error("Run with --project=. or --project=docs")

Pkg.add(Pkg.PackageSpec(name="BMOPFTools", url=source["url"], rev=source["rev"]))
Pkg.instantiate()
info = only(info for info in values(Pkg.dependencies()) if info.name == "BMOPFTools")
info.is_tracking_repo && info.git_revision == source["rev"] ||
    error("Resolved BMOPFTools does not match the requested source commit")
println("BMOPFTools ", info.git_revision, " at ", info.source)
