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

upstream = Pkg.PackageSpec(name="BMOPFTools", url=source["url"], rev=source["rev"])
docs_manifest = joinpath(root, "docs", "Manifest.toml")
if VERSION < v"1.11" && Base.active_project() == project_paths[2] && !isfile(docs_manifest)
    # Julia 1.10 ignores [sources]. Resolve the registered docs dependencies
    # (especially Documenter) together with the pinned upstream before adding
    # the local package. Seeding only the two unregistered packages can freeze
    # a transitive version that Documenter cannot use when Pkg.resolve runs.
    mktempdir() do bootstrap
        try
            bootstrap_project = TOML.parsefile(project_paths[2])
            delete!(bootstrap_project["deps"], "PowerOptLab")
            delete!(bootstrap_project, "sources")
            open(joinpath(bootstrap, "Project.toml"), "w") do io
                TOML.print(io, bootstrap_project)
            end
            Pkg.activate(bootstrap)
            Pkg.add(upstream)
            Pkg.develop(Pkg.PackageSpec(path=root))
            manifest = TOML.parsefile(joinpath(bootstrap, "Manifest.toml"))
            for entries in values(manifest["deps"]), entry in entries
                haskey(entry, "path") || continue
                entry["path"] = abspath(bootstrap, entry["path"])
            end
            open(docs_manifest, "w") do io
                TOML.print(io, manifest)
            end
        finally
            Pkg.activate(joinpath(root, "docs"))
        end
    end
    Pkg.resolve()
end
Pkg.add(upstream)
Pkg.instantiate()
info = only(info for info in values(Pkg.dependencies()) if info.name == "BMOPFTools")
info.is_tracking_repo && info.git_revision == source["rev"] ||
    error("Resolved BMOPFTools does not match the requested source commit")
println("BMOPFTools ", info.git_revision, " at ", info.source)
import BMOPFTools
isdefined(BMOPFTools, :optimize_opf!) ||
    error("The pinned BMOPFTools source must provide optimize_opf!")
