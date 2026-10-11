#!/usr/bin/env julia

# Execute the published code, rather than a second copy of the examples.
# Keep the default pages explicit so removing a page or its marker fails CI.
function run_study_tutorials(; root=dirname(@__DIR__), filenames=[
        "battery_storage_models.md", "constrained_nlls_state_estimation.md",
        "ibr_topology_under_unbalance.md", "ibr_dc_source_and_split_link.md",
        "ibr_ac_harmonics_lcl.md", "helm_vs_nonlinear_power_flow.md",
        "post_opf_operability_study.md"], on_page=nothing)
    isempty(filenames) && error("no study tutorials requested")
    for filename in filenames
        path = joinpath(root, "docs", "src", "tutorials", filename)
        source = read(path, String)
        occursin("<!-- study-executable -->", source) ||
            error("missing executable-study marker: $filename")
        fences = collect(eachmatch(r"(?ms)^```julia[ \t]*\n(.*?)^```[ \t]*$", source))
        isempty(fences) && error("executable study has no Julia blocks: $filename")
        sandbox = Module(gensym(:StudyTutorial))
        Core.eval(sandbox, :(include(path) = Base.include(@__MODULE__, abspath(path))))
        println("Running $filename ($(length(fences)) Julia blocks)")
        cd(root) do
            for fence in fences
                line = countlines(IOBuffer(source[1:prevind(source, fence.offset)])) + 1
                Base.include_string(sandbox, repeat("\n", line) * fence.captures[1], path)
            end
        end
        isnothing(on_page) || Base.invokelatest(on_page, filename, sandbox)
        println("Passed $filename")
    end
    println("Passed $(length(filenames)) executable study tutorials")
    return length(filenames)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        run_study_tutorials()
    elseif length(ARGS) == 2 && ARGS[1] == "--figure-data"
        include(joinpath(@__DIR__, "study_figure_data.jl"))
        run_study_tutorials(on_page=study_figure_exporter(abspath(ARGS[2])))
    else
        error("usage: run_study_tutorials.jl [--figure-data DIRECTORY]")
    end
end
