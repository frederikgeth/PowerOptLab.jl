#!/usr/bin/env julia

# Execute the published code, rather than a second copy of the examples.
# Keep the default pages explicit so removing a page or its marker fails CI.
function run_study_tutorials(; root=dirname(@__DIR__), filenames=[
        "battery_storage_models.md", "constrained_nlls_state_estimation.md"])
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
        println("Passed $filename")
    end
    println("Passed $(length(filenames)) executable study tutorials")
    return length(filenames)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_study_tutorials()
end
