#!/usr/bin/env julia

# Run the exact Julia fences of every DOE page marked doe-executable. Each page
# has a fresh module, matching a fresh REPL with repository-root working directory.
# Runtime failures and assertion failures fail CI; no copied tutorial scripts.
function run_doe_tutorials(; root=dirname(@__DIR__))
    directory = joinpath(root, "docs", "src", "tutorials")
    count = 0
    for filename in sort(readdir(directory))
        endswith(filename, ".md") || continue
        path = joinpath(directory, filename)
        source = read(path, String)
        occursin("<!-- doe-executable -->", source) || continue
        sandbox = Module(gensym(:DOETutorial))
        Core.eval(sandbox, :(include(path) = Base.include(@__MODULE__, abspath(path))))
        fences = collect(eachmatch(r"(?ms)^```julia[ \t]*\n(.*?)^```[ \t]*$", source))
        isempty(fences) && error("executable tutorial has no Julia fences: $filename")
        println("Running $filename ($(length(fences)) Julia blocks)")
        cd(root) do
            for fence in fences
                line = countlines(IOBuffer(source[1:prevind(source, fence.offset)])) + 1
                Base.include_string(sandbox, repeat("\n", line) * fence.captures[1], path)
            end
        end
        count += 1
    end
    count > 0 || error("no executable DOE tutorials found")
    println("Passed $count executable DOE tutorials")
    return count
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_doe_tutorials()
end
