using SetProg
using Documenter, Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

const EXAMPLES = readdir(EXAMPLES_DIR)

for example in EXAMPLES
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR)
    #Literate.notebook(example_filepath, OUTPUT_DIR)
    #Literate.script(example_filepath, OUTPUT_DIR)
end

makedocs(
    sitename = "SetProg",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Index" => "index.md",
        "Tutorials" => map(
            file -> joinpath("generated", file),
            filter(
                file -> endswith(file, ".md"),
                sort(readdir(OUTPUT_DIR)),
            ),
        ),
    ],
)

deploydocs(
    repo   = "github.com/blegat/SetProg.jl.git",
    push_preview = true,
)
