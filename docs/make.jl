using Documenter, BlockHaloArrays

DocMeta.setdocmeta!(BlockHaloArrays, :DocTestSetup, :(using BlockHaloArrays); recursive=true)

makedocs(
    source = "src",
    build = "build",
    clean = true,
    doctest = true,
    repo = "https://github.com/smillerc/BlockHaloArrays.jl",
    highlightsig = true,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    expandfirst = [],
    modules = [BlockHaloArrays],
    pages = [
        "BlockHaloArrays.jl" => "index.md",
        # "Examples" => example_md_files,
        "Reference" => [
            # "MPIHaloArray" => "mpihaloarray.md",
            # "Topology" => "topology.md",
            "API" => "api.md"
        ]
    ],
    sitename = "BlockHaloArrays.jl"
)

deploydocs(repo = "github.com/smillerc/BlockHaloArrays.jl.git",
           branch = "gh-pages",
           devbranch = "main",)
