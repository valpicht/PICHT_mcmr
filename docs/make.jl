using MCMRSimulator
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(MCMRSimulator, :DocTestSetup, :(using MCMRSimulator); recursive=true)
bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), style=:authoryear)

remote = Remotes.GitLab("git.fmrib.ox.ac.uk", "ndcn0236", "MCMRSimulator.jl")

raw_src_files = ["installation.md"]
install_version = "CI_COMMIT_TAG" in keys(ENV) ? ("#" * ENV["CI_COMMIT_TAG"]) : ""
for fn in raw_src_files
    text = read("docs/raw_src/$fn", String)
    text = replace(text, "{install_version}" => install_version)
    write("docs/src/$fn", text)
end

makedocs(;
    modules=[MCMRSimulator],
    authors="Michiel Cottaar <Michiel.Cottaar@ndcn.ox.ac.uk>",
    repo=remote,
    sitename="MCMRSimulator.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        size_threshold_ignore=["api.md"],
        description="Documentation for MCMRSimulator.jl: a Monte Carlo MRI simulator in Julia",
        footer=nothing,
        canonical="https://open.win.ox.ac.uk/pages/ndcn0236/mcmrsimulator.jl/stable/",
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Tutorial (Julia)" => "tutorial_julia.md",
        "Tutorial (CLI)" => "tutorial_cli.md",
        "Geometry" => "geometry.md",
        "MRI/collision properties" => "properties.md",
        "API" => "api.md",
        "References" => "references.md",
    ],
    warnonly=[:missing_docs],
    checkdocs=:exports,
    plugins=[bib],
)
for fn in raw_src_files
    rm("docs/src/$fn")
end


if get(ENV, "CI_COMMIT_REF_NAME", "") == "main" || length(get(ENV, "CI_COMMIT_TAG", "")) > 0
    deploydocs(repo="git.fmrib.ox.ac.uk:ndcn0236/mcmrsimulator.jl.git", branch="pages", devbranch="main")
else
    println("Skipping deployment, because we are local or on a secondary branch.")
end
