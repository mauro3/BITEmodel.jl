using Documenter, BITEModel
# This does not work:
# include("../src/plotting.jl")

makedocs(
    format = :html,
    sitename = "BITE-Model",
    # pages also make the side-bar
    pages = Any[
        "index.md",
        "data-loading.md",
        "forward.md",
        "inverse.md",
        "postproc.md",
        "api.md",
        "index-index.md"]
)
