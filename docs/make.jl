using PseudoSpectral
using Documenter
using DocumenterMarkdown

makedocs(
    modules = [PseudoSpectral],
    format = Markdown(),
    pages = ["README.md"],
)

cp(
    joinpath(@__DIR__, "build", "README.md"),
    joinpath(@__DIR__, "..", "README.md"),
    force=true,
)