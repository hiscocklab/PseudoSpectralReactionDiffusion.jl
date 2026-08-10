using PseudoSpectralReactionDiffusion
using Documenter
using DocumenterMarkdown

makedocs(
    modules = [PseudoSpectralReactionDiffusion],
    format = Markdown(),
    pages = ["README.md"],
)