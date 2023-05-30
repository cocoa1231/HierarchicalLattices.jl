using HierarchicalLattices 
using Documenter

DocMeta.setdocmeta!(HierarchicalLattices, :DocTestSetup, :(using HierarchicalLattices); recursive=true)

makedocs(;
    modules=[HierarchicalLattices],
    authors="Cocoa <cocoathepenguin@protonmail.com> and contributors",
    repo="https://github.com/cocoa1231/HierarchicalLattices.jl/blob/{commit}{path}#{line}",
    sitename="Hierarchical Lattices",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cocoa1231.github.io/HierarchicalLattices.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => "apiref.md"
    ],
)

deploydocs(;
    repo="github.com/cocoa1231/HierarchicalLattices.jl",
    devbranch="main",
)
