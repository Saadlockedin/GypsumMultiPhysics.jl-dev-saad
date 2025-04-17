using GypsumMultiPhysics
using Documenter

DocMeta.setdocmeta!(GypsumMultiPhysics, :DocTestSetup, :(using GypsumMultiPhysics); recursive=true)

makedocs(;
    modules=[GypsumMultiPhysics],
    authors="Hayri Sezer <sezerh24@gmail.com> and contributors",
    sitename="GypsumMultiPhysics.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
